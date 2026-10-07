"""End-to-end verification of the three Python-drivable sizing-field routes:
  (A) tetgen (pyvista/tetgen wheel) with a *background mesh* + per-node
      'target_size' scalar  ->  TetGen -m
  (B) tetgen with a global max volume  -a
  (C) gmsh with a Field stack (already verified in gmsh_bgmesh_demo.py)
Run with the bundled DSH python.
"""
import time
import numpy as np
import scipy.spatial

REG = (6.0 * np.sqrt(2.0)) ** (1.0 / 3.0)
R = 1.0


def ball_surface(n=26):
    import pyvista as pv
    return pv.Sphere(radius=R, theta_resolution=2 * n, phi_resolution=n,
                     start_theta=0, end_theta=360, start_phi=0, end_phi=180
                     ).triangulate()


def radial_bg_mesh(npts=1500, seed=3):
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(npts, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    r = rng.uniform(0, 1, (npts, 1)) ** (1 / 3)
    pts = np.vstack([d * r, np.zeros((1, 3))])
    tri = scipy.spatial.Delaunay(pts)
    keep = np.linalg.norm(pts[tri.simplices].mean(1), axis=1) < 1.0 - 1e-9
    tets = tri.simplices[keep]
    # h(r) = 0.05 + (0.28-0.05) * r      (linear grading centre -> rim)
    h = 0.05 + (0.28 - 0.05) * np.clip(np.linalg.norm(pts, axis=1), 0, 1)
    return pts.astype(np.float64), tets.astype(np.int32), h.astype(np.float64)


def stats(nodes, elems, label, secs):
    P = np.asarray(nodes, float)[np.asarray(elems, int)]
    vol = np.abs(np.einsum('ij,ij->i', P[:, 1] - P[:, 0],
                           np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))) / 6.0
    h = REG * vol ** (1 / 3)
    rad = np.linalg.norm(P.mean(1), axis=1)
    # quality: min dihedral angle
    print(f"--- {label}   [{secs:.2f} s]")
    print(f"    nodes={len(nodes)} tets={len(elems)} vol={vol.sum():.5f} "
          f"(exact {4/3*np.pi:.5f}, err {abs(vol.sum()-4/3*np.pi)/(4/3*np.pi)*100:.3f}%)")
    for lo, hi, lbl in [(0, .25, 'core  r<.25'), (.25, .5, 'mid  .25-.5'),
                        (.5, .75, 'mid  .5-.75'), (.75, 1.02, 'shell r>.75')]:
        m = (rad >= lo) & (rad < hi)
        if m.sum():
            print(f"      {lbl}: n={m.sum():5d} h_eff med={np.median(h[m]):.4f}")
    print(f"      h_eff global: min={h.min():.4f} med={np.median(h):.4f} max={h.max():.4f}")
    return h, vol


def min_dihedral(P):
    """min dihedral angle per tet (4 faces -> 4 dihedral angles)."""
    idx = [(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 3, 1), (1, 2, 3, 0)]
    out = []
    for a, b, c, d in idx:
        n1 = np.cross(P[:, b] - P[:, a], P[:, c] - P[:, a])
        n2 = np.cross(P[:, b] - P[:, a], P[:, d] - P[:, a])
        n1 /= np.linalg.norm(n1, axis=1, keepdims=True) + 1e-300
        n2 /= np.linalg.norm(n2, axis=1, keepdims=True) + 1e-300
        ang = np.degrees(np.arccos(np.clip(np.abs(np.einsum('ij,ij->i', n1, n2)), -1, 1)))
        # dihedral between two faces sharing edge (a,b) is pi - angle(n1,n2)
        out.append(180.0 - ang)
    return np.min(np.stack(out, 1), axis=1)


def main():
    import pyvista as pv
    import tetgen

    surf = ball_surface()
    print(f"surface: {surf.n_points} pts {surf.n_cells} tris\n")

    # ---------------- (A) global max volume, no sizing field -------------
    t0 = time.time()
    tg = tetgen.TetGen(surf)
    n, e, *_ = tg.tetrahedralize(maxvolume=0.004, minratio=1.5, quiet=True)
    stats(n, e, "A) tetgen maxvolume=0.004 (uniform, no field)", time.time() - t0)

    # ---------------- (B) background mesh + per-node target_size ---------
    bp, bt, bh = radial_bg_mesh()
    # pyvista UnstructuredGrid wants the VTK "legacy" cell array:
    #   [nverts, v0, v1, v2, v3, nverts, ...]
    cells = np.hstack([np.full((len(bt), 1), 4, np.int64), bt.astype(np.int64)]).ravel()
    bg = pv.UnstructuredGrid(cells, np.full(len(bt), pv.CellType.TETRA, np.uint8), bp)
    bg.point_data["target_size"] = bh
    t0 = time.time()
    tg2 = tetgen.TetGen(surf)
    try:
        n2, e2, *_ = tg2.tetrahedralize(bgmesh=bg, minratio=1.5, quiet=True)
        h2, v2 = stats(n2, e2, "B) tetgen bgmesh + target_size (graded .mtr)",
                       time.time() - t0)
        P = np.asarray(n2, float)[np.asarray(e2, int)]
        d = min_dihedral(P)
        print(f"      min dihedral: min={d.min():.2f} deg  p1={np.percentile(d,1):.2f} "
              f"med={np.median(d):.2f}  (tets < 5 deg: {(d<5).sum()})")
    except Exception as ex:
        print("B) FAILED:", type(ex).__name__, ex)

    # ---------------- (C) metric + background mesh, quality -------------
    t0 = time.time()
    tg3 = tetgen.TetGen(surf)
    try:
        n3, e3, *_ = tg3.tetrahedralize(bgmesh=bg, quality=True, minratio=1.4,
                                        mindihedral=10.0, quiet=True)
        h3, v3 = stats(n3, e3, "C) bgmesh + quality (minratio 1.4, mindihedral 10)",
                       time.time() - t0)
        P = np.asarray(n3, float)[np.asarray(e3, int)]
        d = min_dihedral(P)
        print(f"      min dihedral: min={d.min():.2f} deg  p1={np.percentile(d,1):.2f} "
              f"med={np.median(d):.2f}  (tets < 5 deg: {(d<5).sum()})")
    except Exception as ex:
        print("C) FAILED:", type(ex).__name__, ex)


if __name__ == "__main__":
    main()
