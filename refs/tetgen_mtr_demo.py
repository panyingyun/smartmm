"""TetGen sizing-field routes, avoiding the crashing `bgmesh=` path of the
pyvista `tetgen` wrapper.  Uses `switches` + `bgmeshfilename` so that the
*actual* TetGen -m + background-mesh code path in the C++ library is exercised.
"""
import os
import time
import tempfile
import numpy as np
import scipy.spatial

REG = (6.0 * np.sqrt(2.0)) ** (1.0 / 3.0)
R = 1.0
W = os.path.join(tempfile.gettempdir(), "dsh_tetgen_sizing")
os.makedirs(W, exist_ok=True)


def ball_surface(n=26):
    import pyvista as pv
    return pv.Sphere(radius=R, theta_resolution=2 * n, phi_resolution=n)


def radial_bg(npts=2000, seed=3):
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(npts, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    r = rng.uniform(0, 1, (npts, 1)) ** (1 / 3)
    pts = np.vstack([d * r, np.zeros((1, 3))])
    tri = scipy.spatial.Delaunay(pts)
    keep = np.linalg.norm(pts[tri.simplices].mean(1), axis=1) < 1.0 - 1e-9
    tets = tri.simplices[keep]
    rad = np.clip(np.linalg.norm(pts, axis=1), 0, 1)
    h = 0.05 + (0.30 - 0.05) * rad           # linear grading 0.05 -> 0.30
    return pts, tets, h


def write_bg(prefix, pts, tets, h):
    """TetGen background mesh: <prefix>.node / .ele / .mtr"""
    with open(prefix + ".node", "w") as f:
        f.write(f"{len(pts)} 3 0 0\n")
        for i, p in enumerate(pts):
            f.write(f"{i} {p[0]:.10g} {p[1]:.10g} {p[2]:.10g}\n")
    with open(prefix + ".ele", "w") as f:
        f.write(f"{len(tets)} 4 0\n")
        for i, t in enumerate(tets):
            f.write(f"{i} {t[0]} {t[1]} {t[2]} {t[3]}\n")
    with open(prefix + ".mtr", "w") as f:
        f.write(f"{len(pts)} 1\n")
        for s in h:
            f.write(f"{s:.10g}\n")


def min_dihedral(P):
    idx = [(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 3, 1), (1, 2, 3, 0)]
    out = []
    for a, b, c, d in idx:
        n1 = np.cross(P[:, b] - P[:, a], P[:, c] - P[:, a])
        n2 = np.cross(P[:, b] - P[:, a], P[:, d] - P[:, a])
        n1 /= np.linalg.norm(n1, axis=1, keepdims=True) + 1e-300
        n2 /= np.linalg.norm(n2, axis=1, keepdims=True) + 1e-300
        ang = np.degrees(np.arccos(np.clip(np.abs(np.einsum('ij,ij->i', n1, n2)), -1, 1)))
        out.append(180.0 - ang)
    return np.min(np.stack(out, 1), axis=1)


def report(nodes, elems, label, secs):
    P = np.asarray(nodes, float)[np.asarray(elems, int)]
    vol = np.abs(np.einsum('ij,ij->i', P[:, 1] - P[:, 0],
                           np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))) / 6.0
    h = REG * vol ** (1 / 3)
    rad = np.linalg.norm(P.mean(1), axis=1)
    d = min_dihedral(P)
    print(f"--- {label}   [{secs:.2f} s]")
    print(f"    nodes={len(nodes):6d} tets={len(elems):7d} vol={vol.sum():.5f} "
          f"(exact {4/3*np.pi:.5f}; err {abs(vol.sum()-4/3*np.pi)/(4/3*np.pi)*100:.3f}%)")
    for lo, hi, lbl in [(0, .25, 'core r<.25'), (.25, .5, 'mid .25-.5'),
                        (.5, .75, 'mid .5-.75'), (.75, 1.02, 'shell r>.75')]:
        m = (rad >= lo) & (rad < hi)
        if m.sum():
            print(f"      {lbl}: n={m.sum():6d}  h_eff med={np.median(h[m]):.4f}"
                  f"   min-dihedral med={np.median(d[m]):.1f} deg")
    print(f"      h_eff: min={h.min():.4f} med={np.median(h):.4f} max={h.max():.4f}")
    print(f"      min dihedral: min={d.min():.2f} p1={np.percentile(d,1):.2f} "
          f"med={np.median(d):.2f} deg; tets<5deg: {(d<5).sum()} "
          f"({(d<5).mean()*100:.3f}%)")
    return h, d


def main():
    import warnings
    warnings.filterwarnings("ignore")
    import tetgen

    surf = ball_surface()
    print(f"surface: {surf.n_points} pts {surf.n_cells} tris\n")

    # (A) uniform, global max volume only
    t0 = time.time()
    n, e, *_ = tetgen.TetGen(surf).tetrahedralize(
        switches="pq1.4a0.002", quiet=True)
    report(n, e, "A) -pq1.4a0.002  (global max volume, uniform)", time.time() - t0)

    # (B) background mesh + .mtr  -> -m
    pts, tets, h = radial_bg()
    pre = os.path.join(W, "ball.b")
    write_bg(pre, pts, tets, h)
    print(f"\nbg mesh: {len(pts)} nodes {len(tets)} tets, h in "
          f"[{h.min():.3f},{h.max():.3f}] -> {pre}.node/.ele/.mtr\n")
    t0 = time.time()
    try:
        tg = tetgen.TetGen(surf)
        n2, e2, *_ = tg.tetrahedralize(switches="pq1.4m",
                                       bgmeshfilename=pre, quiet=True)
        report(n2, e2, "B) -pq1.4m + background mesh/.mtr (graded 0.05->0.30)",
               time.time() - t0)
    except Exception as ex:
        print("B) FAILED:", type(ex).__name__, ex)

    # (C) coarsen + field (the "adaptive remesh" mode)
    t0 = time.time()
    try:
        tg = tetgen.TetGen(surf)
        n3, e3, *_ = tg.tetrahedralize(switches="pq1.4mR", bgmeshfilename=pre,
                                       quiet=True)
        report(n3, e3, "C) -pq1.4mR (field + coarsening)", time.time() - t0)
    except Exception as ex:
        print("C) FAILED:", type(ex).__name__, ex)


if __name__ == "__main__":
    main()
