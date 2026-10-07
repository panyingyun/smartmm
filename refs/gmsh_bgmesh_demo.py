"""Definitive experiment: drive a 3D tet mesh of a ball with BOTH
 (A) a Gmsh Field stack (Box -> Threshold -> MathEval), and
 (B) a genuine background mesh + per-node size scalar (the "background mesh"
     sizing-field route, the one we need for a volume-preserving map),
and measure the achieved element size in different regions.

Element "size" is reported as the median tet edge length multiplied by
(6*sqrt(2))**(1/3) ~= 1.7818... wait, for a REGULAR tet of edge L:
    V = L^3/(6*sqrt(2)),  so L = (6*sqrt(2)*V)^(1/3) = 2.0396 * V^(1/3).
We report  h_eff = (6*sqrt(2)*V_tet)^(1/3), i.e. the edge length the tet
*would* have if it were a regular tet of the same volume.  This is the
standard "volume-equivalent" sizing measure and is directly comparable to
Gmsh's MeshSize / TetGen's -a and .mtr entries.
"""
import sys
import numpy as np
import gmsh

R = 1.0
REG = (6.0 * np.sqrt(2.0)) ** (1.0 / 3.0)  # 2.0396


def extract(apply_box=None, apply_threshold=None,
            bg_nodes=None, bg_tets=None, bg_sizes=None,
            hmin=0.05, hmax=0.30, algorithm3d=1):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ball")
    gmsh.model.occ.addSphere(0, 0, 0, R, tag=1)
    gmsh.model.occ.synchronize()

    if apply_box is not None:
        gmsh.model.mesh.field.add("Box", 1)
        for k, v in apply_box.items():
            gmsh.model.mesh.field.setNumber(1, k, v)
        last = 1
        if apply_threshold is not None:
            gmsh.model.mesh.field.add("Threshold", 2)
            gmsh.model.mesh.field.setNumber(2, "InField", 1)
            for k, v in apply_threshold.items():
                gmsh.model.mesh.field.setNumber(2, k, v)
            last = 2
        gmsh.model.mesh.field.setAsBackgroundMesh(last)

    if bg_nodes is not None:
        # ---- genuine background-mesh sizing field -----------------------
        # Gmsh REFUSES a view built from the current mesh as a background mesh
        # ("Cannot use view based on current mesh for background mesh: you
        #  might want to use a list-based view (.pos file) instead").
        # So we must write a *list-based* .pos (SP = scalar point) view and
        # read it back.  This is exactly the route we need for a
        # volume-preserving map: the map's sample points + local target size.
        import os
        import tempfile
        pos = os.path.join(tempfile.gettempdir(), "dsh_bg_size.pos")
        with open(pos, "w") as f:
            f.write('View "bg_size" {\n')
            for (x, y, z), s in zip(bg_nodes, bg_sizes):
                f.write(f"SP({x:.10g},{y:.10g},{z:.10g}){{{s:.10g}}};\n")
            f.write("};\n")
        tag = gmsh.view.add("bg_size")   # ensure at least one view exists
        gmsh.merge(pos)                  # .pos -> new list-based view
        tags = gmsh.view.getTags()
        tag = tags[-1]
        print(f"    [bgmesh] wrote {pos}; views={tags}; using tag={tag}")
        gmsh.model.mesh.field.add("PostView", 10)
        gmsh.model.mesh.field.setNumber(10, "ViewTag", tag)
        gmsh.model.mesh.field.setAsBackgroundMesh(10)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", algorithm3d)
    gmsh.model.mesh.generate(3)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.asarray(ncoords).reshape(-1, 3)
    tag2i = {int(t): i for i, t in enumerate(ntags)}
    etypes, etags, enodes = gmsh.model.mesh.getElements(3)
    tl = [np.asarray(e) for t, e in zip(etypes, enodes) if t == 4]
    tets = np.vectorize(tag2i.get)(tl[0]).reshape(-1, 4) if tl else np.zeros((0, 4), int)
    gmsh.finalize()

    P = nodes[tets]
    vol = np.abs(np.einsum('ij,ij->i', P[:, 1] - P[:, 0],
                           np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))) / 6.0
    h_eff = REG * vol ** (1.0 / 3.0)
    rad = np.linalg.norm(nodes[tets].mean(1), axis=1)
    return nodes, tets, vol, h_eff, rad


def report(name, nodes, tets, vol, h, rad):
    print(f"--- {name}")
    print(f"    nodes={len(nodes):7d}  tets={len(tets):8d}  vol={vol.sum():.5f}"
          f"  (exact {4/3*np.pi:.5f}, err {abs(vol.sum()-4/3*np.pi)/(4/3*np.pi)*100:.3f}%)")
    for lo, hi, lbl in [(0.0, 0.25, "core  r<0.25"), (0.25, 0.5, "mid  .25-.5"),
                        (0.5, 0.75, "mid  .5-.75"), (0.75, 1.01, "shell r>0.75")]:
        m = (rad >= lo) & (rad < hi)
        if m.sum():
            print(f"      {lbl}: n={m.sum():6d}  h_eff med={np.median(h[m]):.4f}"
                  f"  p10={np.percentile(h[m],10):.4f} p90={np.percentile(h[m],90):.4f}")


if __name__ == "__main__":
    box = dict(VIn=0.05, VOut=0.30,
               XMin=-0.4, XMax=0.4, YMin=-0.4, YMax=0.4, ZMin=-0.4, ZMax=0.4)
    thr = dict(SizeMin=0.05, SizeMax=0.30, DistMin=0.1, DistMax=0.7)
    report("A) Box(VIn=.05 VOut=.30) -> Threshold",
           *extract(apply_box=box, apply_threshold=thr))

    # B) background mesh route: a coarse tet mesh of the whole ball whose nodes
    #    carry a radial size prescription.
    import scipy.spatial
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(900, 3))
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    pts *= R * rng.uniform(0, 1, size=(len(pts), 1)) ** (1 / 3)
    pts = np.vstack([pts, np.zeros((1, 3))])
    hull = scipy.spatial.Delaunay(pts)
    cents = pts[hull.simplices].mean(1)
    keep = np.linalg.norm(cents, axis=1) < R
    bg_tets = hull.simplices[keep]
    rc = np.linalg.norm(pts, axis=1)
    bg_sizes = 0.05 + (0.30 - 0.05) * np.clip(rc / R, 0, 1) ** 1.0
    nodes, tets, vol, h, rad = extract(bg_nodes=pts, bg_tets=bg_tets,
                                       bg_sizes=bg_sizes)
    report("B) background mesh (radial size 0.05 -> 0.30)", nodes, tets, vol, h, rad)
    sys.exit(0)
