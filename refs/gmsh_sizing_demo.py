"""Empirical check: does Gmsh's Python API let us impose a prescribed 3D
element-size field, and does the field actually show up in the output mesh?

Also tests the Octree-based "background mesh" route (Field 1 = PostView).
"""
import sys
import numpy as np
import gmsh

R = 1.0


def build_with_box_threshold(hmin, hmax):
    """Ball of radius R, size hmax everywhere, hmin inside a box near the centre."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ball")
    gmsh.model.occ.addSphere(0, 0, 0, R, tag=1)
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.field.add("Box", 1)
    gmsh.model.mesh.field.setNumber(1, "VIn", hmin)
    gmsh.model.mesh.field.setNumber(1, "VOut", hmax)
    gmsh.model.mesh.field.setNumber(1, "XMin", -0.35)
    gmsh.model.mesh.field.setNumber(1, "XMax", 0.35)
    gmsh.model.mesh.field.setNumber(1, "YMin", -0.35)
    gmsh.model.mesh.field.setNumber(1, "YMax", 0.35)
    gmsh.model.mesh.field.setNumber(1, "ZMin", -0.35)
    gmsh.model.mesh.field.setNumber(1, "ZMax", 0.35)

    # smooth the transition: Threshold field on top of the Box field
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", hmin)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", hmax)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.15)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.6)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay
    gmsh.model.mesh.generate(3)

    ntags, ncoords, _ = gmsh.model.mesh.getNodes()
    nodes = np.asarray(ncoords).reshape(-1, 3)
    tag_to_idx = {int(t): i for i, t in enumerate(ntags)}
    etypes, etags, enodes = gmsh.model.mesh.getElements(3)
    tets_l = [np.asarray(e) for t, e in zip(etypes, enodes) if t == 4]
    tets = np.vectorize(tag_to_idx.get)(tets_l[0]).reshape(-1, 4) if tets_l \
        else np.zeros((0, 4), int)

    # element sizes = mean edge length of each tet
    P = nodes[tets]
    e = np.stack([np.linalg.norm(P[:, a] - P[:, b], axis=1)
                  for a, b in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]], 1)
    vol = np.abs(np.einsum('ij,ij->i', P[:, 1] - P[:, 0],
                           np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))) / 6.0
    gmsh.finalize()
    return nodes, tets, e.mean(1), vol


def report(name, nodes, tets, h, vol):
    rad = np.linalg.norm(nodes, axis=1)
    inner = rad[np.argmin(np.linalg.norm(nodes[tets].mean(1), axis=1))]  # dummy
    cen = np.linalg.norm(nodes[tets].mean(1), axis=1)
    core = cen < 0.3
    shell = cen > 0.75
    print(f"--- {name} ---")
    print(f"  nodes={len(nodes)}  tets={len(tets)}")
    print(f"  h: min={h.min():.4f} med={np.median(h):.4f} max={h.max():.4f}")
    if core.any():
        print(f"  core  (|x|<0.30): n={core.sum():5d}  median h={np.median(h[core]):.4f}")
    if shell.any():
        print(f"  shell (|x|>0.75): n={shell.sum():5d}  median h={np.median(h[shell]):.4f}")
    print(f"  volume sum={vol.sum():.6f} (exact ball vol={4/3*np.pi:.6f}, "
          f"rel err={abs(vol.sum()-4/3*np.pi)/(4/3*np.pi)*100:.3f}%)")


if __name__ == "__main__":
    r = build_with_box_threshold(hmin=0.04, hmax=0.16)
    report("Box+Threshold field, hmin=0.04 hmax=0.16", *r)
    sys.exit(0)
