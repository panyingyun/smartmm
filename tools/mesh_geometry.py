"""Tetrahedralise a CAD geometry with Gmsh, with quality controls, and report quality.

    python tools/mesh_geometry.py model.stp --mesh-size 8 --out out/wheel
    python tools/mesh_geometry.py model.stp --mesh-size 8 --algo3d 4 --netgen --out out/wheel

Why the extra switches exist
----------------------------
The default Gmsh settings are tuned for "normal" solids.  On a thin-walled part with many
small features (fillet, chamfer, thread) they produce a large number of flat / sliver
tetrahedra, for two independent reasons:

1. **Thin walls.**  If the element size ``h`` is not clearly smaller than the wall
   thickness ``t``, a tetrahedron cannot fit inside the wall and becomes a flat sliver.
   Rule of thumb: ``h <= t / 2``.  ``--mesh-size`` must be chosen from a *thickness*
   measurement, not from the bounding box (this script reports the estimated thickness).
2. **Tiny features.**  Gmsh refines the surface at every small fillet and hole, and the
   volume mesh then has to blend those tiny elements into the interior -> slivers.  Either
   simplify the CAD first, or disable the automatic curvature refinement
   (``--no-curvature``) and let only ``--mesh-size`` drive the sizing.

``--netgen`` enables Gmsh's Netgen tetrahedron optimiser, which directly targets the
minimum dihedral angle and is usually the single biggest quality win.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tetparam.mesh import TetMesh  # noqa: E402

ALGO3D = {
    "delaunay": 1,
    "frontal": 4,
    "frontal-del": 5,
    "hxt": 10,
}


# --------------------------------------------------------------------------------------
def quality_report(m: TetMesh, label: str = "") -> dict:
    """Standard tetrahedron quality statistics + a wall-thickness estimate."""
    V = np.abs(m.signed_volumes())
    P = m.vertices[m.tets]
    E = np.array(
        [
            [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3],
        ]
    )
    edges = np.linalg.norm(P[:, E[:, 0]] - P[:, E[:, 1]], axis=2)
    lmin, lmax = edges.min(axis=1), edges.max(axis=1)

    bnd = m.boundary()
    tri = m.vertices[bnd["faces"]]
    area = 0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1
    )
    A = float(area.sum())
    vol = m.total_volume()
    thickness = 2.0 * vol / A if A > 0 else float("nan")

    q = m.quality()
    dih = np.asarray(q["min_dihedral"])
    # radius-edge ratio: the Delaunay-refinement quality measure (regular tet ~ 0.612)
    from tetparam.mesh import circumcenters

    R = np.linalg.norm(circumcenters(m.vertices, m.tets) - P[:, 0], axis=1)
    re_ratio = R / np.maximum(lmin, 1e-300)

    rep = {
        "label": label,
        "n_vertices": m.n_vertices,
        "n_tets": m.n_tets,
        "boundary_vertex_fraction": float(m.is_boundary_vertex.mean()),
        "volume": vol,
        "boundary_area": A,
        "estimated_thickness": thickness,
        "min_dihedral_mean": float(dih.mean()),
        "min_dihedral_p01": float(np.percentile(dih, 1)),
        "min_dihedral_min": float(dih.min()),
        "frac_dihedral_lt_10": float((dih < 10).mean()),
        "frac_dihedral_lt_5": float((dih < 5).mean()),
        "radius_edge_mean": float(re_ratio.mean()),
        "radius_edge_p99": float(np.percentile(re_ratio, 99)),
        "max_edge_over_thickness": float(np.median(lmax) / thickness)
        if thickness == thickness and thickness > 0
        else float("nan"),
    }
    print(f"\n--- {label} ---")
    print(f"  {m.n_vertices} vertices / {m.n_tets} tets  "
          f"(boundary vertices {rep['boundary_vertex_fraction']:.1%})")
    print(f"  volume {vol:.4g}   boundary area {A:.4g}   "
          f"estimated wall thickness t ~ 2V/A = {thickness:.4g}")
    print(f"  min dihedral: mean {rep['min_dihedral_mean']:.2f} deg  "
          f"p01 {rep['min_dihedral_p01']:.2f}  min {rep['min_dihedral_min']:.3f}")
    print(f"  frac(min dihedral < 10 deg) = {rep['frac_dihedral_lt_10']:.2%}   "
          f"(< 5 deg) = {rep['frac_dihedral_lt_5']:.2%}")
    print(f"  radius-edge ratio: mean {rep['radius_edge_mean']:.3f}  "
          f"p99 {rep['radius_edge_p99']:.3f}   (regular tet ~ 0.612)")
    print(f"  median longest edge / thickness = {rep['max_edge_over_thickness']:.3f}"
          f"   -> keep --mesh-size <= t/2 = {thickness / 2:.4g}")
    return rep


# --------------------------------------------------------------------------------------
def mesh_geometry(
    path: Path,
    mesh_size: float | None,
    out: Path,
    algo3d: str = "delaunay",
    optimize: int = 1,
    netgen: bool = True,
    curvature: float | None = None,
    extend_from_boundary: float = 0.0,
    size_min_factor: float = 0.25,
    threads: int = 1,
    heal: float = 0.0,
) -> Path:
    import gmsh

    out.mkdir(parents=True, exist_ok=True)
    msh = out / (path.stem + "_mesh.msh")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    if threads > 1:
        gmsh.option.setNumber("General.NumThreads", threads)
    gmsh.option.setNumber("Geometry.OCCImportLabels", 1)
    gmsh.option.setNumber("Geometry.Tolerance", 1e-6)

    # ---- heal the CAD *before* meshing: merge tiny edges/faces and sew shells.
    #      This is the only real cure for "thousands of tiny features force tiny
    #      boundary triangles", which is what makes the radius-edge ratio explode.
    if heal > 0:
        gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 1)
        gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 1)
        gmsh.option.setNumber("Geometry.OCCSewFaces", 1)
        gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
        gmsh.option.setNumber("Geometry.OCCMakeSolids", 1)

    gmsh.open(str(path))
    if heal > 0:
        try:
            gmsh.model.occ.healShapes(tolerance=float(heal))
            gmsh.model.occ.synchronize()
        except Exception as exc:  # pragma: no cover
            print(f"[heal] skipped: {exc}")

    # ---- sizing
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", float(extend_from_boundary))
    if curvature is None:
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    else:
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", float(curvature))
    if mesh_size:
        gmsh.option.setNumber("Mesh.MeshSizeMax", float(mesh_size))
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(mesh_size) * size_min_factor)

    # ---- algorithm & optimisation
    gmsh.option.setNumber("Mesh.Algorithm3D", ALGO3D[algo3d])
    gmsh.option.setNumber("Mesh.Optimize", int(optimize))
    gmsh.option.setNumber("Mesh.OptimizeNetgen", int(netgen))
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.Binary", 0)

    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh))
    gmsh.finalize()
    print(f"[mesh] wrote {msh}")
    return msh


def estimate_thickness(geometry: Path, diag: float, out: Path) -> float:
    """Coarse-mesh the geometry and estimate the wall thickness via ``t ~ 2V/A``.

    For a plate of thickness ``t``, ``V = A_lat * t`` and ``A = 2 * A_lat``, so ``t = 2V/A``
    is a good global estimate for thin-walled parts.  The element size must be well below
    it (``h <= t/2``), otherwise no tetrahedron fits inside the wall and the mesh fills up
    with flat slivers.
    """
    import gmsh

    import meshio

    from tetparam.mesh import TetMesh

    probe = out / "_probe"
    probe.mkdir(parents=True, exist_ok=True)
    msh = mesh_geometry(geometry, diag / 25.0, probe, netgen=False, size_min_factor=0.5)
    md = meshio.read(str(msh))
    tets = md.cells_dict.get("tetra")
    if tets is None or len(tets) == 0:
        raise SystemExit("coarse probe mesh produced no tetrahedra")
    m = TetMesh(np.asarray(md.points, float), np.asarray(tets, np.int64)).fix_orientation()
    bnd = m.boundary()
    tri = m.vertices[bnd["faces"]]
    A = float(0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum())
    t = 2.0 * m.total_volume() / A
    del gmsh
    print(f"[auto] coarse probe: t ~ 2V/A = {t:.4g}  ->  --mesh-size {t / 2:.4g}")
    return t


# --------------------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("geometry", type=Path)
    ap.add_argument("--mesh-size", default=None,
                    help="target element size, should be <= 1/2 of the wall thickness; "
                         "use 'auto' to measure the thickness from a coarse probe mesh first")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--algo3d", choices=sorted(ALGO3D), default="delaunay")
    ap.add_argument("--optimize", type=int, default=1, choices=[0, 1])
    ap.add_argument("--netgen", action="store_true", default=True)
    ap.add_argument("--no-netgen", dest="netgen", action="store_false")
    ap.add_argument("--curvature", type=float, default=None,
                    help="Mesh.MeshSizeFromCurvature (elements per 2*pi); "
                         "omit/0 to disable automatic curvature refinement")
    ap.add_argument("--extend-from-boundary", type=float, default=0.0)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--heal", type=float, default=0.0,
                    help="heal the CAD first with this tolerance (model units); merges "
                         "small edges/faces -- use ~mesh-size/10 for feature-rich parts")
    ap.add_argument("--auto-factor", type=float, default=1.0,
                    help="with --mesh-size auto: h = auto_factor * t / 2 (default 1.0)")
    ap.add_argument("--no-report", action="store_true")
    args = ap.parse_args(argv)

    if not args.geometry.exists():
        raise SystemExit(f"not found: {args.geometry}")
    out = args.out or Path("out") / f"mesh_{args.geometry.stem}"

    mesh_size = args.mesh_size
    if isinstance(mesh_size, str) and mesh_size.lower() == "auto":
        import gmsh

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(args.geometry))
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(-1, -1)
        gmsh.finalize()
        diag = float(np.linalg.norm([x1 - x0, y1 - y0, z1 - z0]))
        t = estimate_thickness(args.geometry, diag, out)
        mesh_size = t / 2.0 * float(args.auto_factor)
    else:
        mesh_size = float(mesh_size) if mesh_size is not None else None

    msh = mesh_geometry(
        args.geometry, mesh_size, out,
        algo3d=args.algo3d, optimize=args.optimize, netgen=args.netgen,
        curvature=args.curvature, extend_from_boundary=args.extend_from_boundary,
        threads=args.threads, heal=args.heal,
    )

    if not args.no_report:
        import meshio

        md = meshio.read(str(msh))
        tets = md.cells_dict.get("tetra")
        if tets is None:
            raise SystemExit("the mesher produced no tetrahedra")
        m = TetMesh(np.asarray(md.points, float), np.asarray(tets, np.int64))
        m.fix_orientation()
        m.write_vtk(out / "mesh.vtk")
        quality_report(m, f"{args.geometry.name}  algo={args.algo3d} netgen={args.netgen}")
        print(f"[out] {out / 'mesh.vtk'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
