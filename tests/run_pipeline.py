"""Run the whole conformal + optimal-transport pipeline on a mesh you supply.

Examples
--------
    # a TetGen tetrahedral mesh
    python tests/run_pipeline.py part.node --ele part.ele --out out/part

    # a volume mesh (VTK / VTU / Medit / Gmsh .msh) with tetrahedra
    python tests/run_pipeline.py part.vtk --out out/part

    # a *surface* geometry: tetrahedralise it first with Gmsh, then run the pipeline
    python tests/run_pipeline.py part.stl  --mesh-size 0.05 --out out/part
    python tests/run_pipeline.py part.step --mesh-size 0.02 --out out/part

Outputs (in --out, default ``out/pipeline``)
--------------------------------------------
    <out>/input.vtk           the input tetrahedral mesh
    <out>/ball.vtk            stage 1: harmonic ("conformal") map to the unit ball
    <out>/volpres.vtk         stage 2: volume-preserving map to the unit ball
    <out>/remesh_refine1.vtk  stage 3: remeshed solid (template = refined image mesh)
    <out>/remesh_uniform.vtk  stage 3: remeshed solid (quasi-uniform ball template)
    <out>/pipeline.png        distortion figures
    <out>/report.json         every measured quantity
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from tetparam.harmonic import ball_map  # noqa: E402
from tetparam.mesh import TetMesh  # noqa: E402
from tetparam.metrics import distortion_report  # noqa: E402
from tetparam.remesh import pullback_remesh, refine, uniform_ball_template  # noqa: E402
from tetparam.volume_ot import vsem  # noqa: E402

GEOMETRY_SUFFIXES = {".step", ".stp", ".iges", ".igs", ".brep", ".geo"}
SURFACE_SUFFIXES = {".stl", ".off", ".ply", ".obj"}


# --------------------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------------------
def tetrahedralize_geometry(path: Path, mesh_size: float | None, out: Path) -> Path:
    """Turn a CAD model or a closed surface mesh into a tetrahedral mesh with Gmsh.

    * CAD (STEP/IGES/BREP): imported as a solid, then ``generate(3)`` fills it.
    * Closed surface mesh (STL/OFF/PLY/OBJ): Gmsh has to *reconstruct* a volume first --
      classify the facets, rebuild the geometry, close it with a surface loop and add a
      volume -- otherwise ``generate(3)`` would silently produce no tetrahedra.
    """
    try:
        import gmsh
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            f"Gmsh is required to mesh {path.name}. Install it with `pip install gmsh`, "
            f"or pre-tetrahedralise the geometry yourself, e.g. "
            f"`tetgen -pq1.4a0.01 {path.stem}.stl`."
        ) from exc

    msh = out / (path.stem + ".msh")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Geometry.OCCImportLabels", 1)
    gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
    if mesh_size:
        gmsh.option.setNumber("Mesh.MeshSizeMax", float(mesh_size))
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(mesh_size) * 0.4)

    if path.suffix.lower() in SURFACE_SUFFIXES:
        gmsh.merge(str(path))
        gmsh.model.mesh.classifySurfaces(
            np.pi * 40.0 / 180.0, True, True, np.pi
        )
        gmsh.model.mesh.createGeometry()
        surfaces = gmsh.model.getEntities(2)
        if not surfaces:
            gmsh.finalize()
            raise SystemExit(f"{path.name}: no surfaces recovered from the facet soup")
        loop = gmsh.model.geo.addSurfaceLoop([s[1] for s in surfaces])
        gmsh.model.geo.addVolume([loop])
        gmsh.model.geo.synchronize()
    else:
        gmsh.open(str(path))

    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh))
    gmsh.finalize()
    print(f"[0] tetrahedralised {path.name} -> {msh.name}")
    return msh


def load_mesh(path: Path, ele: Path | None, mesh_size: float | None, out: Path) -> TetMesh:
    suffix = path.suffix.lower()
    if suffix in GEOMETRY_SUFFIXES or suffix in SURFACE_SUFFIXES:
        path = tetrahedralize_geometry(path, mesh_size, out)
        suffix = path.suffix.lower()

    if suffix == ".node":
        ele = ele or path.with_suffix(".ele")
        if not Path(ele).exists():
            raise SystemExit(f"no .ele file found next to {path}")
        m = TetMesh.from_tetgen(path, ele)
    elif suffix == ".ele":
        node = path.with_suffix(".node")
        if not node.exists():
            raise SystemExit(f"no .node file found next to {path}")
        m = TetMesh.from_tetgen(node, path)
    else:
        import meshio

        md = meshio.read(str(path))
        tets = md.cells_dict.get("tetra")
        if tets is None or len(tets) == 0:
            raise SystemExit(
                f"{path.name} contains no tetrahedra. If it is a surface mesh, pass it as "
                f".stl/.step so Gmsh can tetrahedralise it first."
            )
        m = TetMesh(np.asarray(md.points, dtype=float), np.asarray(tets, dtype=np.int64))

    n_bad = m.check_orientation()
    if n_bad:
        print(f"[0] fixing {n_bad} inverted tetrahedra")
        m.fix_orientation()
    print(f"[0] {m.summary()}")
    return m


# --------------------------------------------------------------------------------------
# pipeline
# --------------------------------------------------------------------------------------
def run(m: TetMesh, out: Path, template_res: int = 6, remesh: bool = True) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    rec: dict = {
        "input": {"n_vertices": m.n_vertices, "n_tets": m.n_tets,
                  "volume": m.total_volume()},
    }
    m.write_vtk(out / "input.vtk")

    # ---- stage 1: harmonic / conformal map to the unit ball
    t0 = time.perf_counter()
    b = ball_map(m, method="auto")
    y_ball = b["vertices"]
    r0 = distortion_report(m, m.with_vertices(y_ball))
    rec["stage1_harmonic"] = dict(method=b["method"], seconds=time.perf_counter() - t0,
                                  **r0.as_dict())
    print(f"[1] conformal ball map ({b['method']}, {rec['stage1_harmonic']['seconds']:.2f} s)")
    print("    " + r0.summary().replace("\n", "\n    "))

    # ---- stage 2: volume correction (volumetric stretch energy minimisation)
    t0 = time.perf_counter()
    v = vsem(m, y_ball, max_iter=2000)
    y_vol = v["vertices"]
    r1 = distortion_report(m, m.with_vertices(y_vol))
    rec["stage2_volume"] = dict(iterations=v["iterations"], seconds=time.perf_counter() - t0,
                                **r1.as_dict())
    print(f"[2] volume-preserving map ({v['iterations']} it, "
          f"{rec['stage2_volume']['seconds']:.2f} s)")
    print("    " + r1.summary().replace("\n", "\n    "))
    m.with_vertices(y_ball).write_vtk(out / "ball.vtk")
    m.with_vertices(y_vol).write_vtk(out / "volpres.vtk")

    # ---- stage 3: tetrahedral remeshing by pulling templates back through f^-1
    if remesh:
        image = m.with_vertices(y_vol)
        templates = {
            "refine1": refine(image, 1),
            "uniform": uniform_ball_template(image, n_interior=max(200, m.n_vertices)),
        }
        for tag, tmpl in templates.items():
            t0 = time.perf_counter()
            new = pullback_remesh(m, y_vol, template=tmpl)
            q = new.quality()
            rec[f"remesh_{tag}"] = dict(
                n_vertices=new.n_vertices, n_tets=new.n_tets,
                volume_error=new.total_volume() / m.total_volume() - 1.0,
                n_flipped=new.check_orientation(),
                radius_ratio_mean=q["radius_ratio_mean"],
                min_dihedral_min=float(np.min(q["min_dihedral"])),
                seconds=time.perf_counter() - t0,
            )
            rr = rec[f"remesh_{tag}"]
            print(f"[3] remesh ({tag:8s} {new.n_vertices} v / {new.n_tets} t): "
                  f"volume error {rr['volume_error']:+.2e}, flips {rr['n_flipped']}, "
                  f"radius-ratio {rr['radius_ratio_mean']:.2f}, "
                  f"min dihedral {rr['min_dihedral_min']:.2f} deg")
            new.write_vtk(out / f"remesh_{tag}.vtk")

    # ---- figures
    try:
        make_figure(m, y_ball, y_vol, out / "pipeline.png")
        print(f"[fig] {out / 'pipeline.png'}")
    except Exception as exc:  # pragma: no cover
        print(f"[fig] skipped: {exc}")

    with (out / "report.json").open("w", encoding="utf-8") as f:
        json.dump(rec, f, indent=2, default=float)
    print(f"[out] {out / 'report.json'}")
    del template_res
    return rec


def make_figure(m: TetMesh, y_ball: np.ndarray, y_vol: np.ndarray, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def ratios(y):
        v0 = np.abs(m.signed_volumes())
        return m.with_vertices(y).signed_volumes() / np.where(v0 > 0, v0, 1.0)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for a, (tag, y, c) in zip(
        ax[:2],
        [("stage 1: conformal / harmonic", y_ball, "tab:orange"),
         ("stage 2: volume-preserving", y_vol, "tab:blue")],
    ):
        a.hist(np.clip(ratios(y), 0, 4), bins=60, color=c, alpha=0.85)
        a.axvline(1.0, color="k", ls="--", lw=1)
        r = distortion_report(m, m.with_vertices(y))
        a.set_title(f"{tag}\nE_V excess {r.stretch_excess:.2%}, flips {r.n_flipped}")
        a.set_xlabel(r"volume ratio $|f(\tau)|/|\tau|$")
        a.set_ylabel("tetrahedra")
    q = m.with_vertices(y_vol).quality()
    ax[2].hist(q["min_dihedral"], bins=40, color="tab:green", alpha=0.85)
    ax[2].set_xlabel("min dihedral angle (deg)")
    ax[2].set_ylabel("tetrahedra")
    ax[2].set_title("element shape after stage 2")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="conformal + optimal-transport volume parameterisation of a tetrahedral mesh",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("mesh", type=Path,
                    help=".node (+--ele), or a volume mesh (.vtk/.vtu/.msh/.mesh), "
                         "or a surface/CAD geometry (.stl/.step/.iges/.brep) to tetrahedralise")
    ap.add_argument("--ele", type=Path, default=None, help="TetGen .ele file")
    ap.add_argument("--out", type=Path, default=None, help="output directory")
    ap.add_argument("--mesh-size", type=float, default=None,
                    help="target element size when tetrahedralising a geometry")
    ap.add_argument("--no-remesh", action="store_true", help="skip stage 3")
    args = ap.parse_args(argv)

    if not args.mesh.exists():
        raise SystemExit(f"not found: {args.mesh}")
    out = args.out or Path("out") / f"pipeline_{args.mesh.stem}"
    out.mkdir(parents=True, exist_ok=True)

    m = load_mesh(args.mesh, args.ele, args.mesh_size, out)
    run(m, out, remesh=not args.no_remesh)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
