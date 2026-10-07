"""End-to-end demo + figures for the conformal + optimal-transport pipeline.

    python tests/run_demo.py            # all shapes, writes out/report.json + figures

Stages
------
1. harmonic ("conformal") map of the solid onto the unit ball  -> harmonic.ball_map
2. volume correction                                          -> volume_ot.vsem
   (optionally the explicit OT route                            -> volume_ot.ot_vertex_map)
3. tetrahedral remeshing by pulling a template ball mesh back  -> remesh.pullback_remesh
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from tetparam.generate import ball_cube_mesh, demo_mesh  # noqa: E402
from tetparam.harmonic import ball_map, fem_stiffness  # noqa: E402
from tetparam.metrics import distortion_report  # noqa: E402
from tetparam.remesh import pullback_remesh, refine, uniform_ball_template  # noqa: E402
from tetparam.volume_ot import vsem  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out"


def figure(name: str, m, y_ball, y_vol) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # --- volume ratio histograms
    for ax, (tag, y, col) in zip(
        axes[0, :2],
        [("stage 1: harmonic / conformal", y_ball, "tab:orange"),
         ("stage 2: volume-preserving", y_vol, "tab:blue")],
    ):
        r = distortion_report(m, m.with_vertices(y))
        ax.hist(np.clip(_ratios(m, m.with_vertices(y)), 0.0, 4.0), bins=60,
                color=col, alpha=0.85)
        ax.axvline(1.0, color="k", ls="--", lw=1)
        ax.set_title(f"{tag}\nE_V excess {r.stretch_excess:.2%}, flips {r.n_flipped}")
        ax.set_xlabel(r"volume ratio $|f(\tau)|/|\tau|$")
        ax.set_ylabel("tetrahedra")

    # --- tetrahedron quality
    ax = axes[0, 2]
    q0 = m.quality()
    q1 = m.with_vertices(y_vol).quality()
    ax.hist(q0["min_dihedral"], bins=40, alpha=0.6, label="input", color="gray")
    ax.hist(q1["min_dihedral"], bins=40, alpha=0.6, label="volume-preserving map", color="tab:green")
    ax.set_xlabel("min dihedral angle (deg)")
    ax.set_ylabel("tetrahedra")
    ax.legend()
    ax.set_title("element shape (unchanged connectivity)")

    # --- cross sections through the ball
    for ax, (tag, y) in zip(
        axes[1, :2],
        [("ball coordinates, stage 1", y_ball), ("ball coordinates, stage 2", y_vol)],
    ):
        ax.scatter(y[:, 0], y[:, 2], c=np.linalg.norm(y, axis=1), s=4, cmap="viridis")
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(np.cos(th), np.sin(th), "k-", lw=0.6)
        ax.set_aspect("equal")
        ax.set_title(tag)

    # --- radius in the ball, stage 1 vs stage 2 (both per vertex)
    ax = axes[1, 2]
    r1 = np.linalg.norm(y_ball, axis=1)
    r2 = np.linalg.norm(y_vol, axis=1)
    ax.scatter(r1, r2, s=6, alpha=0.7)
    lim = float(max(r1.max(), r2.max()))
    ax.plot([0, lim], [0, lim], "k--", lw=1)
    ax.set_xlabel("radius in the ball after stage 1")
    ax.set_ylabel("radius in the ball after stage 2")
    ax.set_title("how far the volume correction moves points")

    fig.suptitle(f"{name}: conformal (harmonic) map + optimal-transport volume correction")
    fig.tight_layout()
    fig.savefig(OUT / f"pipeline_{name}.png", dpi=120)
    plt.close(fig)


def _ratios(m, target):
    v0 = np.abs(m.signed_volumes())
    v1 = target.signed_volumes()
    return v1 / np.where(v0 > 0, v0, 1.0)


def main(shapes: list[str] | None = None, res: int = 6, ot: bool = False) -> int:
    shapes = shapes or ["pear", "taper", "bent"]
    OUT.mkdir(parents=True, exist_ok=True)
    records = []

    for name in shapes:
        m = demo_mesh(name)
        rec: dict = {
            "shape": name,
            "n_vertices": m.n_vertices,
            "n_tets": m.n_tets,
            "volume": m.total_volume(),
        }
        print(f"\n=== {name}: {m.n_vertices} v, {m.n_tets} t, |M| = {m.total_volume():.5f} ===")

        t0 = time.perf_counter()
        b = ball_map(m, method="auto")
        y_ball = b["vertices"]
        r0 = distortion_report(m, m.with_vertices(y_ball))
        rec["stage1_harmonic"] = dict(
            method=b["method"], seconds=time.perf_counter() - t0, **r0.as_dict()
        )
        print(f"  [1] harmonic ball map ({b['method']}, {rec['stage1_harmonic']['seconds']:.2f} s)")
        print("      " + r0.summary().replace("\n", "\n      "))

        t0 = time.perf_counter()
        v = vsem(m, y_ball, max_iter=1000)
        y_vol = v["vertices"]
        r1 = distortion_report(m, m.with_vertices(y_vol))
        rec["stage2_vsem"] = dict(
            iterations=v["iterations"], seconds=time.perf_counter() - t0, **r1.as_dict()
        )
        print(f"  [2] volumetric stretch energy minimisation ({v['iterations']} it, "
              f"{rec['stage2_vsem']['seconds']:.2f} s)")
        print("      " + r1.summary().replace("\n", "\n      "))

        if ot:
            from tetparam.volume_ot import ot_vertex_map

            t0 = time.perf_counter()
            try:
                o = ot_vertex_map(m, y_ball, neighbours=48, max_iter=25, tol=1e-6)
                r2 = distortion_report(m, m.with_vertices(o["vertices"]))
                rec["stage2_ot"] = dict(
                    n_sites=o["n_sites"], mass_error_rel=o["mass_error_rel"],
                    n_empty=o["diagram"].n_empty, iterations=o["iterations"],
                    seconds=time.perf_counter() - t0, **r2.as_dict(),
                )
                print(f"  [2'] OT power diagram: {o['n_sites']} sites, {o['iterations']} it, "
                      f"rel. mass error {o['mass_error_rel']:.2e}, "
                      f"empty {o['diagram'].n_empty}, E_V excess {r2.stretch_excess:.2%}")
            except Exception as exc:  # pragma: no cover
                rec["stage2_ot"] = {"error": repr(exc)}

        # ---- remeshing by pulling template ball meshes back through f^{-1}
        image = m.with_vertices(y_vol)
        templates = {
            "uniform": uniform_ball_template(image, n_interior=500),
            "refine1": refine(image, 1),
        }
        if res != 6:
            templates["coarse"] = ball_cube_mesh(max(2, res - 2))
        for tag, tmpl in templates.items():
            t0 = time.perf_counter()
            new = pullback_remesh(m, y_vol, template=tmpl)
            q = new.quality()
            rec[f"remesh_{tag}"] = dict(
                n_vertices=new.n_vertices, n_tets=new.n_tets,
                total_volume=new.total_volume(),
                volume_error=new.total_volume() / m.total_volume() - 1.0,
                n_flipped=new.check_orientation(),
                radius_ratio_mean=q["radius_ratio_mean"],
                min_dihedral_mean=q["min_dihedral_mean"],
                min_dihedral_min=float(np.min(q["min_dihedral"])),
                seconds=time.perf_counter() - t0,
            )
            rr = rec[f"remesh_{tag}"]
            print(f"  [3] pullback remesh ({tag:8s}, {tmpl.n_vertices} v / {tmpl.n_tets} t): "
                  f"volume error {rr['volume_error']:+.2e}, flips {rr['n_flipped']}, "
                  f"radius-ratio {rr['radius_ratio_mean']:.2f}, "
                  f"min dihedral {rr['min_dihedral_min']:.2f} deg  ({rr['seconds']:.2f} s)")
            new.write_vtk(OUT / f"remesh_{name}_{tag}.vtk")

        m.write_vtk(OUT / f"input_{name}.vtk")
        m.with_vertices(y_ball).write_vtk(OUT / f"ball_{name}.vtk")
        m.with_vertices(y_vol).write_vtk(OUT / f"volpres_{name}.vtk")
        try:
            figure(name, m, y_ball, y_vol)
            print(f"  [fig] {OUT / f'pipeline_{name}.png'}")
        except Exception as exc:  # pragma: no cover
            print(f"  [fig] failed: {exc}")
        records.append(rec)

    with (OUT / "report.json").open("w", encoding="utf-8") as f:
        json.dump(records, f, indent=2, default=float)
    print(f"\nwrote {OUT / 'report.json'}")
    return 0


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    raise SystemExit(main(shapes=args or None, ot="--ot" in sys.argv))
