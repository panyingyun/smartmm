"""Compare initialisations for the ball map: star (bijective) vs harmonic vs relaxed."""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from tetparam.generate import DEMOS  # noqa: E402
from tetparam.harmonic import (  # noqa: E402
    boundary_map_sphere,
    fem_stiffness,
    fold_free_relax,
    solve_harmonic,
    star_ball_map,
)
from tetparam.metrics import distortion_report  # noqa: E402


def main() -> int:
    for name in ["ball", "pear", "twisted", "bump", "bent", "taper"]:
        m = DEMOS[name]()
        L = fem_stiffness(m)
        print(f"\n=== {name} ({m.n_vertices}v {m.n_tets}t) ===")

        t0 = time.perf_counter()
        ys = star_ball_map(m)
        ts = time.perf_counter() - t0
        r = distortion_report(m, m.with_vertices(ys))
        print(
            f"  star       flips {r.n_flipped:4d}  dist mean {r.distortion_mean:.4f} "
            f"max {r.distortion_max:9.3f}  excess {r.stretch_excess:8.3%}   ({ts * 1e3:.0f} ms)"
        )

        bv = boundary_map_sphere(m, method="laplacian")
        yh, _ = solve_harmonic(m, bv, L=L)
        r = distortion_report(m, m.with_vertices(yh))
        print(
            f"  harmonic   flips {r.n_flipped:4d}  dist mean {r.distortion_mean:.4f} "
            f"max {r.distortion_max:9.3f}  excess {r.stretch_excess:8.3%}"
        )

        t0 = time.perf_counter()
        yr, info = fold_free_relax(m, ys, L=L, max_iter=300)
        tr = time.perf_counter() - t0
        r = distortion_report(m, m.with_vertices(yr))
        print(
            f"  relaxed    flips {r.n_flipped:4d}  dist mean {r.distortion_mean:.4f} "
            f"max {r.distortion_max:9.3f}  excess {r.stretch_excess:8.3%}   "
            f"({tr:.2f} s, {info['iterations']} it, min-ratio {info['min_ratio']:.3g})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
