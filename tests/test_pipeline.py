"""End-to-end stage 1 + 2 demo: ball map then volume correction."""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tetparam.generate import DEMOS  # noqa: E402
from tetparam.harmonic import ball_map  # noqa: E402
from tetparam.metrics import distortion_report  # noqa: E402
from tetparam.volume_ot import solve_ball_ot, vsem  # noqa: E402


def main() -> int:
    for name in ["ball", "pear", "bent", "taper", "twisted", "bump"]:
        m = DEMOS[name]()
        print(f"\n=== {name} ({m.n_vertices}v {m.n_tets}t, |M|={m.total_volume():.5f}) ===")
        t0 = time.perf_counter()
        b = ball_map(m, method="auto")
        r0 = distortion_report(m, m.with_vertices(b["vertices"]))
        print(
            f"  stage 1 [{b['method']:16s}] flips {b['n_flipped']:3d}  "
            f"ratio [{r0.ratio_min:8.4f}, {r0.ratio_max:8.4f}]  "
            f"dist mean {r0.distortion_mean:8.4f} p99 {r0.distortion_p99:8.4f}  "
            f"E_V excess {r0.stretch_excess:9.3%}"
        )

        t1 = time.perf_counter()
        v = vsem(m, b["vertices"], max_iter=800)
        r1 = distortion_report(m, m.with_vertices(v["vertices"]))
        print(
            f"  stage 2 [vsem {v['iterations']:4d} it    ] flips {v['n_flipped']:3d}  "
            f"ratio [{r1.ratio_min:8.4f}, {r1.ratio_max:8.4f}]  "
            f"dist mean {r1.distortion_mean:8.4f} p99 {r1.distortion_p99:8.4f}  "
            f"E_V excess {r1.stretch_excess:9.3%}   ({time.perf_counter() - t1:.2f} s)"
        )

        t2 = time.perf_counter()
        ot = solve_ball_ot(m, v["vertices"], neighbours=64)
        dt = time.perf_counter() - t2
        print(
            f"  stage 2 [OT power diagram]  sites {ot['n_sites']:5d}  "
            f"mass err {ot['mass_error_rel']:.3e}  |P| {ot['domain_volume']:.6f}  "
            f"domain mismatch {ot['domain_mismatch']:+.3e}  "
            f"empty {ot['ot'].diagram.n_empty}  ({dt:.2f} s)"
        )
        print(f"  total {time.perf_counter() - t0:.2f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
