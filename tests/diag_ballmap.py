"""Compare ball-map strategies (star / harmonic / relax / auto) on the demo shapes."""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tetparam.generate import DEMOS  # noqa: E402
from tetparam.harmonic import ball_map  # noqa: E402
from tetparam.metrics import distortion_report  # noqa: E402


def main() -> int:
    for name in ["ball", "pear", "twisted", "bump", "bent", "taper"]:
        m = DEMOS[name]()
        print(f"\n=== {name} ({m.n_vertices}v {m.n_tets}t) ===")
        for method in ("star", "harmonic", "auto"):
            t0 = time.perf_counter()
            r = ball_map(m, method=method, relax=False, verbose=False)
            dt = time.perf_counter() - t0
            rep = distortion_report(m, m.with_vertices(r["vertices"]))
            print(
                f"  {method:9s} -> {r['method']:16s} flips {r['n_flipped']:4d}  "
                f"dist mean {rep.distortion_mean:8.4f}  p99 {rep.distortion_p99:9.4f}  "
                f"max {rep.distortion_max:11.3f}  E_V excess {rep.stretch_excess:9.3%}  "
                f"({dt * 1e3:.0f} ms)"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
