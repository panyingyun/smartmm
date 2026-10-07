"""Smoke test / validation for stage 1 (harmonic ball map) and the mesh core."""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from tetparam.generate import demo_mesh  # noqa: E402
from tetparam.harmonic import harmonic_ball_map  # noqa: E402
from tetparam.metrics import distortion_report  # noqa: E402


def main() -> int:
    for name in ("ball", "pear", "twist", "bump", "bent", "taper", "cube"):
        m = demo_mesh(name)
        print(f"\n=== {name} ===")
        print(m.summary())
        t0 = time.perf_counter()
        res = harmonic_ball_map(m, boundary_method="laplacian")
        dt = time.perf_counter() - t0
        print(f"{res.summary()}   ({dt * 1e3:.0f} ms)")
        rep = distortion_report(m, res.target)
        print(rep.summary())
        b = res.target.is_boundary_vertex
        r = np.linalg.norm(res.target.vertices[b], axis=1)
        print(f"boundary radius: min {r.min():.6f} max {r.max():.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
