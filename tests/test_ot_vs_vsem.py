"""Compare the two volume-correction routes on a small mesh (fast).

Route OT  : the paper's composition  F = phi^{-1} o psi,  F(v_i) = centroid(W_i).
Route VSEM: minimise the volumetric stretch energy with the boundary fixed on the sphere.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from tetparam.generate import ball_cube_mesh, pear, shear_bend, taper, twist  # noqa: E402
from tetparam.harmonic import ball_map  # noqa: E402
from tetparam.metrics import distortion_report  # noqa: E402
from tetparam.volume_ot import ot_vertex_map, vsem  # noqa: E402


def line(tag, m, y):
    r = distortion_report(m, m.with_vertices(y))
    print(
        f"    {tag:26s} flips {r.n_flipped:4d}  ratio [{r.ratio_min:9.5f}, {r.ratio_max:9.5f}]  "
        f"mean {r.ratio_mean:8.5f} std {r.ratio_std:9.2e}  dist mean {r.distortion_mean:8.5f} "
        f"p99 {r.distortion_p99:8.5f}  E_V excess {r.stretch_excess:9.3%}"
    )
    return r


def main() -> int:
    base = ball_cube_mesh(5)  # 750 tets
    shapes = {
        "ball": base,
        "pear": pear(base, 0.45),
        "taper": taper(base, 0.55),
        "twist": twist(base, 0.7),
        "bent": shear_bend(base, 0.30),
    }
    for name, m in shapes.items():
        print(f"\n=== {name}: {m.n_vertices}v {m.n_tets}t ===")
        b = ball_map(m, method="star")
        line(f"stage1 [{b['method']}]", m, b["vertices"])

        t0 = time.perf_counter()
        ot = ot_vertex_map(m, b["vertices"], neighbours=48, max_iter=25, tol=1e-6,
                           snap_boundary=False)
        line(f"OT raw cell centroid ({ot['iterations']}it)", m, ot["vertices"])
        ot2 = ot_vertex_map(m, b["vertices"], neighbours=48, max_iter=25, tol=1e-6,
                            snap_boundary=True)
        line("OT + sphere snap", m, ot2["vertices"])
        print(
            f"      OT diagnostics: sites {ot['n_sites']}, "
            f"rel mass err {ot['mass_error_rel']:.3e}, empty {ot['diagram'].n_empty}, "
            f"domain mismatch {ot['domain_mismatch']:+.2e}, "
            f"({time.perf_counter() - t0:.2f} s for both)"
        )

        t0 = time.perf_counter()
        v = vsem(m, ot2["vertices"], max_iter=800, barrier=1.0)
        line(f"VSEM after OT ({v['iterations']}it)", m, v["vertices"])
        print(f"      VSEM {time.perf_counter() - t0:.2f} s")
        del v
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
