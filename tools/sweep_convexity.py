#!/usr/bin/env python
"""How non-convex can the solid be before the pipeline breaks?

Sweeps the "peanut" amplitude  r(theta) = 1 + a*cos(2 theta)  and reports, for each a:
  - whether the RADIAL-PROJECTION boundary map is a bijection to the sphere
    (test: sum of signed spherical triangle areas == 4*pi; a multi-cover gives 8*pi);
  - the volume-distortion statistics of the plain P1 HARMONIC map to the ball;
  - the same after VSEM (fixed boundary);
  - the VSEM energy deficit E_V(f) - 1.5*V(f).

This is the practical answer to "when does the classical harmonic-to-ball
construction fold, and does VSEM save it?".

Run:  C:\\Python312\\python.exe tools\\sweep_convexity.py
"""
from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, r"E:\panyingyun\smartmm\tools")
import tet_ball_ref as R  # noqa: E402
import demo_failure_modes as D  # noqa: E402


def blob(sub: int, maxvol: float, kind: str, amp: float):
    import gpytoolbox as gpy
    import tetgen as tg

    V, F = gpy.icosphere(sub)
    V = np.asarray(V, dtype=np.float64)
    x, y, z = V[:, 0], V[:, 1], V[:, 2]
    theta = np.arccos(np.clip(z, -1, 1))
    phi = np.arctan2(y, x)
    if kind == "peanut":
        r = 1.0 + amp * np.cos(2 * theta)
    elif kind == "wedge":
        r = 1.0 + amp * np.abs(np.sin(3 * phi) * np.sin(theta))
    elif kind == "crinkly":
        r = 1.0 + amp * np.sin(5 * theta) * np.sin(4 * phi)
    else:
        raise ValueError(kind)
    r = np.maximum(r, 0.15)
    V = V * r[:, None]
    tgen = tg.TetGen(V, np.asarray(F, dtype=np.int32))
    tgen.tetrahedralize(switches=f"pq1.414a{maxvol}Y")
    return np.asarray(tgen.node, dtype=np.float64), np.asarray(tgen.elem, dtype=np.int64)


def main() -> int:
    sub, maxvol = 2, 0.06
    print(f"mesh: icosphere sub={sub}, tetgen maxvol={maxvol}\n")
    hdr = (f"{'shape':8s} {'amp':>5s} {'sphArea/4pi':>11s} {'bndOK':>6s} "
           f"{'A:ratioStd':>10s} {'A:J<0.5':>8s} {'A:flip':>7s} "
           f"{'B:ratioStd':>10s} {'B:min':>8s} {'B:flip':>7s} {'B:gap':>10s}")
    print(hdr)
    print("-" * len(hdr))
    for kind, amps in (("peanut", [0.0, 0.10, 0.20, 0.30, 0.40, 0.45]),
                       ("crinkly", [0.0, 0.10, 0.20, 0.25]),
                       ("wedge", [0.0, 0.10, 0.20, 0.35])):
        for a in amps:
            V, T = blob(sub, maxvol, kind, a)
            n = len(V)
            b, bf = R.boundary_vertices(T, n, V)
            mu = R.tet_volumes(V, T)
            S = V[b] - V[b].mean(axis=0)
            S = S / np.linalg.norm(S, axis=1, keepdims=True)
            ar = D.spherical_tri_areas(S, bf)
            ratio = ar.sum() / (4 * math.pi)
            bnd_ok = abs(ratio - 1.0) < 0.02

            U, _, _ = R.harmonic_ball_map(V, T, S)
            ru = R.volume_distortion_ratios(V, T, U)

            f, hist = R.vsem_solve(V, T, S, n_iter=250, tol=1e-13)
            rf = R.volume_distortion_ratios(V, T, f)
            gap = R.vsem_energy(V, T, f) - 1.5 * float(R.tet_volumes(f, T).sum())

            print(f"{kind:8s} {a:5.2f} {ratio:11.4f} {str(bnd_ok):>6s} "
                  f"{ru.std():10.4f} {float((ru < 0.5).mean()):8.3%} {int((ru <= 0).sum()):7d} "
                  f"{rf.std():10.4f} {rf.min():8.5f} {int((rf <= 0).sum()):7d} {gap:10.3e}")
    print("\nLegend: A = plain P1 harmonic map to the ball; "
          "B = VSEM fixed point with the same (radial-projection) boundary map.")
    print("sphArea/4pi = 1 means the boundary map is a bijection to the sphere; "
          "2 means it double-covers the sphere.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
