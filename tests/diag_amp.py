"""Find deformation amplitudes for which the star-shaped map is fold-free."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from tetparam.generate import ball_tet_mesh, bump, pear, shear_bend, taper, twist  # noqa: E402
from tetparam.harmonic import star_ball_map, _tet_volumes  # noqa: E402


def flips(m, y):
    v = _tet_volumes(y, m.tets)
    return int((v <= 0).sum()), float((v / np.abs(m.signed_volumes())).min())


def main() -> int:
    base = ball_tet_mesh(162, 500)
    trials = {
        "pear": [(pear, a) for a in (0.15, 0.25, 0.35, 0.5)],
        "taper": [(taper, a) for a in (0.75, 0.6, 0.45)],
        "twist": [(twist, a) for a in (0.3, 0.5, 0.7, 1.0)],
        "bump": [(bump, a) for a in (0.05, 0.10, 0.18)],
        "bend": [(shear_bend, a) for a in (0.1, 0.2, 0.3)],
    }
    for name, ts in trials.items():
        for f, a in ts:
            m = f(base, a)
            y = star_ball_map(m)
            nf, mn = flips(m, y)
            print(f"  {name:6s} a={a:<5} flips {nf:5d}/{m.n_tets}  min-ratio {mn:9.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
