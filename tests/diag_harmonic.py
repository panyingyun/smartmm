"""Diagnostics for the harmonic ball map: where do folds come from?"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from tetparam.generate import DEMOS, demo_mesh  # noqa: E402
from tetparam.harmonic import boundary_map_sphere, fem_stiffness, solve_harmonic  # noqa: E402
from tetparam.mesh import TetMesh  # noqa: E402


def boundary_folds(mesh: TetMesh, y: np.ndarray) -> tuple[int, int]:
    """Count boundary triangles whose image is flipped (normal points inward)."""
    b = mesh.boundary()
    faces = b["faces"]
    tri = y[faces]
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    c = tri.mean(axis=1)
    # outward = away from origin (the sphere is centred at the origin)
    out = np.einsum("ij,ij->i", nrm, c)
    return int((out <= 0).sum()), len(faces)


def main() -> int:
    for name in ["ball", "pear", "twisted", "bump", "bent", "taper"]:
        m = DEMOS[name]()
        print(f"\n=== {name}: {m.n_vertices}v {m.n_tets}t ===")
        for method in ("radial", "laplacian"):
            bv = boundary_map_sphere(m, method=method)
            nf, ntot = boundary_folds(m, bv)
            L = fem_stiffness(m)
            y, res = solve_harmonic(m, bv, L=L)
            tgt = m.with_vertices(y)
            sv = tgt.signed_volumes()
            print(
                f"  {method:10s}: boundary folds {nf:4d}/{ntot:5d}   "
                f"interior flips {(sv <= 0).sum():4d}/{m.n_tets}   residual {res:.2e}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
