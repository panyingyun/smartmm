"""Validate the semi-discrete OT solver on small, fully controlled problems."""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from tetparam.generate import fibonacci_sphere  # noqa: E402
from tetparam.laguerre import ConvexPolyhedron, solve_ot, verify_neighbour_cutoff  # noqa: E402


def sites_in_own_cells(sites, weights):
    """Number of sites that are inside their own power cell (must be all of them)."""
    good = 0
    for i in range(len(sites)):
        d = np.sum((sites - sites[i]) ** 2, axis=1) - weights + weights[i]
        good += int(np.argmin(d) == i)
    return good


def main() -> int:
    # convex polyhedron: convex hull of quasi-uniform points on the sphere
    P = fibonacci_sphere(200)
    dom = ConvexPolyhedron.from_points(P)
    print(f"domain: |P| = {dom.volume():.6f}  ({len(dom.planes)} planes)")

    rng = np.random.default_rng(7)
    for n in (12, 40, 120):
        idx = rng.choice(len(P), n, replace=False)
        sites = P[idx] * 0.9
        r = rng.random(n) + 0.3
        masses = r / r.sum() * dom.volume()

        t0 = time.perf_counter()
        res = solve_ot(sites, masses, dom, neighbours=min(n - 1, 48), max_iter=60, tol=1e-8)
        dt = time.perf_counter() - t0
        rel = np.abs(res.mass_error) / masses.mean()
        print(
            f"n={n:4d}: {res.iterations:3d} it, max rel mass err {rel.max():.3e}, "
            f"sum vol {res.diagram.volumes.sum():.6f}/{dom.volume():.6f}, "
            f"empty {res.diagram.n_empty}, sites-in-own-cell "
            f"{sites_in_own_cells(sites, res.weights)}/{n}, "
            f"planes {res.diagram.n_planes_used}  ({dt:.2f} s)"
        )
        print(f"        history tail: {['%.2e' % h for h in res.history[-4:]]}")

    n = 150
    sites = P[rng.choice(len(P), n, replace=False)] * 0.85
    w = np.zeros(n)
    for k in (12, 24, 32, 48):
        err = verify_neighbour_cutoff(sites, w, dom, k, n_samples=8)
        print(f"neighbour cutoff k={k:3d}: max relative volume error {err:.3e}")

    # two sites: an exactly solvable check (a plane at the right offset)
    sites = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])
    m = np.array([0.25, 0.75]) * dom.volume()
    res = solve_ot(sites, m, dom, neighbours=1, max_iter=300, tol=1e-12, use_newton=False)
    print(
        f"2-site check: vol {res.diagram.volumes[0]:.8f} (want {m[0]:.8f}) / "
        f"{res.diagram.volumes[1]:.8f} (want {m[1]:.8f})  it {res.iterations}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
