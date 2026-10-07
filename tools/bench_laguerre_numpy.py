#!/usr/bin/env python
"""Benchmark: pure-Python/numpy power (Laguerre) cell volumes in 3D.

Strategy (the practical fallback for a 3D OMT solver):
  For each site i, clip the convex domain (here: unit cube, and also a
  conforming tet mesh) against the bisector half-spaces
  { x : 2(pj-pi).x <= |pj|^2 - |pi|^2 - wj + wi }  for all neighbours j,
  then integrate volume/centroid of the resulting convex polyhedron by
  decomposing it into tetrahedra from an interior point.

Only *nearby* sites can produce a facet, so we use a uniform spatial grid
to pre-select candidate bisectors: each Laguerre cell is contained in the
ball B(pi, sqrt(wi - wmin)) so a small stencil suffices.

Run:  python tools/bench_laguerre_numpy.py
"""
from __future__ import annotations

import time
import numpy as np

# ----------------------------------------------------------------------------
# Half-space clipping of a convex polyhedron in 3D.
# The polyhedron is kept as a vertex array V (n,3); we do NOT maintain faces,
# we only need volume + centroid, which we get from a tetrahedral fan around
# the (computed) interior point / a fixed reference point using the divergence
# theorem.
# ----------------------------------------------------------------------------


def clip_polyhedron(V: np.ndarray, E: np.ndarray) -> np.ndarray:
    """Clip convex hull of V by half-spaces  E[:, :3] . x <= E[:, 3].

    Clipping one half-space at a time: a vertex is kept if inside, and each
    edge crossing the plane contributes a new vertex.  This needs the edge
    (convex-hull facet) structure -> we recompute it with the convex hull of
    the current vertices at every step (cheap because counts are small).

    Returns the new vertex array (may be empty).
    """
    from scipy.spatial import ConvexHull

    V = np.asarray(V, dtype=np.float64)
    if len(V) < 4:
        return V
    for e in E:
        n = e[:3]
        d = e[3]
        val = V @ n - d
        inside = val <= 1e-12
        if inside.all():
            continue
        if not inside.any():
            return np.zeros((0, 3))
        # new vertices: intersections of the plane with hull edges
        try:
            hull = ConvexHull(V)
        except Exception:
            return np.zeros((0, 3))
        new_pts = []
        for simplex in hull.simplices:
            idx = list(simplex)
            for a in range(3):
                for b in range(a + 1, 3):
                    ia, ib = idx[a], idx[b]
                    va, vb = val[ia], val[ib]
                    if (va > 0) != (vb > 0):
                        t = va / (va - vb)
                        new_pts.append(V[ia] + t * (V[ib] - V[ia]))
        keep = V[inside]
        if new_pts:
            V = np.vstack([keep, np.asarray(new_pts)])
        else:
            V = keep
        if len(V) < 4:
            return np.zeros((0, 3))
    return V


def poly_volume_centroid(V: np.ndarray):
    """Volume + centroid of convex hull of V via tetrahedral fan from centroid
    of vertices.  Exact for a convex polyhedron: pick an interior reference r,
    sum signed volumes of tetra (r, V[f0],V[f1],V[f2]) over hull facets."""
    from scipy.spatial import ConvexHull

    if len(V) < 4:
        return 0.0, np.zeros(3)
    try:
        hull = ConvexHull(V)
    except Exception:
        return 0.0, np.zeros(3)
    r = V.mean(axis=0)
    vol = 0.0
    cen = np.zeros(3)
    for simplex in hull.simplices:
        a, b, c = V[simplex[0]], V[simplex[1]], V[simplex[2]]
        # outward normal of the facet
        nrm = np.cross(b - a, c - a)
        # orient so that the tetra (r,a,b,c) has positive volume
        v6 = np.dot(nrm, r - a)
        if v6 > 0:  # r is on the positive side -> flip
            nrm = -nrm
        v6 = abs(np.dot(np.cross(b - a, c - a), r - a))
        vol_t = v6 / 6.0
        vol += vol_t
        cen += vol_t * (r + a + b + c) / 4.0
    if vol <= 0:
        return 0.0, np.zeros(3)
    return vol, cen / vol


def power_cell_volumes(P: np.ndarray, W: np.ndarray, domain_min, domain_max,
                       stencil=2, verbose=False):
    """Volumes of the power cells Pow_W(p_i) clipped to an axis-aligned box.

    P : (N,3) sites, W : (N,) weights.  Returns (vol, cen) arrays.
    """
    N = len(P)
    lo = np.asarray(domain_min, dtype=np.float64)
    hi = np.asarray(domain_max, dtype=np.float64)
    box = np.array([
        [1, 0, 0, hi[0]], [-1, 0, 0, -lo[0]],
        [0, 1, 0, hi[1]], [0, -1, 0, -lo[1]],
        [0, 0, 1, hi[2]], [0, 0, -1, -lo[2]],
    ], dtype=np.float64)
    boxV = np.array([[lo[0], lo[1], lo[2]], [hi[0], lo[1], lo[2]],
                     [lo[0], hi[1], lo[2]], [lo[0], lo[1], hi[2]],
                     [hi[0], hi[1], lo[2]], [hi[0], lo[1], hi[2]],
                     [lo[0], hi[1], hi[2]], [hi[0], hi[1], hi[2]]])

    # --- uniform grid for candidate neighbours -----------------------------
    extent = hi - lo
    cell = extent.max() / max(int(round(N ** (1.0 / 3.0))) * 1.0, 1)
    cell = max(cell, 1e-9)
    ncell = np.maximum(np.ceil(extent / cell).astype(int), 1)
    # bucket sites
    bidx = np.clip(((P - lo) / cell).astype(int), 0, ncell - 1)
    from collections import defaultdict
    buckets = defaultdict(list)
    for i in range(N):
        buckets[tuple(bidx[i])].append(i)
    # 27-neighbourhood stencil, widened by `stencil`
    offs = [(dx, dy, dz)
            for dx in range(-stencil, stencil + 1)
            for dy in range(-stencil, stencil + 1)
            for dz in range(-stencil, stencil + 1)]

    vol = np.zeros(N)
    cen = np.zeros((N, 3))
    for i in range(N):
        ci = tuple(bidx[i])
        cand = []
        for o in offs:
            key = (ci[0] + o[0], ci[1] + o[1], ci[2] + o[2])
            cand.extend(buckets.get(key, ()))
        cand = np.asarray(cand, dtype=int)
        cand = cand[cand != i]
        if len(cand) == 0:
            v, c = poly_volume_centroid(boxV)
            vol[i], cen[i] = v, c
            continue
        pi = P[i]
        # bisector half-spaces:  2(pj-pi).x <= |pj|^2-|pi|^2 - wj + wi
        A = 2.0 * (P[cand] - pi)
        b = (np.einsum('ij,ij->i', P[cand], P[cand])
             - pi @ pi - W[cand] + W[i])
        E = np.hstack([A, b[:, None]])
        V = boxV
        # cheap pre-filter: skip planes whose distance is far outside the box
        V = clip_polyhedron(V, E)
        v, c = poly_volume_centroid(V)
        vol[i], cen[i] = v, c
    return vol, cen


def main():
    rng = np.random.default_rng(0)
    print(f"{'N':>9} {'t_total(s)':>11} {'t/cell(ms)':>11} "
          f"{'sum_vol':>10} {'vol_err':>10}")
    for N in (100, 1000, 10000):
        P = rng.random((N, 3))
        W = np.zeros(N)
        # make weights roughly "balanced" so cells are comparable in size
        lo = P.min(axis=0) - 0.05
        hi = P.max(axis=0) + 0.05
        t0 = time.perf_counter()
        vol, cen = power_cell_volumes(P, W, lo, hi, stencil=2)
        t1 = time.perf_counter()
        tot = vol.sum()
        boxvol = float(np.prod(hi - lo))
        print(f"{N:>9} {t1-t0:>11.3f} {1000*(t1-t0)/N:>11.3f} "
              f"{tot:>10.4f} {(tot-boxvol)/boxvol:>10.2e}")


if __name__ == "__main__":
    main()
