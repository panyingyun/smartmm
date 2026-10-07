#!/usr/bin/env python
"""Benchmark pure-Python/numpy Laguerre (power) cell volumes in 3D.

Uses:
  * a uniform spatial grid to select candidate bisector half-spaces
    (only nearby sites can share a facet);
  * half-space clipping of a convex polyhedron (convex-hull edge walk);
  * volume + centroid by a tetrahedral fan from an interior reference point.

Run:  python tools/bench2.py
"""
from __future__ import annotations

import time
from collections import defaultdict

import numpy as np
from scipy.spatial import ConvexHull


def clip(V, E):
    """Clip convex hull of V by half-spaces E[:, :3] @ x <= E[:, 3]."""
    for row in E:
        n, d = row[:3], row[3]
        if np.linalg.norm(n) < 1e-300:
            continue
        val = V @ n - d
        ins = val <= 1e-12
        if ins.all():
            continue
        if not ins.any():
            return None
        hull = ConvexHull(V)
        newp = []
        sp = hull.simplices
        for s in sp:
            for a in range(3):
                ia, ib = s[a], s[a + 1] if a + 1 < 3 else s[0]
                va, vb = val[ia], val[ib]
                if (va > 0) != (vb > 0):
                    t = va / (va - vb)
                    newp.append(V[ia] + t * (V[ib] - V[ia]))
        V = np.vstack([V[ins], np.asarray(newp)]) if newp else V[ins]
        if len(V) < 4:
            return None
    return V


def vol_cen(V, r):
    """Volume + centroid of convex hull of V; tetra fan from interior point r."""
    hull = ConvexHull(V)
    sp = hull.simplices
    a = V[sp[:, 0]]
    nrm = np.cross(V[sp[:, 1]] - a, V[sp[:, 2]] - a)
    v6 = np.einsum('ij,ij->i', nrm, r - a)
    flip = v6 > 0
    v6 = np.abs(v6)
    # ensure the fan is taken from inside: use abs and orient consistently
    vol = v6.sum() / 6.0
    cen = (v6[:, None] * (r + V[sp[:, 0]] + V[sp[:, 1]] + V[sp[:, 2]]) / 4.0).sum(0)
    return vol, cen / vol


def laguerre_volumes(P, W, lo, hi, stencil=2):
    N = len(P)
    ext = hi - lo
    ng = max(int(round(N ** (1 / 3))), 1)
    h = ext / ng
    bidx = np.clip(((P - lo) / h).astype(int), 0, ng - 1)
    b = defaultdict(list)
    for i in range(N):
        b[tuple(bidx[i])].append(i)
    offs = [(x, y, z) for x in range(-stencil, stencil + 1)
            for y in range(-stencil, stencil + 1)
            for z in range(-stencil, stencil + 1)]
    box = np.array([[1., 0, 0, hi[0]], [-1, 0, 0, -lo[0]],
                    [0, 1., 0, hi[1]], [0, -1, 0, -lo[1]],
                    [0, 0, 1., hi[2]], [0, 0, -1, -lo[2]]])
    bv = np.array([[lo[0], lo[1], lo[2]], [hi[0], lo[1], lo[2]],
                   [lo[0], hi[1], lo[2]], [lo[0], lo[1], hi[2]],
                   [hi[0], hi[1], lo[2]], [hi[0], lo[1], hi[2]],
                   [lo[0], hi[1], hi[2]], [hi[0], hi[1], hi[2]]])
    rc = (lo + hi) / 2
    vol = np.zeros(N)
    cen = np.zeros((N, 3))
    nclip = 0
    for i in range(N):
        ci = bidx[i]
        cand = []
        for o in offs:
            cand += b.get((ci[0] + o[0], ci[1] + o[1], ci[2] + o[2]), [])
        cand = np.array([c for c in cand if c != i], dtype=int)
        if len(cand) == 0:
            continue
        A = 2.0 * (P[cand] - P[i])
        rhs = (np.einsum('ij,ij->i', P[cand], P[cand]) - P[i] @ P[i]
               - W[cand] + W[i])
        # prune half-spaces that cannot cut the current domain: keep only
        # planes within the bounding sphere radius of the box
        R = 0.5 * np.linalg.norm(ext)
        dist = (A @ rc - rhs) / np.linalg.norm(A, axis=1)
        keep = dist < R * 0.999
        E = np.hstack([A[keep], rhs[keep, None]])
        nclip += len(E)
        V = clip(bv.copy(), E)
        if V is None or len(V) < 4:
            continue
        try:
            vol[i], cen[i] = vol_cen(V, rc)
        except Exception:
            continue
    return vol, cen, nclip


def main():
    rng = np.random.default_rng(0)
    print(f"{'N':>7} {'t_tot(s)':>9} {'ms/cell':>9} {'sumvol':>9} "
          f"{'err':>9} {'clips/cell':>10}")
    for N in (1000, 10000):
        P = rng.random((N, 3))
        W = np.zeros(N)
        lo = np.zeros(3)
        hi = np.ones(3)
        t0 = time.perf_counter()
        vol, cen, nclip = laguerre_volumes(P, W, lo, hi)
        t1 = time.perf_counter()
        print(f"{N:>7} {t1-t0:>9.2f} {1000*(t1-t0)/N:>9.3f} {vol.sum():>9.5f} "
              f"{vol.sum()-1.0:>9.2e} {nclip/N:>10.1f}")


if __name__ == "__main__":
    main()
