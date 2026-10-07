#!/usr/bin/env python
"""Failure-mode evidence: harmonic / Dirichlet map of a NON-convex 3-ball onto the
unit ball, and the VSEM recovery.

The pipeline's real input is a topologically-3-ball solid that is NOT geometrically
a ball.  This script builds such solids and measures
  [A] radial-projection boundary map + P1 harmonic interior
      -> volume distortion + flipped tets        (Goal A failure mode),
  [B] VSEM fixed point from that same boundary   (Goal B recovery),
  [C] VSEM from an authalic (equal-spherical-area) boundary, to show how much the
      choice of boundary map matters.

Run:  C:\\Python312\\python.exe tools\\demo_failure_modes.py
"""
from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, r"E:\panyingyun\smartmm\tools")
import tet_ball_ref as R  # noqa: E402


def make_blob_tet_mesh(subdivisions: int = 3, max_volume: float = 0.03, kind: str = "peanut"):
    """A non-convex solid whose boundary is a genus-0 sphere.

    kind = 'peanut'  : r(theta) = 1 + 0.45*cos(2 theta)     (dumbbell, concave waist)
    kind = 'wedge'   : r = 1 + 0.35*|sin(3 phi) sin(theta)|
    kind = 'crinkly' : r = 1 + 0.25*sin(5 theta) sin(4 phi)
    """
    import gpytoolbox as gpy
    import tetgen as tg

    V, F = gpy.icosphere(subdivisions)
    V = np.asarray(V, dtype=np.float64)
    x, y, z = V[:, 0], V[:, 1], V[:, 2]
    theta = np.arccos(np.clip(z, -1, 1))
    phi = np.arctan2(y, x)
    if kind == "peanut":
        r = 1.0 + 0.45 * np.cos(2 * theta)
    elif kind == "wedge":
        r = 1.0 + 0.35 * np.abs(np.sin(3 * phi) * np.sin(theta))
    elif kind == "crinkly":
        r = 1.0 + 0.25 * np.sin(5 * theta) * np.sin(4 * phi)
    else:
        raise ValueError(kind)
    r = np.maximum(r, 0.15)
    V = V * r[:, None]

    tgen = tg.TetGen(V, np.asarray(F, dtype=np.int32))
    tgen.tetrahedralize(switches=f"pq1.414a{max_volume}Y")
    T = np.asarray(tgen.elem, dtype=np.int64)
    Vt = np.asarray(tgen.node, dtype=np.float64)
    return Vt, T


def spherical_tri_areas(S: np.ndarray, bf: np.ndarray) -> np.ndarray:
    """SIGNED spherical (solid-angle) area of the spherical triangles S[bf].

    Van Oosterom & Strackee (IEEE TBMed 30(2):125-126, 1983) general form:
        tan(Omega/2) = [a b c] / (|a||b||c| + a.b + b.c + c.a),
    with |a|=|b|=|c|=1 on the unit sphere:
        Omega = 2 * atan2( a.(b x c), 1 + a.b + b.c + c.a ).

    Omega > 0 for counter-clockwise (outward) winding and < 0 for clockwise, with
    |Omega| <= 2*pi.  A consistently outward-wound triangulation of the whole
    unit sphere sums to exactly 4*pi; a double cover sums to 8*pi.  The signed
    form is essential -- taking |Omega| hides double covers and folded
    triangles, which is precisely what you need to detect.
    """
    a = S[bf[:, 0]]
    b = S[bf[:, 1]]
    c = S[bf[:, 2]]
    det = np.einsum("ij,ij->i", a, np.cross(b, c))
    den = (1.0
           + np.einsum("ij,ij->i", a, b)
           + np.einsum("ij,ij->i", b, c)
           + np.einsum("ij,ij->i", c, a))
    return 2.0 * np.arctan2(det, den)


def authalic_sphere_boundary(
    V: np.ndarray, T: np.ndarray, n_iter: int = 120, step: float = 0.35
) -> np.ndarray:
    """Cheap spherical equal-area boundary map (stand-in for SEM / SAEM).

    Alternates (i) moving each spherical vertex toward the area-weighted centroid
    of its spherical neighbours (tangentially, then re-normalised to the sphere)
    and (ii) a density-equalising relocation driven by the ratio of the achieved
    spherical area to the target.  Not a substitute for the published SEM/SAEM
    solvers, but it produces a much lower-distortion boundary map than naive
    radial projection.
    """
    n = len(V)
    b, bf = R.boundary_vertices(T, n, V)
    pos = -np.ones(n, dtype=int)
    pos[b] = np.arange(len(b))

    S = V[b] - V[b].mean(axis=0)
    S = S / np.linalg.norm(S, axis=1, keepdims=True)

    # boundary edge list
    e = np.vstack([bf[:, [0, 1]], bf[:, [1, 2]], bf[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    li, lj = pos[e[:, 0]], pos[e[:, 1]]
    keep = (li >= 0) & (lj >= 0)
    li, lj = li[keep], lj[keep]

    nb = len(b)
    for it in range(n_iter):
        # target = equal spherical area per vertex
        area = spherical_tri_areas(S, bf)
        w = np.zeros(nb)
        np.add.at(w, pos[bf[:, 0]], area / 3.0)
        np.add.at(w, pos[bf[:, 1]], area / 3.0)
        np.add.at(w, pos[bf[:, 2]], area / 3.0)
        target = 4.0 * math.pi / nb
        dens = np.maximum(w, 1e-12) / target
        # gradient of density on the sphere -> tangential displacement field
        g = np.zeros((nb, 3))
        np.add.at(g, li, (dens[lj] - dens[li])[:, None] * (S[lj] - S[li]))
        np.add.at(g, lj, (dens[li] - dens[lj])[:, None] * (S[li] - S[lj]))
        # project onto tangent plane and step
        g = g - np.einsum("ij,ij->i", g, S)[:, None] * S
        S = S - step * g
        S = S / np.linalg.norm(S, axis=1, keepdims=True)
    return S


def main() -> int:
    for kind in ("peanut", "wedge", "crinkly"):
        print("=" * 80)
        print(f"shape = {kind}")
        V, T = make_blob_tet_mesh(3, 0.03, kind)
        n = len(V)
        b, bf = R.boundary_vertices(T, n, V)
        mu = R.tet_volumes(V, T)
        print(f"  n_verts={n} (bnd {len(b)}, int {n - len(b)})  n_tets={len(T)}  vol={mu.sum():.5f}")

        # --- naive boundary: radial projection onto the sphere
        Srad = V[b] - V[b].mean(axis=0)
        Srad = Srad / np.linalg.norm(Srad, axis=1, keepdims=True)
        ar = spherical_tri_areas(Srad, bf)
        print(f"  radial BC: spherical tri area min={ar.min():.3e} max={ar.max():.3e} "
              f"std/mean={ar.std() / ar.mean():.3f}  sum={ar.sum():.5f} (4pi={4 * math.pi:.5f})")

        # [A] harmonic interior with the naive boundary
        U, _, _ = R.harmonic_ball_map(V, T, Srad)
        r = R.volume_distortion_ratios(V, T, U)
        print(f"  [A] harmonic + radial BC : ratio mean={r.mean():.4f} std={r.std():.4f} "
              f"min={r.min():.4e} max={r.max():.4f}")
        print(f"      flipped={int((r <= 0).sum())}/{len(T)} ({float((r <= 0).mean()):.3%})"
              f"   J<0.5={float((r < 0.5).mean()):.3%}")

        # [B] VSEM from the naive boundary
        f, hist = R.vsem_solve(V, T, Srad, n_iter=80, tol=1e-13)
        r2 = R.volume_distortion_ratios(V, T, f)
        gap = R.vsem_energy(V, T, f) - 1.5 * float(R.tet_volumes(f, T).sum())
        print(f"  [B] VSEM (radial BC)     : ratio mean={r2.mean():.4f} std={r2.std():.4f} "
              f"min={r2.min():.6f} max={r2.max():.6f}")
        print(f"      flipped={int((r2 <= 0).sum())}   E_V-1.5V(f)={gap:.3e}  iters={len(hist)}")

        # [C] VSEM from a *published* authalic/conformal spherical boundary map is
        # NOT attempted here: the naive density-equalising relaxation that earlier
        # lived in this file is numerically unstable (it collapses the spherical
        # triangles) and is not a substitute for the real solvers.  Use instead
        #   - MATLAB: SphericalAEM(F,V)      (SAEM, Liu & Yueh, doi:10.1137/25M1736979)
        #   - MATLAB: spherical_conformal_map(v,f)  (Choi-Lam-Lui, doi:10.1137/130950008)
        #   - C++:    RiemannMapper / Gu group spherical harmonic map
        print("  [C] authalic/conformal spherical BC: NOT attempted here "
              "(use SphericalAEM.m or spherical_conformal_map.m, see report sec. 2)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
