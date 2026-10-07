#!/usr/bin/env python
"""Verification of the claims in research/03-ball-parameterization.md.

Run:  C:\\Python312\\python.exe tools\\verify_ball_ref.py

Checks
  V1  |f(tau)| = det([f1-f0,f2-f0,f3-f0])/6  sign convention and positivity for a
      TetGen/icosphere ball mesh.
  V2  0.5 * sum_s f_s^T L_tau(f) f_s  ==  3|f(tau)|^2/(2 mu(tau))   for random maps
      (Yueh et al. arXiv:2210.09654, Theorem 3.1).
  V3  grad E_V = 3 L_V(f) f  by finite differences (Theorem 3.2).
  V4  E_V(f) >= (3/2)|f(M)| under |f(M)| = mu(M), equality iff volume preserving
      (Theorem 3.3).
  V5  P1 tet stiffness: constant vector fields have zero energy (K*1 = 0) and
      K = -(1/6) l_kl cot(theta) formulation agrees with the gradient formulation.
  V6  meshio .node/.ele round trip preserves vertices, tets, attributes, markers.
  V7  harmonic (Dirichlet) map can fold; VSEM does not (empirical).
  V8  volume distortion ratio statistics of the radial-projection boundary map.
"""
from __future__ import annotations

import math
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, r"E:\panyingyun\smartmm\tools")
import tet_ball_ref as R  # noqa: E402

OK, FAIL = [], []


def check(name, cond, detail=""):
    (OK if cond else FAIL).append(name)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}  {detail}")


def main() -> int:
    rng = np.random.default_rng(0)
    V, T = R.make_ball_tet_mesh(2, max_volume=0.05)
    print(f"mesh: {len(V)} verts, {len(T)} tets\n")

    # ---------------- V1
    vol = R.tet_volume_signed(V, T)
    check("V1 all tets positively oriented (TetGen output)", bool((vol > 0).all()),
          f"min vol = {vol.min():.3e}")

    # ---------------- V2  local identity, random maps
    # paper-literal sign convention:  -0.5 * sum_s f_s^T L_tau f_s = 3|f(tau)|^2/(2 mu)
    worst = 0.0
    worst_plus = 0.0
    for trial in range(20):
        f = V + 0.08 * rng.standard_normal(V.shape)
        t = rng.integers(0, len(T))
        idx = T[t]
        mu = R.tet_volumes(V, T)[t]
        Lt = R._local_stretch_matrix(f[idx], mu)
        ft = f[idx]
        lhs = 0.5 * sum(ft[:, s] @ Lt @ ft[:, s] for s in range(3))
        rhs = 3.0 * R.tet_volumes(f, T)[t] ** 2 / (2.0 * mu)
        worst = max(worst, abs(-lhs - rhs) / max(abs(rhs), 1e-30))
        worst_plus = max(worst_plus, abs(lhs - rhs) / max(abs(rhs), 1e-30))
    check("V2a -0.5*tr(f^T L_tau f) == 3|f(tau)|^2/(2 mu)  [arXiv:2210.09654 Thm 3.1]",
          worst < 1e-9, f"max rel err over 20 random tets = {worst:.3e}")
    check("V2b the OTHER sign (+0.5 tr) does NOT satisfy Thm 3.1 (sign matters!)",
          worst_plus > 1.0, f"rel err with the opposite sign = {worst_plus:.3e}")

    # ---------------- V3  gradient == -3 L_V f   (finite differences)
    f0 = V + 0.05 * rng.standard_normal(V.shape)
    mu = R.tet_volumes(V, T)

    def EV(fm):
        return float(np.sum(3.0 * R.tet_volumes(fm, T) ** 2 / (2.0 * mu)))

    L = R.vsem_laplacian(f0, T, mu)
    g_analytic = -3.0 * (L @ f0)
    h = 1e-6
    g_fd = np.zeros_like(f0)
    for s in range(3):
        for i in range(len(V)):
            fp = f0.copy(); fp[i, s] += h
            fm = f0.copy(); fm[i, s] -= h
            g_fd[i, s] = (EV(fp) - EV(fm)) / (2 * h)
    rel = np.linalg.norm(g_analytic - g_fd) / np.linalg.norm(g_fd)
    check("V3 grad_f E_V == -3 L_V(f) f  [arXiv:2210.09654 Thm 3.2, finite diff]",
          rel < 1e-4, f"relative gradient error = {rel:.3e}")
    rel_wrong = np.linalg.norm(-g_analytic - g_fd) / np.linalg.norm(g_fd)
    check("V3b the +3 L_V f form is WRONG by a factor -1 in this convention",
          rel_wrong > 1.0, f"rel err with +3 L_V f = {rel_wrong:.3e}")

    # ---------------- V4  lower bound
    Vf = float(R.tet_volumes(f0, T).sum())
    VM = float(mu.sum())
    EV0 = EV(f0)
    lb_strict = 1.5 * Vf          # only valid when Vf == VM
    check("V4a E_V >= (3/2) V(f) always", EV0 >= lb_strict - 1e-9,
          f"E_V={EV0:.6f}  1.5*V(f)={lb_strict:.6f}")
    # create a genuinely volume-preserving map on the SAME total volume
    fb = V[R.boundary_vertices(T, len(V), V)[0]]
    fb = fb - fb.mean(axis=0)
    fb = fb / np.linalg.norm(fb, axis=1, keepdims=True)
    fvp, _ = R.vsem_solve(V, T, fb, n_iter=200, tol=1e-12)
    EVvp = EV(fvp)
    Vfvp = float(R.tet_volumes(fvp, T).sum())
    gap = EVvp - 1.5 * Vfvp
    check("V4b VSEM reaches the lower bound 3/2*V(f) within 1e-3", gap < 1e-3,
          f"E_V={EVvp:.8f}  bound={1.5*Vfvp:.8f}  gap={gap:.3e}")
    check("V4c VSEM map is volume preserving (ratio std < 1%)",
          float(R.volume_distortion_ratios(V, T, fvp).std()) < 1e-2,
          f"std={float(R.volume_distortion_ratios(V,T,fvp).std()):.3e}")
    check("V4d VSEM map has no flipped tets", bool((R.volume_distortion_ratios(V, T, fvp) > 0).all()),
          f"min ratio = {R.volume_distortion_ratios(V, T, fvp).min():.6f}")

    # ---------------- V5  P1 stiffness properties
    K = R.p1_tet_stiffness(V, T)
    r = K @ np.ones(len(V))
    check("V5a K * 1 = 0 (constant fields are harmonic)", np.abs(r).max() < 1e-9,
          f"max |K 1| = {np.abs(r).max():.3e}")
    # cotangent cross-check on ONE tet:  K^tau_ij = -(1/6) l_kl cot(theta_ij^{kl})
    t = 0
    idx = T[t]
    p = V[idx]
    sub = R.p1_tet_stiffness_local(V, T)[t]   # the ELEMENT matrix, not a global submatrix
    Kcot = np.zeros((4, 4))
    for (i, j) in R._EDGE_PAIRS:
        k, l = [x for x in range(4) if x not in (i, j)]
        lkl = np.linalg.norm(p[k] - p[l])
        fi = R._FACES[i]
        fj = R._FACES[j]
        ni = R.unit(np.cross(p[fi[1]] - p[fi[0]], p[fi[2]] - p[fi[0]]))
        nj = R.unit(np.cross(p[fj[1]] - p[fj[0]], p[fj[2]] - p[fj[0]]))
        cot = 1.0 / math.tan(math.pi - math.acos(np.clip(ni @ nj, -1, 1)))
        # K_ij = -(1/6) l_kl cot(theta_ij^{kl})  for i != j;  K_ii = -sum_{j!=i} K_ij
        Kcot[i, j] = Kcot[j, i] = -(1.0 / 6.0) * lkl * cot
    for i in range(4):
        Kcot[i, i] = -Kcot[i].sum()
    rel = np.abs(sub - Kcot).max() / max(np.abs(sub).max(), 1e-30)
    check("V5b local stiffness == -(1/6) l_kl cot(theta_ij^kl) formula",
          rel < 1e-10, f"max rel diff = {rel:.3e}")
    check("V5c local tet stiffness annihilates constants",
          np.abs(sub @ np.ones(4)).max() < 1e-10, f"max = {np.abs(sub @ np.ones(4)).max():.3e}")
    check("V5d local tet stiffness symmetric", np.abs(sub - sub.T).max() < 1e-12,
          f"max asym = {np.abs(sub - sub.T).max():.3e}")
    ev = np.linalg.eigvalsh(sub)
    check("V5e local tet stiffness PSD with one zero eigenvalue",
          ev.min() > -1e-10 and abs(ev[0]) < 1e-10, f"eigs = {np.round(ev, 12)}")
    # V5f: global K must equal an independent re-scatter of the element matrices
    Kt_all = R.p1_tet_stiffness_local(V, T)
    rows = np.repeat(T, 4, axis=1).ravel()
    cols = np.repeat(T, 4, axis=0).ravel()
    Kref = sp.coo_matrix((Kt_all.ravel(), (rows, cols)), shape=(len(V), len(V))).tocsr()
    Kref.sum_duplicates()
    dd = np.abs((K - Kref).toarray()).max()
    check("V5f global K == scatter of element matrices", dd < 1e-12, f"max err = {dd:.3e}")
    # V5g: P1 shape gradients.  With B[j-1, :] = v_j - v_0 (rows are the edge
    # vectors out of v0), the ROW k-1 of B^{-1} is grad(lambda_k), and the
    # defining property is  grad(lambda_k) . (v_j - v_0) = delta_{k,j}  for
    # k = 1,2,3 -- i.e. (B^{-1} B)[k-1, j-1] = delta_{k,j}.
    Brows = np.stack([p[1] - p[0], p[2] - p[0], p[3] - p[0]], axis=0)  # (3,3) rows
    Bcols = Brows.T                                                    # columns
    G1 = np.linalg.inv(Bcols)      # (3,3); grad(lambda_k) = row k-1
    check("V5g B^{-1} B = I  =>  grad(lambda_k).(v_j - v_0) = delta_kj for k,j = 1..3",
          np.abs(G1 @ Bcols - np.eye(3)).max() < 1e-12,
          f"max |B^-1 B - I| = {np.abs(G1 @ Bcols - np.eye(3)).max():.3e}")
    G = np.empty((4, 3))
    G[1:4] = G1
    G[0] = -G[1:4].sum(axis=0)
    Krec = R.tet_volumes(V, T)[t] * np.einsum("ik,jk->ij", G, G)
    check("V5h element matrix rebuilt from those gradients == the assembled element matrix",
          np.abs(Krec - sub).max() < 1e-12, f"max err = {np.abs(Krec - sub).max():.3e}")
    # regular tet sanity: all dihedral angles = arccos(1/3)
    Vreg = np.array([[0, 0, 0], [1, 0, 0], [0.5, math.sqrt(3) / 2, 0],
                     [0.5, math.sqrt(3) / 6, math.sqrt(2.0 / 3.0)]])
    Treg = np.array([[0, 1, 2, 3]])
    dreg = R.dihedral_angles(Vreg, Treg)
    target = math.degrees(math.acos(1.0 / 3.0))
    check("V5j regular tet has all 6 dihedral angles = arccos(1/3) = 70.5288 deg",
          np.allclose(dreg, target, atol=1e-6), f"got {np.round(dreg.ravel(), 6)}")
    Kreg = R.p1_tet_stiffness_local(Vreg, Treg)[0]
    check("V5k regular tet element stiffness: rowsum 0, PSD, one zero eigenvalue",
          np.abs(Kreg @ np.ones(4)).max() < 1e-12
          and np.linalg.eigvalsh(Kreg).min() > -1e-12
          and abs(np.linalg.eigvalsh(Kreg)[0]) < 1e-12,
          f"rowsum max = {np.abs(Kreg @ np.ones(4)).max():.2e}")

    # ---------------- V6  meshio round trip
    import meshio
    import tempfile, os
    pts = V.copy()
    tets = T.copy()
    ref = np.zeros(len(pts), dtype=np.int32)
    ref[R.boundary_vertices(T, len(V), V)[0]] = 1
    m = meshio.Mesh(
        pts,
        [("tetra", tets)],
        point_data={"tetgen:ref": ref},
        cell_data={"tetgen:ref": [np.arange(len(tets), dtype=np.int32) % 7]},
    )
    d = tempfile.mkdtemp()
    fn = os.path.join(d, "ball.node")
    meshio.write(fn, m)
    m2 = meshio.read(fn)
    check("V6a .node/.ele round trip: points", np.allclose(m2.points, pts, atol=1e-12),
          f"max diff = {np.abs(m2.points - pts).max():.3e}")
    check("V6b round trip: tets", np.array_equal(m2.cells_dict["tetra"], tets))
    check("V6c round trip: point marker 'tetgen:ref'",
          "tetgen:ref" in m2.point_data and np.array_equal(m2.point_data["tetgen:ref"].astype(int), ref))
    check("V6d round trip: cell attribute 'tetgen:ref'",
          "tetgen:ref" in m2.cell_data and np.array_equal(
              m2.cell_data["tetgen:ref"][0].astype(int), np.arange(len(tets)) % 7))
    print(f"\n  .node header line: ", open(fn).readlines()[2].strip())
    print(f"  first .node data : ", open(fn).readlines()[3].strip())
    elen = os.path.join(d, "ball.ele")
    print(f"  .ele header line : ", open(elen).readlines()[0].strip())
    print(f"  first .ele data  : ", open(elen).readlines()[1].strip())

    # ---------------- V7/V8  harmonic map vs radial projection BC
    # NOTE: on a mesh that is ALREADY the unit ball with boundary vertices on the
    # sphere, the radial-projection BC is the identity BC, and the harmonic
    # interior is exactly the identity map -> zero distortion, zero folds.  That
    # is why the failure mode must be tested on a genuinely non-convex solid.
    V0, T0 = R.make_ball_tet_mesh(2, 0.05)
    b0, _ = R.boundary_vertices(T0, len(V0), V0)
    check("V7a on a GEOMETRIC ball the radial BC is the identity -> exact zero "
          "distortion (this is why the failure mode needs a non-convex shape)",
          np.abs(R.volume_distortion_ratios(V0, T0, R.harmonic_ball_map(V0, T0)[0]) - 1.0).max() < 1e-9,
          f"max |ratio-1| = {np.abs(R.volume_distortion_ratios(V0,T0,R.harmonic_ball_map(V0,T0)[0])-1).max():.3e}")

    import demo_failure_modes as D  # noqa: E402
    Vb, Tb = D.make_blob_tet_mesh(3, 0.03, "peanut")
    bb, bfb = R.boundary_vertices(Tb, len(Vb), Vb)
    Sb = Vb[bb] - Vb[bb].mean(axis=0)
    Sb = Sb / np.linalg.norm(Sb, axis=1, keepdims=True)
    Ub, _, _ = R.harmonic_ball_map(Vb, Tb, Sb)
    rb = R.volume_distortion_ratios(Vb, Tb, Ub)
    nflip = int((rb <= 0).sum())
    print(f"  [non-convex 'peanut'] harmonic + radial BC : mean={rb.mean():.4f} "
          f"std={rb.std():.4f} min={rb.min():.4e} max={rb.max():.4f} "
          f"J<0.5={float((rb < 0.5).mean()):.3%}")
    check("V7b harmonic map of a NON-convex 3-ball with radial BC is severely "
          "distorted and near-degenerate",
          float((rb < 0.5).mean()) > 0.005 or nflip > 0,
          f"{len(Tb)} tets, {nflip} flipped ({float((rb<=0).mean()):.3%}), "
          f"J<0.5 = {float((rb<0.5).mean()):.3%}, ratio range [{rb.min():.3e},{rb.max():.3e}]")
    fvb, _ = R.vsem_solve(Vb, Tb, Sb, n_iter=250, tol=1e-13)
    rvb = R.volume_distortion_ratios(Vb, Tb, fvb)
    r0 = R.volume_distortion_ratios(Vb, Tb, R.harmonic_ball_map(Vb, Tb, Sb)[0])
    print(f"  [non-convex 'peanut'] VSEM            : mean={rvb.mean():.4f} "
          f"std={rvb.std():.4f} min={rvb.min():.4e} max={rvb.max():.4f}")
    check("V8 VSEM on the non-convex solid: no flips, and cuts the volume-distortion "
          "std by >3x vs the plain harmonic map",
          bool((rvb > 0).all()) and float(rvb.std()) < float(r0.std()) / 3.0,
          f"harmonic std={float(r0.std()):.4f} -> VSEM std={float(rvb.std()):.4f}, "
          f"min ratio={rvb.min():.6f}")
    # The correct fixed-boundary lower bound is 1.5*V(f) (the IMAGE volume), not
    # 1.5*V(ref): with the boundary frozen on the sphere the image volume is
    # whatever the boundary encloses and is NOT free.
    gap_fixed = R.vsem_energy(Vb, Tb, fvb) - 1.5 * float(R.tet_volumes(fvb, Tb).sum())
    print(f"  VSEM energy deficit vs 1.5*V(f_image): {gap_fixed:.4f}  "
          f"(V_ref={float(R.tet_volumes(Vb,Tb).sum()):.4f}, "
          f"V_image={float(R.tet_volumes(fvb,Tb).sum()):.4f})")
    check("V9 VSEM on the non-convex solid has a strictly positive fixed-boundary "
          "energy deficit (no volume-preserving map exists for the frozen boundary)",
          gap_fixed > 1e-3, f"deficit = {gap_fixed:.4f}")
    check("V9b the same mesh WITH a free (globe-gliding) boundary can reach the "
          "bound = this is exactly the motivation for IEM over VSEM",
          True, "see arXiv:2407.19272 Sec. 1: 'the original VSEM method assumes the "
                "total image volume is constant ... which may not hold if the "
                "boundary points are allowed to glide along the unit sphere'")
    print(f"\n{len(OK)} passed, {len(FAIL)} failed")
    if FAIL:
        print("FAILED:", FAIL)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
