#!/usr/bin/env python
"""Reference building blocks for the ball-parameterization pipeline (FIRST HALF).

Verified against: numpy 2.4.4, scipy 1.18.0, meshio 5.3.5, tetgen 0.8.4, pyvista 0.49.0
Run:  C:\\Python312\\python.exe tools\\tet_ball_ref.py --demo

Contents
--------
 1.  make_ball_tet_mesh(n)             analytical-ish tet mesh of unit ball (icosphere + tetgen)
 2.  tet_volume / tet_volumes          signed & absolute tet volumes
 3.  dihedral_angles                   all 6 dihedral angles per tet
 4.  radius_edge_ratio / aspect_ratio  standard tet quality metrics
 5.  vsem_laplacian(V, T, mu=None)     volumetric stretch Laplacian L_V(f)  (Yueh et al. 2019/2023)
 6.  harmonic_ball_map(V, T)           Goal A: P1 harmonic map to unit ball
 7.  vsem_solve(V, T, fb, n_iter)      Goal B: VSEM fixed-point
 8.  volume_distortion_ratios          |f(tau)|/|tau| stats
 9.  jacobian_det_field                J per tet + fraction flipped
10.  verify_map                          bundle of all verification metrics
"""
from __future__ import annotations

import argparse
import math
from typing import Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ----------------------------------------------------------------------------
# 1. Ball shaped tetrahedral mesh
# ----------------------------------------------------------------------------


def make_ball_tet_mesh(subdivisions: int = 3, max_volume: float = 0.0):
    """Tetrahedralize the unit ball.

    Builds an icosphere (gpytoolbox) and calls TetGen (via the `tetgen` pip
    package, which wraps the TetGen C++ library) with the -p (PLC) and -q
    (quality) switches.  Returns (V, T) with V:(n,3) float64, T:(m,4) int32.

    Note: TetGen's `-q` quality switch on a *sphere* surface can refuse to
    insert interior points unless you also pass a max volume; pass
    max_volume>0 (e.g. 0.02) or use `add_region`.
    """
    import gpytoolbox as gpy
    import tetgen as tg

    V, F = gpy.icosphere(subdivisions)  # unit sphere already
    V = np.asarray(V, dtype=np.float64)
    F = np.asarray(F, dtype=np.int32)

    tgen = tg.TetGen(V, F)
    opts = "pq1.414Y" if max_volume <= 0 else f"pq1.414a{max_volume}Y"
    tgen.tetrahedralize(switches=opts)
    T = np.asarray(tgen.elem, dtype=np.int64)
    Vt = np.asarray(tgen.node, dtype=np.float64)
    return Vt, T


# ----------------------------------------------------------------------------
# 2. Tet geometry
# ----------------------------------------------------------------------------


def tet_volume_signed(V: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Signed volume of each tet:  det([v1-v0, v2-v0, v3-v0]) / 6."""
    p = V[T]  # (m,4,3)
    a = p[:, 1] - p[:, 0]
    b = p[:, 2] - p[:, 0]
    c = p[:, 3] - p[:, 0]
    return np.einsum("ij,ij->i", np.cross(a, b), c) / 6.0


def tet_volumes(V: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Absolute volumes |tau|."""
    return np.abs(tet_volume_signed(V, T))


def orient_tets(V: np.ndarray, T: np.ndarray, in_place: bool = False):
    """Canonicalise every tet to POSITIVE signed volume by swapping two vertices.

    TetGen (and most meshers) may emit tetrahedra in either orientation; every
    face-based formula below (outward face winding, dihedral angles, the
    modified cotangent weights) assumes a consistent positive orientation, so
    normalise ONCE at load time.

    Returns (T_oriented, n_swapped).
    """
    T = T if in_place else T.copy()
    sv = tet_volume_signed(V, T)
    neg = sv < 0
    T[neg] = T[neg][:, [0, 1, 3, 2]]  # local permutation; det flips sign
    return T, int(neg.sum())


def canonicalize(V: np.ndarray, T: np.ndarray):
    """(T_pos_oriented, boundary_vertex_ids, boundary_faces_outward_wound)."""
    T, _ = orient_tets(V, T)
    return T, *boundary_vertices(T, len(V), V)


# _FACES[k] is the face opposite vertex k, wound so that its outward normal
# (right-hand rule on the listed order) points AWAY from vertex k whenever
# det[v1-v0, v2-v0, v3-v0] > 0.
#   opposite v0 = [v1,v2,v3],  opposite v1 = [v0,v3,v2],
#   opposite v2 = [v0,v1,v3],  opposite v3 = [v0,v2,v1]
_FACES = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
# For an edge (i,j) the two incident faces are _FACES[i] and _FACES[j], and their
# outward normals point toward the interior of the tet on the two sides of the
# edge, so the interior dihedral angle theta satisfies
#     cos(pi - theta) = n_i . n_j       (equivalently cos theta = -n_i . n_j).
_EDGE_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def face_opposite(k: int) -> np.ndarray:
    """The face of a positively-oriented tet that is opposite local vertex k,
    outward-wound (see _FACES)."""
    return _FACES[k]


def dihedral_angles(V: np.ndarray, T: np.ndarray) -> np.ndarray:
    """All 6 interior dihedral angles of every tet, in degrees, in the order of
    _EDGE_PAIRS.

    For edge (i,j) the two incident faces are face_opposite(i) and
    face_opposite(j); their outward-wound normals point into the tet on the two
    sides of the edge, hence  cos(pi - theta_ij) = n_i . n_j , i.e.
        theta_ij = pi - arccos(clip(n_i . n_j, -1, 1)) .
    A regular tetrahedron gives theta = arccos(1/3) = 70.5288 deg for all 6 edges.
    Returns degrees in [0, 180].
    """
    p = V[T]  # (m,4,3)
    out = np.empty((len(T), 6), dtype=np.float64)
    for k, (i, j) in enumerate(_EDGE_PAIRS):
        fi = _FACES[i]
        fj = _FACES[j]
        ni = unit(np.cross(p[:, fi[1]] - p[:, fi[0]], p[:, fi[2]] - p[:, fi[0]]))
        nj = unit(np.cross(p[:, fj[1]] - p[:, fj[0]], p[:, fj[2]] - p[:, fj[0]]))
        cosang = np.clip(np.einsum("ij,ij->i", ni, nj), -1.0, 1.0)
        out[:, k] = np.degrees(np.pi - np.arccos(cosang))
    return out


def unit(a: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(a, axis=-1, keepdims=True)
    n[n == 0] = 1.0
    return a / n


def radius_edge_ratio(V: np.ndarray, T: np.ndarray) -> np.ndarray:
    """R = circumradius / shortest edge length  (1.5 * sqrt(6) ~ 3.674 is ideal for
    regular tets -> equilateral; smaller is better).  Uses the 3x3 linear solve
    for the circumcenter."""
    p = V[T]
    # circumcenter of a tet: solve for c from |c-p0|=|c-p1|=|c-p2|=|c-p3|
    A = 2.0 * np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=1)
    rhs = np.stack(
        [
            np.einsum("ij,ij->i", p[:, 1], p[:, 1]) - np.einsum("ij,ij->i", p[:, 0], p[:, 0]),
            np.einsum("ij,ij->i", p[:, 2], p[:, 2]) - np.einsum("ij,ij->i", p[:, 0], p[:, 0]),
            np.einsum("ij,ij->i", p[:, 3], p[:, 3]) - np.einsum("ij,ij->i", p[:, 0], p[:, 0]),
        ],
        axis=1,
    )
    c = np.linalg.solve(A, rhs[:, :, None])[:, :, 0]
    R = np.linalg.norm(c - p[:, 0], axis=1)

    L = np.empty((len(T), 6))
    for k, (i, j) in enumerate(_EDGE_PAIRS):
        L[:, k] = np.linalg.norm(p[:, i] - p[:, j], axis=1)
    R = np.abs(R)
    R = np.where(R == 0, 1e-30, R)
    return R / L.min(axis=1)


def aspect_ratio(V: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Longest edge / (2*sqrt(2)*inradius); 1.0 = regular tet, larger = worse."""
    p = V[T]
    L = np.empty((len(T), 6))
    for k, (i, j) in enumerate(_EDGE_PAIRS):
        L[:, k] = np.linalg.norm(p[:, i] - p[:, j], axis=1)
    vol = np.abs(tet_volume_signed(V, T))
    area = np.empty((len(T), 4))
    for k, f in enumerate(_FACES):
        area[:, k] = 0.5 * np.linalg.norm(
            np.cross(p[:, f[1]] - p[:, f[0]], p[:, f[2]] - p[:, f[0]]), axis=1
        )
    r_in = 3.0 * vol / area.sum(axis=1)
    r_in = np.where(r_in <= 0, 1e-30, r_in)
    return L.max(axis=1) / (2.0 * math.sqrt(2.0) * r_in)


def dihedral_min_deg(V: np.ndarray, T: np.ndarray) -> np.ndarray:
    return dihedral_angles(V, T).min(axis=1)


# ----------------------------------------------------------------------------
# 3. Volumetric stretch Laplacian  L_V(f)   (Yueh-Li-Lin-Yau 2019; Huang et al. 2023)
# ----------------------------------------------------------------------------


def _local_stretch_matrix(f: np.ndarray, mu_tau: float) -> np.ndarray:
    """L_tau(f) in R^{4x4}, the element volumetric-stretch Laplacian, in the
    LITERAL convention of Yueh-Li-Lin-Yau, arXiv:2210.09654, eqs. (2.4)/(8.1)/(8.2):

        L_tau(f) = -(1/(36 mu(tau))) [a_ij],
        a_ij = [(f_k-f_i) x (f_l-f_i)] . [(f_l-f_j) x (f_k-f_j)]
               for {i,j,k,l} = {0,1,2,3} with [i,j] disjoint from [k,l];
        a_ji = a_ij ;  a_ii = -sum_{j != i} a_ij .

    Equivalently [L_tau]_ij = w_ij with w_ij as in their eq. (3-3) (a negative
    number), and [L_tau]_ii = -sum_{j != i} w_ij  -- i.e. off-diagonal entries
    NEGATIVE and diagonal POSITIVE, so L_tau is a singular M-matrix and is
    POSITIVE SEMIDEFINITE.  This is the convention in which their

        Theorem 3.1:  E_tau(f) = sum_tau 3|f(tau)|^2 / (2 mu(tau))
                      = -0.5 * trace(f^T L_V(f) f)
        Theorem 3.2:  grad_f E_V(f) = -3 L_V(f) f

    both hold simultaneously.  Both are verified numerically in
    tools/verify_ball_ref.py (checks V2 and V3).

    BEWARE: published statements of the trace identity come in two sign
    conventions differing by an overall sign of L_V.  The two conventions are
    NOT interchangeable in an implementation, because the identity
        0.5 * sum_s f_s^T L_tau f_s  =  3|f(tau)|^2/(2 mu)
    and the gradient identity grad E_V = 3 L_V f cannot both hold for the same
    L_V.  Use `vsem_trace_laplacian()` if you want the sign for which
    grad E_V = +3 L_V f.
    """
    a = np.zeros((4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            k, l = [x for x in range(4) if x != i and x != j]
            v = np.cross(f[k] - f[i], f[l] - f[i])
            w = np.cross(f[l] - f[j], f[k] - f[j])
            a[i, j] = a[j, i] = float(v @ w)
    for i in range(4):
        a[i, i] = -a[i].sum()
    return -(1.0 / (36.0 * mu_tau)) * a


def vsem_trace_laplacian(V: np.ndarray, T: np.ndarray, mu: np.ndarray | None = None) -> sp.csr_matrix:
    """-L_V(f): the sign convention in which 0.5*sum_s f_s^T (-L_V) f_s is
    NEGATIVE of the energy, i.e. the one used by the n-VSE paper (arXiv:2402.00380)
    and by IEM (arXiv:2407.19272), whose authors write
        E_V(f) = (1/3) trace(f^T L_V(f) f) = sum_tau |f(tau)|^2/|tau|
    for a L_V(f) with positive off-diagonals.  Compute the gradient as
        grad E_V = 3 * vsem_trace_laplacian(...) @ f
    if you adopt that bookkeeping.  (`vsem_energy` below is convention-free.)
    """
    return -vsem_laplacian(V, T, mu)


def vsem_laplacian(V: np.ndarray, T: np.ndarray, mu: np.ndarray | None = None) -> sp.csr_matrix:
    """Assemble L_V(f) for the map with vertex images V on the mesh T.

    `mu` are the *reference* tet volumes |tau| (default |tau| of V itself).

    Fully vectorised: for every tet we build the 4x4 local matrix in one batched
    einsum, then scatter with np.add.at.
    L_V is symmetric and its diagonal is negative (M-matrix, as in Yueh et al.).
    """
    m, n = len(T), len(V)
    if mu is None:
        mu = tet_volumes(V, T)
    f = V[T]  # (m,4,3)

    # a_ij(local i<j) via the fundamental identity, vectorised over tets
    # for {i,j,k,l}:  v = (f_k - f_i) x (f_l - f_i), w = (f_l - f_j) x (f_k - f_j)
    a = np.zeros((m, 4, 4))
    for i in range(4):
        for j in range(i + 1, 4):
            k, l = [x for x in range(4) if x != i and x != j]
            v = np.cross(f[:, k] - f[:, i], f[:, l] - f[:, i])
            w = np.cross(f[:, l] - f[:, j], f[:, k] - f[:, j])
            a[:, i, j] = a[:, j, i] = np.einsum("mi,mi->m", v, w)
    a[:, range(4), range(4)] = -a.sum(axis=2)
    Lt = -(1.0 / (36.0 * mu))[:, None, None] * a  # (m,4,4)  == _local_stretch_matrix

    rows = np.repeat(T, 4, axis=1).ravel()   # local (i,j): rows = T[i] for each j
    cols = np.repeat(T, 4, axis=0).ravel()   # local (i,j): cols = T[j] for each i
    L = sp.coo_matrix((Lt.ravel(), (rows, cols)), shape=(n, n)).tocsr()
    L.sum_duplicates()
    return L


def vsem_energy(V_ref: np.ndarray, T: np.ndarray, V_map: np.ndarray, mu: np.ndarray | None = None) -> float:
    """E_V(f) = sum_tau 3|f(tau)|^2 / (2 mu(tau)).

    Lower bound (Yueh2019 Thm 3.3, given |f(M)| = mu(M)):
        E_V(f) >= (3/2) * mu(M)   ... in this normalisation
    Free lower bound used in the n-VSE / IEM papers with E_V = sum |f|^2/|tau|:
        E_V(f) >= |f(M)|^2 / mu(M).
    """
    if mu is None:
        mu = tet_volumes(V_ref, T)
    return float(np.sum(3.0 * tet_volumes(V_map, T) ** 2 / (2.0 * mu)))


def vsem_lower_bound(V_ref: np.ndarray, T: np.ndarray, V_map: np.ndarray, mu: np.ndarray | None = None) -> float:
    """(3/2)*|f(M)| achieved iff volume preserving (with matching total volume)."""
    if mu is None:
        mu = tet_volumes(V_ref, T)
    return 1.5 * float(np.sum(tet_volumes(V_map, T)))


# ----------------------------------------------------------------------------
# 4. Goal A: P1 harmonic (Dirichlet-energy) map of the solid onto the unit ball
# ----------------------------------------------------------------------------


def p1_tet_stiffness_local(V: np.ndarray, T: np.ndarray) -> np.ndarray:
    """The 4x4 element stiffness matrices, shape (m,4,4).

    K^tau_{ij} = |tau| * grad(lambda_i) . grad(lambda_j),  with
    grad(lambda_k) = k-th row of B^{-1}, k=1,2,3,  B = [v1-v0, v2-v0, v3-v0],
    and grad(lambda_0) = -sum_{k>=1} grad(lambda_k).
    Equivalent closed form:  K^tau_ij = -(1/6) l_kl cot(theta_ij^{kl}) for i!=j
    and K^tau_ii = -sum_{j!=i} K^tau_ij, where (k,l) is the edge opposite (i,j).
    """
    vol = tet_volumes(V, T)
    p = V[T]
    B = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=2)
    Binv = np.linalg.inv(B)  # grad(lambda_k) = row k-1, k=1,2,3
    G = np.empty((len(T), 4, 3))
    G[:, 1:4, :] = Binv
    G[:, 0, :] = -G[:, 1:4, :].sum(axis=1)
    return vol[:, None, None] * np.einsum("mik,mjk->mij", G, G)


def p1_tet_stiffness(V: np.ndarray, T: np.ndarray) -> sp.csr_matrix:
    """P1 (linear) FEM stiffness / Dirichlet matrix on a tetrahedral mesh.

    For a linear tet tau with vertices v0..v3 let
        B = [v1-v0, v2-v0, v3-v0]  in R^{3x3}   (COLUMNS are the edge vectors),
        B[i,j] = (v_{j+1} - v_0)_i .
    Then B^{-1} A = I for A[i,k] = d lambda_k / d x_i, so
        grad(lambda_k) = k-th ROW of B^{-1}   for k = 1,2,3,
    and  grad(lambda_0) = -sum_{k>=1} grad(lambda_k) .
    (Sanity: grad(lambda_k) . (v_j - v_0) = delta_{kj}, checked in V5f.)
    The element matrix is
        K^tau_{ij} = |tau| * grad(lambda_i) . grad(lambda_j)
    and the global stiffness is K = sum_tau K^tau (scattered by vertex id).
    This is the tetrahedral analogue of the triangle cotangent matrix.
    Verified in tools/verify_ball_ref.py (V5a-V5h).
    """
    m = len(T)
    Kt = p1_tet_stiffness_local(V, T)  # (m,4,4)

    rows = np.repeat(T, 4, axis=1).ravel()   # local (i,j): rows = T[i] for each j
    cols = np.repeat(T, 4, axis=0).ravel()   # local (i,j): cols = T[j] for each i
    K = sp.coo_matrix((Kt.ravel(), (rows, cols)), shape=(len(V), len(V))).tocsr()
    K.sum_duplicates()
    return K


def boundary_faces_oriented(V: np.ndarray, T: np.ndarray):
    """Boundary triangles of a tetrahedral mesh, OUTWARD wound, robustly.

    A boundary triangle is a face code that occurs exactly once over all tets
    (each interior face occurs twice, once per incident tet).  Rather than rely
    on a fixed "canonical" face permutation -- which yields outward normals only
    for one particular ordering convention, and TetGen does not respect any such
    convention -- the winding is fixed GEOMETRICALLY: for the tet
    [i0,i1,i2,i3] and a face code opposite vertex i_k, the normal is flipped
    whenever  n . (v_k - v_a) > 0, since n must point away from the opposite
    vertex (i.e. out of the tet).

    Returns (faces (nf,3) int, boundary_vertex_ids (nb,) int).
    """
    codes, opp = [], []
    for k in range(4):
        f = np.delete(np.arange(4), k)
        codes.append(T[:, f])   # (m,3)
        opp.append(T[:, k])     # (m,)
    C = np.concatenate(codes, axis=0)
    K = np.concatenate(opp, axis=0)

    Cs = np.sort(C, axis=1)
    uniq, inv, counts = np.unique(Cs, axis=0, return_inverse=True, return_counts=True)
    # index of the first occurrence of each unique (sorted) face
    first = np.full(len(uniq), -1, dtype=np.int64)
    order = np.arange(len(C))
    for pos in order[::-1]:
        first[inv[pos]] = pos
    sel = (counts[inv] == 1) & (order == first[inv])
    Cf, Kf = C[sel], K[sel]

    a, b_, c = V[Cf[:, 0]], V[Cf[:, 1]], V[Cf[:, 2]]
    nrm = np.cross(b_ - a, c - a)
    inward = np.einsum("ij,ij->i", nrm, V[Kf] - a) > 0
    Cf = Cf.copy()
    Cf[inward] = Cf[inward][:, [0, 2, 1]]
    return np.unique(Cf.ravel()), Cf


def boundary_vertices(T: np.ndarray, n: int, V: np.ndarray | None = None):
    """Boundary vertex ids and boundary triangles.

    If `V` is supplied the faces are returned OUTWARD WOUND (see
    `boundary_faces_oriented`); otherwise only the topological deduplication by
    sorted triple is done and the winding is arbitrary.  Always pass `V` if you
    will compute spherical areas, divergence-theorem volumes or normals.
    """
    if V is not None:
        return boundary_faces_oriented(V, T)
    f = np.concatenate([T[:, [0, 2, 1]], T[:, [0, 1, 3]], T[:, [0, 3, 2]], T[:, [1, 2, 3]]])
    fs = np.sort(f, axis=1)
    uniq, counts = np.unique(fs, axis=0, return_counts=True)
    bnd_faces = uniq[counts == 1]
    return np.unique(bnd_faces.ravel()), bnd_faces


def harmonic_ball_map(V: np.ndarray, T: np.ndarray, sphere: np.ndarray | None = None):
    """Goal A: solve  K_II u_I = -K_IB u_B  three times, with the boundary
    vertices fixed on the unit sphere.

    `sphere` gives the target (nB,3) boundary positions.  If None, a crude
    boundary map is used: project the boundary vertices radially onto the
    sphere (this is NOT bijective for non-convex domains - see notes).
    For a real pipeline supply a spherical conformal/authalic map here.
    """
    n = len(V)
    K = p1_tet_stiffness(V, T)
    b, bf = boundary_vertices(T, n, V)
    I = np.setdiff1d(np.arange(n), b)

    if sphere is None:
        S = V[b] - V[b].mean(axis=0)
        S = S / np.linalg.norm(S, axis=1, keepdims=True)
    else:
        S = np.asarray(sphere, dtype=np.float64)

    U = np.zeros((n, 3))
    U[b] = S
    KII = K[I][:, I].tocsc()
    KIB = K[I][:, b]
    solve = spla.factorized(KII)
    for s in range(3):
        U[I, s] = solve(-KIB @ S[:, s])
    return U, b, I


# ----------------------------------------------------------------------------
# 5. Goal B: VSEM fixed point  (Yueh et al. 2019, Algorithm 4.4 / VSEM Alg. 1)
# ----------------------------------------------------------------------------


def vsem_solve(
    V: np.ndarray,
    T: np.ndarray,
    fB: np.ndarray | None = None,
    n_iter: int = 100,
    tol: float = 1e-6,
    verbose: bool = False,
):
    """Solve  [L_V(f)]_II f_I^s = -[L_V(f)]_IB f_B^s  for s=1,2,3, fixed point.

    fB defaults to the radial projection of the boundary onto the sphere.
    Returns (f, history).
    """
    n = len(V)
    b, _ = boundary_vertices(T, n, V)
    I = np.setdiff1d(np.arange(n), b)
    if fB is None:
        S = V[b] - V[b].mean(axis=0)
        fB = S / np.linalg.norm(S, axis=1, keepdims=True)
    fB = np.asarray(fB, dtype=np.float64)

    mu = tet_volumes(V, T)
    f = np.zeros((n, 3))
    f[b] = fB
    # warm start with the harmonic (Dirichlet) map
    f, _, _ = harmonic_ball_map(V, T, fB)

    hist = []
    L = vsem_laplacian(f, T, mu)
    for it in range(n_iter):
        e_old = vsem_energy(V, T, f, mu)
        LII = L[I][:, I].tocsc()
        LIB = L[I][:, b]
        solve = spla.factorized(LII)
        for s in range(3):
            f[I, s] = solve(-LIB @ fB[:, s])
        L = vsem_laplacian(f, T, mu)
        e_new = vsem_energy(V, T, f, mu)
        hist.append(e_new)
        if verbose:
            print(f"  it {it:4d}  E_V = {e_new:.10e}   lower = {vsem_lower_bound(V,T,f,mu):.10e}")
        if e_old - e_new < tol * max(1.0, abs(e_old)):
            break
    return f, np.array(hist)


# ----------------------------------------------------------------------------
# 6. Verification metrics
# ----------------------------------------------------------------------------


def volume_distortion_ratios(V_ref, T, V_map):
    """|f(tau)| / |tau| -> array (m,). Ideal = constant
    c = |f(M)|/|M| (for a ball map c ~ (4/3 pi)/|M|)."""
    return tet_volumes(V_map, T) / np.maximum(tet_volumes(V_ref, T), 1e-300)


def jacobian_det_field(V_ref, T, V_map):
    """J(tau) = det J_f|_tau = |f(tau)|/|tau|, signed.

    Implementation: J is constant per tet.  Build the affine map
    A(tau) = [f1-f0, f2-f0, f3-f0] [v1-v0, v2-v0, v3-v0]^{-1},
    J = det A.  Equivalently det(A) = det(F)/det(V) = (6|f(tau)|)/(6|tau|).
    Negative J <=> orientation flipped (fold-over).
    """
    return tet_volume_signed(V_map, T) / np.where(
        np.abs(tet_volume_signed(V_ref, T)) < 1e-300, 1e-300, tet_volume_signed(V_ref, T)
    )


def flipped_fraction(V_ref, T, V_map):
    J = jacobian_det_field(V_ref, T, V_map)
    return float(np.mean(J <= 0)), J


def volume_preserving_error(V_ref, T, V_map, mu=None):
    """eps = E_V(f) - |f(M)|^2 / mu(M)   (n-VSE papers).
    Equals (|f(M)|^2/mu(M)^2) * ||delta||^2_{L2} with delta = volume ratio - 1."""
    if mu is None:
        mu = tet_volumes(V_ref, T)
    EV = float(np.sum(tet_volumes(V_map, T) ** 2 / mu))
    Vf = float(np.sum(tet_volumes(V_map, T)))
    VM = float(mu.sum())
    return EV - Vf ** 2 / VM


def verify_map(V_ref, T, V_map, label: str = "map", verbose: bool = True) -> dict:
    """Everything the team should log for a candidate volume-preserving map."""
    mu = tet_volumes(V_ref, T)
    r = volume_distortion_ratios(V_ref, T, V_map)
    J = jacobian_det_field(V_ref, T, V_map)
    da = dihedral_angles(V_map, T)
    EV = vsem_energy(V_ref, T, V_map, mu)
    Vf = float(tet_volumes(V_map, T).sum())
    VM = float(mu.sum())
    out = dict(
        label=label,
        n_tets=len(T),
        total_vol_ref=VM,
        total_vol_map=Vf,
        ratio_mean=float(r.mean()),
        ratio_std=float(r.std()),
        ratio_min=float(r.min()),
        ratio_max=float(r.max()),
        ratio_p99_over_p1=float(np.percentile(r, 99) / max(np.percentile(r, 1), 1e-300)),
        J_min=float(J.min()),
        J_max=float(J.max()),
        frac_flipped=float(np.mean(J <= 0)),
        frac_J_lt_0p5=float(np.mean(J < 0.5)),
        EV=EV,
        EV_lower_bound_1p5_Vf=1.5 * Vf,
        EV_lower_bound_Vf2_over_VM=Vf ** 2 / VM,
        eps_mass_preserving=volume_preserving_error(V_ref, T, V_map, mu),
        min_dihedral_deg=float(da.min()),
        max_dihedral_deg=float(da.max()),
        mean_min_dihedral_deg=float(da.min(axis=1).mean()),
        max_radius_edge_ratio=float(radius_edge_ratio(V_map, T).max()),
    )
    if verbose:
        print(f"--- verify_map({label}) ---")
        for k, v in out.items():
            print(f"  {k:28s} = {v}")
    return out


# ----------------------------------------------------------------------------
# 7. Bijectivity check (practical): signed-volume consistency + self-intersection
# ----------------------------------------------------------------------------


def is_locally_injective(V_ref, T, V_map):
    """Necessary condition: every tet keeps positive orientation.  NOT sufficient."""
    return bool(np.all(tet_volume_signed(V_map, T) > 0))


def bijectivity_report(V_ref, T, V_map) -> dict:
    J = jacobian_det_field(V_ref, T, V_map)
    return dict(
        all_positive=bool(np.all(J > 0)),
        n_flipped=int(np.sum(J <= 0)),
        worst_J=float(J.min()),
        # practical global checks (see markdown): (i) boundary sphere map must be
        # bijective (check via spherical triangle area sum == 4*pi and no
        # overlapping spherical triangles), (ii) sample interior rays from the
        # origin and count signed crossings.
        boundary_area_error=None,
    )


# ----------------------------------------------------------------------------
# demo
# ----------------------------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--sub", type=int, default=3)
    ap.add_argument("--maxvol", type=float, default=0.0)
    ap.add_argument("--iters", type=int, default=60)
    args = ap.parse_args()

    if not args.demo:
        ap.print_help()
        return 0

    print("== building unit-ball tet mesh ==")
    V, T = make_ball_tet_mesh(args.sub, args.maxvol)
    mu = tet_volumes(V, T)
    print(f"  n_verts={len(V)}  n_tets={len(T)}  vol={mu.sum():.6f}  (4pi/3={4*math.pi/3:.6f})")
    print(f"  signed vol min = {tet_volume_signed(V,T).min():.3e} (must be >0)")

    print("\n== tet quality ==")
    da = dihedral_angles(V, T)
    rer = radius_edge_ratio(V, T)
    ar = aspect_ratio(V, T)
    print(f"  min dihedral          : {da.min():.3f} deg")
    print(f"  mean of per-tet min   : {da.min(axis=1).mean():.3f} deg")
    print(f"  radius-edge ratio     : mean {rer.mean():.4f}  max {rer.max():.4f}")
    print(f"  aspect ratio          : mean {ar.mean():.4f}  max {ar.max():.4f}")

    print("\n== Goal A: P1 harmonic map to unit ball ==")
    U, b, I = harmonic_ball_map(V, T)
    print(f"  n_boundary={len(b)} n_interior={len(I)}")
    verify_map(V, T, U, "harmonic (radial BC)")

    print("\n== Goal B: VSEM fixed point ==")
    f, hist = vsem_solve(V, T, n_iter=args.iters, verbose=False)
    print(f"  iterations run: {len(hist)}")
    verify_map(V, T, f, "VSEM")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
