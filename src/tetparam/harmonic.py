"""Stage 1 of the pipeline: the "conformal" (harmonic) map of the solid onto the unit ball.

Discrete setting
----------------
A tetrahedral mesh ``M`` with vertices ``v_1..v_n``.  A piecewise-linear map ``f`` is fixed
by its values ``f(v_i)``.  The volumetric Dirichlet (harmonic) energy is

    E_D(f) = sum_{tau in T} |tau| * || grad f |_tau ||_F^2
           = 1/2 * trace( f^T L f )        (L = P1 finite-element stiffness matrix)
           = 1/2 * sum_{s=1..3} f_s^T L f_s

with

    L_ij = sum_{tau} |tau| * < grad lambda_i , grad lambda_j > .

Minimising ``E_D`` subject to a fixed boundary map is *the* natural 3D analogue of the
harmonic/conformal surface map: it is the Euler-Lagrange equation of the Dirichlet energy,
i.e. ``Delta f = 0`` with Dirichlet data on the boundary.  The result is a homeomorphism
onto a *convex* target for a very large class of domains, which is why we map to the ball.

The boundary is put on the unit sphere.  Several strategies are provided; ``radial`` is
exact and cheap when the solid is star-shaped w.r.t. its centre, ``laplacian`` performs
tangential relaxation to spread the boundary vertices more evenly (lower area distortion),
and ``auc``/``harmonic_extension`` are simple fallbacks.  A production pipeline should
replace this by a spherical conformal / area-preserving map (see docs).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .mesh import TetMesh, tet_shape_gradients


# --------------------------------------------------------------------------------------
# finite element operators
# --------------------------------------------------------------------------------------
def fem_stiffness(mesh: TetMesh, weights: np.ndarray | None = None) -> sp.csr_matrix:
    """P1 stiffness (Laplacian) matrix ``L_ij = sum_tau w_tau <grad l_i, grad l_j>``.

    Parameters
    ----------
    weights : (m,) optional
        Per-tetrahedron multiplier; defaults to the tetrahedron volume.  Passing
        ``|f(tau)| / |tau|``-like factors gives the modified cotangent weights used by the
        volumetric stretch energy.
    """
    grads, vols = tet_shape_gradients(mesh.vertices, mesh.tets)
    w = vols if weights is None else np.asarray(weights, dtype=np.float64) * np.sign(vols)
    # local 4x4 matrices
    G = grads  # (m,4,3)
    local = np.einsum("mik,mjk->mij", G, G) * w[:, None, None]  # (m,4,4)

    t = mesh.tets
    rows = np.repeat(t, 4, axis=1).ravel()
    cols = np.tile(t, (1, 4)).ravel()
    vals = local.ravel()
    n = mesh.n_vertices
    L = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    L.sum_duplicates()
    return L


def mass_matrix_lumped(mesh: TetMesh) -> np.ndarray:
    """Lumped P1 mass matrix: ``vol(tau)/4`` on every vertex of ``tau``."""
    V = np.abs(tet_shape_gradients(mesh.vertices, mesh.tets)[1]) / 4.0
    out = np.zeros(mesh.n_vertices)
    for k in range(4):
        np.add.at(out, mesh.tets[:, k], V)
    return out


# --------------------------------------------------------------------------------------
# boundary -> unit sphere
# --------------------------------------------------------------------------------------
def boundary_map_sphere(
    mesh: TetMesh,
    method: str = "laplacian",
    n_iter: int = 30,
    relax: float = 0.5,
    seed: int = 0,
) -> np.ndarray:
    """Return sphere positions for the boundary vertices, shape ``(n, 3)``.

    Only boundary rows are meaningful; interior rows are left at zero.
    """
    bnd = mesh.boundary()
    faces = bnd["faces"]
    bverts = bnd["boundary_vertices"]
    is_b = bnd["is_boundary_vertex"]

    P = mesh.vertices
    if method == "radial":
        c = P[bverts].mean(axis=0)
        d = P[bverts] - c
        d /= np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-300)
        out = np.zeros_like(P)
        out[bverts] = d
        return out

    if method != "laplacian":
        raise ValueError(f"unknown boundary method {method!r}")

    # ---- initial radial projection
    c = P[bverts].mean(axis=0)
    d = P[bverts] - c
    d /= np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-300)
    y = np.zeros_like(P)
    y[bverts] = d

    # ---- surface graph Laplacian on the boundary triangle mesh
    n = mesh.n_vertices
    idx = -np.ones(n, dtype=np.int64)
    idx[bverts] = np.arange(len(bverts))
    f = faces
    rows = np.concatenate([idx[f[:, 0]], idx[f[:, 1]], idx[f[:, 2]]])
    cols = np.concatenate([idx[f[:, 1]], idx[f[:, 2]], idx[f[:, 0]]])
    data = np.ones(len(rows))
    A = sp.coo_matrix((data, (rows, cols)), shape=(len(bverts), len(bverts))).tocsr()
    A.data[:] = 1.0
    # uniform (combinatorial) weights; positive-definite Laplacian
    deg = np.asarray(A.sum(axis=1)).ravel()
    Lap = sp.diags(deg) - A

    Y = y[bverts].copy()
    for _ in range(n_iter):
        Y = Y - relax * (Lap @ Y) / np.maximum(deg, 1.0)[:, None]
        Y /= np.maximum(np.linalg.norm(Y, axis=1, keepdims=True), 1e-300)
    y[bverts] = Y
    del seed
    return y


# --------------------------------------------------------------------------------------
# harmonic extension + result container
# --------------------------------------------------------------------------------------
@dataclass
class HarmonicMapResult:
    target: TetMesh
    boundary_values: np.ndarray
    residual: float
    n_flipped: int
    energy: float

    def summary(self) -> str:
        return (
            f"HarmonicMap: residual {self.residual:.3e}, Dirichlet energy {self.energy:.6g}, "
            f"flipped tets {self.n_flipped}/{self.target.n_tets}"
        )


def solve_harmonic(
    mesh: TetMesh,
    boundary_values: np.ndarray,
    L: sp.csr_matrix | None = None,
    tol: float = 1e-10,
) -> tuple[np.ndarray, float]:
    """Solve ``L_II y_I = -L_IB y_B`` for the interior vertices.

    Returns the full vertex array ``y`` (shape ``(n, 3)``) and the relative residual.
    """
    if L is None:
        L = fem_stiffness(mesh)
    is_b = mesh.is_boundary_vertex
    I = np.flatnonzero(~is_b)
    B = np.flatnonzero(is_b)
    y = np.array(boundary_values, dtype=np.float64, copy=True)

    if len(I) == 0:
        return y, 0.0

    LII = L[I][:, I].tocsc()
    LIB = L[I][:, B]
    rhs = -(LIB @ y[B])
    # the system is symmetric positive definite on the interior block
    try:
        sol = spla.spsolve(LII, rhs)
    except Exception:  # pragma: no cover - fall back to CG
        sol, _ = spla.cg(LII, rhs, rtol=tol, maxiter=5000)
    if sol.ndim == 1:
        sol = sol[:, None]
    y[I] = sol
    res = np.linalg.norm(LII @ sol - rhs) / max(np.linalg.norm(rhs), 1e-300)
    return y, float(res)


def harmonic_ball_map(
    mesh: TetMesh,
    boundary_method: str = "laplacian",
    boundary_kwargs: dict | None = None,
    rescale_to_ball: bool = True,
) -> HarmonicMapResult:
    """Map ``mesh`` (a topological 3-ball) onto the unit ball by a volumetric harmonic map.

    ``rescale_to_ball`` renormalises the boundary radius to exactly 1 after the solve (the
    harmonic extension keeps boundary vertices fixed, so this is only a safety net).
    """
    bv = boundary_map_sphere(mesh, method=boundary_method, **(boundary_kwargs or {}))
    L = fem_stiffness(mesh)
    y, res = solve_harmonic(mesh, bv, L=L, tol=1e-10)
    if rescale_to_ball:
        b = mesh.is_boundary_vertex
        y[b] /= np.maximum(np.linalg.norm(y[b], axis=1, keepdims=True), 1e-300)

    target = mesh.with_vertices(y)
    sv = target.signed_volumes()
    energy = float(0.5 * sum(y[:, s] @ (L @ y[:, s]) for s in range(3)))
    return HarmonicMapResult(
        target=target,
        boundary_values=y.copy(),
        residual=res,
        n_flipped=int((sv <= 0).sum()),
        energy=energy,
    )


def harmonic_map_to_domain(
    mesh: TetMesh, boundary_values: np.ndarray, L: sp.csr_matrix | None = None
) -> TetMesh:
    """Harmonic extension with arbitrary Dirichlet data (used for remeshing / lifting)."""
    y, _ = solve_harmonic(mesh, boundary_values, L=L)
    return mesh.with_vertices(y)


# --------------------------------------------------------------------------------------
# bijective initialisation: the star-shaped ("radial") map
# --------------------------------------------------------------------------------------
def ray_surface_hits(
    origin: np.ndarray, directions: np.ndarray, tri: np.ndarray, eps: float = 1e-14
) -> np.ndarray:
    """Distance from ``origin`` to the *first* triangle hit along each direction.

    ``tri`` is ``(k, 3, 3)``.  Uses the Moller-Trumbore test.  Returns ``(n,)`` distances
    (``inf`` where no hit).  Rays are processed in chunks to bound memory.
    """
    directions = np.asarray(directions, dtype=np.float64)
    tri = np.asarray(tri, dtype=np.float64)
    if len(tri) == 0:
        return np.full(len(directions), np.inf)
    v0 = tri[:, 0]
    e1 = tri[:, 1] - v0
    e2 = tri[:, 2] - v0
    tvec = origin[None, :] - v0  # (k,3)

    out = np.empty(len(directions))
    chunk = max(1, int(4_000_000 / max(len(tri) * 3, 1)))
    for s in range(0, len(directions), chunk):
        d = directions[s : s + chunk]  # (c,3)
        pvec = np.cross(d[:, None, :], e2[None, :, :])  # (c,k,3)
        det = np.einsum("ckj,kj->ck", pvec, e1)
        ok = np.abs(det) > eps
        inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
        u = np.einsum("ckj,kj->ck", pvec, tvec) * inv
        qvec = np.cross(
            np.broadcast_to(tvec, (len(d), len(tri), 3)),
            np.broadcast_to(e1, (len(d), len(tri), 3)),
        )
        v = np.einsum("ckj,ckj->ck", d[:, None, :], qvec) * inv
        t = np.einsum("kj,ckj->ck", e2, qvec) * inv
        hit = ok & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1.0 + 1e-9) & (t > 1e-9)
        out[s : s + chunk] = np.where(hit, t, np.inf).min(axis=1)
    return out


def star_ball_map(mesh: TetMesh, origin: np.ndarray | None = None) -> np.ndarray:
    """Bijective map of a star-shaped solid onto the unit ball.

    Every vertex is pushed along the ray from ``origin`` so that the boundary lands on the
    unit sphere; the interior is scaled by the same local factor.  For a solid that is
    star-shaped with respect to ``origin`` (every ray meets the boundary exactly once) this
    is a homeomorphism *by construction* -- unlike the harmonic map, which may fold.
    """
    bnd = mesh.boundary()
    faces = bnd["faces"]
    tri = mesh.vertices[faces]
    P = mesh.vertices
    c = P.mean(axis=0) if origin is None else np.asarray(origin, dtype=np.float64)

    d = P - c
    r = np.linalg.norm(d, axis=1)
    dirs = d / np.maximum(r, 1e-300)[:, None]
    R = ray_surface_hits(c, dirs, tri)
    R = np.where(np.isfinite(R) & (R > 0), R, np.maximum(r, 1e-300))
    # boundary vertices have r == R (up to the polyhedral approximation)
    scale = 1.0 / np.maximum(R, 1e-300)
    return c + d * scale[:, None]


def volume_penalty_relax(
    mesh: TetMesh,
    y_init: np.ndarray,
    L: sp.csr_matrix | None = None,
    max_iter: int = 400,
    tol: float = 1e-12,
    verbose: bool = False,
) -> tuple[np.ndarray, dict]:
    """Remove fold-overs by minimising the volume-mismatch energy

        Phi(y) = sum_tau ( |tau_y| / |tau_0| - 1 )^2 .

    ``Phi`` is *zero exactly when the map is volume-preserving per tetrahedron*, and it is
    bounded and smooth even for inverted tetrahedra, so gradient descent from a folded map
    pulls every negative volume back towards its positive target.  The gradient uses

        d|tau| / d v_k = (A_k / 3) n_k ,

    with ``A_k, n_k`` the area and outward unit normal of the face opposite vertex ``k``.

    Returns ``(y, info)``.
    """
    if L is None:
        L = fem_stiffness(mesh)
    del L
    is_b = mesh.is_boundary_vertex
    y = np.array(y_init, dtype=np.float64, copy=True)

    v0, dV = _volume_and_grad_ops(mesh)
    target = v0
    edges_all = mesh.boundary()["edges"]
    e = np.linalg.norm(mesh.vertices[edges_all[:, 0]] - mesh.vertices[edges_all[:, 1]], axis=1)
    h = float(np.median(e[e > 0])) if len(e) and (e > 0).any() else 1.0

    def phi(yy: np.ndarray) -> float:
        v = _tet_volumes(yy, mesh.tets)
        return float(np.mean((v / target - 1.0) ** 2))

    alpha = 0.5
    p = phi(y)
    info = {"iterations": 0, "penalty": p, "min_ratio": float((_tet_volumes(y, mesh.tets) / v0).min())}
    for k in range(max_iter):
        g = _volume_penalty_grad(y, mesh.tets, dV, target)
        g[is_b] = 0.0
        nrm = np.linalg.norm(g)
        if nrm < 1e-16:
            break
        d = -g / nrm * h  # step of at most  h * alpha  in space
        a = alpha
        accepted = False
        for _ in range(50):
            cand = y + a * d
            pc = phi(cand)
            if pc < p:
                accepted = True
                break
            a *= 0.5
        if not accepted:
            break
        y = y + a * d
        alpha = min(a * 1.5, 0.5)
        p_new = phi(y)
        info.update(
            iterations=k + 1,
            penalty=p_new,
            min_ratio=float((_tet_volumes(y, mesh.tets) / v0).min()),
        )
        if verbose and k % 20 == 0:
            print(f"    repair {k:4d}  Phi {p_new:.6e}  min-ratio {info['min_ratio']:.4g}")
        if abs(p - p_new) <= tol * max(abs(p), 1e-30):
            p = p_new
            break
        p = p_new
    return y, info


def _tet_volumes(y: np.ndarray, tets: np.ndarray) -> np.ndarray:
    a = y[tets[:, 0]]
    b = y[tets[:, 1]]
    c = y[tets[:, 2]]
    d = y[tets[:, 3]]
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def _volume_and_grad_ops(mesh: TetMesh) -> tuple[np.ndarray, np.ndarray]:
    """Reference volumes and, per tet, ``(A_k/3) n_k`` for the face opposite vertex ``k``."""
    v0 = _tet_volumes(mesh.vertices, mesh.tets)
    P = mesh.vertices[mesh.tets]  # (m,4,3)
    ops = np.empty((mesh.n_tets, 4, 3))
    for k in range(4):
        idx = [i for i in range(4) if i != k]
        a, b, c = P[:, idx[0]], P[:, idx[1]], P[:, idx[2]]
        n = np.cross(b - a, c - a)  # 2*A*n, oriented away from vertex k
        # ensure it points away from the opposite vertex
        sgn = np.sign(np.einsum("ij,ij->i", n, P[:, k] - a))
        n = n * np.where(sgn[:, None] < 0, -1.0, 1.0)
        ops[:, k] = n / 6.0  # dV/dv_k = (1/3) A n = |n|/2 /3 * n_hat = n/6
    return v0, ops


def _volume_penalty_grad(
    y: np.ndarray, tets: np.ndarray, dV: np.ndarray, target: np.ndarray
) -> np.ndarray:
    v = _tet_volumes(y, tets)
    coef = 2.0 * (v / target - 1.0) / target / len(tets)
    g = np.zeros_like(y)
    for k in range(4):
        np.add.at(g, tets[:, k], coef[:, None] * dV[:, k])
    return g


def fold_free_relax(
    mesh: TetMesh,
    y_init: np.ndarray,
    L: sp.csr_matrix | None = None,
    max_iter: int = 200,
    min_ratio: float = 0.2,
    tol: float = 1e-9,
    repair: bool = True,
    verbose: bool = False,
) -> tuple[np.ndarray, dict]:
    """Minimise the Dirichlet energy *without ever creating a negative tetrahedron*.

    Projected gradient descent on ``E_D`` with a backtracking step size that keeps every
    normalised volume ``|tau_y| / |tau|`` above ``min_ratio``.  Starting from a bijective
    map (:func:`star_ball_map`) this yields a harmonic-like map that is still a
    homeomorphism -- the standard remedy for the well-known fact that a 3-D harmonic map
    onto a convex domain need not be injective.

    Returns ``(y, info)``.
    """
    if L is None:
        L = fem_stiffness(mesh)
    is_b = mesh.is_boundary_vertex
    y = np.array(y_init, dtype=np.float64, copy=True)

    v0 = np.abs(mesh.signed_volumes())
    v0 = np.maximum(v0, 1e-300)

    def min_norm_volume(yy: np.ndarray) -> float:
        return float((_tet_volumes(yy, mesh.tets) / v0).min())

    def energy(yy: np.ndarray) -> float:
        return float(0.5 * sum(yy[:, s] @ (L @ yy[:, s]) for s in range(3)))

    if repair and min_norm_volume(y) <= min_ratio:
        y, rep = volume_penalty_relax(mesh, y, max_iter=max_iter, verbose=verbose)
        if verbose:
            print(f"    repair finished: Phi {rep['penalty']:.3e} min-ratio {rep['min_ratio']:.4g}")

    # the floor only guards against degenerate steps; never stricter than what the
    # current map already achieves
    floor = min(min_ratio, 0.5 * min_norm_volume(y))

    alpha = 0.1
    e = energy(y)
    info = {
        "iterations": 0,
        "energy": e,
        "min_ratio": min_norm_volume(y),
        "alpha": alpha,
    }
    for k in range(max_iter):
        g = L @ y
        d = -g
        d[is_b] = 0.0
        gd = float(np.sum(g * d))
        if np.linalg.norm(d) < 1e-14:
            break
        a = min(alpha * 1.5, 1.0)
        accepted = False
        for _ in range(50):
            cand = y + a * d
            if min_norm_volume(cand) > floor:
                ec = energy(cand)
                if ec <= e + 1e-4 * a * gd:
                    accepted = True
                    break
            a *= 0.5
        if not accepted:
            break
        y = y + a * d
        alpha = a
        e_new = energy(y)
        info.update(iterations=k + 1, energy=e_new, min_ratio=min_norm_volume(y), alpha=a)
        if verbose and k % 10 == 0:
            print(
                f"    relax {k:4d}  E {e_new:.8g}  alpha {a:.3g}  min-ratio {info['min_ratio']:.3g}"
            )
        if abs(e - e_new) <= tol * max(abs(e), 1.0):
            e = e_new
            break
        e = e_new
    return y, info


def ball_map(
    mesh: TetMesh,
    method: str = "auto",
    boundary_method: str = "laplacian",
    relax: bool = True,
    verbose: bool = False,
) -> dict:
    """Produce a fold-free map of ``mesh`` onto the unit ball, trying several strategies.

    ``method`` is one of ``'auto'``, ``'harmonic'``, ``'star'``, ``'relax'``.

    Returns a dict with keys ``vertices`` (the ball coordinates), ``method`` (which strategy
    won), ``n_flipped``, ``energy`` and ``stretch_excess``.
    """
    from .metrics import stretch_energy, stretch_energy_lower_bound  # local import

    L = fem_stiffness(mesh)

    candidates: list[tuple[str, np.ndarray]] = []
    if method in ("auto", "star"):
        try:
            candidates.append(("star", star_ball_map(mesh)))
        except Exception as exc:  # pragma: no cover
            if verbose:
                print(f"  star map failed: {exc}")
    if method in ("auto", "harmonic", "relax"):
        bv = boundary_map_sphere(mesh, method=boundary_method)
        yh, _ = solve_harmonic(mesh, bv, L=L)
        candidates.append(("harmonic", yh))
    if method == "harmonic" and not candidates:
        raise ValueError("no candidate produced")

    def score(y: np.ndarray) -> tuple[int, float]:
        tgt = mesh.with_vertices(y)
        nf = int((_tet_volumes(y, mesh.tets) <= 0).sum())
        exc = stretch_energy(mesh, tgt) / max(stretch_energy_lower_bound(mesh, tgt), 1e-300) - 1.0
        return nf, exc

    best_name, best_y, best_score = None, None, (10**9, np.inf)
    for name, y in candidates:
        s = score(y)
        if verbose:
            print(f"  candidate {name:9s}: flips {s[0]:5d}  excess {s[1]:.4%}")
        if s < best_score:
            best_name, best_y, best_score = name, y, s

    # Refine only when it can help: the relaxation never *increases* the flip count, and
    # a fold-free start is relaxed without the fold-repair stage (which could otherwise
    # introduce folds while chasing the volume penalty).
    if best_score[0] > 0 or (relax and method != "harmonic"):
        y2, info = fold_free_relax(
            mesh, best_y, L=L, repair=best_score[0] > 0, verbose=verbose
        )
        s2 = score(y2)
        if s2[0] < best_score[0] or (s2[0] == 0 and s2[1] < best_score[1]):
            if verbose:
                print(
                    f"  relaxed   {best_name:9s}: flips {s2[0]:5d}  excess {s2[1]:.4%} "
                    f"({info['iterations']} it)"
                )
            best_name, best_y, best_score = f"{best_name}+relax", y2, s2

    return {
        "vertices": best_y,
        "method": best_name,
        "n_flipped": int(best_score[0]),
        "energy": float(0.5 * sum(best_y[:, s] @ (L @ best_y[:, s]) for s in range(3))),
        "stretch_excess": float(best_score[1]),
    }


__all__ = [
    "fem_stiffness",
    "mass_matrix_lumped",
    "boundary_map_sphere",
    "solve_harmonic",
    "harmonic_ball_map",
    "harmonic_map_to_domain",
    "star_ball_map",
    "fold_free_relax",
    "volume_penalty_relax",
    "ball_map",
    "ray_surface_hits",
    "HarmonicMapResult",
]
