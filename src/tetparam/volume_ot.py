"""Stage 2/3: from the ball map to a *volume-preserving* map.

Two routes, and they are the same mathematics seen from two sides.

Route A -- explicit optimal mass transport
------------------------------------------
The harmonic map ``h`` pushes the volume measure of the solid forward to a measure ``mu`` on
the ball whose density is piecewise constant on the image mesh::

    mu|_{h(tau_i)} = |tau_i| / |h(tau_i)| .

Finding the map that rearranges the *uniform* Lebesgue measure of the ball into ``mu`` is a
Brenier / L2 optimal transport problem.  Its semi-discrete discretisation is a
capacity-constrained Laguerre (power) tessellation: find weights ``f`` such that the power
cells of the sites tile the ball with ``|Lag_i(f)| = |tau_i| * |B| / |M|``.  That is exactly
:func:`tetparam.laguerre.solve_ot`, and the resulting assignment ``cell_i <-> tau_i`` is the
optimal transport plan.  ``T(x) = p_i`` on ``Lag_i`` is the optimal (Brenier) map from the
uniform measure to the atomic approximation of ``mu``.

Route B -- the variational form of the same problem (no explicit OT solve)
-------------------------------------------------------------------------
Let ``r_tau = |f(tau)| / |tau|`` be the volume ratio.  The volumetric stretch energy

    E_V(f) = sum_tau |f(tau)|^2 / |tau| = sum_tau |tau| r_tau^2

satisfies, because the boundary is pinned on the unit sphere (so ``sum_tau |f(tau)| = |B|``
is a constant of the problem, see the note below),

    E_V(f) >= (sum_tau |f(tau)|)^2 / sum_tau |tau| = V(f)^2 / V(M),

with **equality if and only if every ``r_tau`` is the same constant**, i.e. iff ``f`` is
volume-preserving.  Minimising ``E_V`` therefore *is* the optimal transport problem in
variational form -- the statement is the lower-bound theorem of Yueh et al. / Huang et al.
and is the basis of VSEM and n-VSE.  This route needs no power diagram at all, and it is
what we use to actually produce the map; Route A is used to expose the transport plan, to
control a prescribed measure and to drive the tetrahedral remeshing.

The constant-volume property: with every boundary vertex on the unit sphere and no folded
tetrahedron, the image of the map is the polyhedral ball bounded by the image of the
boundary, whose volume is fixed.  Hence ``sum_tau |f(tau)|`` is automatically constant on
the fold-free set, and minimising ``E_V`` drives ``r_tau`` to its mean.

(Both normalisations of the volumetric stretch energy appear in the literature; see
``metrics.stretch_energy`` (this one) and ``stiffness_stretch_energy`` (the modified
cotangent-weight form, equivalent up to a factor).)
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize
from scipy.spatial import ConvexHull, cKDTree

from .harmonic import _tet_volumes, _volume_and_grad_ops, fem_stiffness
from .laguerre import ConvexPolyhedron, laguerre_cells, solve_ot
from .mesh import TetMesh


# --------------------------------------------------------------------------------------
# Route A: the explicit optimal transport problem on the ball
# --------------------------------------------------------------------------------------
def ball_domain(y: np.ndarray, boundary_mask: np.ndarray) -> ConvexPolyhedron:
    """The convex polyhedron that the ball map maps the solid onto.

    It is the convex hull of the images of the boundary vertices.  For all the meshes built
    by :mod:`tetparam.generate` (and for any mesh whose boundary is a convex triangulated
    sphere) this is exactly the image region; the residual mismatch is reported by
    :func:`domain_mismatch`.
    """
    return ConvexPolyhedron.from_points(y[boundary_mask])


def polyhedron_volume_from_faces(vertices: np.ndarray, faces: np.ndarray) -> float:
    """Volume of the region bounded by a closed triangle mesh (divergence theorem)."""
    a = vertices[faces[:, 0]]
    b = vertices[faces[:, 1]]
    c = vertices[faces[:, 2]]
    return float(np.einsum("ij,ij->i", np.cross(a, b), c).sum() / 6.0)


def domain_mismatch(mesh: TetMesh, y: np.ndarray) -> float:
    """Relative volume difference between the image mesh and the convex hull of its boundary."""
    bnd = mesh.boundary()
    v_mesh = abs(polyhedron_volume_from_faces(y, bnd["faces"]))
    dom = ball_domain(y, bnd["is_boundary_vertex"])
    v_hull = dom.volume()
    return float(v_mesh / v_hull - 1.0)


def solve_ball_ot(
    mesh: TetMesh,
    y: np.ndarray,
    sites: str = "tet_centroids",
    neighbours: int | list = 64,
    max_iter: int = 60,
    tol: float = 1e-8,
    verbose: bool = False,
):
    """Solve the capacity-constrained Laguerre tessellation of the image ball.

    ``sites='tet_centroids'`` gives one site per image tetrahedron with mass ``|tau_i|`` --
    the discretisation of ``mu`` used by the transport construction.
    ``sites='vertices'`` gives one site per mesh vertex with the barycentric dual-cell mass
    ``sum_{tau ni v} |tau| / 4`` -- the discretisation used for measure-driven remeshing.

    Returns a dict with the solved :class:`~tetparam.laguerre.OTResult`, the domain and the
    masses.
    """
    bnd = mesh.boundary()
    dom = ball_domain(y, bnd["is_boundary_vertex"])

    if sites == "tet_centroids":
        pts = y[mesh.tets].mean(axis=1)
        masses = np.abs(mesh.signed_volumes())
    elif sites == "vertices":
        pts = y
        v = np.abs(mesh.signed_volumes()) / 4.0
        masses = np.zeros(mesh.n_vertices)
        for k in range(4):
            np.add.at(masses, mesh.tets[:, k], v)
    else:
        raise ValueError("sites must be 'tet_centroids' or 'vertices'")

    keep = masses > 0
    res = solve_ot(
        pts[keep],
        masses[keep],
        dom,
        neighbours=neighbours,
        max_iter=max_iter,
        tol=tol,
        verbose=verbose,
    )
    return {
        "ot": res,
        "domain": dom,
        "sites": pts[keep],
        "masses": masses[keep],
        "mass_error_rel": float(
            np.abs(res.mass_error).max() / max(masses[keep].mean(), 1e-300)
        ),
        "domain_volume": dom.volume(),
        "domain_mismatch": domain_mismatch(mesh, y),
        "n_sites": int(keep.sum()),
    }


# --------------------------------------------------------------------------------------
# Route A':  the composition used by the paper --  F = phi^{-1} o psi
# --------------------------------------------------------------------------------------
def cell_centroids(diagram) -> np.ndarray:
    """Volume centroid (mass centre) of every Laguerre cell.

    For a tetrahedron with apex ``g`` and vertices ``a,b,c`` the centroid is
    ``(g+a+b+c)/4``, so the cell centroid is the volume-weighted average of its
    tetrahedral fan through the cell's vertex centroid.
    """
    out = np.zeros((len(diagram.cells), 3))
    for i, verts in enumerate(diagram.cells):
        if len(verts) < 4:
            out[i] = diagram.sites[i]
            continue
        g = verts.mean(axis=0)
        simp = diagram.simplices[i]
        a = verts[simp]  # (f,3,3)
        vol = np.abs(
            np.einsum("ij,ij->i", np.cross(a[:, 0] - g, a[:, 1] - g), a[:, 2] - g)
        ) / 6.0
        cent = 0.25 * (g[None, :] + a.sum(axis=1))
        out[i] = (vol[:, None] * cent).sum(axis=0) / max(vol.sum(), 1e-300)
    return out


def ot_vertex_map(
    mesh: TetMesh,
    y_ball: np.ndarray,
    neighbours: int | list = 48,
    max_iter: int = 40,
    tol: float = 1e-6,
    snap_boundary: bool = True,
    verbose: bool = False,
) -> dict:
    """The optimal-transport composition of the paper: ``F = phi^{-1} o psi``.

    ``psi`` is the harmonic (ball) map, so the sites are the images of the *mesh vertices*;
    the mass attached to vertex ``i`` is its barycentric dual volume
    ``nu_i = sum_{tau ni v_i} |tau| / 4`` (this is the measure used by the authors -- it was
    verified against the ``weight`` attribute of their released test meshes).  Solving the
    semi-discrete OT problem gives the power diagram; the inverse transport map
    ``phi^{-1}`` sends each site to the mass centre of its power cell, hence

        F(v_i) = centroid( W_i ) .

    ``snap_boundary`` re-projects the boundary vertices onto the unit sphere afterwards,
    which keeps the image region equal to the ball (the raw cell centroids of boundary
    vertices lie slightly inside).
    """
    bnd = mesh.boundary()
    is_b = bnd["is_boundary_vertex"]
    dom = ball_domain(y_ball, is_b)

    v4 = np.abs(mesh.signed_volumes()) / 4.0
    masses = np.zeros(mesh.n_vertices)
    for k in range(4):
        np.add.at(masses, mesh.tets[:, k], v4)

    keep = masses > 0
    res = solve_ot(
        y_ball[keep],
        masses[keep],
        dom,
        neighbours=neighbours,
        max_iter=max_iter,
        tol=tol,
        verbose=verbose,
    )
    cent = cell_centroids(res.diagram)
    y_new = np.array(y_ball, dtype=np.float64, copy=True)
    y_new[keep] = cent
    if snap_boundary:
        b = np.linalg.norm(y_new[is_b], axis=1, keepdims=True)
        y_new[is_b] = y_new[is_b] / np.maximum(b, 1e-300)

    return {
        "vertices": y_new,
        "weights": res.weights,
        "diagram": res.diagram,
        "iterations": res.iterations,
        "mass_error_rel": float(np.abs(res.mass_error).max() / max(masses[keep].mean(), 1e-300)),
        "domain": dom,
        "domain_mismatch": domain_mismatch(mesh, y_ball),
        "n_sites": int(keep.sum()),
        "snap_boundary": snap_boundary,
    }


# --------------------------------------------------------------------------------------
# Route B: volumetric stretch energy minimisation (the variational form)
# --------------------------------------------------------------------------------------def stiffness_stretch_energy(mesh: TetMesh, y: np.ndarray) -> float:
    """``E_V(f) = 1/2 trace(f^T L_V(f) f)`` with modified cotangent weights (n-VSE form)."""
    v = np.abs(_tet_volumes(y, mesh.tets))
    L = fem_stiffness(mesh, weights=v)
    return float(0.5 * sum(y[:, s] @ (L @ y[:, s]) for s in range(3)))


def _energy_and_grad(
    x: np.ndarray,
    mesh: TetMesh,
    y_full: np.ndarray,
    free: np.ndarray,
    dV: np.ndarray,
    v0: np.ndarray,
    barrier: float,
) -> tuple[float, np.ndarray]:
    """``E_V`` (this module's normalisation) and its gradient w.r.t. the free vertices."""
    y = y_full.copy()
    y[free] = x.reshape(-1, 3)
    v = _tet_volumes(y, mesh.tets)
    r = v / v0
    E = float(np.sum(r * v))  # sum |tau_y|^2 / |tau_0| = sum (v/v0) * v
    coef = 2.0 * r
    g = np.zeros_like(y)
    for k in range(4):
        np.add.at(g, mesh.tets[:, k], coef[:, None] * dV[:, k])
    if barrier > 0:
        neg = v < 0
        if neg.any():
            E += barrier * float(np.sum((v[neg] / v0[neg]) ** 2))
            coef_b = 2.0 * barrier * (v[neg] / v0[neg] ** 2)
            gb = np.zeros_like(y)
            for k in range(4):
                np.add.at(gb, mesh.tets[neg, k], coef_b[:, None] * dV[neg, k])
            g = g + gb
    return E, g[free].ravel()


def vsem(
    mesh: TetMesh,
    y_init: np.ndarray,
    max_iter: int = 500,
    tol: float = 1e-14,
    barrier: float = 1.0,
    fix_boundary: bool = True,
    verbose: bool = False,
) -> dict:
    """Minimise the volumetric stretch energy with the boundary fixed on the sphere.

    This is the "no explicit optimal transport" route to a volume-preserving map: by the
    lower-bound theorem the minimiser has a constant volume ratio, i.e. it *is* an optimal
    transport map (the L2-optimal one among the maps that agree on the boundary).

    ``barrier`` adds ``lambda * sum max(0, -r_tau)^2`` so a fold cannot be created silently;
    a strictly fold-free start plus this term keeps the minimiser injective in practice.

    Returns a dict with ``vertices``, ``energy``, ``iterations``, ``n_flipped`` and the
    convergence history.
    """
    bnd = mesh.boundary()
    is_b = bnd["is_boundary_vertex"] if fix_boundary else np.zeros(mesh.n_vertices, bool)
    free = np.flatnonzero(~is_b)
    v0, dV = _volume_and_grad_ops(mesh)
    v0 = np.where(np.abs(v0) > 0, v0, 1e-300)

    x0 = y_init[free].ravel()
    history: list[float] = []

    def fun(x):
        E, g = _energy_and_grad(x, mesh, y_init, free, dV, v0, barrier)
        history.append(E)
        return E, g

    res = minimize(
        fun,
        x0,
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": max_iter, "ftol": tol, "gtol": 1e-12, "maxls": 50},
    )
    y = y_init.copy()
    y[free] = res.x.reshape(-1, 3)
    if fix_boundary:
        # keep the boundary exactly on the unit sphere
        b = np.linalg.norm(y[is_b], axis=1, keepdims=True)
        y[is_b] = y[is_b] / np.maximum(b, 1e-300)
    v = _tet_volumes(y, mesh.tets)
    if verbose:
        print(f"  VSEM: {res.nit} it, E_V {res.fun:.8g}, flips {int((v <= 0).sum())}")
    return {
        "vertices": y,
        "energy": float(res.fun),
        "iterations": int(res.nit),
        "n_flipped": int((v <= 0).sum()),
        "history": history,
        "success": bool(res.success),
        "message": str(res.message),
    }


# --------------------------------------------------------------------------------------
# IEM: let the boundary glide on the sphere
# --------------------------------------------------------------------------------------
def boundary_relax_sphere(
    mesh: TetMesh,
    y: np.ndarray,
    dV: np.ndarray,
    v0: np.ndarray,
    iters: int = 20,
    step: float = 0.2,
    verbose: bool = False,
) -> np.ndarray:
    """Projected gradient descent of ``E_V`` on the *boundary*, constrained to the sphere.

    This is the essential difference between VSEM and IEM: VSEM freezes the boundary
    mapping, which makes the lower bound of the volumetric stretch energy unreachable
    (the boundary alone fixes part of the volume distribution).  IEM lets the boundary
    points glide along the unit sphere, which is exactly a Riemannian gradient step:

        y_b  <-  y_b - alpha * (I - y_b y_b^T) grad E_V ,   then renormalise.

    The interior is kept at its fixed-boundary optimum between the boundary steps.
    """
    is_b = mesh.is_boundary_vertex
    free = np.flatnonzero(~is_b)
    y = y.copy()
    for _ in range(iters):
        v = _tet_volumes(y, mesh.tets)
        E0 = float(np.sum(v * v / v0))
        gb = _boundary_gradient(mesh, y, dV, v0)
        yb = y[is_b]
        # tangential (Riemannian) component of the gradient on the sphere
        gt = gb - np.sum(gb * yb, axis=1, keepdims=True) * yb
        if np.linalg.norm(gt) < 1e-14:
            break
        alpha = step
        accepted = False
        E1 = E0
        for _ in range(25):
            cand = y.copy()
            yb_new = yb - alpha * gt
            yb_new /= np.maximum(np.linalg.norm(yb_new, axis=1, keepdims=True), 1e-300)
            cand[is_b] = yb_new
            vv = _tet_volumes(cand, mesh.tets)
            if (vv > 0).all():
                E1 = float(np.sum(vv * vv / v0))
                if E1 < E0:
                    y = cand
                    accepted = True
                    break
            alpha *= 0.5
        if not accepted:
            break
        if verbose:
            print(f"    boundary relax: E_V {E0:.8g} -> {E1:.8g}  alpha {alpha:.3g}")
    del free
    return y


def _boundary_gradient(mesh: TetMesh, y: np.ndarray, dV: np.ndarray, v0: np.ndarray):
    """Gradient of ``E_V`` with respect to the boundary vertices (3-D Cartesian)."""
    v = _tet_volumes(y, mesh.tets)
    coef = 2.0 * v / v0
    g = np.zeros_like(y)
    for k in range(4):
        np.add.at(g, mesh.tets[:, k], coef[:, None] * dV[:, k])
    return g[mesh.is_boundary_vertex]


def iem(
    mesh: TetMesh,
    y_init: np.ndarray,
    outer: int = 10,
    inner_max_iter: int = 500,
    boundary_iters: int = 30,
    step: float = 0.2,
    verbose: bool = False,
) -> dict:
    """Isovolumetric-energy-style minimisation: alternate interior solve and boundary glide.

    Returns the same dict as :func:`vsem`.

    Reference: Liu, Huang, Lin, Yueh, "Isovolumetric Energy Minimization for Ball-Shaped
    Volume-Preserving Parameterizations of 3-Manifolds", arXiv:2407.19272 -- which shows
    that freeing the boundary on the sphere removes the fixed-boundary deficit of VSEM.
    """
    y = np.array(y_init, dtype=np.float64, copy=True)
    v0, dV = _volume_and_grad_ops(mesh)
    v0 = np.where(np.abs(v0) > 0, v0, 1e-300)
    E = float(np.sum(_tet_volumes(y, mesh.tets) ** 2 / v0))
    it_total = 0
    for o in range(outer):
        r = vsem(mesh, y, max_iter=inner_max_iter)
        y = r["vertices"]
        it_total += r["iterations"]
        y = boundary_relax_sphere(mesh, y, dV, v0, iters=boundary_iters, step=step,
                                  verbose=verbose)
        v = _tet_volumes(y, mesh.tets)
        E_new = float(np.sum(v * v / v0))
        if verbose:
            print(f"  IEM outer {o}: E_V {E:.8g} -> {E_new:.8g}")
        if abs(E - E_new) <= 1e-12 * max(abs(E), 1.0):
            E = E_new
            break
        E = E_new
    y = vsem(mesh, y, max_iter=inner_max_iter)["vertices"]
    v = _tet_volumes(y, mesh.tets)
    return {
        "vertices": y,
        "energy": float(np.sum(v * v / v0)),
        "iterations": it_total,
        "n_flipped": int((v <= 0).sum()),
        "success": True,
        "message": "iem",
    }


# --------------------------------------------------------------------------------------
# full stage
# --------------------------------------------------------------------------------------
def volume_preserving_map(
    mesh: TetMesh,
    ball: dict | None = None,
    method: str = "vsem",
    vse_iter: int = 500,
    verbose: bool = False,
) -> dict:
    """Run stage 1 + stage 2 and return the volume-preserving map to the ball.

    ``ball`` is the output of :func:`tetparam.harmonic.ball_map`; it is computed if omitted.
    Returns a dict with ``vertices``, the intermediate ball map, the VSEM result and the
    optimal-transport diagnostics.
    """
    from .harmonic import ball_map as _ball_map

    if ball is None:
        ball = _ball_map(mesh, method="auto", verbose=verbose)
    y0 = ball["vertices"]

    out: dict = {"ball": ball, "vertices": y0, "method": "identity"}
    if method in ("vsem", "both"):
        r = vsem(mesh, y0, max_iter=vse_iter, verbose=verbose)
        out["vsem"] = r
        out["vertices"] = r["vertices"]
        out["method"] = "vsem"

    if method in ("ot", "both"):
        try:
            ot = solve_ball_ot(mesh, out["vertices"], verbose=verbose)
            out["ot"] = ot
            if out["method"] == "identity":
                out["method"] = "ot"
        except Exception as exc:  # pragma: no cover
            out["ot_error"] = repr(exc)
            if verbose:
                print(f"  OT solve failed: {exc}")
    return out


__all__ = [
    "ball_domain",
    "polyhedron_volume_from_faces",
    "domain_mismatch",
    "solve_ball_ot",
    "cell_centroids",
    "ot_vertex_map",
    "vsem",
    "stiffness_stretch_energy",
    "volume_preserving_map",
]
