"""Discrete optimal mass transport on a convex polyhedron, via Laguerre / power diagrams.

This is the numerical heart of the "optimal transport" half of the pipeline.

Setting
-------
Given

* a convex polyhedral domain ``P`` in R^3 (half-space representation),
* ``N`` sites ``p_1..p_N`` in ``P`` with prescribed masses ``m_i > 0``, ``sum m_i = |P|``,

find weights ``f_1..f_N`` such that the **Laguerre cells**

    Lag_i(f) = { x in P : |x - p_i|^2 - f_i  <=  |x - p_j|^2 - f_j   for all j }

have volumes ``|Lag_i(f)| = m_i``.  This is exactly the semi-discrete optimal transport
problem between the Lebesgue measure on ``P`` and the atomic measure ``sum_i m_i
delta_{p_i}``: the map ``T(x) = p_i`` for ``x in Lag_i(f)`` is the (Brenier) optimal
transport map, and its inverse is the "cell assignment".

The duality used here is the standard one (Aurenhammer, Merigot, Levy):

    maximize  D(f) = sum_i m_i f_i + integral_P min_i ( |x - p_i|^2 - f_i ) dx

    dD/df_i         = m_i - |Lag_i(f)|                       (gradient)
    -d^2 D/df_i df_j = w_ij = |F_ij| / (2 |p_i - p_j|)        (a graph Laplacian)

with ``|F_ij|`` the area of the face shared by cells ``i`` and ``j``.  ``D`` is concave,
so a damped Newton method with an Armijo line search converges globally; the Hessian is a
weighted graph Laplacian, which is singular only by the constant vector (adding a constant
to every weight does not change the diagram), so we take the minimum-norm Newton step.

Everything is exact polyhedral geometry -- no quadrature, no Monte-Carlo.  Convex cell
volumes come from ``qhull`` through :mod:`scipy.spatial`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, HalfspaceIntersection, cKDTree

# --------------------------------------------------------------------------------------
# convex polyhedron
# --------------------------------------------------------------------------------------


@dataclass
class ConvexPolyhedron:
    """Convex polyhedron as an intersection of half-spaces ``A x + b <= 0``.

    ``points`` (the vertices, any point of the convex hull) is kept for convenience.
    """

    planes: np.ndarray  # (k, 4): [nx, ny, nz, b]  meaning  n . x + b <= 0
    points: np.ndarray  # (v, 3) points whose convex hull is the polyhedron
    tol: float = 1e-9

    def __post_init__(self) -> None:
        self.planes = np.asarray(self.planes, dtype=np.float64)
        self.points = np.asarray(self.points, dtype=np.float64)

    # ---------------------------------------------------------------- constructors
    @classmethod
    def from_points(cls, points: np.ndarray, tol: float = 1e-9) -> "ConvexPolyhedron":
        pts = np.asarray(points, dtype=np.float64)
        hull = ConvexHull(pts)
        eq = hull.equations  # normalised so that  eq[:, :3] . x + eq[:, 3] <= 0
        return cls(planes=eq.copy(), points=pts, tol=tol)

    @classmethod
    def from_halfspaces(
        cls, planes: np.ndarray, interior_point: np.ndarray, tol: float = 1e-9
    ) -> "ConvexPolyhedron":
        planes = np.asarray(planes, dtype=np.float64)
        hs = HalfspaceIntersection(planes, np.asarray(interior_point, dtype=np.float64))
        return cls(planes=planes, points=np.asarray(hs.intersections), tol=tol)

    # ---------------------------------------------------------------- queries
    def volume(self) -> float:
        try:
            return float(ConvexHull(self.points).volume)
        except Exception:  # pragma: no cover
            return 0.0

    def sample_interior(self, rng: np.random.Generator | None = None) -> np.ndarray:
        rng = rng or np.random.default_rng(0)
        return self.points.mean(axis=0) + 1e-3 * rng.normal(size=3)

    def scaled(self, factor: float) -> "ConvexPolyhedron":
        """Scale about the origin (used to fit a target region)."""
        return ConvexPolyhedron.from_points(self.points * factor, tol=self.tol)

    def __repr__(self) -> str:  # pragma: no cover
        return f"ConvexPolyhedron({len(self.planes)} planes, {len(self.points)} vertices, |P|={self.volume():.6g})"


# --------------------------------------------------------------------------------------
# power / Laguerre diagram
# --------------------------------------------------------------------------------------


@dataclass
class PowerDiagram:
    sites: np.ndarray
    weights: np.ndarray
    cells: list[np.ndarray]  # vertices of each cell (possibly empty)
    volumes: np.ndarray
    neighbour_face_areas: list[dict[int, float]]  # cell i -> {j: area of shared face}
    n_empty: int = 0
    simplices: list[np.ndarray] = field(default_factory=list)  # hull triangles per cell
    n_planes_used: int = 0  # total half-spaces actually handed to qhull (profiling)

    def gradient(self, masses: np.ndarray) -> np.ndarray:
        """``dD/df = m_i - |Lag_i|``."""
        return masses - self.volumes

    def hessian(self, sites: np.ndarray, n: int | None = None) -> sp.csr_matrix:
        """The (positive semi-definite) graph Laplacian ``Lambda = -Hessian(D)``."""
        n = n or len(self.sites)
        rows, cols, vals = [], [], []
        for i, nb in enumerate(self.neighbour_face_areas):
            for j, area in nb.items():
                if j <= i:
                    continue
                d = float(np.linalg.norm(sites[i] - sites[j]))
                if d <= 0:
                    continue
                w = area / (2.0 * d)
                rows += [i, j, i, j]
                cols += [i, j, j, i]
                vals += [w, w, -w, -w]
        A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
        diag = np.asarray(-A.sum(axis=1)).ravel()
        return (A + sp.diags(diag)).tocsr()


def _domain_plane_subset(dom: np.ndarray, point: np.ndarray, radius: float) -> np.ndarray:
    """Domain planes that can possibly cut a cell contained in ``B(point, radius)``.

    ``dom`` rows are ``[n, b]`` with ``|n| = 1`` and the plane at signed distance
    ``n . x + b`` from ``x``.  A plane further than ``radius`` from the point cannot touch
    the cell.
    """
    dist = np.abs(dom[:, :3] @ point + dom[:, 3])
    keep = dist <= radius
    return dom[keep] if keep.any() else dom


def laguerre_cells(
    sites: np.ndarray,
    weights: np.ndarray,
    domain: ConvexPolyhedron,
    neighbours: list[np.ndarray] | int = 48,
    tol: float = 1e-9,
    want_faces: bool = False,
    domain_pad: float = 1.25,
) -> PowerDiagram:
    """Compute the Laguerre tessellation of ``domain`` for the given sites and weights.

    ``neighbours`` is either an explicit neighbour list or the number ``k`` of nearest
    sites considered for each cell (a superset of the true power-diagram neighbours in
    practice; see :func:`verify_neighbour_cutoff`).

    ``domain_pad`` enlarges the radius used to discard domain planes that cannot cut the
    cell; this is a pure speed optimisation whose result is exact whenever the cell really
    is contained in that ball (it always is for a reasonable padding factor).
    """
    sites = np.asarray(sites, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    n = len(sites)

    if isinstance(neighbours, (int, np.integer)):
        k = min(int(neighbours), n - 1)
        tree = cKDTree(sites)
        dist, idx = tree.query(sites, k=k + 1)
        nbr_list = [np.asarray(row[1:], dtype=np.int64) for row in idx]
        radius = domain_pad * np.maximum(dist[:, -1], 1e-12)
    else:
        nbr_list = [np.asarray(x, dtype=np.int64) for x in neighbours]
        radius = np.full(n, np.inf)

    dom = domain.planes
    rng = np.random.default_rng(12345)

    cells: list[np.ndarray] = []
    simplices: list[np.ndarray] = []
    vols = np.empty(n)
    faces: list[dict[int, float]] = []
    n_empty = 0
    n_planes = 0

    for i in range(n):
        js = nbr_list[i]
        # bisector with j:  2 (p_j - p_i) . x  +  (|p_i|^2 - |p_j|^2 + f_j - f_i) <= 0
        d = sites[js] - sites[i]
        b = (
            float(sites[i] @ sites[i])
            - np.einsum("ij,ij->i", sites[js], sites[js])
            + weights[js]
            - weights[i]
        )
        dsub = _domain_plane_subset(dom, sites[i], radius[i])
        hs = np.empty((len(js) + len(dsub), 4))
        hs[: len(js), :3] = 2.0 * d
        hs[: len(js), 3] = b
        hs[len(js) :] = dsub
        n_planes += len(hs)

        verts = _cell_vertices(hs, sites[i], rng)
        if verts is None or len(verts) < 4:
            cells.append(np.zeros((0, 3)))
            simplices.append(np.zeros((0, 3), dtype=np.int64))
            vols[i] = 0.0
            faces.append({})
            n_empty += 1
            continue

        try:
            hull = ConvexHull(verts)
            vol = float(hull.volume)
        except Exception:
            cells.append(np.zeros((0, 3)))
            simplices.append(np.zeros((0, 3), dtype=np.int64))
            vols[i] = 0.0
            faces.append({})
            n_empty += 1
            continue

        cells.append(verts)
        simplices.append(hull.simplices)
        vols[i] = vol
        faces.append(
            _face_areas(hull, verts, i, js, sites, weights, dom, tol) if want_faces else {}
        )

    return PowerDiagram(
        sites, weights, cells, vols, faces, n_empty, simplices, n_planes
    )


def _cell_vertices(hs: np.ndarray, site: np.ndarray, rng: np.random.Generator) -> np.ndarray | None:
    """Vertices of ``{x : hs @ [x,1] <= 0}``; ``None`` when empty / degenerate."""
    for interior in _interior_candidates(hs, site, rng):
        try:
            out = HalfspaceIntersection(hs, interior, qhull_options="QJ")
            v = np.asarray(out.intersections, dtype=np.float64)
            if len(v) >= 4:
                return v
        except Exception:
            continue
    return None


def _interior_candidates(hs: np.ndarray, site: np.ndarray, rng: np.random.Generator):
    """Yield candidate strictly-interior points, best first."""
    # 1. the site itself (inside its own cell for a valid diagram)
    yield site
    # 2. shifted site
    for s in (1e-7, 1e-5, 1e-3):
        yield site + s * rng.normal(size=3)
    # 3. Chebyshev centre of the half-space intersection (linear program)
    #    maximize t  s.t.  A x + b + t <= 0,  t <= 1
    A = hs[:, :3]
    b = hs[:, 3]
    c = np.array([0.0, 0.0, 0.0, -1.0])
    A_ub = np.hstack([A, np.ones((len(A), 1))])
    bounds = [(None, None)] * 3 + [(0.0, 1.0)]
    try:
        res = linprog(c, A_ub=A_ub, b_ub=-b, bounds=bounds, method="highs")
        if res.success and res.x[3] > 0:
            yield res.x[:3]
    except Exception:
        return


def _face_areas(
    hull: ConvexHull,
    verts: np.ndarray,
    i: int,
    js: np.ndarray,
    sites: np.ndarray,
    weights: np.ndarray,
    dom_planes: np.ndarray,
    tol: float,
) -> dict[int, float]:
    """Area of the face shared with each neighbouring cell (keyed by site index).

    Fully vectorised over the facets of the cell.
    """
    tri = verts[hull.simplices]  # (f,3,3)
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = 0.5 * np.linalg.norm(cross, axis=1)
    if not np.any(area > 0):
        return {}
    cent = tri.mean(axis=1)  # (f,3)

    d = sites[js] - sites[i]  # (k,3)
    b = (
        float(sites[i] @ sites[i])
        - np.einsum("ij,ij->i", sites[js], sites[js])
        + weights[js]
        - weights[i]
    )
    # value of each bisector plane at every facet centroid:  (k, f)
    val = 2.0 * (d @ cent.T) + b[:, None]
    scale = np.maximum(np.linalg.norm(2.0 * d, axis=1), 1.0)  # (k,)
    rel = np.abs(val) / scale[:, None]
    best = np.argmin(rel, axis=0)
    ok = rel[best, np.arange(len(cent))] < 1e-6

    out: dict[int, float] = {}
    for f_idx in np.flatnonzero(ok):
        jj = int(js[best[f_idx]])
        out[jj] = out.get(jj, 0.0) + float(area[f_idx])
    del dom_planes, tol, cross
    return out


def verify_neighbour_cutoff(
    sites: np.ndarray,
    weights: np.ndarray,
    domain: ConvexPolyhedron,
    k: int,
    n_samples: int = 12,
    seed: int = 0,
) -> float:
    """Relative volume error introduced by the kNN neighbour cutoff (0 is exact).

    Recomputes a random sample of cells twice -- once with only the ``k`` nearest
    neighbours, once with *every* site -- and returns the largest relative volume
    difference.  Use it to pick ``k`` for a given point set.
    """
    rng = np.random.default_rng(seed)
    n = len(sites)
    full = laguerre_cells(sites, weights, domain, neighbours=max(n - 1, 1))
    idx = rng.choice(n, size=min(n_samples, n), replace=False)
    err = 0.0
    kk = min(k, n - 1)
    tree = cKDTree(sites)
    for i in idx:
        _, nn = tree.query(sites[i], k=kk + 1)
        nbrs = {int(x) for x in nn[1:]}
        dom = domain.planes
        d = sites[list(nbrs)] - sites[i]
        b = (
            float(sites[i] @ sites[i])
            - np.einsum("ij,ij->i", sites[list(nbrs)], sites[list(nbrs)])
            + weights[list(nbrs)]
            - weights[i]
        )
        hs = np.vstack([np.hstack([2.0 * d, b[:, None]]), dom])
        v = _cell_vertices(hs, sites[i], rng)
        if v is None:
            err = 1.0
            break
        try:
            vol = float(ConvexHull(v).volume)
        except Exception:
            vol = 0.0
        err = max(err, abs(vol - full.volumes[i]) / max(full.volumes[i], 1e-300))
    return float(err)


# --------------------------------------------------------------------------------------
# solver
# --------------------------------------------------------------------------------------


@dataclass
class OTResult:
    weights: np.ndarray
    diagram: PowerDiagram
    mass_error: np.ndarray
    iterations: int
    history: list[float] = field(default_factory=list)

    def summary(self) -> str:
        e = self.mass_error
        rel = e / max(self.diagram.volumes.mean(), 1e-300)
        return (
            f"OTResult: {self.iterations} iterations, "
            f"max |vol - m| / mean(vol) = {np.abs(rel).max():.3e}, "
            f"empty cells {self.diagram.n_empty}, dual {self.history[-1]:.8g}"
        )


def dual_value(diagram: PowerDiagram, masses: np.ndarray) -> float:
    """The Kantorovich dual ``D(f) = sum m_i f_i + integral_P min_i(|x-p_i|^2 - f_i)``.

    The integral is exact: each cell is split into tetrahedra through its vertex centroid,
    reusing the hull connectivity already computed for the cell volume.
    """
    total = float(masses @ diagram.weights)
    for i, verts in enumerate(diagram.cells):
        if len(verts) < 4:
            continue
        g = verts.mean(axis=0)
        simp = diagram.simplices[i]
        _, mom = _tet_fan_moments(g, verts, simp, diagram.sites[i])
        total += mom - diagram.weights[i] * diagram.volumes[i]
    return total


def _tet_fan_moments(apex: np.ndarray, verts: np.ndarray, simplices: np.ndarray, p: np.ndarray):
    """``sum_T |T| * E_{x~T}|x-p|^2`` over the tetrahedra ``(apex, triangle)``."""
    a = verts[simplices]  # (k,3,3)
    e1 = a[:, 0] - apex
    e2 = a[:, 1] - apex
    e3 = a[:, 2] - apex
    vol = np.abs(np.einsum("ij,ij->i", np.cross(e1, e2), e3)) / 6.0
    q = apex - p
    # barycentric moments on a tetrahedron: E[s_i]=1/4, E[s_i^2]=1/10, E[s_i s_j]=1/20
    qq = float(np.dot(q, q))
    qe = (e1 @ q) + (e2 @ q) + (e3 @ q)
    ee = (
        np.einsum("ij,ij->i", e1, e1)
        + np.einsum("ij,ij->i", e2, e2)
        + np.einsum("ij,ij->i", e3, e3)
    )
    cross = (
        np.einsum("ij,ij->i", e1, e2)
        + np.einsum("ij,ij->i", e1, e3)
        + np.einsum("ij,ij->i", e2, e3)
    )
    val = qq + 0.5 * qe + 0.1 * ee + 0.1 * cross
    total_vol = float(vol.sum())
    return total_vol, float(np.sum(vol * val))


def solve_ot(
    sites: np.ndarray,
    masses: np.ndarray,
    domain: ConvexPolyhedron,
    weights0: np.ndarray | None = None,
    neighbours: list[np.ndarray] | int = 48,
    max_iter: int = 40,
    tol: float = 1e-6,
    use_newton: bool = True,
    verbose: bool = False,
    exact_line_search_value: bool = True,
) -> OTResult:
    """Solve the semi-discrete OT problem: weights with ``|Lag_i| = m_i``.

    Parameters
    ----------
    sites : (N,3)
    masses : (N,)  positive, ``sum(masses) == domain.volume()`` (a warning is issued
        otherwise and the masses are rescaled).
    domain : :class:`ConvexPolyhedron`
    weights0 : (N,) optional initial weights.  The default ``0`` gives the ordinary Voronoi
        diagram of the sites, which has no empty cell when the sites are distinct; the
        alternative ``|p_i|^2`` gives the "conical" diagram used in some OT papers and is a
        poor start here because most cells are then empty.

    Returns
    -------
    :class:`OTResult`
    """
    sites = np.asarray(sites, dtype=np.float64)
    masses = np.asarray(masses, dtype=np.float64)
    n = len(sites)

    target = domain.volume()
    total = masses.sum()
    if abs(total - target) > 1e-8 * max(target, 1.0):
        masses = masses * (target / total)

    if weights0 is None:
        weights = np.zeros(n)
    else:
        weights = np.array(weights0, dtype=np.float64)

    # neighbour lists (fixed during the iteration -- the diagram topology is stable)
    if isinstance(neighbours, (int, np.integer)):
        k = min(int(neighbours), n - 1)
        tree = cKDTree(sites)
        dist, idx = tree.query(sites, k=k + 1)
        nbr_list: list[np.ndarray] | int = [np.asarray(row[1:], dtype=np.int64) for row in idx]
        # a site can only move by a fraction of the local site spacing before the diagram
        # topology changes -- this is the standard step clamp of semi-discrete OT solvers
        spacing = float(np.median(dist[:, 1])) if k >= 1 else 1.0
        step_cap = (4.0 * max(spacing, 1e-12)) ** 2
    else:
        nbr_list = neighbours
        step_cap = np.inf

    history: list[float] = []
    diag = laguerre_cells(sites, weights, domain, nbr_list, want_faces=use_newton)
    err = np.abs(diag.volumes - masses).max() / max(masses.mean(), 1e-300)
    history.append(float(err))
    if verbose:
        print(f"  iter   0  rel. mass error {err:.3e}  empty {diag.n_empty}")

    it = 0
    d0 = dual_value(diag, masses) if exact_line_search_value else 0.0
    for it in range(1, max_iter + 1):
        g = masses - diag.volumes  # gradient of the dual
        if np.abs(g).max() / max(masses.mean(), 1e-300) < tol:
            break

        if use_newton:
            Lam = diag.hessian(sites, n)
            step = spla.lsqr(Lam, g, atol=1e-12, btol=1e-12, iter_lim=500)[0]
        else:
            face = np.array(
                [sum(v.values()) for v in diag.neighbour_face_areas], dtype=np.float64
            )
            face[face <= 0] = face[face > 0].mean() if (face > 0).any() else 1.0
            step = g / face

        if not np.all(np.isfinite(step)) or np.linalg.norm(step) == 0:
            step = g / max(np.abs(g).max(), 1e-300)
        # clamp the step so that the diagram cannot change topology in one move
        mx = np.abs(step).max()
        if mx * 1.0 > step_cap:
            step = step * (step_cap / mx)

        slope = float(g @ step)
        alpha = 1.0
        accepted = False
        for _ in range(14):
            cand = weights + alpha * step
            d1_diag = laguerre_cells(sites, cand, domain, nbr_list, want_faces=use_newton)
            # An empty cell makes the dual unbounded above (the term m_i f_i grows while
            # the integral does not), so a step that empties a cell must never be accepted
            # -- this is the roll-back rule of the damped Newton method.
            if d1_diag.n_empty == 0:
                if exact_line_search_value:
                    d1 = dual_value(d1_diag, masses)
                    ok = d1 >= d0 + 1e-4 * alpha * slope
                else:
                    ok = np.abs(masses - d1_diag.volumes).max() < np.abs(g).max()
                if ok:
                    weights = cand
                    diag = d1_diag
                    if exact_line_search_value:
                        d0 = d1
                    accepted = True
                    break
            alpha *= 0.5
        if not accepted:
            if diag.n_empty > 0:
                # cannot make progress from a diagram with empty cells: stop cleanly
                break
            break

        err = float(np.abs(diag.volumes - masses).max() / max(masses.mean(), 1e-300))
        history.append(err)
        if verbose:
            print(
                f"  iter {it:3d}  rel. mass error {err:.3e}  alpha {alpha:.3g}  "
                f"empty {diag.n_empty}"
            )
        if err < tol:
            break

    return OTResult(
        weights=weights,
        diagram=diag,
        mass_error=diag.volumes - masses,
        iterations=it,
        history=history,
    )


__all__ = [
    "ConvexPolyhedron",
    "PowerDiagram",
    "laguerre_cells",
    "solve_ot",
    "dual_value",
    "verify_neighbour_cutoff",
    "OTResult",
]
