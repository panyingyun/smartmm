"""Stage 4: turn the volume-preserving map into an actual tetrahedral mesh of the solid.

The map ``f : M -> B`` produced by stages 1-2 is a piecewise-linear homeomorphism.  Its
inverse is therefore also piecewise linear, and it can be evaluated cheaply:

1. locate the query point ``q`` in the *image* mesh ``f(M)`` (a few candidate tetrahedra
   from a k-d tree on the image centroids, then a barycentric test),
2. read off the barycentric coordinates ``lambda``,
3. the pre-image is ``sum_k lambda_k v_k`` with ``v_k`` the *original* vertices.

Consequently, any tetrahedral mesh of the ball can be pulled back to a tetrahedral mesh of
``M`` **with the same connectivity**: the result is automatically conforming, has exactly
the same number of elements as the template, and its element sizes are the template sizes
divided by the local Jacobian of ``f``.

That is the "tetrahedral decomposition" step:

* uniform template + volume-preserving ``f``  ->  uniform-quality mesh of ``M``;
* any template density ``rho`` on the ball    ->  mesh of ``M`` with density
  ``rho(f(x)) * |det Df(x)|`` (the change-of-variables formula), which is exactly what the
  measure-controllable variant of the pipeline prescribes.

The inverse-map route avoids a second tetrahedralisation of ``M`` altogether; if the caller
instead wants the boundary of ``M`` re-meshed to a new resolution, the recommended recipe is
to sample points inside the ball, pull them back, and hand those points *together with the
original boundary surface* to a robust mesher (fTetWild/TetGen) -- see ``docs/``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from .mesh import TetMesh, tet_shape_gradients


@dataclass
class PointLocation:
    tet_index: np.ndarray  # (q,) index of the containing tetrahedron, -1 if outside
    barycentric: np.ndarray  # (q, 4) barycentric coordinates
    found: np.ndarray  # (q,) bool


def locate_points(
    image_mesh: TetMesh, queries: np.ndarray, candidates: int = 48, tol: float = 1e-6
) -> PointLocation:
    """Locate points of the ball in the image mesh ``f(M)``.

    A k-d tree over the tetrahedron centroids proposes ``candidates`` tetrahedra; the
    barycentric test decides.  Because the image mesh is a conforming tetrahedralisation,
    exactly one candidate contains the point (up to the mesh boundary and up to the
    floating-point tolerance ``tol`` -- points that lie exactly on a shared face or vertex
    need a tolerance of order ``1e-6`` relative to a unit ball).
    """
    queries = np.asarray(queries, dtype=np.float64)
    cent = image_mesh.centroids()
    tree = cKDTree(cent)
    k = min(candidates, image_mesh.n_tets)
    _, idx = tree.query(queries, k=k)

    grads, _ = tet_shape_gradients(image_mesh.vertices, image_mesh.tets)
    V = image_mesh.vertices[image_mesh.tets]  # (m,4,3)

    tet_index = np.full(len(queries), -1, dtype=np.int64)
    bary = np.zeros((len(queries), 4))
    found = np.zeros(len(queries), dtype=bool)

    for qi, cand in enumerate(idx):
        if np.ndim(cand) == 0:
            cand = np.array([cand])
        q = queries[qi]
        for t in cand:
            t = int(t)
            # lambda_k(x) = delta_{k0} + grad(lambda_k) . (x - v_0)
            d = q - V[t, 0]
            l0 = 1.0 + float(np.dot(grads[t, 0], d))
            l1 = float(np.dot(grads[t, 1], d))
            l2 = float(np.dot(grads[t, 2], d))
            lam4 = np.array([l0, l1, l2, 1.0 - l0 - l1 - l2])
            if (lam4 >= -tol).all() and (lam4 <= 1.0 + tol).all():
                tet_index[qi] = t
                bary[qi] = lam4
                found[qi] = True
                break
    return PointLocation(tet_index, bary, found)


def inverse_map(mesh: TetMesh, f_vertices: np.ndarray, queries: np.ndarray, **kw) -> np.ndarray:
    """``f^{-1}`` evaluated at points of the ball.

    Points that fall outside the image mesh (they can, by a fraction of a percent, because
    the image region is a triangulated sphere while a template may reach the convex hull)
    are projected onto the closest boundary triangle of the image mesh first, which is the
    mathematically correct pre-image for a boundary point.
    """
    image = mesh.with_vertices(f_vertices)
    loc = locate_points(image, queries, **kw)
    orig = mesh.vertices[mesh.tets]
    out = np.zeros((len(queries), 3))
    if loc.found.any():
        t = loc.tet_index[loc.found]
        lam = loc.barycentric[loc.found]
        out[loc.found] = np.einsum("qk,qkj->qj", lam, orig[t])
    if (~loc.found).any():
        q = queries[~loc.found]
        out[~loc.found] = _project_to_boundary(mesh, f_vertices, q, candidates=8)
    return out


def _project_to_boundary(
    mesh: TetMesh, f_vertices: np.ndarray, queries: np.ndarray, candidates: int = 8
) -> np.ndarray:
    """Pre-image of the closest point of the image boundary surface (per query)."""
    bnd = mesh.boundary()
    faces = bnd["faces"]
    tri = f_vertices[faces]  # (k,3,3)
    cent = tri.mean(axis=1)
    tree = cKDTree(cent)
    k = min(candidates, len(faces))
    _, idx = tree.query(queries, k=k)

    src_tri = mesh.vertices[faces]  # (k,3,3)
    out = np.empty((len(queries), 3))
    for qi in range(len(queries)):
        cand = np.atleast_1d(idx[qi])
        q = queries[qi]
        best = np.inf
        best_p = None
        for t in cand:
            t = int(t)
            p, lam = _closest_point_on_triangle(q, tri[t])
            d = np.linalg.norm(p - q)
            if d < best:
                best = d
                best_p = lam @ src_tri[t]
        out[qi] = best_p
    return out


def _closest_point_on_triangle(q: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Closest point of ``q`` on triangle ``t`` (3,3), and its clamped barycentric coords."""
    a, b, c = t
    ab, ac, aq = b - a, c - a, q - a
    d1, d2 = ab @ aq, ac @ aq
    if d1 <= 0 and d2 <= 0:
        return a, np.array([1.0, 0.0, 0.0])
    bp, bq = q - b, c - b
    d3, d4 = ab @ bp, ac @ bp
    if d3 >= 0 and d4 <= d3:
        return b, np.array([0.0, 1.0, 0.0])
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        return a + v * ab, np.array([1 - v, v, 0.0])
    cp, cq = q - c, b - c
    d5, d6 = ab @ cp, ac @ cp
    if d6 >= 0 and d5 <= d6:
        return c, np.array([0.0, 0.0, 1.0])
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        return a + w * ac, np.array([1 - w, 0.0, w])
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b), np.array([0.0, 1 - w, w])
    denom = 1.0 / (va + vb + vc)
    v, w = vb * denom, vc * denom
    u = 1 - v - w
    return u * a + v * b + w * c, np.array([u, v, w])


def pullback_remesh(
    mesh: TetMesh,
    f_vertices: np.ndarray,
    template: TetMesh | None = None,
    template_resolution: int = 6,
    **kw,
) -> TetMesh:
    """Tetrahedral mesh of ``M`` obtained by pulling a template ball mesh through ``f^{-1}``.

    ``template`` defaults to the structured cubed-sphere mesh of the unit ball at
    ``template_resolution``.  The output has the template's connectivity, so it is
    conforming by construction and has exactly the prescribed element count.
    """
    from .generate import ball_cube_mesh

    if template is None:
        template = ball_cube_mesh(template_resolution)
    pts = inverse_map(mesh, f_vertices, template.vertices, **kw)
    return TetMesh(pts, template.tets.copy())


def element_size_field(mesh: TetMesh, f_vertices: np.ndarray,
                       dV: np.ndarray | None = None) -> np.ndarray:
    """Predicted element size of the pulled-back mesh, per image tetrahedron.

    For a template of local size ``h_B`` in the ball, the pulled-back edge length scales
    like ``h_B / |det Df|^{1/3}`` (the standard change-of-variables rule for a metric).
    """
    if dV is None:
        dV = np.linalg.det(
            np.einsum(
                "mka,mkb->mab",
                f_vertices[mesh.tets],
                tet_shape_gradients(mesh.vertices, mesh.tets)[0],
            )
        )
    r = np.abs(dV)
    return r ** (-1.0 / 3.0)


def uniform_ball_template(
    image_mesh: TetMesh, n_interior: int = 400, seed: int = 0, smooth: int = 6
) -> TetMesh:
    """A quasi-uniform tetrahedral mesh of the *same* ball region as ``image_mesh``.

    The point set is ``n_interior`` quasi-uniform points inside the ball plus the boundary
    vertices of ``image_mesh``; its Delaunay tetrahedralisation therefore reproduces the
    image boundary (every boundary vertex of the image mesh is an extreme point of the
    point set), which is what makes the pull-back well defined everywhere.
    """
    from scipy.spatial import Delaunay

    from .generate import _odv_smooth, uniform_ball_points
    from .mesh import TetMesh as _TM

    rng = np.random.default_rng(seed)
    bnd = image_mesh.boundary()
    surface = image_mesh.vertices[bnd["boundary_vertices"]]
    # interior points: keep them well inside so the Delaunay boundary is the given surface
    inter = uniform_ball_points(n_interior, rng, radius=0.95)
    pts = np.vstack([surface, inter])
    frozen = np.arange(len(surface))
    if smooth > 0:
        pts = _odv_smooth(pts, frozen, n_iter=smooth, relax=0.5)
        # the smoothing must not push anything outside the ball nor off the surface
        pts[: len(surface)] = surface
        r = np.linalg.norm(pts, axis=1, keepdims=True)
        pts = np.where(r > 0.995, pts * (0.995 / np.maximum(r, 1e-300)), pts)
    tri = Delaunay(pts)
    t = _TM(pts, tri.simplices.astype(np.int64)).fix_orientation()
    v = t.signed_volumes()
    keep = np.abs(v) > 1e-12 * max(np.abs(v).max(), 1.0)
    if not keep.all():
        used = np.unique(t.tets[keep])
        remap = -np.ones(t.n_vertices, dtype=np.int64)
        remap[used] = np.arange(len(used))
        t = _TM(t.vertices[used], remap[t.tets[keep]])
    return t


def refine(mesh: TetMesh, levels: int = 1) -> TetMesh:
    """Uniform refinement of a tetrahedral mesh (each tet split by its vertex centroid).

    The boundary surface is preserved exactly, so the refined mesh is conforming with the
    original boundary -- the safest template when the ball map is strongly distorted.
    """
    verts = [mesh.vertices]
    tets = mesh.tets
    n = mesh.n_vertices
    for _ in range(levels):
        cent = verts[-1][tets].mean(axis=1)
        cidx = np.arange(n, n + len(tets))
        new_v = np.vstack([verts[-1], cent])
        parts = []
        for k in range(4):
            idx = [i for i in range(4) if i != k]
            parts.append(np.stack([tets[:, idx[0]], tets[:, idx[1]], tets[:, idx[2]], cidx], axis=1))
        tets = np.vstack(parts)
        verts.append(new_v)
        n = len(new_v)
    out = TetMesh(verts[-1], tets)
    from .mesh import signed_tet_volumes

    bad = signed_tet_volumes(out.vertices, out.tets) < 0
    if bad.any():
        tt = out.tets[bad].copy()
        tt[:, [1, 2]] = tt[:, [2, 1]]
        out.tets[bad] = tt
    return out


__all__ = [
    "PointLocation",
    "locate_points",
    "inverse_map",
    "pullback_remesh",
    "uniform_ball_template",
    "refine",
    "element_size_field",
]
