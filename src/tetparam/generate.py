"""Synthetic tetrahedral meshes of a topological 3-ball, for tests and demos.

Everything here is built from scipy only: no external mesh generator needed, so the
whole pipeline stays reproducible.  For production input use :meth:`TetMesh.from_tetgen`.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import Delaunay

from .mesh import TetMesh, signed_tet_volumes


# --------------------------------------------------------------------------------------
# point sets on / inside the unit ball
# --------------------------------------------------------------------------------------
def fibonacci_sphere(n: int) -> np.ndarray:
    """``n`` quasi-uniform points on the unit sphere (Fibonacci / golden angle spiral)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)  # polar angle, uniform in cos
    golden = np.pi * (1.0 + 5.0**0.5)
    theta = golden * i
    return np.stack(
        [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], axis=1
    )


def uniform_ball_points(n: int, rng: np.random.Generator, radius: float = 1.0) -> np.ndarray:
    """``n`` points uniform in the ball of the given radius."""
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    u = rng.random(n) ** (1.0 / 3.0)
    return d * (u * radius)[:, None]


def _odv_smooth(
    pts: np.ndarray,
    frozen: np.ndarray,
    n_iter: int = 8,
    relax: float = 0.7,
) -> np.ndarray:
    """Optimal-Delaunay-vertex smoothing of an interior point cloud.

    Every free point is moved towards the volume-weighted average of the circumcentres of
    its incident tetrahedra, which is the centroid of its Voronoi dual cell.  This is the
    cheapest way to remove the slivers that a raw Delaunay of random points always has --
    and slivers are exactly what makes the later tetrahedron maps fold.
    """
    from .mesh import circumcenters

    pts = pts.copy()
    free = np.setdiff1d(np.arange(len(pts)), frozen)
    for _ in range(n_iter):
        tri = Delaunay(pts)
        tets = tri.simplices.astype(np.int64)
        vols = np.abs(signed_tet_volumes(pts, tets))
        cc = circumcenters(pts, tets)
        acc = np.zeros((len(pts), 3))
        wsum = np.zeros(len(pts))
        for k in range(4):
            np.add.at(acc, tets[:, k], vols[:, None] * cc)
            np.add.at(wsum, tets[:, k], vols)
        tgt = acc / np.maximum(wsum, 1e-300)[:, None]
        old = pts[free].copy()
        pts[free] = old + relax * (tgt[free] - old)
    return pts


def ball_tet_mesh(
    n_surface: int = 162,
    n_interior: int = 400,
    seed: int = 0,
    radii: tuple[float, ...] = (0.35, 0.7),
    smooth: int = 8,
) -> TetMesh:
    """Tetrahedralise a polyhedral approximation of the unit ball.

    The mesh is the Delaunay tetrahedralisation of

    * ``n_surface`` points spread on the unit sphere, and
    * ``n_interior`` points on a few concentric shells,

    so the boundary of the mesh is the convex hull of the surface points and the solid is
    *exactly* a convex polyhedron -- a property the optimal-transport stage relies on.
    ``smooth`` optimal-Delaunay-vertex iterations improve the tetrahedron quality.
    """
    rng = np.random.default_rng(seed)
    pts = [fibonacci_sphere(n_surface)]
    per_shell = max(1, n_interior // max(1, len(radii)))
    for r in radii:
        pts.append(uniform_ball_points(per_shell, rng, radius=r))
    pts.append(np.zeros((1, 3)))
    P = np.vstack(pts)
    frozen = np.arange(n_surface)

    if smooth > 0:
        P = _odv_smooth(P, frozen, n_iter=smooth)
        # keep the boundary exactly on the unit sphere
        P[:n_surface] /= np.linalg.norm(P[:n_surface], axis=1, keepdims=True)

    tri = Delaunay(P)
    tets = tri.simplices.astype(np.int64)
    m = TetMesh(P, tets).fix_orientation()
    m = _drop_degenerate(m)
    return m


def _drop_degenerate(m: TetMesh, tol: float = 1e-12) -> TetMesh:
    v = signed_tet_volumes(m.vertices, m.tets)
    scale = np.abs(v).max()
    keep = np.abs(v) > tol * max(scale, 1.0)
    if keep.all():
        return m
    used = np.unique(m.tets[keep])
    remap = -np.ones(m.n_vertices, dtype=np.int64)
    remap[used] = np.arange(len(used))
    return TetMesh(m.vertices[used], remap[m.tets[keep]])


# --------------------------------------------------------------------------------------
# smooth deformations that make the solid *hard* (non-convex, anisotropic)
# --------------------------------------------------------------------------------------
def twist(mesh: TetMesh, angle: float = 1.2, axis: int = 2) -> TetMesh:
    """Twist the solid about ``axis`` by ``angle`` radians per unit length."""
    V = mesh.vertices.copy()
    others = [i for i in range(3) if i != axis]
    a, b = others
    t = V[:, axis]
    c, s = np.cos(angle * t), np.sin(angle * t)
    x, y = V[:, a].copy(), V[:, b].copy()
    V[:, a] = c * x - s * y
    V[:, b] = s * x + c * y
    return mesh.with_vertices(V)


def taper(mesh: TetMesh, factor: float = 0.45, axis: int = 2) -> TetMesh:
    """Linearly taper the cross-section along ``axis``."""
    V = mesh.vertices.copy()
    t = V[:, axis]
    lo, hi = t.min(), t.max()
    w = 1.0 + (factor - 1.0) * (t - lo) / max(hi - lo, 1e-12)
    for i in range(3):
        if i != axis:
            V[:, i] *= w
    return mesh.with_vertices(V)


def bump(mesh: TetMesh, amp: float = 0.35, k: int = 2, axis: int = 1) -> TetMesh:
    """Radially ripple the solid -- creates a genuinely non-convex, non-star-shaped body."""
    V = mesh.vertices.copy()
    r = np.linalg.norm(V, axis=1)
    theta = np.arctan2(V[:, 1], V[:, 0])
    scale = 1.0 + amp * np.cos(k * theta) * r**2
    V *= scale[:, None]
    if axis is not None and False:
        pass
    return mesh.with_vertices(V)


def shear_bend(mesh: TetMesh, amount: float = 0.5, axis: int = 0) -> TetMesh:
    """A soft S-shaped bend: ``x += amount * sin(pi * y)``."""
    V = mesh.vertices.copy()
    a = (axis + 1) % 3
    V[:, axis] += amount * np.sin(np.pi * V[:, a])
    return mesh.with_vertices(V)


def pear(mesh: TetMesh, amount: float = 0.5) -> TetMesh:
    """Pear-shaped solid: the classic stress test for ball parameterisation."""
    V = mesh.vertices.copy()
    z = V[:, 2]
    s = 1.0 + amount * z * (1.0 - z * z) * 1.2
    V[:, 0] *= s
    V[:, 1] *= s
    return mesh.with_vertices(V)


# --------------------------------------------------------------------------------------
# structured fallbacks
# --------------------------------------------------------------------------------------
def cube_tet_mesh(n: int = 4, extent: float = 1.0) -> TetMesh:
    """Regular ``n x n x n`` hexahedral grid, each hex split into 6 tetrahedra (Kuhn).

    The boundary is a *cube*, not a sphere: useful to test the pipeline on a shape that
    is convex but whose boundary is far from the ball.
    """
    lin = np.linspace(-extent, extent, n + 1)
    X, Y, Z = np.meshgrid(lin, lin, lin, indexing="ij")
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    def vid(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k

    # 6-tet Kuhn decomposition of a unit cube (consistent across cells -> conforming)
    kuhn = [
        [0, 1, 3, 7],
        [0, 1, 7, 5],
        [0, 2, 3, 7],
        [0, 2, 6, 7],
        [0, 4, 5, 7],
        [0, 4, 6, 7],
    ]
    tets = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                corner = [
                    vid(i, j, k),
                    vid(i + 1, j, k),
                    vid(i, j + 1, k),
                    vid(i + 1, j + 1, k),
                    vid(i, j, k + 1),
                    vid(i + 1, j, k + 1),
                    vid(i, j + 1, k + 1),
                    vid(i + 1, j + 1, k + 1),
                ]
                for t in kuhn:
                    tets.append([corner[x] for x in t])
    m = TetMesh(P, np.asarray(tets, dtype=np.int64)).fix_orientation()
    return _drop_degenerate(m)


def cubed_sphere(points: np.ndarray) -> np.ndarray:
    """The COBE / cubed-sphere diffeomorphism from the cube ``[-1,1]^3`` onto the unit ball.

        g_i(x) = x_i * sqrt(1 - (x_j^2 + x_k^2)/2 + x_j^2 x_k^2 / 3)

    It is smooth, bijective, fixes the coordinate axes, and maps the boundary of the cube
    onto the unit sphere -- so a structured cube grid becomes a structured ball mesh whose
    boundary is *convex* and whose tetrahedra keep a good shape.
    """
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    out = np.empty_like(points)
    out[:, 0] = x * np.sqrt(np.maximum(1.0 - (y * y + z * z) / 2.0 + y * y * z * z / 3.0, 0.0))
    out[:, 1] = y * np.sqrt(np.maximum(1.0 - (z * z + x * x) / 2.0 + z * z * x * x / 3.0, 0.0))
    out[:, 2] = z * np.sqrt(np.maximum(1.0 - (x * x + y * y) / 2.0 + x * x * y * y / 3.0, 0.0))
    return out


def ball_from_cube(mesh: TetMesh) -> TetMesh:
    """Map a cube-shaped tetrahedral mesh onto the unit ball by :func:`cubed_sphere`."""
    return mesh.with_vertices(cubed_sphere(mesh.vertices))


def ball_cube_mesh(n: int = 6) -> TetMesh:
    """Structured tetrahedral mesh of the unit ball (cubed-sphere grid).

    ``6 n^3`` tetrahedra, boundary convex, boundary vertices exactly on the unit sphere.
    """
    return ball_from_cube(cube_tet_mesh(n))


DEMOS = {
    "ball": lambda: ball_cube_mesh(6),
    "twisted": lambda: twist(ball_cube_mesh(6), 0.7),
    "pear": lambda: pear(ball_cube_mesh(6), 0.35),
    "bump": lambda: bump(ball_cube_mesh(6), 0.15, 2),
    "bent": lambda: shear_bend(ball_cube_mesh(6), 0.25),
    "taper": lambda: taper(ball_cube_mesh(6), 0.55),
    "squash": lambda: taper(ball_cube_mesh(6), 0.30),
}


def demo_mesh(name: str = "pear", **kwargs) -> TetMesh:
    if name == "cube":
        return cube_tet_mesh(kwargs.get("n", 6))
    if name in DEMOS:
        return DEMOS[name]()
    raise KeyError(f"unknown demo {name!r}; choose from {sorted(DEMOS) + ['cube']}")


__all__ = [
    "fibonacci_sphere",
    "uniform_ball_points",
    "ball_tet_mesh",
    "cube_tet_mesh",
    "ball_from_cube",
    "twist",
    "taper",
    "bump",
    "shear_bend",
    "pear",
    "demo_mesh",
    "DEMOS",
]
