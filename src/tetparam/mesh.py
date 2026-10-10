"""Tetrahedral mesh core: storage, I/O, geometry, quality.

Deliberately dependency-light: only numpy + (optionally) scipy and meshio.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

# --------------------------------------------------------------------------------------
# low level geometry helpers
# --------------------------------------------------------------------------------------

#: the four (local) faces of a tetrahedron, each oriented so that its normal points
#: *outward* for a positively oriented tet ``[v0, v1, v2, v3]`` (i.e. ``det > 0``).
TET_FACES = np.array(
    [
        [0, 2, 1],  # opposite v3
        [0, 1, 3],  # opposite v2
        [0, 2, 3],  # opposite v1
        [1, 2, 3],  # opposite v0
    ],
    dtype=np.int64,
)

#: the six edges of a tetrahedron.
TET_EDGES = np.array(
    [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]], dtype=np.int64
)


def signed_tet_volumes(verts: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Signed volume of every tetrahedron (positive for consistent orientation)."""
    p0 = verts[tets[:, 0]]
    p1 = verts[tets[:, 1]]
    p2 = verts[tets[:, 2]]
    p3 = verts[tets[:, 3]]
    return np.einsum("ij,ij->i", np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0


def tet_shape_gradients(verts: np.ndarray, tets: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-tet gradients of the four barycentric coordinates and signed volumes.

    Returns
    -------
    grads : (m, 4, 3) float64
        ``grads[t, k]`` is ``grad lambda_k`` on tetrahedron ``t``.
    vols : (m,) float64
        signed volumes.
    """
    p0 = verts[tets[:, 0]]
    p1 = verts[tets[:, 1]]
    p2 = verts[tets[:, 2]]
    p3 = verts[tets[:, 3]]
    # edge matrix E = [p1-p0 | p2-p0 | p3-p0], columns
    E = np.stack([p1 - p0, p2 - p0, p3 - p0], axis=2)  # (m,3,3)
    det = np.linalg.det(E)
    inv = np.linalg.inv(E)  # rows of inv^T are grad lambda_1..3 for the apexes
    # lambda_1..3 gradients are the columns of inv^T == rows of inv
    g1, g2, g3 = inv[:, 0, :], inv[:, 1, :], inv[:, 2, :]
    g0 = -(g1 + g2 + g3)
    grads = np.stack([g0, g1, g2, g3], axis=1)
    return grads, det / 6.0


def circumcenters(verts: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Circumcentres of tetrahedra (may be far outside for slivers)."""
    p0 = verts[tets[:, 0]]
    p1 = verts[tets[:, 1]]
    p2 = verts[tets[:, 2]]
    p3 = verts[tets[:, 3]]
    E = np.stack([p1 - p0, p2 - p0, p3 - p0], axis=2)
    rhs = 0.5 * np.stack(
        [
            np.einsum("ij,ij->i", p1 - p0, p1 - p0),
            np.einsum("ij,ij->i", p2 - p0, p2 - p0),
            np.einsum("ij,ij->i", p3 - p0, p3 - p0),
        ],
        axis=1,
    )
    try:
        sol = np.linalg.solve(E, rhs[..., None])[..., 0]
    except np.linalg.LinAlgError:  # pragma: no cover - sliver fallback
        sol = np.einsum("mij,mj->mi", np.linalg.pinv(E), rhs)
    return p0 + sol


# --------------------------------------------------------------------------------------
# mesh container
# --------------------------------------------------------------------------------------


@dataclass
class TetMesh:
    """A tetrahedral mesh of a solid region.

    Attributes
    ----------
    vertices : (n, 3) float64
    tets : (m, 4) int64, positively oriented (signed volume > 0)
    """

    vertices: np.ndarray
    tets: np.ndarray
    _boundary: dict | None = field(default=None, repr=False, compare=False)

    # ---------------------------------------------------------------- construction
    def __post_init__(self) -> None:
        self.vertices = np.ascontiguousarray(self.vertices, dtype=np.float64)
        self.tets = np.ascontiguousarray(self.tets, dtype=np.int64)
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 3:
            raise ValueError("vertices must have shape (n, 3)")
        if self.tets.ndim != 2 or self.tets.shape[1] != 4:
            raise ValueError("tets must have shape (m, 4)")

    @property
    def n_vertices(self) -> int:
        return int(self.vertices.shape[0])

    @property
    def n_tets(self) -> int:
        return int(self.tets.shape[0])

    # ---------------------------------------------------------------- geometry
    def signed_volumes(self) -> np.ndarray:
        return signed_tet_volumes(self.vertices, self.tets)

    def volumes(self) -> np.ndarray:
        return np.abs(self.signed_volumes())

    def total_volume(self) -> float:
        return float(self.volumes().sum())

    def centroids(self) -> np.ndarray:
        return self.vertices[self.tets].mean(axis=1)

    def quality(self) -> dict[str, np.ndarray | float]:
        """Sliver / shape quality indicators for every tetrahedron.

        * ``radius_ratio``: circumradius / (6 * inradius); 3 for a regular tet, larger =
          worse.  This is the classic "radius ratio" of mesh generation.
        * ``min_dihedral`` / ``max_dihedral``: degrees.
        * ``aspect``: longest edge / (sqrt(2) * ... ) -- we report
          ``max_edge / (6*sqrt(2)*inradius)`` restricted to well-shaped tets.
        """
        V = np.abs(self.signed_volumes())
        cc = circumcenters(self.vertices, self.tets)
        R = np.linalg.norm(cc - self.vertices[self.tets[:, 0]], axis=1)
        # inradius r = 3V / S with S the surface area of the tet
        S = self._tet_surface_areas()
        r = np.where(S > 0, 3.0 * V / np.maximum(S, 1e-300), 0.0)
        radius_ratio = np.where(r > 0, R / np.maximum(r, 1e-300), np.inf)

        dih = self._dihedral_angles()
        P = self.vertices[self.tets]
        edges = np.linalg.norm(P[:, TET_EDGES[:, 0]] - P[:, TET_EDGES[:, 1]], axis=2)
        max_edge = edges.max(axis=1)
        return {
            "volume": V,
            "radius_ratio": radius_ratio,
            "min_dihedral": np.degrees(dih.min(axis=1)),
            "max_dihedral": np.degrees(dih.max(axis=1)),
            "max_edge": max_edge,
            "min_dihedral_mean": float(np.degrees(dih.min(axis=1)).mean()),
            "radius_ratio_mean": float(radius_ratio.mean()),
            "radius_ratio_max": float(radius_ratio.max()),
            "n_slivers": int((np.degrees(dih.min(axis=1)) < 5.0).sum()),
        }

    def _tet_surface_areas(self) -> np.ndarray:
        P = self.vertices[self.tets]
        F = P[:, TET_FACES]  # (m, 4, 3, 3)
        a = np.linalg.norm(np.cross(F[:, :, 1] - F[:, :, 0], F[:, :, 2] - F[:, :, 0]), axis=2)
        return 0.5 * a.sum(axis=1)

    def _dihedral_angles(self) -> np.ndarray:
        """All six dihedral angles of every tet, in radians."""
        P = self.vertices[self.tets]
        F = P[:, TET_FACES]  # (m,4,3,3)
        n = np.cross(F[:, :, 1] - F[:, :, 0], F[:, :, 2] - F[:, :, 0])
        n = n / np.maximum(np.linalg.norm(n, axis=2, keepdims=True), 1e-300)
        # face pairs sharing an edge: (0,1) share edge 0-1, (0,2)->0-2, (0,3)->0-3,
        # (1,2)->1-2? use local table
        pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        out = np.empty((self.n_tets, 6))
        for k, (a, b) in enumerate(pairs):
            c = np.einsum("ij,ij->i", n[:, a], n[:, b])
            out[:, k] = np.arccos(np.clip(c, -1.0, 1.0))
        return out

    # ---------------------------------------------------------------- topology
    def boundary(self) -> dict:
        """Boundary triangles and vertex classification.

        Returns a dict with
        ``faces`` (k,3) int, ``face_owners`` (k,) int, ``is_boundary_vertex`` (n,) bool,
        ``boundary_vertices`` (b,) int, ``interior_vertices`` (n-b,) int, ``edges`` (e,2).
        """
        if self._boundary is not None:
            return self._boundary

        all_faces = self.tets[:, TET_FACES].reshape(-1, 3)  # (4m, 3)
        owners = np.repeat(np.arange(self.n_tets), 4)
        key = np.sort(all_faces, axis=1)
        order = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
        ks = key[order]
        same = np.all(ks[1:] == ks[:-1], axis=1)
        # faces that appear only once are on the boundary
        is_dup = np.zeros(len(ks), dtype=bool)
        is_dup[1:] |= same
        is_dup[:-1] |= same
        bnd = order[~is_dup]
        faces = all_faces[bnd]
        face_owners = owners[bnd]

        # orient boundary faces outward w.r.t. their owner tet
        tri = self.vertices[faces]
        nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        # centroid of the owner tet relative to the face centroid
        tet_c = self.vertices[self.tets[face_owners]].mean(axis=1)
        face_c = tri.mean(axis=1)
        flip = np.einsum("ij,ij->i", nrm, face_c - tet_c) < 0
        faces[flip] = faces[flip][:, [0, 2, 1]]
        del nrm

        is_bnd_v = np.zeros(self.n_vertices, dtype=bool)
        is_bnd_v[faces.ravel()] = True
        bverts = np.flatnonzero(is_bnd_v)
        iverts = np.flatnonzero(~is_bnd_v)

        edges = np.vstack(
            [
                np.sort(self.tets[:, [0, 1]], axis=1),
                np.sort(self.tets[:, [0, 2]], axis=1),
                np.sort(self.tets[:, [0, 3]], axis=1),
                np.sort(self.tets[:, [1, 2]], axis=1),
                np.sort(self.tets[:, [1, 3]], axis=1),
                np.sort(self.tets[:, [2, 3]], axis=1),
            ]
        )
        edges = np.unique(edges, axis=0)

        self._boundary = {
            "faces": faces,
            "face_owners": face_owners,
            "is_boundary_vertex": is_bnd_v,
            "boundary_vertices": bverts,
            "interior_vertices": iverts,
            "edges": edges,
        }
        return self._boundary

    @property
    def boundary_faces(self) -> np.ndarray:
        return self.boundary()["faces"]

    @property
    def is_boundary_vertex(self) -> np.ndarray:
        return self.boundary()["is_boundary_vertex"]

    def vertex_incidence(self):
        """Return ``tets_of_vertex`` (list of arrays) and ``vertex_tet_counts``."""
        counts = np.bincount(self.tets.ravel(), minlength=self.n_vertices)
        order = np.argsort(self.tets.ravel(), kind="stable")
        flat = self.tets.ravel()[order]
        starts = np.searchsorted(flat, np.arange(self.n_vertices))
        ends = np.searchsorted(flat, np.arange(self.n_vertices), side="right")
        owners = np.repeat(np.arange(self.n_tets), 4)[order]
        return [owners[s:e] for s, e in zip(starts, ends)], counts

    def dual_vertex_volumes(self) -> np.ndarray:
        """Barycentric dual cell volume of each vertex: sum of |tet| / 4."""
        V = self.volumes() / 4.0
        out = np.zeros(self.n_vertices)
        for k in range(4):
            np.add.at(out, self.tets[:, k], V)
        return out

    def check_orientation(self) -> int:
        """Number of negatively oriented tetrahedra (should be 0)."""
        return int((self.signed_volumes() <= 0).sum())

    def boundary_topology(self) -> dict:
        """Euler characteristic and genus of the boundary surface.

        The pipeline needs the solid to be a **topological 3-ball**: the boundary must be a
        single closed surface of genus 0 (a sphere).  A wheel with a central bore or with
        lightening holes has genus >= 1 and cannot be mapped to a ball as it is -- it must
        first be cut into simply-connected pieces, or handled as a solid torus / with the
        genus-one variant of the algorithm.

        Returns ``{'components', 'genus', 'euler', 'n_boundary_vertices',
        'n_boundary_edges', 'n_boundary_faces', 'is_ball'}``.
        """
        import scipy.sparse as sp
        from scipy.sparse.csgraph import connected_components

        faces = self.boundary_faces
        if len(faces) == 0:
            return {
                "components": 0, "genus": 0.0, "euler": 0, "n_boundary_vertices": 0,
                "n_boundary_edges": 0, "n_boundary_faces": 0, "is_ball": False,
            }
        edges = np.unique(
            np.vstack(
                [
                    np.sort(faces[:, [0, 1]], axis=1),
                    np.sort(faces[:, [1, 2]], axis=1),
                    np.sort(faces[:, [2, 0]], axis=1),
                ]
            ),
            axis=0,
        )
        n = self.n_vertices
        rows = np.concatenate([edges[:, 0], edges[:, 1]])
        cols = np.concatenate([edges[:, 1], edges[:, 0]])
        A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)).tocsr()
        _, labels = connected_components(A, directed=False)
        used = np.unique(faces)
        components = int(len(np.unique(labels[used])))
        V, E, F = int(len(used)), int(len(edges)), int(len(faces))
        euler = V - E + F
        genus = (2 * components - euler) / 2.0
        return {
            "components": components,
            "genus": float(genus),
            "euler": int(euler),
            "n_boundary_vertices": V,
            "n_boundary_edges": E,
            "n_boundary_faces": F,
            "is_ball": bool(components == 1 and abs(genus) < 1e-9),
        }

    def fix_orientation(self) -> "TetMesh":
        s = self.signed_volumes()
        bad = s < 0
        if bad.any():
            t = self.tets[bad].copy()
            t[:, [1, 2]] = t[:, [2, 1]]
            self.tets[bad] = t
        return self

    def copy(self) -> "TetMesh":
        return TetMesh(self.vertices.copy(), self.tets.copy())

    def with_vertices(self, vertices: np.ndarray) -> "TetMesh":
        return TetMesh(np.asarray(vertices, dtype=np.float64), self.tets.copy())

    # ---------------------------------------------------------------- I/O
    @classmethod
    def from_tetgen(cls, node: str | Path, ele: str | Path) -> "TetMesh":
        node = Path(node)
        ele = Path(ele)
        with node.open("r", encoding="utf-8", errors="replace") as f:
            header = f.readline().split()
            n, dim = int(header[0]), int(header[1])
            n_attr = int(header[2]) if len(header) > 2 else 0
            n_bnd = int(header[3]) if len(header) > 3 else 0
            verts = np.empty((n, 3))
            for i in range(n):
                parts = f.readline().split()
                verts[i] = [float(x) for x in parts[1 : 1 + min(dim, 3)]]
                if dim < 3:
                    verts[i, dim:] = 0.0
            del n_attr, n_bnd
        with ele.open("r", encoding="utf-8", errors="replace") as f:
            header = f.readline().split()
            m, k = int(header[0]), int(header[1])
            tets = np.empty((m, 4), dtype=np.int64)
            for i in range(m):
                parts = f.readline().split()
                tets[i] = [int(x) - 1 for x in parts[1:5]]
            assert k >= 4
        return cls(verts, tets)

    def write_tetgen(self, node: str | Path, ele: str | Path) -> None:
        node = Path(node)
        ele = Path(ele)
        with node.open("w", encoding="utf-8") as f:
            f.write(f"{self.n_vertices} 3 0 0\n")
            for i, v in enumerate(self.vertices):
                f.write(f"{i + 1} {v[0]:.17g} {v[1]:.17g} {v[2]:.17g}\n")
        with ele.open("w", encoding="utf-8") as f:
            f.write(f"{self.n_tets} 4 0\n")
            for i, t in enumerate(self.tets):
                f.write(f"{i + 1} {t[0] + 1} {t[1] + 1} {t[2] + 1} {t[3] + 1}\n")

    @classmethod
    def from_vtk(cls, path: str | Path) -> "TetMesh":
        import meshio

        m = meshio.read(str(path))
        tets = m.cells_dict.get("tetra")
        if tets is None:
            raise ValueError(f"no tetra cells in {path}")
        return cls(np.asarray(m.points, dtype=np.float64), np.asarray(tets, dtype=np.int64))

    def write_vtk(self, path: str | Path, point_data: dict | None = None,
                  cell_data: dict | None = None) -> None:
        import meshio

        cells = [("tetra", self.tets)]
        m = meshio.Mesh(self.vertices, cells, point_data=point_data, cell_data=cell_data)
        m.write(str(path))

    # ---------------------------------------------------------------- reporting
    def summary(self) -> str:
        q = self.quality()
        return (
            f"TetMesh: {self.n_vertices} vertices, {self.n_tets} tets, "
            f"boundary verts {int(self.is_boundary_vertex.sum())}, "
            f"volume {self.total_volume():.6g}, "
            f"radius-ratio mean {q['radius_ratio_mean']:.3f} max {q['radius_ratio_max']:.3f}, "
            f"min-dihedral mean {q['min_dihedral_mean']:.2f} deg, slivers {q['n_slivers']}"
        )


__all__ = [
    "TetMesh",
    "TET_FACES",
    "TET_EDGES",
    "signed_tet_volumes",
    "tet_shape_gradients",
    "circumcenters",
]
