"""Distortion metrics and acceptance tests for volumetric maps.

All metrics are defined so that a *volume-preserving* map has

    ratio == 1 everywhere,  jacobian > 0 everywhere,  E_V == V(f) == lower bound.

Definitions
-----------
For a piecewise-linear map ``f`` on a tetrahedral mesh ``M`` and a tetrahedron ``tau``:

* signed volume ratio      ``r_tau = |f(tau)| / |tau|``
* volume distortion        ``D_tau = max(r_tau, 1/r_tau) - 1``        (0 iff preserved)
* volumetric stretch energy ``E_V(f) = sum_tau |f(tau)|^2 / |tau|``
  with the theorem (Yueh et al.; Huang et al.) ::

      E_V(f) >= V(f)^2 / V(M),      equality  <=>  f is volume-preserving,

  where ``V(f) = sum_tau |f(tau)|``.  Normalising the image gives the standard form
  ``E_V(f) >= V(M)`` for ``V(f) = V(M)``.
* isovolumetric energy     ``E_I(f) = V(e)/V(f) * E_V(f) - V(f)`` (arXiv:2407.19272)
* Jacobian determinant     ``J_tau = |f(tau)| / |tau|`` (signed), flips detected by ``J <= 0``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .mesh import TetMesh, tet_shape_gradients


def volume_ratios(src: TetMesh, dst: TetMesh) -> np.ndarray:
    """Signed ``|f(tau)| / |tau|`` for the map ``src -> dst`` (same connectivity)."""
    v0 = src.signed_volumes()
    v1 = dst.signed_volumes()
    return np.where(np.abs(v0) > 1e-300, v1 / v0, np.nan)


def jacobian_determinants(src: TetMesh, dst: TetMesh) -> np.ndarray:
    """Per-tetrahedron Jacobian determinant of the affine map ``tau -> f(tau)``.

    ``J_tau = det( D f )`` = signed volume ratio.  This is the field plotted as the
    standard "volume distortion" colour map.
    """
    return volume_ratios(src, dst)


def stretch_energy(src: TetMesh, dst: TetMesh) -> float:
    """Volumetric stretch energy ``sum_tau |f(tau)|^2 / |tau|``."""
    v0 = np.abs(src.signed_volumes())
    v1 = dst.signed_volumes()
    return float(np.sum(np.where(v0 > 1e-300, v1**2 / v0, 0.0)))


def isovolumetric_energy(src: TetMesh, dst: TetMesh) -> float:
    V_e = src.total_volume()
    V_f = dst.total_volume()
    if V_f <= 0:
        return float("inf")
    return float(V_e / V_f * stretch_energy(src, dst) - V_f)


def stretch_energy_lower_bound(src: TetMesh, dst: TetMesh) -> float:
    """``V(f)^2 / V(M)`` -- the theoretical minimum of :func:`stretch_energy`."""
    V_f = dst.total_volume()
    return float(V_f**2 / src.total_volume())


def barycentric_gradient_field(src: TetMesh, dst: TetMesh) -> np.ndarray:
    """Per-tet Jacobian matrices ``D f`` of the affine map, shape ``(m, 3, 3)``.

    ``Df = sum_k f(v_k) (grad lambda_k)^T``, so ``Df[a, b] = sum_k f_k[a] grad_k[b]``.
    """
    g0, _ = tet_shape_gradients(src.vertices, src.tets)
    fv = dst.vertices[src.tets]  # (m,4,3)
    return np.einsum("mka,mkb->mab", fv, g0)


@dataclass
class DistortionReport:
    ratio_min: float
    ratio_max: float
    ratio_mean: float
    ratio_std: float
    distortion_mean: float
    distortion_max: float
    distortion_p99: float
    n_flipped: int
    total_volume_ratio: float
    stretch_energy: float
    stretch_lower_bound: float
    stretch_excess: float
    isovolumetric_energy: float
    max_singular_value: float
    min_singular_value: float

    def as_dict(self) -> dict:
        return dict(self.__dict__)

    def summary(self) -> str:
        return (
            "DistortionReport(\n"
            f"  volume ratio  min {self.ratio_min:.6f}  max {self.ratio_max:.6f}  "
            f"mean {self.ratio_mean:.6f}  std {self.ratio_std:.3e}\n"
            f"  distortion    mean {self.distortion_mean:.4f}  p99 {self.distortion_p99:.4f}  "
            f"max {self.distortion_max:.4f}   (0 = perfect)\n"
            f"  flipped tets  {self.n_flipped}\n"
            f"  total volume ratio {self.total_volume_ratio:.6f}\n"
            f"  E_V {self.stretch_energy:.8g}   lower bound {self.stretch_lower_bound:.8g}   "
            f"excess {self.stretch_excess:.4%}\n"
            f"  E_I {self.isovolumetric_energy:.6g}\n"
            f"  singular values of Df: [{self.min_singular_value:.4g}, {self.max_singular_value:.4g}]\n"
            ")"
        )


def distortion_report(src: TetMesh, dst: TetMesh) -> DistortionReport:
    r = volume_ratios(src, dst)
    r = r[np.isfinite(r)]
    r = r[np.abs(r) > 0]
    d = np.maximum(np.abs(r), 1.0 / np.abs(r)) - 1.0
    Df = barycentric_gradient_field(src, dst)
    sv = np.linalg.svd(Df, compute_uv=False)
    E = stretch_energy(src, dst)
    lb = stretch_energy_lower_bound(src, dst)
    # "flipped" counts *image* tetrahedra with non-positive orientation.  Note that the
    # ratio-based test can report a "negative ratio" only when the *source* mesh itself has
    # an inverted tetrahedron (a bad input), which must be reported separately.
    n_flip = int((dst.signed_volumes() <= 0).sum())
    return DistortionReport(
        ratio_min=float(r.min()),
        ratio_max=float(r.max()),
        ratio_mean=float(r.mean()),
        ratio_std=float(r.std()),
        distortion_mean=float(d.mean()),
        distortion_max=float(d.max()),
        distortion_p99=float(np.percentile(d, 99)),
        n_flipped=n_flip,
        total_volume_ratio=float(dst.total_volume() / src.total_volume()),
        stretch_energy=E,
        stretch_lower_bound=lb,
        stretch_excess=float(E / lb - 1.0) if lb > 0 else float("inf"),
        isovolumetric_energy=isovolumetric_energy(src, dst),
        max_singular_value=float(sv.max()),
        min_singular_value=float(sv.min()),
    )


def map_is_injective_heuristic(mesh: TetMesh) -> bool:
    """Cheap necessary condition for injectivity: no negative tetrahedron."""
    return bool((mesh.signed_volumes() > 0).all())


__all__ = [
    "volume_ratios",
    "jacobian_determinants",
    "stretch_energy",
    "isovolumetric_energy",
    "stretch_energy_lower_bound",
    "barycentric_gradient_field",
    "distortion_report",
    "DistortionReport",
    "map_is_injective_heuristic",
]
