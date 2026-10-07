"""tetparam -- conformal + optimal-transport volume parameterisation of tetrahedral meshes.

Modules
-------
mesh        tetrahedral mesh container, TetGen/VTK I/O, quality metrics
generate    synthetic ball/cube tetrahedral meshes and smooth deformations
harmonic    stage 1: P1 finite-element harmonic (volumetric "conformal") map to the ball
laguerre    stage 2: discrete optimal mass transport via Laguerre / power diagrams
volume_ot   stage 3: composition into a volume-preserving map + tetrahedral remeshing
metrics     distortion / volume-preservation / injectivity measurements
"""
from __future__ import annotations

__version__ = "0.1.0"

from .mesh import TetMesh  # noqa: F401
from .laguerre import ConvexPolyhedron, laguerre_cells, solve_ot  # noqa: F401
from .harmonic import fem_stiffness, harmonic_ball_map  # noqa: F401
from .metrics import distortion_report  # noqa: F401

__all__ = [
    "TetMesh",
    "ConvexPolyhedron",
    "laguerre_cells",
    "solve_ot",
    "fem_stiffness",
    "harmonic_ball_map",
    "distortion_report",
]
