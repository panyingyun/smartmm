"""Tetrahedralise a geometry (or load a mesh) and report ONLY its boundary topology.

    python tools/check_topology.py <geometry|mesh> [mesh-size]

Tells you whether the solid is a topological 3-ball, which the pipeline requires.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("run_pipeline", ROOT / "tests" / "run_pipeline.py")
rp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rp)  # type: ignore[union-attr]

path = Path(sys.argv[1]).resolve()
size = float(sys.argv[2]) if len(sys.argv) > 2 else None
m = rp.load_mesh(path, None, size, ROOT / "out" / "_topocheck")

t = m.boundary_topology()
print()
print(f"boundary components : {t['components']}")
print(f"genus               : {t['genus']:.0f}")
print(f"Euler characteristic: {t['euler']}   (sphere = 2)")
print(f"boundary V/E/F      : {t['n_boundary_vertices']} / "
      f"{t['n_boundary_edges']} / {t['n_boundary_faces']}")
print(f"is topological ball : {'YES' if t['is_ball'] else 'NO'}")
