"""Report tetrahedron quality for an existing mesh (VTK / Medit / TetGen).

    python tools/mesh_quality.py <mesh.vtk|mesh.msh|nodes.node --ele t.ele> [label]
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("mg", ROOT / "tools" / "mesh_geometry.py")
mg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mg)  # type: ignore[union-attr]

from tetparam.mesh import TetMesh  # noqa: E402

path = Path(sys.argv[1]).resolve()
label = sys.argv[2] if len(sys.argv) > 2 else path.name

if path.suffix.lower() == ".node":
    ele = path.with_suffix(".ele")
    m = TetMesh.from_tetgen(path, ele)
else:
    import meshio

    md = meshio.read(str(path))
    tets = md.cells_dict.get("tetra")
    if tets is None:
        raise SystemExit(f"{path.name}: no tetrahedra")
    m = TetMesh(np.asarray(md.points, float), np.asarray(tets, np.int64))

m.fix_orientation()
mg.quality_report(m, label)
