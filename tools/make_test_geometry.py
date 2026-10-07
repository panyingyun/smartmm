"""Build a small test geometry (STEP + STL) so the meshing commands can be verified."""
from __future__ import annotations

import sys
from pathlib import Path

OUT = Path(r"E:\panyingyun\smartmm\out\geomtest")
OUT.mkdir(parents=True, exist_ok=True)

import gmsh  # noqa: E402

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("demo")
# a box minus a sphere: non-trivial geometry with a curved internal cavity
gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
gmsh.model.occ.addSphere(0.5, 0.5, 0.5, 0.32, 2)
gmsh.model.occ.cut([(3, 1)], [(3, 2)])
gmsh.model.occ.synchronize()

step = OUT / "demo.step"
gmsh.write(str(step))
gmsh.finalize()
print("wrote", step, step.stat().st_size, "bytes")

# a plain triangulated surface for the TetGen route
import meshio  # noqa: E402
import numpy as np  # noqa: E402

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("sphere")
gmsh.model.occ.addSphere(0, 0, 0, 1.0)
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.25)
gmsh.model.mesh.generate(2)
pts, tris = gmsh.model.mesh.getNodes()[1].reshape(-1, 3), None
tris = gmsh.model.mesh.getElements(2)[2][0].reshape(-1, 3) - 1
gmsh.finalize()
stl = OUT / "sphere.stl"
meshio.write(str(stl), meshio.Mesh(pts, [("triangle", tris)]))
print("wrote", stl, stl.stat().st_size, "bytes,", len(pts), "vertices", len(tris), "triangles")
del sys, np
