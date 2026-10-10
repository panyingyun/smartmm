"""Inspect a CAD geometry without meshing it: bounding box, solids, shells, faces.

Usage:
    python tools/inspect_geometry.py <model.step|.stp|.iges|.brep>
"""
from __future__ import annotations

import sys
from pathlib import Path

import gmsh

path = Path(sys.argv[1] if len(sys.argv) > 1 else "").resolve()
if not path.exists():
    raise SystemExit(f"not found: {path}")

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(str(path))
gmsh.model.occ.synchronize()

xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
print(f"file      : {path.name}  ({path.stat().st_size / 1e6:.1f} MB)")
print(f"bbox      : x [{xmin:.4g}, {xmax:.4g}]  y [{ymin:.4g}, {ymax:.4g}]  z [{zmin:.4g}, {zmax:.4g}]")
print(f"extent    : {dx:.4g} x {dy:.4g} x {dz:.4g}   diagonal {((dx**2 + dy**2 + dz**2) ** 0.5):.4g}")
print(f"entities  : volumes {len(gmsh.model.getEntities(3))}  "
      f"surfaces {len(gmsh.model.getEntities(2))}  "
      f"curves {len(gmsh.model.getEntities(1))}  "
      f"points {len(gmsh.model.getEntities(0))}")

diag = (dx**2 + dy**2 + dz**2) ** 0.5
for frac in (0.02, 0.05, 0.1):
    h = diag * frac
    print(f"  --mesh-size {h:<10.4g} ~ {frac * 100:.0f}% of diagonal")
print(f"\n建议：--mesh-size 取对角线的 2%~5%，即 {diag * 0.02:.4g} ~ {diag * 0.05:.4g}")

gmsh.finalize()
