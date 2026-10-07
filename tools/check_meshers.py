import importlib, sys

for m in ("gmsh", "tetgen", "meshio", "pyvista", "pytetwild", "mmgpy", "trimesh", "numpy", "scipy"):
    try:
        mod = importlib.import_module(m)
        v = getattr(mod, "__version__", "?")
        print(f"{m:12s} OK   {v}")
    except Exception as e:
        print(f"{m:12s} --   {type(e).__name__}")
print("python:", sys.executable)
