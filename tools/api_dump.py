import inspect
import sys

sys.path.insert(0, r"E:\panyingyun\smartmm\src")

import tetparam
from tetparam import harmonic, laguerre, mesh, metrics, remesh, volume_ot

for mod in (mesh, harmonic, laguerre, volume_ot, remesh, metrics):
    print(f"\n===== {mod.__name__} =====")
    fns = [
        (n, o) for n, o in vars(mod).items()
        if (inspect.isfunction(o) or inspect.isclass(o))
        and not n.startswith("_")
        and getattr(o, "__module__", "") == mod.__name__
    ]
    for n, o in fns:
        kind = "class" if inspect.isclass(o) else "def  "
        doc = (inspect.getdoc(o) or "").split("\n")[0][:78]
        print(f"  {kind} {n:28s} {doc}")
