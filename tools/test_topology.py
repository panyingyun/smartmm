import sys

sys.path.insert(0, r"E:\panyingyun\smartmm\src")
from tetparam.generate import ball_cube_mesh, cube_tet_mesh  # noqa: E402

for nm, m in (("ball(n=6)", ball_cube_mesh(6)), ("cube(n=4)", cube_tet_mesh(4))):
    t = m.boundary_topology()
    print(
        f"{nm:12s} chi={t['euler']:5d}  genus={t['genus']:.0f}  "
        f"components={t['components']}  is_ball={t['is_ball']}"
    )
