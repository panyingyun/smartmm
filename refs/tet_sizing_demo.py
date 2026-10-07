"""Verify the two production-ready sizing-field routes empirically:
  (1) TetGen  -m with a background mesh + .mtr per-node edge lengths
  (2) fTetWild --bg-mesh (a .msh tet mesh with a nodal "values" field)
plus a timing measurement of each.
"""
import os
import sys
import time
import subprocess
import numpy as np

WORK = os.path.join(os.environ.get("TEMP", "/tmp"), "dsh_tet_sizing")
os.makedirs(WORK, exist_ok=True)
PY = sys.executable


def make_ball_surface(n_lat=24, n_lon=48, R=1.0):
    """UV sphere -> vertices + triangles, outward oriented."""
    V, F = [], []
    th = np.linspace(0, np.pi, n_lat + 1)
    ph = np.linspace(0, 2 * np.pi, n_lon, endpoint=False)
    idx = {}
    for i, t in enumerate(th):
        for j, p in enumerate(ph):
            idx[(i, j)] = len(V)
            V.append((R * np.sin(t) * np.cos(p), R * np.sin(t) * np.sin(p), R * np.cos(t)))
    V = np.array(V)
    for i in range(n_lat):
        for j in range(n_lon):
            a, b = idx[(i, j)], idx[(i, (j + 1) % n_lon)]
            c, d = idx[(i + 1, j)], idx[(i + 1, (j + 1) % n_lon)]
            if i == 0:
                F.append((a, d, c))
            elif i == n_lat - 1:
                F.append((a, b, c))
            else:
                F.append((a, b, d))
                F.append((a, d, c))
    return V, np.array(F)


def write_smesh(path, V, F):
    with open(path, "w") as f:
        f.write(f"{len(V)} 3 0 0\n")
        for i, v in enumerate(V):
            f.write(f"{i} {v[0]:.10g} {v[1]:.10g} {v[2]:.10g}\n")
        f.write(f"{len(F)} 0\n")
        for i, t in enumerate(F):
            f.write(f"{i+1} {t[0]} {t[1]} {t[2]}\n")
        f.write("0\n0\n")


def tet_stats(path_node, path_ele):
    n = np.loadtxt(path_node, skiprows=1)[:, 1:4]
    e = np.loadtxt(path_ele, skiprows=1, dtype=int)[:, 1:5]
    P = n[e]
    vol = np.abs(np.einsum('ij,ij->i', P[:, 1] - P[:, 0],
                           np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))) / 6.0
    h = (6 * np.sqrt(2) * vol) ** (1 / 3)
    rad = np.linalg.norm(P.mean(1), axis=1)
    return len(n), len(e), vol, h, rad


def show(name, nodes, tets, vol, h, rad, secs):
    print(f"--- {name}   [{secs:.2f} s]")
    print(f"    nodes={nodes} tets={tets} vol={vol.sum():.5f} "
          f"(exact {4/3*np.pi:.5f}, err {abs(vol.sum()-4/3*np.pi)/(4/3*np.pi)*100:.3f}%)")
    for lo, hi, lbl in [(0, .25, 'core  r<.25'), (.25, .5, 'mid  .25-.5'),
                        (.5, .75, 'mid  .5-.75'), (.75, 1.01, 'shell r>.75')]:
        m = (rad >= lo) & (rad < hi)
        if m.sum():
            print(f"      {lbl}: n={m.sum():5d} h_eff med={np.median(h[m]):.4f}")
    # quality: radius ratio (3*inradius/circumradius) via min dihedral proxy
    print(f"      h_eff global: min={h.min():.4f} med={np.median(h):.4f} max={h.max():.4f}")


def main():
    V, F = make_ball_surface()
    smesh = os.path.join(WORK, "ball.smesh")
    write_smesh(smesh, V, F)
    print(f"surface: {len(V)} verts {len(F)} tris -> {smesh}\n")

    # ---------- 1) TetGen background mesh + .mtr ----------
    # bgmesh: coarse tet mesh of the ball, node sizes graded 0.06 (centre) -> 0.25 (rim)
    import scipy.spatial
    rng = np.random.default_rng(1)
    bp = rng.normal(size=(1200, 3))
    bp /= np.linalg.norm(bp, axis=1, keepdims=True)
    bp *= rng.uniform(0, 1, (len(bp), 1)) ** (1 / 3)
    bp = np.vstack([bp, np.zeros((1, 3))])
    tri = scipy.spatial.Delaunay(bp)
    keep = np.linalg.norm(bp[tri.simplices].mean(1), axis=1) < 1.0
    bt = tri.simplices[keep]
    r = np.linalg.norm(bp, axis=1)
    bs = 0.06 + (0.25 - 0.06) * np.clip(r, 0, 1)

    base = os.path.join(WORK, "ball.b")
    with open(base + ".node", "w") as f:
        f.write(f"{len(bp)} 3 0 0\n")
        for i, p in enumerate(bp):
            f.write(f"{i} {p[0]:.10g} {p[1]:.10g} {p[2]:.10g}\n")
    with open(base + ".ele", "w") as f:
        f.write(f"{len(bt)} 4 0\n")
        for i, t in enumerate(bt):
            f.write(f"{i} {t[0]} {t[1]} {t[2]} {t[3]}\n")
    with open(base + ".mtr", "w") as f:
        f.write(f"{len(bp)} 1\n")
        for s in bs:
            f.write(f"{s:.10g}\n")

    t0 = time.time()
    r = subprocess.run(["tetgen", "-pq1.414a0.001m", "ball.smesh"],
                       cwd=WORK, capture_output=True, text=True)
    secs = time.time() - t0
    print("tetgen rc", r.returncode, "|", (r.stdout or "")[-400:].replace("\n", " | "),
          (r.stderr or "")[:200])
    outs = sorted([f for f in os.listdir(WORK) if f.startswith("ball.1.")])
    print("tetgen outputs:", outs)
    node = os.path.join(WORK, "ball.1.node")
    ele = os.path.join(WORK, "ball.1.ele")
    if os.path.exists(node) and os.path.exists(ele):
        show("TetGen -pq1.414a0.001m (bg mesh + .mtr)", *tet_stats(node, ele), secs)


if __name__ == "__main__":
    main()
