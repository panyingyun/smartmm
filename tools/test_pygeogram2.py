import numpy as np
import geogram

geogram.initialize()
rng = np.random.default_rng(0)
P = rng.random((10, 3))
W = np.zeros(10)
verts = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                  [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]])
tets = np.array([[0, 1, 2, 5], [0, 2, 3, 6], [0, 1, 5, 4], [0, 4, 5, 7],
                 [0, 5, 6, 7], [0, 2, 6, 7], [0, 3, 6, 7], [0, 1, 4, 7],
                 [0, 1, 7, 5]], dtype=np.uint32)
v = geogram.Voronoi(P, W, verts, tets)
q = np.asarray(v.q); t = np.asarray(v.t); ts = np.asarray(v.tseed)
print("q rows", q.shape, "t rows", t.shape)
print("do q[0..9] equal the sites P?",
      np.allclose(q[:10], P))
print("q[10] =", q[10])
# Are the 4 vertices of each 'tet' coplanar-equidistant -> check power equality
i = 0
a, b, c, d = (q[t[i][k]] for k in range(4))
print("vol6 of t[0]:", np.dot(np.cross(b - a, c - a), d - a))
# try: cells = convex hull of {q[j] : j in t[ts==i]} clipped to domain
from scipy.spatial import ConvexHull
for i in (0, 6):
    idx = np.unique(t[ts == i].ravel())
    pts = q[idx]
    print(f"cell {i}: {len(idx)} dual vertices, hull?", end=" ")
    try:
        h = ConvexHull(pts)
        print("volume=", h.volume)
    except Exception as e:
        print("degenerate", e)
print("total volume via hull of all cell dual verts:",
      sum(ConvexHull(q[np.unique(t[ts == i].ravel())]).volume
          for i in range(10) if len(np.unique(t[ts == i].ravel())) >= 4))
