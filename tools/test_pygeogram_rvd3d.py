import numpy as np
import geogram

geogram.initialize()

# 10 random sites in the unit cube, zero weights -> ordinary Voronoi
rng = np.random.default_rng(0)
P = rng.random((10, 3))
W = np.zeros(10)

# domain = unit cube as a tet mesh (6 tets around the main diagonal)
verts = np.array([
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
    [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1],
], dtype=np.float64)
tets = np.array([
    [0, 1, 2, 3], [1, 2, 3, 7], [1, 2, 4, 7],
    [1, 3, 5, 7], [2, 3, 6, 7], [1, 2, 3, 7],
], dtype=np.uint32)
# use a clean 5-tet decomposition instead
tets = np.array([
    [0, 1, 2, 5], [0, 2, 3, 6], [0, 1, 5, 4],
    [0, 4, 5, 7], [0, 5, 6, 7], [0, 2, 6, 7],
    [0, 3, 6, 7], [0, 1, 4, 7], [0, 1, 7, 5],
], dtype=np.uint32)

v = geogram.Voronoi(P, W, verts, tets)
print("dimension:", v.dimension)
print("seeds shape:", np.asarray(v.seeds).shape)
print("t shape:", np.asarray(v.t).shape, " tadj:", np.asarray(v.tadj).shape,
      " tseed:", np.asarray(v.tseed).shape)
print("q (dual vertices) shape:", np.asarray(v.q).shape)
print("first few t:", np.asarray(v.t)[:4])
print("first few tseed:", np.asarray(v.tseed)[:8])

# --- volume of each Laguerre cell = sum over tets t with tseed t and whose
# 4 dual vertices are all finite, of |det|/6
q = np.asarray(v.q)
t = np.asarray(v.t)
tseed = np.asarray(v.tseed)
print("max tseed:", tseed.max(), "nb sites:", len(P))
vol = np.zeros(len(P))
for i, tet in enumerate(t):
    # tet: 4 global indices into q (dual vertices)
    if (tet < 0).any():
        continue
    a, b, c, d = q[tet[0]], q[tet[1]], q[tet[2]], q[tet[3]]
    vol6 = abs(np.dot(np.cross(b - a, c - a), d - a))
    vol[tseed[i]] += vol6 / 6.0
print("total volume:", vol.sum(), "(domain volume = 1.0)")
print("per-cell volumes:", np.round(vol, 5))
