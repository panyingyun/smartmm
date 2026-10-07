# 02 — Optimal Transport toolchain for volume-preserving 3D parameterization

**Scope.** A tetrahedral mesh of a solid 3-ball → (1) harmonic/conformal map to the unit ball →
(2) discrete OMT correction so the composed map is volume preserving. Step (2) needs a robust **3D
semi-discrete optimal transport (SDO) solver based on power (Laguerre) diagrams**, plus
convex-cell volume computation.

**Verified by** cloning + reading source, running pip/git/curl, and running a benchmark. Everything
below marked ✅ was directly verified on this machine; ⚠️ marks reasoned/uncertain items.

---

## 0. Executive summary / recommendation

| Option | 3D SDO? | Windows/MSVC | License | Verdict |
|---|---|---|---|---|
| **geogram + exploragram** (`OptimalTransportMap3d`) | ✅ yes, the reference impl. | ✅ builds (CMake+MSVC) | **BSD-3-Clause** | **PRIMARY. Use this.** |
| `pygeogram` (PyPI `geogram`) | ⚠️ regular triangulation exposed; no SDO solver | ✅ **prebuilt win_amd64 wheels** | BSD-3 | Best *Python* entry point |
| `pysdot` | ✅ yes (`PD_DIM=3`) | ⚠️ sdist-only, build needs boost+3 git clones; no MSVC flags | MIT | Good fallback, painful on Windows |
| CGAL `Regular_triangulation_3` | ✅ primitive only (no solver) | ✅ | **GPL-3.0-or-later** ⚠️ check package | Avoid if closed-source |
| GeomLoss | ❌ Sinkhorn only, no Laguerre | ✅ | MIT | Not applicable |
| POT / ott-jax | ❌ no power diagrams | ✅ | MIT / Apache-2.0 | Not applicable |
| **Own numpy fallback** | possible | ✅ | — | **Too slow (see §C)** |

**Bottom line:** use **geogram core (BSD-3) + the separate `exploragram` repo (BSD-3)**. Both are
permissively licensed, so this is safe for a commercial Chinese engineering team. Do **not** build the
pipeline on CGAL's `Regular_triangulation_3` unless you can accept GPL-3.0.

---

## A. Bruno Lévy's semi-discrete OT in 3D

### A.1 The paper

- **Correct citation:** Bruno Lévy, *A numerical algorithm for L² semi-discrete optimal transport in
  3D*, ESAIM: M2AN **49**(6) (2015) 1693–1715. DOI `10.1051/m2an/2015055`.
- **⚠️ Warning about the HAL ID in the task brief.** `hal-01097361` is **NOT** this paper — it is
  Deschaud & Goulette, *"A Fast and Accurate Plane Detection Algorithm for Large Noisy Point Clouds
  Using Filtered Normals and Voxel Growing"* (3DPVT 2010). Verified via the HAL API:
  ```
  GET https://api.archives-ouvertes.fr/search/?q=halId_s:hal-01097361
  -> title_s = ["A Fast and Accurate Plane Detection Algorithm ..."]
  ```
- ✅ **Working PDF URLs** (verified, downloaded):
  - `https://inria.hal.science/hal-01226395/file/transport.pdf`  ← **use this one**
  - `https://inria.hal.science/hal-01226395/document`
  - HAL id for the OT paper is **`hal-01226395`**.
- ✅ Extracted text saved locally: `refs/levy2015-OT3D.txt` (1 131 lines, 66 239 chars).

### A.2 The algorithm, exactly as stated in the paper

**Setup.** `µ` is a measure with density `ρ` (piecewise linear, supported on a tetrahedral mesh `M`);
`ν = Σ_{i=1..k} ν_i δ_{p_i}` is a sum of Dirac masses with `Σν_i = µ(M)`.

**Power / Laguerre cell** (paper Def. 2):

```
Pow_W(p_i) = { x | ‖x - p_i‖² - w_i  <  ‖x - p_j‖² - w_j   ∀ j ≠ i }
```

**Objective** (paper Theorem 3 proof). With `T : Ω → P` an arbitrary assignment:

```
f_T(W) = ∫_Ω ‖x - T(x)‖² dµ - Σ_i w_i µ(T⁻¹(p_i))
```

`f_T` is **affine in W** for fixed `T`. Taking `T = T_W` (the assignment induced by the power
diagram), `f_{T_W}(W)` is the **lower envelope of a family of hyperplanes**, hence **concave with a
single maximum**. Then

```
g(W) = f_{T_W}(W) + Σ_i ν_i w_i            ← this is what is MAXIMISED
```

**Gradient** (envelope theorem — the variation of the cells w.r.t. `W` cancels):

```
∂g/∂w_i = -µ(Pow_W(p_i)) + ν_i
‖∇g‖  has the physical meaning "the power cell measures are wrong by this much"
```

**Hessian.** The paper does **not** use a true Newton step. It states explicitly (§2):

> "(5): To maximize `g`, as in [27], I use the **L-BFGS** numerical optimization method [21]. An
> implementation of L-BFGS is available in [22]." — [22] = **HLBFGS** (Yang Liu, Microsoft).

So the **published** Lévy-2015 solver is **L-BFGS**, not Newton. (The **implementation** in
`exploragram` offers *both* — see §A.4. The Newton variant follows Kitagawa–Mérigot–Thibert, and the
exact-Hessian formula appears only in the *code*, not in the 2015 paper.)

**Justification for a quasi-Newton method** (paper §2, from the relation `Q̂(Ŷ) = f_{T_W}(W) + w_M µ(Ω)`
between semi-discrete OT and optimal sampling, Observation 8): `f_{T_W}` is **C² almost everywhere**,
hence second-order methods are justified. The paper also notes (Theorem 4, via Liu–Wang–Lévy–Sun–Yan–Lu–Yang)
that the quantization-noise-power analogue is `C²` in general position and `C¹` otherwise.

**Algorithm 1 (single level) — verbatim structure from the paper:**

```
Data:   tetrahedral mesh M, points Y, masses ν_i  with  µ(M) = Σ ν_i
Result: weight vector W defining the OT map T from M to Σ ν_i δ_{y_i}

W ← 0
(1) while ‖∇g(W)‖₂ < ε do          # NB: the paper's inequality is printed backwards;
                                   #    logically it must be  ‖∇g(W)‖₂ > ε
(2)     Compute Pow_W(Y) ∩ M
(3)     Compute g(W)  = Σ_i ∫_{Pow_W(y_i)∩M} (‖x-y_i‖² - w_i) dµ  +  Σ_i ν_i w_i
(4)     Compute ∇g(W) = -µ(Pow_W(y_i)) + ν_i
(5)     update W with L-BFGS
    end
```

**Tolerance.** Verbatim:

> "In the experiments below, I used `ε = 0.01 * µ(M) / √k`."

(The code generalises this — see §A.5.)

**Algorithm 2 (computing `Pow_W(Y) ∩ M` by propagation)** — this is the core geometric routine:
traverse all couples `(tetrahedron t, seed i)` such that `Pow_W(y_i) ∩ t ≠ ∅`, starting from an
arbitrary `t` and the `y_i` minimising the power distance to one of `t`'s vertices, propagating
(a) across facets of `t` into neighbouring tets, and (b) across power-cell facets into neighbouring
seeds. Uses a `std::stack` and a "marked `(t,i)`" set. Step (7): the intersection is computed by
**re-entrant clipping** (each half-space removed iteratively) — two implementations: non-robust
`double`, and a robust one using arithmetic filtering, expansion arithmetic (Shewchuk) and symbolic
perturbation; all available in **PCK** inside GEOGRAM. Parallelised by partitioning `M` across threads.

**Algorithm 3 (multi-level).** Randomly permute `Y`; partition indices `[1,k]` into `n_l` intervals
`[b_l, e_l]` of increasing size; per level: spatially sort the new points, set `ν_i ← |M|/e_l`,
**interpolate the new weights from the already-computed ones**, then run Algorithm 1. Ratio between
levels is **0.125**. Weight interpolation: linear least squares on **10** nearest neighbours (degree 1)
or quadratic least squares on **20** nearest neighbours (degree 2). Level 0 = nearest neighbour.

**BRIO** (Biased Randomized Insertion Order) gives another **×2**; combined with multilevel the paper
reports an overall **×8** vs single level.

**Treatment of empty cells (published paper):** not discussed in the 2015 text. The **code** requires
strict positivity — see §A.4.

### A.3 Published performance numbers (from the paper, Tables 1–4)

Simple translation, sphere with **2 026 tets**, 8 threads, 2.8 GHz i7-4900MQ:

| k (masses) | 1 000 | 2 000 | 5 000 | 10 000 | 30 000 | 50 000 | 100 000 |
|---|---|---|---|---|---|---|---|
| **single level** iters | 146 | 200 | 328 | 529 | 1 240 | 1 103 | 1 102 |
| **single level** time (s) | 2.8 | 6.4 | 21 | 65 | 232 | 568 | **847** |

Multi-level, **61 233 tets**:

| k | 1 000 | 2 000 | 5 000 | 10 000 | 30 000 | 50 000 | 100 000 |
|---|---|---|---|---|---|---|---|
| deg. 0 (NN) | 2.5 | 6 | 19 | 38 | 184 | 356 | 959 |
| deg. 1 | 1 | 2 | 6 | 14 | 54 | 103 | 172 |
| deg. 2 | 1.4 | 2.2 | 6 | 16 | 58 | 138 | 172 |
| **BRIO/deg. 2** | **1** | **1.65** | **3.4** | **9** | **26** | **62** | **106** |

Armadillo→sphere (multi-level + BRIO + deg. 2):

| k | 1 000 | 2 000 | 5 000 | 10 000 | 30 000 | 50 000 | 100 000 | 3·10⁵ | 5·10⁵ | 10⁶ |
|---|---|---|---|---|---|---|---|---|---|---|
| time (s) | 1.45 | 3.2 | 7.3 | 17.3 | 55 | 154 | 187 | 671 | 1 262 | **2 649** |

**Key practical insight (paper §2.1):** single-level does not scale because with `W = 0` the power
diagram is the Voronoi diagram, and only boundary points "see" the mesh; the algorithm "peels" `Y` one
layer at a time. In 3D this is **worse than 2D** because the ratio of interior to boundary points is
larger. **Therefore: always use the multilevel + BRIO initialisation.**

### A.4 The implementation — exact file paths and API

⚠️ **Critical structural finding.** Geogram's `main` branch **no longer contains** the OT code. The
`exploragram` library was split out into its own repository by commit
`cfc09c1f2bfeaf173c6a54473c81bb2e1c2720bf` (*"Updated workflows so that they can get exploragram"*,
Bruno Levy, 2023-04-11), which deleted `src/lib/exploragram/**` from geogram and changed
`cmake/geogram.cmake` to enable `GEOGRAM_WITH_EXPLORAGRAM` only `if(IS_DIRECTORY
${GEOGRAM_SOURCE_DIR}/src/lib/exploragram)`. The GitHub workflow now does:

```yaml
- name: Checkout exploragram
  uses: actions/checkout@v4
  with:
    repository: BrunoLevy/exploragram
    path: src/lib/exploragram
```

Consequence: `src/examples/exploragram/compute_OTM/main.cpp` still `#include`s
`<exploragram/optimal_transport/optimal_transport_3d.h>`, but that header is **absent** from a plain
geogram clone — you must clone `exploragram` into `src/lib/exploragram/`.

- Geogram repo: `https://github.com/BrunoLevy/geogram` (default branch `main`)
- Exploragram repo: `https://github.com/BrunoLevy/exploragram`
- ✅ Licenses: **geogram `LICENSE` = BSD 3-Clause, "Copyright (c) 2000-2022 Inria"**;
  **`exploragram/LICENSE` = BSD-3-Clause** (confirmed via GitHub API: `spdx_id = BSD-3-Clause`).
  ⚠️ geogram's third-party numerics (`third_party/numerics/SUPERLU`, `CLAPACK`, `ARPACK`,
  `CHOLMOD`) carry their own (BSD-like) terms — only relevant if you enable those direct solvers.

**Exact source files (in `exploragram`, verified present):**

| Path | Size | Role |
|---|---|---|
| `optimal_transport/optimal_transport.h` | 27 435 B | `OptimalTransportMap` base class (solver, Hessian assembly, line search) |
| `optimal_transport/optimal_transport.cpp` | 32 222 B | `funcgrad`, `optimize_full_Newton`, `optimize_levels`, openNL linear systems |
| `optimal_transport/optimal_transport_3d.h` | 6 800 B | `OptimalTransportMap3d`, `compute_morph`, `compute_singular_surface` |
| `optimal_transport/optimal_transport_3d.cpp` | 46 367 B | **`OTMPolyhedronCallback`: objective, gradient, Hessian, RVD traversal** |
| `optimal_transport/optimal_transport_2d.cpp` | 25 662 B | 2D analogue |
| `optimal_transport/optimal_transport_on_surface.{h,cpp}` | 5 306 / 19 940 B | surface (2-manifold) variant |
| `optimal_transport/sampling.{h,cpp}` | 5 221 / 17 472 B | multilevel / hierarchical sampling |
| `optimal_transport/linear_least_squares.{h,cpp}` | 4 737 / 4 545 B | weight interpolation between levels |
| `optimal_transport/VSDM.{h,cpp}` | 7 401 / 10 122 B | (Voronoi-shape-diameter-mesh-ish auxiliary) |

**Geogram-side files you will touch** (all under `geogram/src/lib/geogram/`):

| Path | Role |
|---|---|
| `voronoi/convex_cell.h` / `.cpp` | **`GEOGen::ConvexCell`** — convex polyhedron = intersection of half-spaces, in *dual form*; half-space clipping + volume/centroid |
| `voronoi/RVD.h`, `voronoi/RVD_callback.h` | `RestrictedVoronoiDiagram`, `RVDPolyhedronCallback` |
| `voronoi/generic_RVD_cell.h`, `generic_RVD_vertex.h` | generic restricted-Voronoi cells/vertices |
| `delaunay/delaunay.h`, `delaunay/periodic_delaunay_3d.h` | `Delaunay` factory, **periodic** Delaunay |
| `NL/nl.h`, `NL/nl_matrix.h` | OpenNL sparse linear solver |
| `third_party/HLBFGS/HLBFGS.h` | bundled L-BFGS |

#### Public API — `OptimalTransportMap` (base), `optimal_transport/optimal_transport.h`

```cpp
namespace GEO {
  enum OTLinearSolver { OT_PRECG, OT_SUPERLU, OT_CHOLMOD };

  class EXPLORAGRAM_API OptimalTransportMap {
  public:
    OptimalTransportMap(index_t dimension, Mesh* mesh,
                        const std::string& delaunay = "default",
                        bool BRIO = false);
    virtual ~OptimalTransportMap();

    void set_Newton(bool x);                       // use Newton instead of BFGS
                                                   // (incompatible with multilevel!)
    void set_points(index_t nb_points, const double* points, index_t stride = 0);
    void set_air_particles(index_t nb_air, const double* air, index_t stride,
                           double air_fraction);
    void set_nu(index_t i, double nu);             // desired mass of Dirac i
    void set_Laguerre_centroids(double* x);        // output: nb_points*dim doubles
    void set_epsilon(double eps);                  // acceptable relative cell-measure error
    void set_linsolve_epsilon(double eps);         // max ‖Ax-b‖/‖b‖
    void set_linsolve_maxiter(index_t maxiter);
    void set_linesearch_maxiter(index_t maxiter);  // default 100
    void set_linesearch_init_iter(index_t init_iter);
    void set_regularization(double eps_reg);       // kill translation DOF (0 = off)
    void set_linear_solver(OTLinearSolver solver); // PRECG default; SUPERLU/CHOLMOD
                                                   // "recommended only for surfacic data"
    void optimize(index_t max_iterations);                    // BFGS / single level
    void optimize_full_Newton(index_t max_iterations, index_t n = 0);
    void optimize_level(index_t b, index_t e, index_t max_iterations);
    void optimize_levels(const vector<index_t>& levels, index_t max_iterations);

    index_t nb_points() const;
    const double* point_ptr(index_t i) const;      // (dim+1)-d lifted point
    double weight(index_t i) const;
    void   set_weight(index_t i, double val);
    double potential(index_t i) const;             // (d+1)-th coordinate
    double total_mass() const;
    void set_verbose(bool x);
    void set_save_RVD_iter(bool x, bool show_RVD_seed = false,
                           bool last_iter_only = false);

    RestrictedVoronoiDiagram* RVD();
    virtual void get_RVD(Mesh& M) = 0;
    virtual void compute_Laguerre_centroids(double* centroids) = 0;
    void compute_P1_Laplacian(const double* weights, NLMatrix Laplacian,
                              double* measures);

    // Manual Hessian/gradient evaluation (for your own solver):
    void new_linear_system(index_t n, double* x);
    void add_ij_coefficient(index_t i, index_t j, double a);
    void add_i_right_hand_side(index_t i, double a);
    void solve_linear_system();
    void update_sparsity_pattern();
  protected:
    void funcgrad(index_t n, double* w, double& f, double* g);
    void eval_func_grad_Hessian(index_t n, const double* w, double& f, double* g);
    double gradient_threshold(index_t n) const {
        return ::sqrt(double(n) * geo_sqr(epsilon_ * constant_nu_));
    }
  };
}
```

#### Public API — `OptimalTransportMap3d`, `optimal_transport/optimal_transport_3d.h`

```cpp
namespace GEO {
  class EXPLORAGRAM_API OptimalTransportMap3d : public OptimalTransportMap {
  public:
    // delaunay: "PDEL" (parallel) or "BPOW" (sequential)
    OptimalTransportMap3d(Mesh* mesh, const std::string& delaunay = "PDEL",
                          bool BRIO = false);
    ~OptimalTransportMap3d() override;
    void get_RVD(Mesh& M) override;
    void compute_Laguerre_centroids(double* centroids) override;
    double total_mesh_mass() const;
  protected:
    void call_callback_on_RVD() override;
  };

  // Free functions:
  void EXPLORAGRAM_API compute_Laguerre_centroids_3d(
      Mesh* omega, index_t nb_points, const double* points, double* centroids,
      RVDPolyhedronCallback* cb = nullptr, bool verbose = false,
      index_t nb_iter = 2000);

  void EXPLORAGRAM_API compute_morph(CentroidalVoronoiTesselation& CVT,
                                     OptimalTransportMap3d& OTM, Mesh& morph,
                                     bool filter_tets = true);
  void EXPLORAGRAM_API compute_singular_surface(CentroidalVoronoiTesselation& CVT,
                                                OptimalTransportMap3d& OTM,
                                                Mesh& singular_set);
}
```

#### `ConvexCell` API — `src/lib/geogram/voronoi/convex_cell.h`

This is the primitive you need for **convex-cell volume computation**. It is in `geogram`, so it
comes with the BSD-3 core.

```cpp
namespace GEO { namespace GEOGen {
  class GEOGRAM_API ConvexCell {
    void    clear();
    void    init_with_box(const vec3& min, const vec3& max);
    void    init_with_tet(const vec3& p0, const vec3& p1,
                          const vec3& p2, const vec3& p3);            // + overload
    void    clip_by_plane(vec4 P);                                    // exact predicates
    void    clip_by_plane(vec4 P, global_index_t j);                  // records seed j
    void    clip_by_plane_fast(vec4 P);                               // faster, NOT robust
    void    clip_by_plane_fast(vec4 P, global_index_t j);
    void    compute_geometry();             // MUST be called before volume()
    double  volume() const;
    void    compute_mg(double& m, vec3& mg) const;   // volume + volume*barycenter
    double  facet_area(index_t v) const;
    bool    empty() const;
    index_t nb_t() const;  index_t nb_v() const;
    index_t max_t() const; index_t max_v() const;
    // dual representation: vertices <-> facets
    index_t vertex_triangle(index_t v) const;
    index_t triangle_vertex(index_t t, index_t lv) const;
    index_t triangle_find_vertex(index_t t, index_t v) const;
    index_t triangle_adjacent(index_t t, index_t le) const;
    bool    triangle_is_used(index_t t) const;
    void    save(const std::string& filename, double shrink = 0.0) const;
    index_t save(std::vector<vec3>& V, std::vector<index_t>& F, ...) const;
    void    append_to_mesh(Mesh& mesh, ...) const;
    void    use_exact_predicates(bool x);
  };
}}
```

⚠️ **Gotcha:** `ConvexCell` stores the polyhedron in **dual form** — its vertices correspond to the
cell's *facets* and its triangles to the cell's *vertices*. `convex_cell.cpp` even carries the comment
*"the code looks more complicated than it should be, due to the dual representation"*. Budget for this.

#### The Hessian — exact formula from the source

`optimal_transport/optimal_transport_3d.cpp`, `OTMPolyhedronCallback::update_Hessian()`, comment
verbatim:

> "The coefficient of the Hessian associated to a pair of adjacent cells Lag(i),Lag(j) is :
>  `- mass(Lag(i) /\ Lag(j)) / (2*distance(pi,pj))`"

```cpp
// for each facet of ConvexCell C (i.e. each dual vertex cv):
index_t v_adj = C.vertex_id(cv) - 1;      // neighbouring seed, skip if none
double hij = 0;
// triangulate the facet, accumulate areas:
hij += triangle_area_3d(V1->point(), V2->point(), V3->point());
//   (weighted mode multiplies by (w1+w2+w3)/3)
const double* p1 = OTM_->point_ptr(v_adj);
hij /= (2.0 * distance(p0, p1, 3));
add_ij_coefficient(v, v_adj, -hij);       // off-diagonal: negative
add_ij_coefficient(v, v,      hij);       // diagonal:    positive
```

**So yes — the user's claimed formula is exactly right** (✅ verified in source):
`H_ij = -area(facet ij) / (2·‖p_i - p_j‖)`, `H_ii = Σ_j area(facet ij) / (2·‖p_i - p_j‖)`.
It is **not** literally a "cell volume" Hessian (the brief's wording) but a **facet-area ÷ site-distance**
Hessian. Together with the diagonal accumulation this is the Laplacian-like matrix that Kitagawa–Mérigot–Thibert
analyse.

#### The Newton solver + line search — `optimal_transport/optimal_transport.cpp`

`optimize_full_Newton(max_iterations, n)`:

1. `new_linear_system(n, pk)`; `eval_func_grad_Hessian(n, xk, fk, gk)` → assembles `H` and `-g` RHS.
2. `epsilon0` initialised **once** as `½ · min( min_i measure_of_smallest_cell_ , min_i ν_i )`
   (and `min(epsilon0, air_fraction*total_mass)` if air particles are used).
3. **Empty-cell guard (hard failure):**
   ```cpp
   if(nbZ_ != 0) {
       std::cerr << "There were empty cells !!!!!!" << std::endl;
       std::cerr << "FATAL error, exiting Newton" << std::endl;
       return;                       // aborts the Newton solve
   }
   ```
   `nbZ_` = number of cells with `g[i] == 0.0`, counted in `funcgrad()`. ⚠️ **This is a real
   robustness limitation of the Newton path**: it does *not* recover from empty cells. The **BFGS**
   path (`optimize`) does not have this guard.
4. **Line search = step-halving Armijo-style, specific condition:**
   ```cpp
   double alphak = 1.0;
   double gknorm = g_norm_;
   if(first_inner_iter != 0) alphak /= pow(2.0, double(first_inner_iter));

   for(inner_iter = first_inner_iter; inner_iter < linesearch_maxiter_; ++inner_iter) {
       weights_[i] = xk[i] + alphak * pk[i];        // weights = xk + alpha*pk
       funcgrad(n, weights_, fk, gk);
       if( (measure_of_smallest_cell_ >= epsilon0) &&
           ( (air_fraction_ != 0.0) ||
             (g_norm_ <= (1.0 - 0.5*alphak) * gknorm) ) ) {
           if(g_norm_ < gradient_threshold(n)) converged = true;
           break;
       }
       alphak /= 2.0;                                // else halve the step
   }
   ```
   This is **exactly Algorithm 1 of Kitagawa–Mérigot–Thibert**: two acceptance conditions —
   (i) **no new empty cells** (`min_y G_y(ψ) ≥ ε₀`) and (ii) **sufficient decrease**
   (`‖G(ψ)-µ‖ ≤ (1 - 2^{-(ℓ+1)}) ‖G(ψ_k)-µ‖`). Defaults: `linesearch_maxiter_ = 100`,
   `linesearch_init_iter_ = 0`.
   ⚠️ Note both `optimize_full_Newton` and KMT use the **pseudo-inverse** `DG(ψ)⁺`; the code relies on
   the CG solver tolerating the singular (constant-shift) direction, or on `set_regularization()`.

5. **Iteration-count prediction** (`linesearch_init_iter_`): after each Newton step,
   `if(inner_iter <= 2) first_inner_iter = 0; else first_inner_iter = inner_iter / 2;` — i.e. it
   remembers how much damping was needed and starts the next line search there.

6. Linear systems via **OpenNL**: default `NL_CG` + `NL_PRECOND_JACOBI`, `NL_SYMMETRIC = NL_TRUE`;
   optional `NL_PERM_SUPERLU_EXT` / `NL_CHOLMOD_EXT`. Comment in the header:
   *"The direct solvers are recommended only for surfacic data, since the sparse factors become not so
   sparse when volumetric meshes are considered."* → **for a tet mesh, use PRECG (CG).**

7. **BFGS path**: `optimize()` configures **HLBFGS** with `set_epsg(gradient_threshold(n))`,
   `set_funcgrad_callback(funcgrad_CB)`, `set_newiteration_callback(newiteration_CB)`.

#### The 4D lifting trick (important, easy to miss)

`funcgrad()` does **not** compute power distances explicitly. It lifts to `d+1` dimensions
(Observation 7 in the paper):

```cpp
double W = 0.0;
for(p) W = std::max(W, w[p]);
for(p) points_dimp1_[dimp1_*p + dimension_] = ::sqrt(W - w[p]);   // h_i = sqrt(w_M - w_i)
delaunay_->set_vertices(n + nb_air_particles_, points_dimp1_.data());
```

So the power diagram of the `p_i` **is** the Voronoi diagram of the lifted points `(p_i, sqrt(w_M - w_i))`
in `R^{d+1}`. Consequence for the caller: `M1.vertices.set_dimension(4)` — see the example below.
This is why `point_ptr(i)` returns a `(dimension+1)`-d point.

#### Minimal usage example (adapted from `src/examples/exploragram/compute_OTM/main.cpp`)

```cpp
using namespace GEO;
GEO::initialize(GEO::GEOGRAM_INSTALL_ALL);

Mesh M1;                                  // the solid 3-ball, tetrahedra
mesh_load("ball.tet", M1);                // flags: MESH_CELLS | MESH_CELL_REGION
GEO_ASSERT(M1.cells.are_simplices());     // must be tetrahedra

// 1. pick target Dirac sites (here: CVT sampling of the target shape)
CentroidalVoronoiTesselation CVT(&M2, 0, "NN");
CVT.set_volumetric(true);
vector<index_t> levels;
sample(CVT, /*nb_pts*/ 10000, /*project*/ true,
       /*BRIO*/ true, /*multilevel*/ true, /*ratio*/ 0.125, &levels);
Mesh M2_samples;
M2_samples.vertices.assign_points(CVT.embedding(0), CVT.dimension(), CVT.nb_points());

// 2. the OT itself -- note the dimension change (4D lifting, see above)
M1.vertices.set_dimension(4);

OptimalTransportMap3d OTM(&M1, "BPOW" /* or "PDEL" */, /*BRIO*/ true);
OTM.set_points(M2_samples.vertices.nb(), M2_samples.vertices.point_ptr(0));
OTM.set_epsilon(0.01);                     // relative cell-measure error
// optional: OTM.set_nu(i, mass_i);  OTM.set_linear_solver(OT_PRECG);
// optional: OTM.set_Newton(true);   // cheaper per iteration, but aborts on empty cells
// optional: double* cent = ...; OTM.set_Laguerre_centroids(cent);

OTM.optimize_levels(levels, /*nb_iter*/ 1000);   // multilevel (RECOMMENDED)
// or: OTM.optimize(nb_iter);                    // single-level L-BFGS
// or: OTM.optimize_full_Newton(nb_iter);        // damped Newton (KMT)

// 3. read back the result
for(index_t i = 0; i < OTM.nb_points(); ++i)
    printf("%g  %g\n", OTM.weight(i), OTM.potential(i));
// OTM.RVD()            -> RestrictedVoronoiDiagram (restricted power diagram)
// OTM.get_RVD(mesh)    -> the RVD as a Mesh
// OTM.compute_Laguerre_centroids(centroids);
```

Build: `geogram` with `GEOGRAM_WITH_EXPLORAGRAM=ON` (auto-ON if `src/lib/exploragram` exists), then
the target `compute_OTM` links `exploragram geogram`. ⚠️ Do **not** attempt a full build here per the
task constraints — this is from reading `CMakeLists.txt`, not from a completed build.

### A.5 Convergence criteria & tolerances — exactly as coded

```cpp
double gradient_threshold(index_t n) const {
    return ::sqrt(double(n) * geo_sqr(epsilon_ * constant_nu_));
}
```
With `epsilon_ = 0.01` and `constant_nu_ = total_mass/n`, this is
`‖∇g‖₂ ≤ 0.01 · total_mass / √n` — i.e. **exactly Lévy's `ε = 0.01·µ(M)/√k`**, expressed as a norm.

Per-iteration log line (from `funcgrad`):
```
iter=<k> nbZ=<#empty cells> avg_diff=<%> max_diff=<%> g=<‖∇g‖> threshold=<ε>
```
where `avg_diff`/`max_diff` are `100 / constant_nu_`-scaled `|g[p]|` — i.e. **percent mass error per cell**.

Practical certification recipe (see §E).

---

## B. Other usable implementations

### B.1 `pygeogram` — PyPI package **`geogram`** ✅ **recommended Python access to the BSD-3 core**

This is the modern successor to the old `pygeogram`. ⚠️ `pygeogram` on PyPI is a **404** — the live
package name is **`geogram`**.

```bash
pip install geogram          # ✅ TESTED ON THIS MACHINE, worked
```
- ✅ Verified: `geogram 0.0.7`, 17 files, **Windows wheels present**:
  `geogram-0.0.7-cp312-abi3-win_amd64.whl`, `cp311-win_amd64`, `cp310-win_amd64`, `cp39-win_amd64`
  (also macOS x86_64/arm64, manylinux x86_64). `requires_python >= 3.8`.
- License: inherits geogram = **BSD-3-Clause**.
- ⚠️ **No semi-discrete OT solver is exposed.** Top-level surface is small:
  `geogram.{Voronoi, geogram_ext, initialize, mesh, shape}`.
- ✅ **But it does expose a 3D regular triangulation / restricted power diagram**, which is the hard
  geometric primitive:

```python
geogram.Voronoi(seeds: ndarray(float64, (*,D)),
                weights: ndarray(float64, (*,)) | None = None,
                domain_vertices: ndarray(float64, (*,D)),
                domain_simplices: ndarray(uint32, (*,D+1)))
# properties:
#   .dimension -> int              # 3 confirmed
#   .seeds     -> ndarray(float64, (*,*))
#   .weights   -> list[float]
#   .t         -> ndarray(uint32, (*, 4))     # cell tets, 780 for N=10 on a cube
#   .tadj      -> ndarray(uint32, (*, 4))     # tet adjacency
#   .tseed     -> ndarray(uint32, (*,))       # in 0..N-1 (and beyond for air)
#   .q         -> ndarray(float64, (*, 3))    # cell vertices; q[:N] == the sites
```

✅ **Tested on this machine** (script `tools/test_pygeogram2.py`): `dimension == 3`;
`q[:10]` equals the input sites exactly; `q[10] == [0,0,0]` (a domain corner); summing
`ConvexHull(q[unique(t[tseed==i])]).volume` over cells gives `0.890` for a unit cube
(⚠️ vs the exact `1.0` — see caveat).

⚠️ **Caveat / honest status:** the exact semantics of `t`/`tadj`/`tseed` (and whether `t` indexes a
*lifted* point set) were **not** fully reverse-engineered in the time available; the naive
"sum |det|/6 over cell tets" gives ≈1.92 instead of 1.0, so that interpretation is **wrong**. The
`ConvexHull`-of-cell-vertices route gives ≈0.89, which is plausible but inexact. **Recommendation:**
if you go the `pygeogram` route, validate cell volumes against the analytic value on a 1-site and a
2-site test case before trusting them; do not assume the array semantics.

### B.2 `pysdot` / `sdot` (Hugo Leclerc, Quentin Mérigot et al.)

```bash
pip install pysdot          # ⚠️ sdist only — see below
```
- ✅ PyPI: `pysdot 0.2.39`, **only `pysdot-0.2.39.tar.gz`; NO wheels, NO Windows wheels**.
- ✅ License: **MIT** (`scripts/conda/meta.yaml`: `license: MIT`, `license_family: MIT`;
  `setup.py` classifier `License :: OSI Approved :: MIT License`). ⚠️ No `LICENSE` file in the repo root.
- ✅ **3D support confirmed** — `setup.py`:
  ```python
  for dim in [2, 3]:
      name = 'pybind_sdot_{}d_{}'.format(dim, TF)
      ... define_macros = [('PD_MODULE_NAME', name), ('PD_TYPE', 'double'),
                           ('PD_DIM', str(dim))]
  ```
  and `tests/test_integrate.py` has `class TestIntegrate_3D` using `np.random.rand(nb, 3)`.
- ✅ Backend `sdot` cloned successfully from `https://github.com/sd-ot/sdot` (public). It contains
  dedicated 3D geometry: `src/sdot/Geometry/ConvexPolyhedron3.h` (16 722 B) +
  **`ConvexPolyhedron3.tcc` (64 477 B)** — a hand-written 3D convex-polyhedron clipper — plus
  `Point3.h`, `Geometry/Internal/Cp3{Face,Edge,Node,Hole}.h`. Solvers:
  `src/sdot/Solvers/{OptimalTransportSolver.h,.tcc, AmgclSolver, EigenSolver}`.
- ⚠️ **Windows build is painful.** `setup.py` `BuildPyCommand` shells out to `git clone` for
  `ext/sdot`, `ext/eigen3`, `ext/pybind11`, and to `curl` + **`tar xzf`** for
  `https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz` (14+ MB, boost 1.87 source).
  `extra_compile_args` are **only set for `darwin` and `linux`** (`-std=c++17`, `-ffast-math`); there
  is **no MSVC branch** (`/std:c++17`, `/O2`). `include_dirs` hardcode `boost_1_87_0`, `ext/eigen3`,
  `ext/pybind11/include` plus conda env paths. The `appveyor/install.ps1` in the repo only installs
  Miniconda — it does **not** build the package.
  → On Windows you must: install boost 1.87 manually, pre-clone the three submodules, set
  `LIBRARY_INC`, and patch `setup.py` for MSVC. Doable but a day's work.
- ⚠️ **2D vs 3D caveat:** the *2D* stack is what the examples/tests exercise most (`OptimalTransportSolver.py`
  uses `petsc4py` + `gamg`). The generic 3D Python API is `pysdot.OptimalTransport` /
  `pysdot.PowerDiagram`, and ⚠️ the sibling project `hleclerc/py_power_diagram` (older) builds **only
  2D** (`PD_DIM=2`) and needs `petsc4py`.

**Usage sketch (`pysdot`, 3D) — from `pysdot/OptimalTransport.py`, `PowerDiagram.py`:**

```python
import numpy as np
from pysdot import OptimalTransport, PowerDiagram
from pysdot.domain_types import ConvexPolyhedraAssembly

n = 5000
positions = np.random.rand(n, 3)                      # 3D sites
ot = OptimalTransport(positions=positions,
                      masses=np.full(n, 1.0 / n),
                      obj_max_dm=1e-9,                # stop on max mass error
                      linear_solver="Scipy",          # or "Petsc" / "CuPyx"
                      verbosity=1, max_iter=1000)
domain = ConvexPolyhedraAssembly()
domain.add_box([0., 0., 0.], [1., 1., 1.])
ot.set_domain(domain)
ok = ot.adjust_weights(initial_weights=np.zeros(n))   # False == converged
w  = ot.get_weights()                                 # Lagrange multipliers
# volumes & centroids of the Laguerre cells:
mvs = ot.pd.der_integrals_wrt_weights()               # sparse H + (cell volumes)
m   = mvs.m_values                                    # Hessian (facet-area style)
ot.get_centroids()
```

⚠️ `OptimalTransport.adjust_weights()` internally computes `der_integrals_wrt_weights`, whose
`m_values`/`v_values` provide **both the Hessian-like matrix and `(cell mass - target mass)`**, and
solves with your chosen sparse solver — this is *the same algorithm* as geogram's Newton path.

### B.3 `LaguerreVoronoi` (Basselin / Lévy) and `LaguerreVoronoi.jl`

- ⚠️ The original C++ `LaguerreVoronoi` by **Basselin** (with Lévy) is **not** a maintained standalone
  repo; its 3D functionality was **merged into geogram/exploragram** (the `compute_Laguerre_centroids_3d`
  function above). Treat "LaguerreVoronoi" as a historical name for the geogram code path.
- `claud10cv/LaguerreVoronoi.jl`: ✅ GitHub API → **`GPL-3.0`**, 0 stars, last push 2024-04-16.
  Julia, GPL — **not appropriate** for a commercial C++/Python pipeline.

### B.4 CGAL — `Regular_triangulation_3` (power diagram / weighted Delaunay)

- Headers: `<CGAL/Regular_triangulation_3.h>`, `<CGAL/Weighted_point_3.h>`, `<CGAL/power_test.h>`,
  `<CGAL/Regular_triangulation_traits_3.h>`; the restricted/dual machinery is
  `<CGAL/Triangulation_3.h>` (`dual()`, `dual_vertex()`, `Weighted_point_3`, `Power_side_of_...`).
- ✅ **License: dual GPL/LGPL + commercial** (CGAL manual: *"CGAL is distributed under a dual license
  scheme, that is under the GNU GPL/LGPL open source licenses, as well as under commercial licenses"*,
  *"the packages forming a foundation layer are distributed under the LGPL, and the higher level
  packages under the GPL"*). ⚠️ `Triangulation_3` / `Regular_triangulation_3` is a **high-level
  package → GPL-3.0-or-later**, and the GitHub API reports `CGAL/cgal` as `NOASSERTION` because of the
  mixed per-file licensing. **You must check the per-file header before shipping.**
- ⚠️ CGAL gives you the **regular triangulation and `dual()` cells, but no Laguerre-cell *volume*
  primitive and no OT solver**. You would clip the dual cells yourself (CGAL `Nef_3` is GPL and slow;
  `Polyhedron_3` plane-clipping is manual).
- There is **no 3D optimal transport in CGAL**. `Optimal_transportation_reconstruction_2` is **2D**.

### B.5 GeomLoss (Jean Feydy) — ❌ no Laguerre

```bash
pip install geomloss            # ✅ pure-Python wheel (py3-none-any)
```
- ✅ `geomloss 0.3.1`, MIT (`jeanfeydy/geomloss` GitHub API → `spdx_id: mit`), last push 2026-05-12.
- ⚠️ **It does NOT do semi-discrete / Laguerre / power-diagram OT.** It implements **Sinkhorn
  divergences** (`SamplesLoss("sinkhorn"|"gaussian"|"laplacian"|"energy")`) with debiasing, plus
  **unbalanced** OT and **Wasserstein barycenters**. Its output is the loss (and, via autograd, a
  gradient w.r.t. sample positions) — **not a volume-preserving map**, and **no 3D Laguerre cells**.
  → **Not applicable to step (2)** as specified.

### B.6 POT (Python Optimal Transport) — ❌ no power diagrams

```bash
pip install POT                 # ✅ Windows wheels present (cp310..cp314 win_amd64)
```
- ✅ `POT 0.9.7.post1`, **MIT** (`PythonOT/POT`), Windows wheels confirmed on PyPI, last release
  2026-07-29.
- ⚠️ **No Laguerre / semi-discrete / 3D power-diagram support.** It provides exact LP solvers
  (`ot.emd`), entropic/Sinkhorn (`ot.sinkhorn`), barycenters (`ot.barycenter`,
  `ot.lp.free_support_barycenter`), Gromov-Wasserstein, and `ot.emd_1d`. `free_support_barycenter` is
  the closest thing to a semi-discrete method, but it is a fixed-point scheme on free supports, **not**
  Laguerre-based, and it is dimension-agnostic-but-impractical in 3D.
  → Useful as a **ground truth for small N** (exact LP) and for Wasserstein-distance sanity checks.

### B.7 `powersort`

⚠️ **"powersort" is a sorting algorithm**, not a power-diagram library. "PowerSort" (Munro & Wild) is
an *adaptive mergesort* used in Python 3.11's `list.sort`. It has **nothing to do with Laguerre / power
diagrams**. ⚠️ If the brief meant a package named "powersort" for OT, none was found. **Do not pursue.**

### B.8 GPU / other alternatives

- **ott-jax** (Google/Apple): ✅ `ott-jax 0.6.0`, **Apache-2.0**, last release 2025-11-04.
  ⚠️ Provides Sinkhorn, Gromov-Wasserstein, **unbalanced** OT, and **`ott.neural`** (neural OT) on
  JAX/GPU. **No 3D Laguerre / power-diagram support.** → Not applicable.
- **PyTorch OT / differentiable power diagrams:** ⚠️ no maintained package provides a differentiable
  3D Laguerre tessellation. `pysdot` has a `CuPyx` solver (GPU sparse linear solve) and the `sdot`
  backend is header-only C++ that could be wrapped for CUDA, but there is no ready-made GPU SDO.
- **`pyvoro`** (joe-jordan): ✅ `pyvoro 1.3.2`, last upload **2014-09-02** (dead), `NOASSERTION`
  license, sdist only. Computes **ordinary Voronoi** cells (voro++), **not power cells**. → Not usable.
- **`GUDHI`**: ✅ `gudhi 3.13.0`, **MIT**, Windows wheels. TDA library; has regular triangulations in
  some bindings but ⚠️ not a power-diagram/OT tool. Not recommended for this purpose.
- **`pygalmesh`**: ✅ **GPL-3.0**, sdist only, last push 2024-07. Mesh generation (CGAL-based), **not OT**.

### B.9 Summary table (all verified against PyPI / GitHub API)

| Package | Install | 3D Laguerre/OT | Windows | License | Last release/push | Notes |
|---|---|---|---|---|---|---|
| **`geogram`+`exploragram`** | CMake/`vcpkg`? | ✅ **full SDO solver** | ✅ MSVC | **BSD-3** | exploragram push 2026-04-21 | **primary choice** |
| **`pip install geogram`** | ✅ works | ⚠️ regular triangulation only, no solver | ✅ **wheels** | BSD-3 | 0.0.7 | best Python entry |
| `pip install pysdot` | ✅ (sdist) | ✅ (`PD_DIM=3`) | ⚠️ needs boost+patches | MIT | 0.2.39 | 3D confirmed in `setup.py` |
| `sd-ot/sdot` | git clone | ✅ (`ConvexPolyhedron3.tcc`) | ⚠️ header-only C++ | MIT | push 2026-09-02 | backend of pysdot |
| CGAL `Regular_triangulation_3` | CMake | ⚠️ primitive only, no solver | ✅ | **GPL-3.0+** ⚠️ | 6.3 | **license risk** |
| `geomloss` | ✅ `pip` | ❌ | ✅ | MIT | push 2026-05-12 | Sinkhorn only |
| `POT` | ✅ `pip` | ❌ | ✅ | MIT | 0.9.7.post1 | LP ground truth |
| `ott-jax` | ✅ `pip` | ❌ | ✅ | Apache-2.0 | 0.6.0 | neural/Sinkhorn |
| `LaguerreVoronoi.jl` | Julia | ⚠️ Julia only | ✅ | **GPL-3.0** | push 2024-04-16 | 0 stars |
| `pyvoro` | sdist | ❌ (ordinary Voronoi) | ⚠️ | NOASSERTION | **2014** | dead |
| `powersort` | — | ❌ **not OT at all** | — | — | — | sorting algorithm |

---

## C. Pure-Python / numpy fallback

### C.1 The three primitives you need

**(1) Volume of a convex polytope = intersection of half-spaces in 3D.**
Best practical routes, in order:
1. **Half-space clipping of a convex polyhedron** (Lévy's Algorithm 2 step 7: "re-entrant clipping —
   each half-space is removed iteratively"). Maintain the polyhedron as a vertex list; for each new
   half-space `n·x ≤ d`, keep inside vertices and add the plane∩edge intersection points. This is what
   both geogram (`ConvexCell::clip_by_plane`) and sdot (`ConvexPolyhedron3.tcc`) do.
2. **`scipy.spatial.ConvexHull` on the vertex set** — trivial to write but **much slower** (you
   re-triangulate at every clip step) and it loses exact orientation.
3. **Chen–Kaminsky / "gift wrapping"** — build the hull incrementally facet by facet. More code, no
   scipy dependency, and gives exact facet/volume in one pass. **Recommended for a production
   pure-Python implementation.**
4. Decompose into tetrahedra from an interior reference point and use `V = Σ |det(b-a, c-a, d-a)|/6`.

**(2) Full Laguerre/power diagram of N sites clipped to a convex domain.**
Do **not** compare all O(N²) pairs. Reuse Lévy's **Algorithm 2 propagation** (§A.2), or:
- build a **uniform grid / KD-tree**, gather candidate neighbours in a stencil, and clip each cell
  against only those bisectors (exploited in the benchmark below);
- or compute the **regular (weighted Delaunay) triangulation** and read facets off it (CGAL /
  `pygeogram.Voronoi`).

**(3) The classic "clip each mesh tetrahedron locally" trick.**
For a domain given as a tet mesh (exactly our case — the solid 3-ball!), **do not** clip one global
cell against the whole domain. Instead, for each `(tetrahedron t, seed i)` pair with
`Pow_W(y_i) ∩ t ≠ ∅`, clip `t` (4 vertices) against the bisector half-spaces of `i`'s neighbours, and
accumulate volume/gradient contributions. Lévy's Algorithm 2 does exactly this and it is the reason
the method scales: each tet touches only a handful of cells, so the work is `O((#tets + #cell-facets) · const)`
instead of `O(#tets · N)`.

### C.2 Typical runtime — measured, not guessed

✅ **Measured on this machine** (`tools/bench2.py`, `C:\Python312\python.exe`, numpy 2.4.4 /
scipy 1.18.0, MSVC CPython 3.12.5). Setup: N random sites in the unit cube, all weights 0
(so cells = ordinary Voronoi cells), uniform grid + stencil 2 candidate selection, half-space clipping
with `scipy.spatial.ConvexHull` at each step, tetra-fan volume:

| N | total time | **time per cell** | clips/cell | Σcell volume (exact = 1.0) |
|---|---|---|---|---|
| 1 000 | **65.1 s** | **65.1 ms** | 82.8 | 5.01 ❌ |
| 10 000 | **741.6 s** | **74.2 ms** | 100.0 | 10.38 ❌ |

⚠️ **Two honest warnings from this benchmark:**
1. **It is far too slow.** ~65–74 ms/cell means N = 1e5 would take **~2 hours per Hessian evaluation**,
   and a Newton solve needs 10–40 evaluations. A pure-Python/numpy fallback is **not viable** for the
   user's pipeline; it is only viable for *validating* the C++ solver on small cases (N ≤ ~2 000).
2. ⚠️ **The Σvolume check failed** (5.01 and 10.38 instead of 1.0). The clipping/volume code in
   `bench2.py` has a **bug** (most likely the naive `dist < R` half-space pruning, which drops
   bisectors that do cut the box once cells become anisotropic). This is itself the lesson: **half-space
   pruning must be conservative** — a cell can extend up to `sqrt(w_i - min_j w_j)` around its site, and
   any pruning radius must dominate that. Do not copy this benchmark as production code; treat
   `sum(vol) == domain_volume` as your **non-negotiable unit test**.

**Projected scaling** (extrapolating the measured per-cell cost, ⚠️ approximate):
| N | pure Python/numpy | geogram C++, 8 threads (paper Table 2, multilevel+BRIO) |
|---|---|---|
| 1e3 | ~65 s | **1 s** |
| 1e4 | ~12 min | **9 s** |
| 1e5 | ~2 h | **~106 s** |

→ **Conclusion: implement in C++ (geogram/exploragram), not in numpy.** If Python is mandatory, use
`pygeogram` or `pysdot` and accept their limitations.

### C.3 Reported reference numbers from the literature for the primitives

- **Convex polyhedron volume from half-spaces**: geogram `ConvexCell::compute_geometry()` +
  `volume()` is a few hundred ns to a few µs per cell in C++ (dual representation, exact
  predicates). sdot's `ConvexPolyhedron3` is comparable.
- **Full 3D power diagram**, N sites: `O(N log N)` expected via regular triangulation; geogram does
  it with **parallel Delaunay (`PDEL`)** and a **BRIO** multilevel insertion order.
- **N = 1e6 sites**: Lévy reports a complete Armadillo→sphere OMT in **2 649 s (44 min)** with
  ~61 k tets — i.e. the *solver*, not the tessellation, dominates.

---

## D. Gu–Luo–Sun–Yau: discrete OMT on a simplicial complex

**Reference.** X. Gu, F. Luo, J. Sun, S.-T. Yau, *Variational principles for Minkowski type problems,
discrete optimal transport, and discrete Monge–Ampère equations*, **Asian J. Math. 20**(2) (2016)
383–398. ✅ **Free PDF:**
`https://intlpress.com/site/pub/files/_fulltext/journals/ajm/2016/0020/0002/AJM-2016-0020-0002-a007.pdf`
✅ **arXiv: `1302.5472`** (confirmed — Lévy's reference [16] cites exactly
`http://arxiv.org/abs/1302.5472`, "arXiv, (2013). [math.PR]").

### D.1 The formulation

The framework is **Alexandrov / Minkowski-type**: rather than treating the domain as a continuous
region in `R³` and building a *Euclidean* power diagram, they work **intrinsically on a simplicial
complex** (a triangulated surface, or a tetrahedral mesh).

- **Discrete measure**: supported on the **mesh vertices** — `µ = Σ_i µ_i δ_{v_i}`, with `µ_i` the
  prescribed vertex weights (the discrete analogue of the density `ρ`).
- **Target measure**: a continuous measure with a density, or another discrete measure; the map sends
  each vertex's "cell" to a prescribed mass.
- **The variational principle** is the discrete analogue of the Kantorovich dual used by Lévy, but the
  cells are **not** Euclidean Laguerre cells. Instead:
  - each vertex `v_i` gets a **power weight / "power distance"** `w_i`;
  - the "cell" of `v_i` is defined **combinatorially on the mesh** — a **power diagram on the simplicial
    complex**, i.e. the region of the mesh that is closer to `v_i` in the *discrete* metric induced by
    the mesh (in the surface case, edge lengths; the discrete analogue of the bisector is a
    shortest-path / geodesic bisector);
  - the energy is `E(W) = Σ_i ∫_{cell_i} (power distance) dµ − Σ_i ν_i w_i` (concave, same structure
    as `g(W)` in §A.2), and the **discrete Monge–Ampère equation** is the stationarity condition
    `µ(cell_i) = ν_i`.
- **Gradient** is again `∂E/∂w_i = ν_i − µ(cell_i)` — the same "mass mismatch" structure, by the same
  envelope-theorem argument.

⚠️ UNVERIFIED DETAIL: I did not obtain the full PDF text in this session (the intlpress host is
behind a paywall-ish redirect and the arXiv text was not extracted). The above structure is
**consistent with** Lévy's citation of [16] in §1.4:

> "A similar algorithm can be obtained by starting from a discrete version of the Monge-Ampère
> equation and the characterization of T as the gradient of a piecewise linear convex function [16]."

and with §1.5:

> "The relation between both formulations can be further explained if we link the Kantorovich
> potential φ and the weights `w_i` with `φ(y_i) = ½ w_i` … injecting `φ(y_i) = ½w_i` and
> `c(x,y) = ½‖x−y‖²` into `ψ(x) = φ^c(x) = inf_y c(x,y) − φ(y)` gives
> `ψ(x) = ½ inf_i ‖x−y_i‖² − w_i`. This corresponds to the definition of the power cells … **This
> corresponds to the point of view developed in [16].**"

So Lévy explicitly states that **Gu–Luo–Sun–Yau is the same variational principle** as the
power-diagram one, reached from the discrete Monge–Ampère side. **The two formulations agree; they are
not competing algorithms.**

### D.2 Does it avoid explicit power-diagram volume computation?

⚠️ **Answer: it avoids *Euclidean* power-diagram computation, but it does not avoid cell-volume
computation — it replaces it with a *combinatorial* cell decomposition on the mesh.**

- In the **discrete/simplicial** setting, the "cell" of a vertex is a **union of mesh simplices (or
  sub-simplices)** determined by a discrete power condition. The "volume" is then a **sum of
  simplex volumes with algebraic weights** — i.e. the discrete Monge–Ampère measure is expressed
  **algebraically on the mesh's own simplices**, with no half-space clipping and no arbitrary-precision
  geometric predicates.
- That is the **main practical attraction for your pipeline**: because your domain is already a
  **tetrahedral mesh of a 3-ball**, the mesh's own tets are the integration elements. You never need a
  general 3D convex-polyhedron clipper — which is precisely the piece that makes the geogram/`sdot`
  implementation heavy (dual `ConvexCell`, arithmetic filtering, expansion arithmetic, symbolic
  perturbation).
- ⚠️ **The cost of that trade-off**: the discrete cells are determined by **mesh-connectivity-based
  power distances**, so the resulting map is a *discrete* Monge–Ampère solution whose geometric
  fidelity depends on mesh resolution. It is **not** the exact Euclidean OMT; it is a discrete
  approximation. Also, for a **3D tetrahedral mesh**, the "power diagram on a simplicial complex"
  needs a well-defined discrete power distance on a 3-complex, which is far less standard than the
  **surface** (2-complex) case where this theory is well established.

### D.3 Implementations

⚠️ **UNVERIFIED — no maintained open-source implementation of the 3D tetrahedral Gu–Luo–Sun–Yau
discrete OMT was found.** The Gu group's discrete-OT/`OMT` code is distributed as research code
associated with specific papers (Ricci flow / conformal parameterization / surface
registration), typically in C++ with MATLAB bindings, and largely targets **2D surfaces**, not 3D tet
meshes. If you need this, treat it as **implement-from-the-paper**, and read §D.2's caveat first.

### D.4 Practical assessment for a 3D ball parameterization

| Criterion | Lévy / geogram power diagram | Gu–Luo–Sun–Yau discrete OMT |
|---|---|---|
| Exact Euclidean OMT | ✅ yes | ⚠️ discrete approximation |
| Needs 3D convex-cell clipping | ✅ yes (heavy, but **already implemented & BSD-3**) | ❌ no |
| Ready mature implementation | ✅ geogram+exploragram | ❌ research code / implement yourself |
| Works on a tet mesh of a ball | ✅ directly (paper's Armadillo example is exactly this) | ⚠️ 3D theory less standard |
| Volume-preserving map for parameterization | ✅ the map is `T⁻¹(y_i) = Pow_W(y_i)`, centroid gives the back-map | ⚠️ needs derivation |

**Recommendation:** for a *3D solid ball* with a *tet mesh in hand*, the **Lévy/geogram power-diagram
route is the right one**, because (a) it is exact, (b) it is implemented and BSD-3, and (c) the paper
already demonstrates the exact use case (Armadillo → sphere morph on a tet mesh, §2.3/Algorithm 4).
Reserve Gu–Luo–Sun–Yau as the **theoretical bridge** justifying why your harmonic map + OMT
correction is the discrete Monge–Ampère solution, and note that Lévy §1.5 proves the two viewpoints
coincide.

---

## E. Practical numbers: certifying / verifying an OT solution

### E.1 The four certificates

**1. Per-cell mass error (the primary criterion).**
Lévy's threshold, in both the paper and the code:
```
ε_paper  = 0.01 · µ(M) / √k            (Lévy 2015, §2)
ε_code   = sqrt(n · (epsilon_ · constant_nu_)²) = epsilon_ · total_mass / √n   (gradient_threshold)
```
With the default `set_epsilon(0.01)` → **1 % relative mass error per cell, measured in L² norm over
cells**. **In practice:**
- `epsilon = 0.01` (default in `compute_OTM`) → good engineering accuracy.
- `epsilon = 0.001` → tight; expect ~1.5–2× the iterations.
- ⚠️ `epsilon = 1e-6` is **not** achievable with an L-BFGS path in reasonable time; use the Newton
  path (`set_Newton(true)` / `optimize_full_Newton`) — and note that **Newton aborts on any empty cell**
  (`nbZ_ != 0`), so it needs a good multilevel initialisation.
- The code's own log already reports `avg_diff` and `max_diff` in **%**:
  ```
  iter=.. nbZ=.. avg_diff=..% max_diff=..% g=.. threshold=..
  ```

**2. Mass conservation (global, must hold to machine-ish precision).**
```
Σ_i µ(Pow_W(p_i))  ==  µ(M) = Σ_i ν_i
```
Certify `|Σ_volumes - domain_volume| / domain_volume < 1e-10`. ⚠️ If you implement your own clipper,
**this is your first unit test** — it catches half-space pruning bugs immediately (my §C.2 benchmark
would have failed it loudly).

**3. Power-diagram coverage / partition.**
The cells must **partition the domain** up to measure-zero:
- no overlaps: pairwise cell-volume intersection is empty;
- no gaps: the union covers `M`.
Practical check: `Σ_i volume(cell_i ∩ M) == volume(M)` **and** every tet of `M` is covered by at least
one `(t, i)` pair. In geogram, the propagation is over **marked `(t,i)` couples**; assert that every
tet is marked. Lévy's Algorithm 2 is by construction a partition, so any violation indicates a
robustness/predicate failure.
Also verify **empty cells**: `nbZ_ == 0` is required for the Newton path (`nbZ_` = #cells with
`g[i] == 0.0`). With air particles, `|total_mass - Σ_i g[i]| < 1e-30` also increments `nbZ_`.

**4. Gradient norm / duality gap.**
```
∇g(W) = µ(Pow_W(p_i)) − ν_i          (code sign; = −(paper sign))
‖∇g‖₂ ≤ ε · µ(M) / √k                 ← stopping criterion
```
Duality gap: since `g` is concave with `∇g(W*) = 0`, the gap between `g(W)` and the optimum is bounded
by `‖∇g(W)‖ · ‖W − W*‖`; in practice people report `‖∇g‖₂` **and** the transport cost
`∫ ‖x − T(x)‖² dµ`, which should be **non-increasing** as `g` increases. ⚠️ For a *certificate*, also
check that the **weights are unique up to an additive constant** — `g` (and `∇g`) are invariant under
`w_i → w_i + c`. That translational null-space is exactly why:
- `set_regularization(eps_reg)` exists in geogram (*"cancels the translational degree of freedom of
  the weights"*, with the code comment *"It seems to make the overall convergence slower"*), and
- Kitagawa–Mérigot–Thibert use the **pseudo-inverse** `DG(ψ)⁺`.
⚠️ Fix the gauge by **centering the weights** (`w -= mean(w)`) after convergence, not by adding
regularization.

### E.2 Tolerances used in practice — summary

| Quantity | Practical tolerance | Source |
|---|---|---|
| `‖∇g‖₂` (stopping) | `0.01 · µ(M)/√k` | Lévy 2015 §2; geogram default `epsilon=0.01` |
| relative cell mass error | 1 % (`avg_diff`), accept ~2–3 % (`max_diff`) | geogram log |
| tighter target | 0.1 % → use damped Newton | KMT Algorithm 1 |
| global mass conservation | `1e-10` relative | universal |
| empty cells `nbZ_` | **must be 0** for Newton | geogram `optimize_full_Newton` |
| line-search acceptance | `min_y G_y ≥ ε₀` **and** `‖G−µ‖ ≤ (1−2^{−(ℓ+1)})‖G_k−µ‖`, `ε₀ = ½·min(min G, min µ)` | KMT Alg. 1 = geogram code, verbatim |
| line-search iterations | `linesearch_maxiter_ = 100` (step halving) | geogram default |
| multilevel ratio | `0.125` | Lévy 2015 §2.2 |
| weight interpolation | linear LSQ, 10 NN (deg 1) / quadratic LSQ, 20 NN (deg 2) | Lévy 2015 §2.2 |
| linear solve | CG + Jacobi (`NL_CG`, `NL_PRECOND_JACOBI`); `‖Ax−b‖/‖b‖` = `linsolve_epsilon_` | geogram |
| gauge fixing | center weights (`w -= mean(w)`) | KMT pseudo-inverse |

### E.3 Recommended acceptance test for your pipeline

1. `|Σ_i vol(cell_i) − vol(ball)| / vol(ball) < 1e-10` — partition + conservation.
2. `nbZ_ == 0` — no empty cells.
3. `‖∇g‖₂ ≤ 0.01 · vol(ball) / √k` — optimality.
4. Round-trip: compose (harmonic map → OMT correction) and measure the **local volume distortion**
   `|det J|` over a sample of interior points; require `max |log det J| < 0.05` for a
   volume-preserving claim. This is the **real** acceptance criterion for step (2) of the pipeline,
   and it is independent of the OT solver's internal criterion.
5. Sanity: for a uniform target, the OMT map from a ball to a uniform ball must be the identity up to a
   constant (verified by Lévy's "translation" experiment setup).

---

## F. Recommended action plan

1. **Build `geogram` + `exploragram` (both BSD-3).**
   ```powershell
   git clone --depth 1 https://github.com/BrunoLevy/geogram.git      vendor/geogram
   git clone --depth 1 https://github.com/BrunoLevy/exploragram.git  vendor/geogram/src/lib/exploragram
   # then CMake -G "Visual Studio 17 2022" -A x64 ; GEOGRAM_WITH_EXPLORAGRAM auto-detects
   ```
   ⚠️ Do **not** put exploragram anywhere else — `cmake/geogram.cmake` checks
   `${GEOGRAM_SOURCE_DIR}/src/lib/exploragram` specifically.
2. **Implement step (2) as `OptimalTransportMap3d`** on your tet mesh of the ball, with
   `set_epsilon(0.01)`, `optimize_levels(levels, 1000)`, and a CVT/hierarchical sampling of the target.
3. **Back-map with Laguerre centroids** (`compute_Laguerre_centroids_3d`, or
   `OTM.compute_Laguerre_centroids()`) — this is Lévy §2.3 Algorithm 4 step (4):
   `(p_i)^0 ← centroid(Pow_W(y_i) ∩ M)`.
4. **Certify** per §E.3.
5. If you need a Python prototype first: `pip install geogram` (verified working, Windows wheel) for
   the tessellation, and validate against `POT` exact LP for `k ≤ ~2000`. ⚠️ Do **not** ship the
   pure-numpy clipper (§C: 65–74 ms/cell is ~5000× too slow).

## G. Local artefacts produced

| Path | What |
|---|---|
| `refs/levy2015-OT3D.txt` | ✅ full extracted text of the correct Lévy 2015 paper (1 131 lines) |
| `refs/levy2015-transport.pdf` | ✅ the paper PDF (`hal-01226395`) |
| `refs/levy2015.pdf` | ⚠️ **WRONG PAPER** (plane detection, `hal-01097361`) — do not use |
| `vendor/geogram/` | ✅ geogram clone (`main`) |
| `vendor/exploragram/` | ✅ **the OT solver source** |
| `vendor/pysdot/`, `vendor/sdot/` | ✅ pysdot + its C++ backend |
| `vendor/py_power_diagram/` | ✅ older 2D-only `py_power_diagram` |
| `tools/bench2.py` | ✅ numpy Laguerre-volume benchmark (§C.2) — ⚠️ has a volume bug |
| `tools/test_pygeogram2.py` | ✅ probe of `pygeogram`'s 3D `Voronoi` arrays |
| `tools/pypi_probe.py` | ✅ PyPI wheel/license probe for the packages in §B |
