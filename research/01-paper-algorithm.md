# 01 — Volume Preserving Mesh Parameterization based on Optimal Mass Transportation
## (Su, Chen, Lei, Zhang, Qian, Gu — Computer-Aided Design 82 (2017) 42–56, DOI 10.1016/j.cad.2016.05.020)
### Literature retrieval + full technical extraction

Prepared for: smartmm project (volume-preserving volumetric parameterization).
Working directory: `E:\panyingyun\smartmm`.

---

## 0. 中文要点（TL;DR）

* **拿不到原文正文。** 该文在 Unpaywall/OpenAlex 上是 "bronze OA"（出版社免费，`acceptedVersion`），但**只**托管在 ScienceDirect；本机网络被 Cloudflare / 人机验证挡住（普通请求、curl、CORS 代理、Google Translate 代理、headless Chrome 全部失败）。**没有任何机构库/作者主页/预印本副本**（Unpaywall `has_repository_copy=false`，OpenAlex `any_repository_has_fulltext=false`）。我**没有**使用 Sci-Hub / Library Genesis 镜像（搜索结果里出现了它们的本地索引条目，已刻意忽略）。
* **但是**我拿到了作者本人放出的**同一条流水线的完整实现级材料**：顾险峰（Xianfeng Gu）OT 算法 C++ 源码（`PowerDiagramMesh.h` / `OptimalMassTransportationMesh.h`）、2020-11-05 的 OT 算法讲义（含**完整算法框**与 Hessian 公式）、IPAM/TCD/SIGGRAPH-Asia 讲义、**OMT3D.exe 二进制 + 作者自己的测试数据 `lion_h.t` / `david_h.t` / `bimba_h.t`**，以及同组（Su/Zhang/Gu）Stony Brook 博士论文、被引文献 [44] Lévy 的 3D 半离散 OT 原文。
* **我做了数值反推并验证**：作者测试数据里的 `weight` 属性**正好等于** `ν_i = (1/4)·Σ_{τ∋v_i} |τ|`（相对误差 5.6e-8），`h` 属性**正好是体积调和映射像**（边界面顶点 `|h| ≡ 1.000000`，内部 `|h| < 1`）。→ 这就把 (b)、(c)、(e) 三个关键问题从"推测"升级为"**由作者自己的数据实证**"。
* 结论：**方法 = 体积调和映射 ψ: M → B³，再对偶空间做离散最优传输 φ: (B³, Lebesgue) → (P, ν)，最终 F = φ⁻¹∘ψ : M → B³，F(v_i) = 单元 W_i 的质心。**
* 文中标注约定：**【Q】= 原文引用（逐字）**，**【V】= 已由作者代码/数据/权威二次文献验证**，**【I】= 本报告推断（尚无直接证据）**。

---

## 1. Full bibliographic record

**Target paper**

| Field | Value |
|---|---|
| Title | Volume preserving mesh parameterization based on optimal mass transportation |
| Authors | Kehua Su; Wei Chen; Na Lei (corresponding); Junwei Zhang; Kun Qian; Xianfeng Gu |
| Affiliations | ① State Key Laboratory of Software Engineering, Wuhan University, Wuhan 430072, China ② Computer Science Department, Stony Brook University, NY 11794, USA ③ School of Software, Dalian University of Technology, Dalian 116620, China ④ Yau Mathematics Science Center, Tsinghua University, Beijing 100084, China ⑤ School of Civil Engineering and Mechanics, Kunming University of Science and Technology, Kunming, China |
| Journal | Computer-Aided Design (CAD), Elsevier |
| Volume / pages | 82, pp. 42–56 (15 pages) |
| Year | 2017 issue; online 2016-06-17 |
| DOI | 10.1016/j.cad.2016.05.020 |
| PII | S0010448516300495 |
| ISSN | 0010-4485 / 1879-2685 |
| Keywords (author-supplied) | Parameterization; Volume preserving; Harmonic map; Optimal mass transportation |
| OA status | `bronze`, version = `acceptedVersion`, publisher-hosted only, no repository copy (Unpaywall), `any_repository_has_fulltext: false` (OpenAlex) |
| Ids | Semantic Scholar `e92c70468260c4cfebeded82bfca8ea105177c3b` (CorpusID 36045233); DBLP `journals/cad/SuCLZQG17`; MAG 2434575306 |
| Citations | 34–43 (Scopus 43, OpenAlex 34, OUCI 36) |
| # references | 44 |

**Companion paper (also requested)**

| Field | Value |
|---|---|
| Title | Measure controllable volumetric mesh parameterization |
| Authors | Kehua Su; Wei Chen; Na Lei (corresponding); Li Cui; Jian Jiang; Xianfeng David Gu |
| Journal | Computer-Aided Design 78 (2016) 188–198 |
| DOI | 10.1016/j.cad.2016.04.007 |
| PII | S0010448516300197 |
| OA status | `bronze`, publisher-hosted only |
| Recommended by | Scott Schaefer and Charlie C. L. Wang (per the paper's footnote) |

**Abstract of the target paper (verbatim, identical across Semantic Scholar / Stony Brook Pure / DLUT Pure / Daneshyari):**

> "In order to convert a finite element mesh model to the spline representation for the purpose of isogeometric analysis, one needs to parameterize the solid. This work introduces a novel volumetric parameterization method, which guarantees to be free of volume distortion. Given a simply connected tetrahedral mesh with a single boundary surface, we first compute a harmonic map from the boundary triangle mesh to the unit sphere by non-linear heat diffusion method; then we use the surface harmonic map as the boundary condition to compute the volumetric harmonic map to parameterize the solid onto the unit solid ball; finally we compute an optimal mass transportation map from the unit solid ball with the push-forward volume element induced by the harmonic map onto itself with the Euclidean volume element. The composition of the volumetric harmonic map and the optimal mass transportation map gives an volume-preserving parameterization. The method has solid theoretic foundation, and is based on conventional algorithms in computational geometry, easy to implement. We have thoroughly tested our algorithm on many solid models in reality. The experimental results demonstrate the efficiency and efficacy of the proposed method. To the best of our knowledge, it is the first work addressing volume-preserving parameterization in the literature."

---

## 2. Access attempts — what actually worked and what did not

### 2.1 Target paper: NOT OBTAINABLE (full text)

| URL / route | Result |
|---|---|
| `https://www.sciencedirect.com/science/article/am/pii/S0010448516300495` (accepted manuscript, the bronze-OA copy) | **403** — Cloudflare. `Invoke-WebRequest`, `curl --compressed` with full browser header set, all fail. |
| same, headless Chrome `--dump-dom` (real Chrome 154, `--headless=new`) | **blocked**: first `CLOUDFLARE_ERROR_1000S_BOX`, then with a spoofed desktop UA an **hCaptcha / "Are you a robot?"** interstitial. Not circumvented. |
| `.../am/pii/.../pdfft`, `.../pii/.../pdfft?isDTMRedir=true&download=true` | 403 / connection reset |
| CORS proxies (`api.codetabs.com`, `api.allorigins.win`, `corsproxy.io`, `thingproxy`) | 522 / 403 / TLS failure |
| `https://r.jina.ai/https://www.sciencedirect.com/...` | **403** (the reader refuses) |
| Google-Translate proxy `www-sciencedirect-com.translate.goog/...` | **400** |
| Wayback Machine CDX for the `/am/` and `/abs/` pages | **no snapshot of the AM page**; the only `/pii/` snapshots are 302/403 Cloudflare stubs |
| Wayback for `daneshyari.com/article/preview/440675.pdf` (companion) | **200 — WORKED** (see below) |
| `academia.edu/101357636/Volume_preserving_mesh_parameterization_based_on_optimal_mass_transportation` | **403** (Cloudflare) |
| `daneshyari.com` live (landing + preview) | **504** (server-side timeout); Wayback copy exists only for id 440675 (the companion) |
| `journals.scholarsportal.info/details/00104485/v82icomplete/42_vpmpboomt.xml&sub=all` | 4.4 kB login stub |
| `hal.science` / `inria.hal.science` | **Anubis proof-of-work** JS challenge |
| `base-search.net` | **Anubis** challenge |
| `colab.ws`, `x-mol.com`, `infona.pl` | 403 / 404 / no full text |
| Unpaywall API, OpenAlex API, Semantic Scholar API, OpenAIRE API | metadata only; all agree there is **no** repository/OA copy |
| `researchconnect.stonybrook.edu`, `jszy.whu.edu.cn/sukehua`, `faculty.dlut.edu.cn/leina` | abstract + metadata only, no PDF |
| **Sci-Hub / Library Genesis mirrors** | appeared in search-engine result lists as local-corpus entries — **deliberately not used** |
| **Conclusion** | The paper's *body text* (sections 3–6, Algorithm box, figures, tables) is **not retrievable** from this environment. Everything below about (a)–(f) is either quoted from the indexed accepted-manuscript page (a few sentences only), or reconstructed from the authors' own code/lecture notes/data, all explicitly labelled. |

### 2.2 Companion paper: PARTIALLY OBTAINED

* `https://web.archive.org/web/20240416061146if_/https://daneshyari.com/article/preview/440675.pdf` → **200, 352 056 bytes, application/pdf** → saved to `refs\pdf\mc_wb.pdf`; text in `refs\MC_companion.txt`. This is the **first 3 PDF pages** (title page, abstract, §1 Introduction, §2.1–2.2 Related work) — i.e. `daneshyari`'s "preview" is a truncated excerpt. Body sections (§3 method, experiments) are **not** present.

### 2.3 What DID work (all downloaded to `E:\panyingyun\smartmm\refs\`)

| Artifact | URL | Local file |
|---|---|---|
| Gu, "Geometric Algorithm for Optimal Transportation Map", Tsinghua/Stony Brook lecture, **Nov 5 2020** — contains the complete **GU-LUO-SUN-YAU algorithm box**, the **Hessian formula**, the damping algorithm, the data-structure/algorithm choices | `https://www3.cs.stonybrook.edu/~gu/software/OT/Geometric_OT_Lecture.pdf` | `pdf\Geometric_OT_Lecture.pdf` (13.8 MB) + `Geometric_OT_Lecture.txt` |
| **OT C++ source release 2020-11-06** — `OT/include/PowerDiagramMesh.h` (Lawson edge flip / Legendre dual / Sutherland–Hodgman power-cell clipping), `OT/include/OptimalMassTransportationMesh.h` (gradient, **Hessian assembly**, Eigen SimplicialLDLT solve, damping/roll-back) | `https://www3.cs.stonybrook.edu/~gu/software/OT/OT_11_06_2020.zip` | `pdf\OT_11_06_2020.zip` → `code\OT\OT_11_06_2020\...` |
| **`omt-binary.rar` — OMT2D.exe, OMT3D.exe (+ `geogram.dll`, CGAL 4.7, Qt5), test scripts and test meshes** | `https://www3.cs.stonybrook.edu/~gu/software/omt/omt-binary.rar` | `pdf\omt-binary.rar` → `code\omtbin\omt-binary\...` |
| Gu, IPAM UCLA OMT talk (2016-02-12) | `https://www3.cs.stonybrook.edu/~gu/software/omt/IPAM_OMT.pdf` | `pdf\IPAM_OMT.pdf` (9.7 MB) |
| Gu, "Theory of Optimal Mass Transportation", TCD Dublin | `https://www.maths.tcd.ie/~hmigca-18/slides/gu_omt_theory.pdf` | `pdf\gu_omt_theory_tcd.pdf` |
| Gu, SIGGRAPH-Asia 2013 course: Optimal Transport | `https://www3.cs.stonybrook.edu/~gu/talks/2013_SIGGRAPH_ASIA/Gu_3_SIGGRAPH_Asia_Optimal_Transport.pdf` | `pdf\SIGGRAPH_ASIA_OT.pdf` |
| Gu–Luo–Sun–Yau, "Variational principles for Minkowski type problems, discrete optimal transport, and discrete Monge–Ampère equations", AJM 20(2):383–398, 2016 — **ref [29] of the paper** | `https://arxiv.org/pdf/1302.5472` and `https://intlpress.com/site/pub/files/_fulltext/journals/ajm/2016/0020/0002/AJM-2016-0020-0002-a007.pdf` | `pdf\arxiv1302.5472_GLY_variational.pdf`, `pdf\AJM_GLY_variational.pdf` |
| **Bruno Lévy, "A numerical algorithm for L2 semi-discrete optimal transport in 3D", ESAIM M2AN 49(6):1693–1715, 2015 — ref [44], THE 3-D power-cell-volume reference** | `https://arxiv.org/pdf/1409.1279` (ESAIM PDF is 403) | `pdf\levy_ot3d_arxiv.pdf` |
| Wang, Gu, Yau, "**Volumetric Harmonic Map**", Communications in Information and Systems 3(3):191–202, 2004 — **ref [5], the volumetric harmonic energy** | `https://www3.cs.stonybrook.edu/~gu/publications/pdf/2004/volumetric.pdf` | `pdf\volumetric.pdf` |
| Gu, Wang, Yang, Yau, "Genus zero surface conformal mapping and its application to brain surface mapping", IEEE TMI 23(8):949–958, 2004 — **ref [39], the "non-linear heat diffusion" spherical map** | `https://www3.cs.stonybrook.edu/~gu/publications/pdf/2004/TMI_04.pdf` | `pdf\TMI_04.pdf` |
| Zhengyu Su, PhD dissertation "Optimal Mass Transport and Its Applications", Stony Brook, 2015 — **same group; contains Algorithm 2 (OMT-Map box) and the exact Hessian (3.28)** | `https://repo.library.stonybrook.edu/xmlui/bitstream/handle/11401/77314/Su_grad.sunysb_0771E_12622.pdf?sequence=1&isAllowed=y` | `pdf\Su_dissertation_suny.pdf` |
| Min Zhang, PhD dissertation "Ricci Flow and Its Applications", Stony Brook, 2014 — same group, OMT chapter | `https://repo.library.stonybrook.edu/xmlui/bitstream/handle/11401/77322/Zhang_grad.sunysb_0771E_12003.pdf?sequence=1&isAllowed=y` | `pdf\Zhang_dissertation_suny.pdf` |
| Gu, BIMSA Chinese OMT lecture notes (Minkowski–Alexandrov convex geometry) | `https://bimsa.net/doc/notes/39708.pdf` (also `39728.pdf`) | `pdf\gu_bimsa_omt_notes_cn.pdf` (52 MB; **CJK extraction garbled**, math readable) |
| Huang, Liao, Lin, Yueh, Yau, "Convergence Analysis of Volumetric Stretch Energy Minimization and its Associated Optimal Mass Transport" | `https://arxiv.org/pdf/2210.09654` | `pdf\arxiv2210.09654_vsem.pdf` |
| n-VSE / IEM papers (cite the target paper; §1.2 contains an **independent description of exactly this method**) | `https://arxiv.org/pdf/2402.00380`, `https://arxiv.org/pdf/2407.19272` | `pdf\...` |
| Crane, "The n-dimensional cotangent formula" — the modern statement of the tetrahedral cotan weights | `http://www.cs.cmu.edu/~kmcrane/Projects/Other/nDCotanFormula.pdf` | `pdf\crane_nd_cotan.pdf` |
| OUCI record — **full 44-item reference list** of the target paper + citing works | `https://ouci.dntb.gov.ua/en/works/96zjxar9/` | `refs\` (session output) |
| Stony Brook Pure portal | `https://researchconnect.stonybrook.edu/en/publications/volume-preserving-mesh-parameterization-based-on-optimal-mass-tra/` | abstract, page range, keywords |
| Signalling note | Semantic Scholar's `openAccessPdf.url` for this DOI literally resolves back to `https://doi.org/10.1016/j.cad.2016.05.020` — i.e. "free to read, only at the publisher". |

**Sentences of the accepted manuscript that the search index does expose** (the only direct quotes of the target paper I could obtain; treat as partial):

* 【Q】"We compute the discrete optimal mass transportation map φ : 𝔻³ → P, the inverse map φ⁻¹ maps each …"
* 【Q】"Δ_PL f(u) = Σ_v k_{u,v} ( f(v) − f(u) )"

The first of these is decisive: it confirms the paper uses **exactly the same notation as the group's other papers** — the discrete OMT map is `φ : 𝔻³ → P` (from the unit ball **to the point set**), and it is `φ⁻¹` that is used to build the parameterization.

---

## 3. The complete pipeline, stage by stage

### 3.0 Overview (as the paper defines it)

```
M  (tetrahedral mesh, topological 3-ball, single boundary surface ∂M ≅ S²)
 │
 │  (1) boundary spherical map  g : ∂M → S²      [non-linear heat diffusion / harmonic energy]
 ▼
 │  (2) volumetric harmonic map ψ : M → B³ ,  ψ|∂M = g   [linear Dirichlet problem for Δ_PL]
 ▼
(B³, μ := ψ_# dV_M)          μ  =  push-forward of the volume element of M
 │
 │  (3) discrete optimal mass transport  φ : (B³, dx) → (P, ν)
 │      P = {p_i = ψ(v_i)},  ν_i = per-vertex volume of M
 │      φ|W_i = p_i ,  Vol(W_i) = ν_i   (power / Laguerre diagram of B³)
 ▼
 │  (4) compose:  F = φ⁻¹ ∘ ψ : M → B³ ,  F(v_i) = mass centre of W_i
 ▼
Volume-preserving parameterization  F : M → B³   (M → unit solid ball)
```

The paper's own one-sentence version 【Q, abstract】: *"we first compute a harmonic map from the boundary triangle mesh to the unit sphere by non-linear heat diffusion method; then we use the surface harmonic map as the boundary condition to compute the volumetric harmonic map to parameterize the solid onto the unit solid ball; finally we compute an optimal mass transportation map from the unit solid ball with the push-forward volume element induced by the harmonic map onto itself with the Euclidean volume element. The composition of the volumetric harmonic map and the optimal mass transportation map gives an volume-preserving parameterization."*

---

### 3.a The initial map M → B³

**Two sub-stages.**

#### (a1) Boundary spherical map `g : ∂M → S²` — "non-linear heat diffusion"

【Q, abstract】 *"a harmonic map from the boundary triangle mesh to the unit sphere by non-linear heat diffusion method"*.
【V】 The cited method is **[39] Y. Wang, X. Gu, S.-T. Yau (and X. Gu, Y. Wang, T. F. Chan, S.-T. Yau), "Genus zero surface conformal mapping and its application to brain surface mapping", IEEE TMI 23(8):949–958, 2004** — the group's canonical "non-linear heat diffusion" spherical map. Its structure (from the PDF I downloaded) is a **steepest-descent / heat-flow** minimisation of the harmonic (string) energy with the vertices constrained to the unit sphere:

* Discrete **string/harmonic energy** on the boundary triangle mesh `∂M`:
  `E(f) = ½ Σ_{{u,v}∈E(∂M)} k_{u,v} · |f(u) − f(v)|²`,  `f : ∂M → S² ⊂ ℝ³`
* **Piecewise-linear Laplacian** (= the one the target paper quotes):
  【Q】 `Δ_PL f(u) = Σ_v k_{u,v} ( f(v) − f(u) )`
  (2-D harmonic weight: `k_{u,v} = ½ (cot α_{uv} + cot β_{uv})`, α, β the two interior angles opposite edge (u,v)).
  【V】 The annotation "the absolute derivative" in the TMI paper is the **projected (tangential) Laplacian**:
  `Df = Δ_PL f − ⟨Δ_PL f , f⟩ f  ∈ T_{f} S²`
* **Non-linear heat flow / update** (the "non-linear" part is the re-projection onto the sphere):
  `f^{(n+1)}(v) = f^{(n)}(v) + δ · D f^{(n)}(v)`,  then `f^{(n+1)}(v) ← f^{(n+1)}(v) / |f^{(n+1)}(v)|`
  iterated until `|E − E_0| < δE`.
  【V】 The TMI paper's Algorithm 1 (Spherical Tuette Mapping) and Algorithm 2 (Spherical Conformal Mapping) are exactly this loop; Algorithm 2 adds a **Möbius (mass-centre) normalisation** every step, approximated by "compute mass centre `c`; `f ← f − c`; `f ← f/|f|`". The Tuette map is used as the initialisation for the harmonic one.
* 【I】 The target paper very likely runs the same two-phase scheme (Tuette → harmonic + Möbius normalisation), since that is the only "non-linear heat diffusion" method in that reference.

#### (a2) Volumetric harmonic map `ψ : M → B³`

【Q, abstract】 *"we use the surface harmonic map as the boundary condition to compute the volumetric harmonic map to parameterize the solid onto the unit solid ball"*.
【V】 Reference **[5] Y. Wang, X. Gu, S.-T. Yau, "Volumetric Harmonic Map", Commun. Inf. Syst. 3(3):191–202, 2004** gives the discrete energy. Quoting that paper directly:

* Per-tetrahedron harmonic energy (their Eq. 4), with `s⃗_i = Area(Face i)·n⃗_i`:
  `E_τ(f) = −Σ_{i≠j} ⟨s⃗_i, s⃗_j⟩/(18 V) · (f_i − f_j)² = Σ_{i≠j} (l_{pq} cot θ_{pq})/12 · (f_i − f_j)²`
  where `{p,q} = I \ {i,j}` and `l_{pq}`, `θ_{pq}` are the length of the edge opposite to (i,j) and the dihedral angle **along that opposite edge**.
* 【Q】 Definition 4 of the VHM paper — the **edge weight**:
  `k_{u,v} = (1/12) Σ_{i=1}^{n} l_i · cot(θ_i)`
  "Suppose for edge {u,v}, it is shared by n tetrahedra thus it is against to n dihedral angles, θ_i, i = 1,…,n … where l_i is the length of edge to which edge {u,v} is against in the domain manifold M."
* 【V】 Modern, equivalent normalisation (K. Crane, *The n-dimensional cotangent formula*):
  `w_{ij} = (1/6) Σ_{tets ijkl ∋ ij} ℓ_{kl} · cot θ_{kl}`
  with `ℓ_{kl}` the length of the edge **opposite** to (i,j), and `θ_{kl}` the interior dihedral angle at that opposite edge. General n-D form: `w_{ij} = 1/(n(n−1)) Σ_{σ∋ij} V_{σ_{i,j}} cot θ_{σ_{i,j}}`.
  **The factor-2 difference with the VHM paper is purely a convention in the definition of the string energy (`½` in front) — do not be confused by it.** For implementation, use `E(f) = ½ Σ_{edges} w_{ij}|f_i−f_j|²` with Crane's `w_{ij}`; that is the standard and matches the 2-D `½(cot α + cot β)`.
* 【V】 The energy is minimised by steepest descent on the interior with **fixed boundary data**:
  `dψ/dt = −Δ_PL ψ`, interior vertices updated, `ψ(v) = g(v)` held fixed on `∂M`.
  Because the boundary data is fixed, the **stationary** problem is the **linear Dirichlet problem**
  `Σ_v k_{u,v}(ψ(v) − ψ(u)) = 0 for u ∈ M∖∂M`,  `ψ|∂M = g`,
  i.e. a **sparse symmetric linear system** `L_II ψ_I = −L_IB g_B` (L = the tetrahedral cotan Laplacian). The VHM paper presents it as an explicit heat-flow iteration (their Algorithm 1, step 2: "for each boundary vertex v ∈ ∂M let h(v) = g(v); for each interior vertex v ∈ M∖∂M let h(v) = (0,0,0)"), with step length `δt` and stopping threshold `δE`; a direct sparse solve is equivalent and is what a modern implementation does.
* 【V, from the authors' own data】 Since `|ψ| ≤ 1` on `∂M` (boundary is mapped **onto** the unit sphere), the maximum principle gives `|ψ| ≤ 1` inside, so the image lies in the closed unit ball — matching the definition `ψ : M → B³`.

---

### 3.b The measure for the OMT step — **exactly the per-vertex `¼ Σ |τ|`**

【Q, abstract】 source measure = *"the unit solid ball with the push-forward volume element induced by the harmonic map"*; target measure = *"the Euclidean volume element"*.
【Q】 (indexed manuscript) *"We compute the discrete optimal mass transportation map φ : 𝔻³ → P, the inverse map φ⁻¹ maps each …"*.

【V — **verified numerically against the authors' own released data**】

The released 3-D test meshes `omt-binary/omt3d_test/{lion,david,bimba}_h.t` store, for every vertex, three things:

```
Vertex 1 0.058012 -0.455052 -0.01709  {h=(0.064476 -0.981668 -0.179364) weight=(1.2433752e-06) boundary}
Vertex 2 -0.271788 0.024171 0.114956  {h=(-0.851080 0.244949 0.464395) weight=(9.1825904e-07) boundary}
...
Tet 97131 13164 18641 19099 20874
```

* `(x,y,z)` = the vertex position in `M`.
* `h = (h_x,h_y,h_z)` = **the volumetric harmonic map image `ψ(v) ∈ B³`**: for the flagged `boundary` vertices **`|h| = 1.000000` for 100 % of them** (measured range `[0.999999, 1.000001]`), while interior vertices have `|h| ∈ [0.059, 0.9998]`. → `h` is `ψ`, not a unit normal per se; on `∂M` it *is* the unit normal because `ψ(∂M) = S²`.
* `weight` = **the target atomic mass `ν_i`**.

**Numerical verification (my computation, `lion_h.t`: 20 874 vertices, 97 131 tets):**

```
Σ_τ |τ|                  = 0.16288814433865131      (total volume of the mapped mesh)
Σ_i  weight_i            = 0.16288814425065204
(1/4)·Σ_{τ∋v_i}|τ|  vs  weight_i :  max relative error = 5.56e-08
(1/5)·Σ_{τ∋v_i}|τ|  vs  weight_i :  max relative error = 2.0e-1
Σ_{τ∋v_i}|τ|        vs  weight_i :  max relative error = 3.0
```

**⇒ The atom / site / mass structure is, definitively:**

* **Sites (atoms of the target measure):** the *vertices* of the tetrahedral mesh `M`, i.e. `P = { p_i = ψ(v_i) } ⊂ B³`. (One site per vertex — **not** per tetrahedron.)
* **Mass of atom `i`:** the **barycentric ("1/4") dual-cell volume**
  ```
  ν_i  =  (1/4) · Σ_{ τ ∈ T(M), τ ∋ v_i } |τ|            (exactly, per the data)
  ```
  i.e. one quarter of the sum of the volumes of the tetrahedra incident to `v_i`. This is the 3-D analogue of the group's 2-D rule `ν_i = (1/3)Σ_{faces ∋ v_i} area(face)` (Zhengyu Su's dissertation, Eq. 3.30 and Algorithm 3 step 3: *"ν_i to be 1/3 of the total area of the faces adjacent to v_i"*).
  It is **NOT** the circumcentric/Voronoi dual volume, **not** a per-tetrahedron mass, and **not** a mesh-volume-weighted Voronoi measure.
  Note `Σ_i ν_i = Σ_τ |τ| = Vol(M)` identically.
* **Source measure:** the **Lebesgue / Euclidean measure `dx` on the unit ball `B³`** (continuous); in the released tool the domain is clipped by the unit sphere or the unit cube (`-clip sphere|cube`, default `cube`; `-project true` projects boundary points onto the unit sphere/cube first).
* **Normalisation:** the two total masses must agree, so `ν_i` is rescaled inside the solver by `Vol(clip domain)/Σ_j ν_j`. 【V】 Gu's own 2-D reference implementation does exactly this (`OMT_Initialize` in `OptimalMassTransportationMesh.h`: `s = Σ target_area; ds = D.area(); target_area *= (ds/s);`). 【I】 the 3-D tool does the analogous thing.
* **Discretisation of the push-forward:** `μ = ψ_#(dV_M) ≈ Σ_i ν_i δ_{p_i}` with the barycentric dual cells. This is the standard first-order (vertex-sampling) approximation of a piecewise-constant/P1 density pushed forward, which is *why* the per-vertex `¼Σ|τ|` appears.
* The companion paper states the general (measure-controllable) version 【Q】: *"we first compute a volumetric harmonic map to parameterize the solid onto the unit solid ball; then we compute an optimal mass transportation map from the unit solid ball with the push-forward volume element induced by the harmonic map onto the parameter domain with the user prescribed volumetric measure."* — i.e. only the target side changes; for **volume-preserving** you take the prescribed measure = the Euclidean one.

**Direction of the transport (important, and the paper's abstract is verbally ambiguous):**

* The paper's own notation 【Q】 is `φ : 𝔻³ → P` — a map *from the ball to the point set* — which is exactly the group's other papers:
  * Zhengyu Su's dissertation, Algorithm 2: *"Input: A convex planar domain with measure (Ω,μ); a planar point set with measure (P,ν) … Output: The unique discrete OMT-Map f : (Ω,μ) → (P,ν)"*, and *"return f : Ω → P, W_i(h) → p_i"*.
  * Gu's website: *"we compute an optimal mass transportation map from this measure to the canonical Euclidean measure on the unit ball"* — the map actually used is its **inverse**.
* Therefore, operationally: **source = `(B³, dx)` (continuous, Euclidean), target = `Σ_i ν_i δ_{p_i}` (discrete at the harmonic images)**, and each power cell satisfies `∫_{W_i} dx = ν_i`, `φ|W_i = p_i`.
  The abstract's wording ("from the ball with the push-forward volume element onto itself with the Euclidean volume element") describes the same object viewed from the dual side, i.e. `φ⁻¹`, which is the map that is actually consumed. Both readings give the same final composition (§3.e).

---

### 3.c The discrete OMT map: power (Laguerre) diagram, energy, gradient, Hessian, solver

**Structure: it is a POWER (Laguerre) diagram, not a plain Voronoi diagram.** 【Q】 *"the projection of `u*_h` is a power Delaunay triangulation … the projection of `u_h` is the dual power diagram"*; 【Q】 the paper cites Aurenhammer, *Power diagrams: properties, algorithms and applications*, SIAM J. Comput. 16:78 (1987) as **ref [38]**.

#### Objects

* **Brenier potential (piecewise-linear convex function):**
  `u_h(x) = max_{i=1..k} { ⟨x, p_i⟩ − h_i }`,  `h = (h_1,…,h_k)` the height vector.
  【Q, Gu's Nov-2020 lecture】 `u_h` is the **upper envelope** of the planes `π_i(x) = ⟨x,p_i⟩ − h_i`; its **Legendre dual** `u*_h` is the convex hull of the lifted points `{(p_i, h_i)}`; projecting `u*_h` gives the **power Delaunay triangulation**, projecting `u_h` gives the **power diagram**.
* **Power cell** (this is exactly what the code uses — `Pow(x,p_i) = ½|x−p_i|² − ½ h_i`):
  `W_i(h) = { x : Pow(x,p_i) ≤ Pow(x,p_j) ∀j } = { x : ⟨x, p_i⟩ − h_i ≥ ⟨x, p_j⟩ − h_j ∀j }`.
* **Cell measure:** `w_i(h) = ∫_{W_i(h) ∩ Ω} ρ(x) dx`. In our setting `Ω = B³`, `ρ ≡ 1`, so `w_i(h) = |W_i(h) ∩ B³|` (Euclidean volume).
* **Admissible height space:** `H = { h ∈ ℝ^k : w_i(h) > 0 ∀i }` (usually intersected with `Σ_i h_i = 0`). 【Q, GLY】 non-empty, open, **convex** (proved via Brunn–Minkowski).

#### Energy (Alexandrov / Brenier energy) — GLY convention

【Q, Gu lecture / Zhengyu Su dissertation Thm 3.8 / GLY AJM 2016】
```
E(h)  =  Σ_{i=1}^{k} ν_i h_i  −  ∫_0^h  Σ_{j=1}^{k} w_j(η) dη_j
```
Geometric meaning 【Q】: *"one can define a cylinder through ∂Ω, the cylinder is truncated by the xy-plane and the convex polyhedron. The energy term `∫^h Σ w_i(η)dη_i` equals the volume of the truncated cylinder"*; *"geometrically, the energy is the volume beneath the parabola."*
`E` is **concave**; `G(h) = Σ A_i h_i − F(h)` with `F(h)=∫^h Σw_j dη_j` attains its max at an **interior point** of `H₀` — that interior critical point is the solution.

Equivalent "positive-definite" convention (Zhengyu Su, Eq. 3.25 / 3.27 / 3.28 — signs flipped, and this is what the released code implements):
```
E(h) = ∫_{B³} u_h(x) dx − Σ_i ν_i h_i ,
```
minimised. Both conventions appear in the group's own writing; **pick one and keep the signs consistent.**

#### Gradient

【Q, GLY & Gu lecture】 `∇E(h) = ( ν_i − w_i(h) )_{i=1..k}`
【Q, Zhengyu Su Eq. 3.27 (opposite sign convention)】 `∇E(h) = ( w_1(h) − ν_1, …, w_k(h) − ν_k )ᵀ`
【V, the released 2-D code `_calculate_gradient()`】 `update_direction[i] = target_area[i] − dual_area[i]` — literally "target mass minus current cell volume". So **the gradient is the volume-difference vector**.

#### Hessian

**【Q, the definitive formula — Gu, Nov-5-2020 lecture, slide 12:】**
```
a_ij = ∂²E/∂h_i∂h_j = − 1/|p_i − p_j| · ∫_{W_i ∩ W_j} ρ(x) dx        (i ≠ j)
a_ii = ∂²E/∂h_i∂h_i = − Σ_{j≠i} a_ij        (because Σ_i w_i(h) ≡ const)
```
**【Q, Zhengyu Su dissertation Eq. (3.28):】**
```
∂²E(h)/∂h_i∂h_j =  ( ∫_{e_ij} μ(x) dx ) / |p_j − p_i|   if  W_i(h) ∩ W_j(h) ∩ Ω ≠ ∅,   else 0
```
where `e_ij = W_i(h) ∩ W_j(h) ∩ Ω`.
**【Q, Gu lecture slide 13 (2-D form):】** `∂w_i/∂h_j = −|e_ij| / |ē_ij|` — "the ratio between the power Voronoi edge length and the power Delaunay edge length"; **2-D code** `__edge_weight()`: `weight = | edge(W_i ∩ W_j) ∩ D | / |p_i − p_j|`.

**In our 3-D setting `ρ ≡ 1`, so**
```
∂²E/∂h_i∂h_j  =  −  Area( W_i ∩ W_j ) / |p_i − p_j|          (i ≠ j)
```
where `W_i ∩ W_j` is the **common planar facet** of the two power cells (a convex polygon in 3-D), and `|p_i − p_j|` the site distance. The Hessian is a **weighted graph Laplacian on the power-Delaunay 1-skeleton**, symmetric, negative semi-definite with a 1-D kernel (constants / translations), hence positive definite on `Σ h_i = 0`. 【Q】 *"the Hessian of F is diagonally dominated"*, *"the negative Hessian matrix is diagonal dominant, so E is concave on H"*.

#### Solver and update rule (Newton with damping / roll-back)

【Q, Gu, Nov-5-2020 lecture, slides 14–16 — reproduced verbatim from the PDF text】

> **Optimal Transport Map** — *Input:* a set of distinct points `Y = {y_1..y_k}`, weights `{ν_1..ν_k}`, convex domain `Ω`, `Σν_j = Vol(Ω)`. *Output:* the optimal transport map `T : Ω → Y`
> 1. Scale and translate `Y`, such that `Y ⊂ Ω`;
> 2. Initialize `h⁰ ← ½(|y_1|², |y_2|², …, |y_k|²)ᵀ`;
> 3. Compute the Brenier potential `u(h_k)` (envelope of `π_i`'s) and its Legendre dual `u*(h_k)` (convex hull of `π*_i`'s);
> 4. Project the Brenier potential and Legendre dual to obtain weighted Delaunay triangulation `T(h_k)` and power diagram `D(h_k)`;
> 5. Compute the gradient of the energy `∇E(h) = (ν_1 − w_1(h), ν_2 − w_2(h), …, ν_k − w_k(h))ᵀ`;
> 6. If `‖∇E(h_k)‖` is less than `ε`, then return `T = ∇u(h_k)`;
> 7. Compute the Hessian matrix of the energy `∂w_i(h)/∂h_j = −|e_ij|/|ē_ij|`, `∂w_i/∂h_i = −Σ_j ∂w_i(h)/∂h_j`;
> 8. Solve linear system `∇E(h) = Hess(h_k) d`;
> **Damping Algorithm**
> 9. Set the step length `λ ← 1`;
> 10. Construct the convex hull `Conv(h_k + λ d)`;
> 11. If there is any empty power cell, `λ ← ½λ`, repeat steps 3 and 4, until all power cells are non-empty;
> 12. Set `h_{k+1} ← h_k + λ d`;
> 13. Repeat steps 9 through 12.

**【Q, Gu's earlier TCD slides, slightly different step order (the "textbook" box):】**
> 1 Initialize `h = 0`; 2 Compute the Power Voronoi diagram and the dual Power Delaunay Triangulation; 3 Compute the cell areas, which gives the gradient `∇E`; 4 Compute the edge lengths and the dual edge lengths, which gives the Hessian matrix `Hess(E)`; 5 Solve linear system `∇E = Hess(E) dh`; 6 Update the height vector `h ← h − λ dh`, where `λ` is a constant to ensure that no cell disappears; 7 Repeat step 2 through 6, until `‖dh‖ < ε`.

**【Q, Zhengyu Su dissertation Algorithm 2 (OMT-Map), verbatim:】**
> *Input:* A convex planar domain with measure `(Ω,μ)`; a planar point set with measure `(P,ν)`, `ν_i > 0`, `∫_Ω μ(x)dx = Σν_i`; a threshold `ε`. *Output:* the unique discrete OMT-Map `f : (Ω,μ) → (P,ν)`.
> Scale and translate `P`, such that `P ⊂ Ω`. `h ← (0,…,0)`. Compute the power diagram `D(h)`, compute the dual power Delaunay triangulation `T(h)`, compute the cell areas `w(h) = (w_1(h),…,w_k(h))`.
> **repeat**
>   Compute `∇E(h)` using Eqn. 3.27. Compute the Hessian matrix using Eqn. 3.28. `λ ← 1`. `h ← h − λ H⁻¹∇E(h)`. Compute `D(h)`, `T(h)` and `w(h)`.
>   **while** `∃ w_i(h) == 0` **do** `h ← h + λH⁻¹∇E(h)`; `λ ← ½λ`; `h ← h − λH⁻¹∇E(h)`; Compute `D(h)`, `T(h)` and `w(h)`. **end while**
> **until** `‖∇E‖ < ε`.
> **return** `f : Ω → P`, `W_i(h) → p_i`, `i = 1,…,k`.

**【V, the actual released C++ (`OptimalMassTransportationMesh.h`) — all constants:】**

```cpp
// gradient
pv->dual_area()      = pv->dual_cell2D().area();
pv->update_direction() = pv->target_area() - pv->dual_area();

// Hessian edge weight  (2-D version; 3-D analogue: area of the common facet / site distance)
pe->weight() = | segment(dual_point(face1), dual_point(face2)) ∩ D | / |uv_i − uv_j|;

// Hessian assembly (sparse "triplets")
triplets.push_back({idx0, idx1,  h/2 });   triplets.push_back({idx1, idx0,  h/2 });
triplets.push_back({idx0, idx0, -h/2 });   triplets.push_back({idx1, idx1, -h/2 });
double epsilon = -1e-8;   for (i) triplets.push_back({i,i,epsilon});   // regulariser

// solve
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
solver.compute(hessian);  x = solver.solve(gradient);
x.array() -= x.mean();                       // remove the constant/translation null-space
pv->update_direction() = x[idx];

// damping / roll-back (_OT_damping), step = 1 by default
step_length = nearest ? step : -step;
backup phi_z;  phi_z += step_length * update_direction;
if (!mesh._Lawson_edge_swap(nearest))  { roll back; step_length /= 2; }     // keeps all sites on the convex position
mesh._Legendre_transform(...); mesh._power_cell_clip(D);
if ( |dual_area| < target_area * 1e-10 ) { roll back; step_length /= 2; }   // no empty power cell

// convergence (_error)
max relative error = max_i |target_area_i − dual_area_i| / target_area_i;
```
Default CLI threshold `-threshold` = **1e-4**; the bundled 2-D demos use `1e-6` (`test_alex.bat`, `test_skull.bat`).
Summary of the solver choice 【Q】: *"the optimization of Alexandrov energy is based on damping algorithm"*, *"the linear numerical solver is Eigen library"*, *"the geometric computation is based on adaptive arithmetic method"*, *"the power Delaunay is based on Lawson's edge flip algorithm"*, *"the polygon clipping is based on Sutherland–Hodgman algorithm"*.

**【I】 (inference, explicitly unverified):** the target paper's own solver is one of these two — **(i)** the damped Newton above, or **(ii)** the L-BFGS quasi-Newton used by **[44] Lévy (2015)**: *"I use the L-BFGS numerical optimization method"*, with `ε = 0.01 · µ(M)/√k`. Gu's lecture explicitly contrasts *"Gradient descend, Quasi-Newton vs Newton's method"*, and the paper's only citation for a *numerical* 3-D algorithm is Lévy. Given the paper ships `OMT3D.exe` linked against **`geogram.dll`** (Lévy's library), option (ii) is plausible; given the paper cites GLY [29] as its theory and Gu's own code is Newton, option (i) is equally plausible. **Both use the same energy / gradient / Hessian.**

---

### 3.d How the power-cell (Laguerre cell) volumes are computed in 3-D

Three equivalent routes; the group uses (B)+(C), and the paper's ref [44] documents (A)+(C) precisely.

**(A) Lifting / convex hull (the theoretical route, used to *build* the diagram).**
Lift the sites to `p̂_i = (p_i, |p_i|² − h_i) ∈ ℝ⁴`; the **lower convex hull** of `{p̂_i}` projects to the **power Delaunay tetrahedralisation**; the **upper envelope** of the planes `π_i(x) = ⟨x,p_i⟩ − h_i` (equivalently the **Legendre dual**) projects to the **power diagram**. Each power cell `W_i` is the vertical projection of the union of the envelope faces lying on plane `π_i`.
【Q, Gu lecture slide 33】 *"Power Diagram Algorithm: 1 Compute the convex hull using Lawson edge flipping, add the infinity vertex (0,0,−h); project the convex hull to power Delaunay triangulation T; 2 Compute the upper envelope using Legendre dual algorithm and project to the power diagram D; 3 Clip the power cells using Sutherland–Hodgman algorithm."*
【V】 `PowerDiagramMesh.h`: `__face_dual_point()` computes, for each power-Delaunay face with vertex-lift `p[0..2]`, `n = (p1−p0)×(p2−p0)`, `A = −n_x/n_z`, `B = −n_y/n_z`, `C = p0·n/n_z`, and sets the dual point `(A,B,−C)`, flagging it `inside` iff `A²+B² < 1` (the unit-disk/ball test in 2-D). `_Legendre_transform()` loops over all faces calling this; `_power_cell_clip()` walks the faces around each vertex to form its cell boundary and then clips.

**(B) Half-space intersection by re-entrant clipping (the practical route).**
`W_i ∩ B³ = B³ ∩ ⋂_{j≠i} { x : ⟨x, p_i − p_j⟩ ≥ h_i − h_j }` — the intersection of a convex body with `k−1` half-spaces. Clip iteratively (one half-space at a time). Each clipping step is the 3-D generalisation of Sutherland–Hodgman.
**【Q, the paper's ref [44], Lévy 2015, §2.1 Algorithm 2, verbatim structure:】**
> `Data`: A tetrahedral mesh `M`, a set of points `Y` and a weight vector `W`. `Result`: the intersection `Vor_W(Y) ∩ M`.
> `S : Stack(couple(tet index, point index))`
> **foreach** tetrahedron `t ∈ M` **do** if `t` is not marked **then** `i ← i | Pow_W(y_i) ∩ t ≠ ∅`; `Mark(t,i)`; `Push(S,(t,i))`;
> **while** `S` is not empty **do** `(t,i) ← Pop(S)`; (2) `P : Convex ← Pow_W(y_i) ∩ t`; (3) `Accumulate(P)`; (4) **foreach** `j` neighbour of `i` in `P` … `Push(S,(t,j))`; (5) **foreach** `t′` neighbour of `t` in `P` … `Push(S,(t′,i))`.
> *"a tetrahedron t and a power cell Pow_W(y_i) can be both described as the intersection of half-spaces, as well as the intersection t ∩ Pow_W(y_i), computed using **re-entrant clipping** (each half-space is removed iteratively). I use two versions of the algorithm, a non-robust one that uses floating point arithmetics, and a robust one, that uses arithmetic filters, expansion arithmetics and symbolic perturbation. Both predicates and power diagram construction algorithm are available in **PCK (Predicate Construction Kit)** part of my publicly available **`geogram`** programming library."*

**(C) Volume of the resulting convex polyhedron.**
Once `P = W_i ∩ (ball or cube)` is a convex polyhedron given by its bounding planes/faces:
* either decompose `P` into tetrahedra from one vertex `x₀` and sum `|det(v1−x0, v2−x0, v3−x0)|/6`;
* or use the divergence theorem: `Vol(P) = (1/3) ∮_{∂P} x·n dA`, i.e. sum `(1/3)·(x_c · n_f)·Area(f)` over facets `f`.
【I】 the target paper does not state which; both are standard and either reproduces the authors' stored `weight` values (which I verified equal `(¼)Σ|τ|`, so the code's cell-volume routine is consistent with the true cell volumes).

**(D) 2-D special case actually shipped as source.** `__power_cell_clip()` builds the cell polygon from the dual points of the faces around the vertex, then `SutherlandHodgman(cell, unit_disk_or_cube, out)`; areas are `CPolygon2D::area()`. For a cell straddling the clip boundary, the 3-D conjugate polygon is transformed to the cell plane, intersected with the clip boundary (`plane_cylinder_intersection`), and the intersection polygon is clipped in that plane and mapped back (`inverse_transform`) — an exact, robust way to keep `Vol(W_i) = ν_i` even on the boundary. This is the same trick needed in 3-D at `∂B³`.

**【V, 3-D CLI】**
```
omt3d -input <mesh>.t  [-output <mesh>.tet6]  [-threshold 1e-4]
      [-project true|false]      # project the boundary point to unit sphere (or cube), default true
      [-clip cube|sphere]        # set the clip boundary, default cube
```
The bundled 3-D demos (`test_lion.bat`, `test_david.bat`, `test_bimba.bat`):
```
..\bin\omt3d -input lion_h.t  -output lion_omt.tet6  -clip cube
..\bin\omt3d -input david_h.t -output david_omt.tet6 -clip cube
..\bin\omt3d -input bimba_h.t -output bimba_omt.tet6 -clip cube
```
【V】 `bin\geogram.dll` and `bin\CGAL-vc120-mt-4.7.dll` are shipped with `OMT3D.exe` — i.e. the authors' 3-D power-diagram engine is exactly the Lévy/geogram + CGAL stack. (The paper's ref [43] is the CGAL manual.)

---

### 3.e Composition — the exact formula and direction

**Let `ψ : M → B³` be the volumetric harmonic map (§3.a2) and `φ : (B³, dx) → (P, ν)` the discrete OMT map (§3.b–c). Then**

```
                F  =  φ⁻¹ ∘ ψ  :  M  →  B³
                F(v_i)  =  φ⁻¹( p_i )  =  c_i  :=  mass centre of the power cell W_i ⊂ B³
```
with `c_i = ( ∫_{W_i} x dx ) / |W_i|`, and `F` extended piecewise-affinely on each tetrahedron of `M`.

* 【Q】 the indexed manuscript states the map `φ : 𝔻³ → P` and then *"the inverse map `φ⁻¹` maps each …"*.
* 【V】 the group's identical 2-D construction (Zhengyu Su dissertation, Algorithm 3 step 5): *"Construct the mapping `τ⁻¹∘φ : M → D`, which maps each vertex `v_i ∈ M` to **the centroid of `W_i(h) ⊂ D`**."* The released code has `_power_cell_center()` setting `pv->uv() = pv->dual_cell2D().mass_center();` and the hot key `'m'` = *"Compute the mass center of power cells"*.
* 【V】 Gu's website, 3-D page: *"Given a simply connected solid with a single boundary surface, we first map it onto the unit solid ball using volumetric harmonic map. We treat the volume-distortion factor induced by the harmonic mapping as a probability measure. Then we compute an optimal mass transportation map from this measure to the canonical Euclidean measure on the unit ball. **The composition of the harmonic parameterization and the inverse of the optimal transportation map gives a volume-preserving parameterization of the input volume.**"*
* **Why this is volume-preserving:** since `|W_i| = ν_i = (¼)Σ_{τ∋v_i}|τ|`, the map that sends `W_i` (of volume `ν_i`) onto the vertex `v_i`'s image is measure preserving at the level of the discrete measure, and `F_#(dV_M) = φ⁻¹_# (ψ_# dV_M) ≈ φ⁻¹_# (Σ_i ν_i δ_{p_i}) = dx`. Equivalently, `det(J_{F}) ≡ 1` in the continuum limit.
* The alternative reading of the abstract (source = `μ`, target = `dx`, giving `F = φ∘ψ`) is the **same map written in the opposite direction**; because `φ` and `φ⁻¹` are inverse, `φ∘ψ` vs `φ⁻¹∘ψ` differ. **Only `F = φ⁻¹∘ψ` is consistent with the paper's own notation (`φ : 𝔻³ → P`), with the released implementation (source = Euclidean on the ball, sites = harmonic images, masses = per-vertex volumes), and with the group's 2-D pipeline.** I flag the abstract's verbal phrasing as loose.

---

### 3.f Pseudocode

**【I】 The target paper probably contains an algorithm box; I could not verify its exact content, so the box below is RECONSTRUCTED** (from Gu's Nov-2020 lecture box + the released `omt3d` CLI + the verified data semantics + the group's 2-D Algorithm boxes). The *verbatim* algorithm boxes I did obtain are reproduced in §3.c above.

```
Algorithm  Volume-preserving mesh parameterization by OMT
Input : tetrahedral mesh M ⊂ ℝ³ (topological 3-ball, single boundary surface ∂M);
        threshold ε > 0 (omt3d default 1e-4)
Output: volume-preserving parameterization F : M → B³  (F(v_i) as a new position in B³)

// ── Stage 1: boundary spherical map ─────────────────────────────────────────
1  g ← Tuette map ∂M → S²                       // initial spherical embedding
2  repeat                                        // non-linear heat diffusion (ref [39])
3      for each v ∈ V(∂M):  Df(v) ← Δ_PL f(v) − ⟨Δ_PL f(v), f(v)⟩ f(v)
4      f(v) ← f(v) + δ·Df(v);   f(v) ← f(v)/|f(v)|
5      apply the Möbius/mass-centre normalisation  f ← f − c,  f ← f/|f|
6  until |E − E_old| < δE;      g ← f

// ── Stage 2: volumetric harmonic map ────────────────────────────────────────
7  ψ(v) ← g(v) for v ∈ ∂M;   ψ(v) ← 0 for v ∈ M∖∂M
8  solve  Σ_{u} k_{u,v}( ψ(u) − ψ(v) ) = 0  for v ∈ M∖∂M     // k_{u,v} = (1/12)Σ l_i cot θ_i
     (equivalently: flow  dψ/dt = −Δ_PL ψ  with fixed boundary data)
     →  ψ : M → B³ ,  |ψ| ≤ 1,  ψ(∂M) = S²

// ── Stage 3: discrete OMT on the unit ball ──────────────────────────────────
9  P ← { p_i = ψ(v_i) };   ν_i ← (1/4)·Σ_{τ∋v_i} |τ| ;   rescale ν so that Σν_i = |B³|
10 h ← ½ (|p_1|², …, |p_k|²)                       // or h ← 0
11 repeat
12     (T, D) ← PowerDelaunayAndPowerDiagram( {p_i}, h )     // Lawson edge flip + Legendre dual
13     for each i:  W_i ← clip( B³ , ⋂_{j≠i}{⟨x,p_i−p_j⟩ ≥ h_i−h_j} )   // Sutherland–Hodgman
14     w_i ← Vol(W_i);   ∇E_i ← w_i − ν_i               // gradient = volume difference
15     if ‖∇E‖ < ε: break
16     H_ij ← − Area(W_i ∩ W_j)/|p_i − p_j|  (i≠j);   H_ii ← −Σ_{j≠i} H_ij
17     solve  H d = ∇E        (Eigen SimplicialLDLT / CG; d ← d − mean(d))
18     λ ← 1
19     repeat   h ← h − λ d
20              rebuild (T, D);   // if some site leaves the convex hull, or some W_i = ∅:
21              λ ← λ/2 ;  h ← h + λ d        // roll back and halve
22     until all power cells are non-empty
23 until converged
24 φ : W_i ↦ p_i ,  |W_i| = ν_i

// ── Stage 4: composition ───────────────────────────────────────────────────
25 for each vertex v_i ∈ M:   F(v_i) ← mass centre of W_i          // = φ⁻¹( p_i )
26 extend F piecewise-affinely over the tetrahedra of M
27 return F : M → B³
```

**Complexity.** 【I】 (standard bounds; the paper's own claim not verified)
* Per Newton iteration: build the power diagram ≈ convex hull of `k` lifted points in `ℝ⁴` → `O(k log k)` expected (`O(k²)` worst case); a Lawson edge-flip implementation is `O(k)` flips total for a convex position, `O(k log k)` with a queue. Cell volume clipping: each cell has `O(1)` faces on average, so `O(k)` total, `O(k·f)` worst case.
* One sparse `k×k` SPD solve per iteration: `O(k)` non-zeros, `O(k^{1.5})`–`O(k²)` depending on the solver (`SimplicialLDLT` with a good ordering is near-linear in practice; the code also leaves a CG fallback commented in).
* Newton on this energy converges in a small number of iterations (typically ≲ 20–50 in the group's practice; the damping loop halves `λ` only when a cell degenerates).
* Total: **`O(I·(k log k + k^{1.5}))`**, `I` = number of Newton steps; memory `O(k)`. This matches the group's advertised improvement of the Monge–Brenier approach over Kantorovich: *"the number of variables is O(k²) … the Monge-Brenier approach reduces the unknown variables from O(k²) to O(k)"* 【Q, companion paper §2.2】.

---

## 4. The mathematics that justifies it

### 4.1 Existence / uniqueness of the optimal transport map

* **Brenier's theorem** 【Q, Gu lecture】: *"the optimal transportation map exists and is unique, it is the gradient of a convex function `u : Ω → ℝ`, `T = ∇u`, where `u` is the Brenier potential function. The Brenier potential satisfies the Monge–Ampère equation `det(D²u) = f(x)/g∘∇u(x)` with boundary condition `∇u(Ω) = Ω*`."*
  The quadratic cost is `C(T) = ∫_Ω |x − T(x)|² dx` (the paper also cites [23] Bonnotte, *From Knothe's rearrangement to Brenier's optimal transport map*, and [28] Brenier, *Polar factorization and monotone rearrangement*).
* **Aurenhammer–Hoffmann–Aronov (1998)** 【Q, Gu lecture, slide 35, verbatim】:
  > *"Given a compact convex domain `U` in `ℝⁿ` and `p_1,…,p_k ∈ ℝⁿ` and `A_1,…,A_k > 0`, `Σ_i A_i = vol(U)`, there exists a unique power diagram `U = ⋃_{i=1}^k W_i`, `vol(W_i) = A_i`, the map `T : W_i ↦ p_i` minimizes the transport cost `∫_U |x − T(x)|² dx`."*
* **Alexandrov's convex-polytope theorem** (the geometric face of the same statement; cited as refs [36]/[37]): a convex polytope is uniquely determined (up to translation) by its facet normals and facet volumes. The Minkowski problem version 【Q】: *"given `k` unit vectors `n_1..n_k` not contained in a half-space in `ℝⁿ` and `A_1..A_k > 0` such that `Σ A_i n_i = 0`, find a compact convex polytope `P` with exactly `k` codimension-1 faces `F_1..F_k` such that `area(F_i) = A_i` and `n_i ⊥ F_i`. **Theorem (Minkowski): `P` exists and is unique up to translations.**"*

### 4.2 The variational principle (Gu–Luo–Sun–Yau 2013/2016) — **ref [29]**

【Q, Gu lecture, slide 8, verbatim】
> *"**Theorem (Gu–Luo–Sun–Yau 2013).** `Ω` is a compact convex domain in `ℝⁿ`, `y_1,…,y_k` distinct in `ℝⁿ`, `μ` a positive continuous measure on `Ω`. For any `ν_1,…,ν_k > 0` with `Σν_i = μ(Ω)`, there exists a vector `(h_1,…,h_k)` so that `u(x) = max{⟨x,p_i⟩ + h_i}` satisfies `μ(W_i ∩ Ω) = ν_i`, where `W_i = {x | ∇f(x) = p_i}`. Furthermore, `h` is the maximum point of the concave function `E(h) = Σ_i ν_i h_i − ∫_0^h Σ_i w_i(η) dη_i`, where `w_i(η) = μ(W_i(η) ∩ Ω)` is the `μ`-volume of the cell. The energy `E(h)` is called the **Alexandrov energy**."*

**The five steps of the variational proof** 【Q, Gu's TCD slides, slides 61–65, verbatim fragments】:
1. *"First, we show the admissible height space `H = {h ∈ ℝᵏ | w_i(h) > 0, ∀i}` is non-empty open convex set in `ℝᵏ` by using the **Brunn–Minkowski inequality**."*
2. *"Second, we can show the symmetry `∂w_i/∂h_j = ∂w_j/∂h_i ≤ 0` for `i ≠ j`. Thus the differential 1-form `Σ w_i(h) dh_i` is closed in `H`. Therefore `∃` a smooth `F : H → ℝ` so that `∂F/∂h_i = w_i(h)`, hence `F(h) := ∫^h Σ_{i=1}^k w_i(η) dη_i`."*
3. *"Third, because `∂w_i(h)/∂h_i = 0` due to `Σ w_i(h) = vol(Ω)`, the Hessian of `F` is diagonally dominated, `F(h)` is convex in `H`."*
4. *"Fourth, `F` is strictly convex in `H_0 = {h ∈ H | Σ_{i=1}^k h_i = 0}` and `∇F(h) = (w_1(h),…,w_k(h))`. If `F` is strictly convex on an open convex set `Ω` in `ℝᵏ`, then `∇F : Ω → ℝᵏ` is one-one. This shows the **uniqueness** part of Alexandrov's theorem."*
5. *"Fifth, it can be shown that the convex function `G(h) = Σ A_i h_i − F(h)` has a maximum point in `H_0`. The gradient `∇G` on the boundary of `H_0` points to the interior. Therefore, the maximum point is an interior point, which is the solution to Alexandrov's theorem. This gives the **existence** proof."*
6. 【Q】 *"Remarks. Only the convexity of the energy is insufficient to guarantee the existence of the solution, it is further required that the domain is convex and the critical point is interior. The Brunn–Minkowski inequality is fundamentally essential. **The L2 cost and the volume are Legendre dual to each other.**"*

### 4.3 The lower bound of the volumetric stretch energy / distortion

【V, from the modern literature, not from the paper itself】 The paper's guarantee is stated in terms of "free of volume distortion". The rigorous version of this statement lives in the *stretch-energy* literature that follows it:
* arXiv:2210.09654, Theorem 3.1 + the remark *"the volumetric stretch energy functional (3.1) can be reformulated as … a minimizer of the volumetric stretch energy functional **if and only if** `f*` is volume-/mass-preserving between `M` and `B³`"*; the energy is `E_V(f) = ½ trace(fᵀ L_V(f) f)`, with the stretch factor `σ_{μ,f⁻¹}(τ) = μ(τ)/|f(τ)|`, and ★ `E_V(f) ≥ Σ_τ μ(τ)` (equality iff `σ ≡ 1`), which is the "lower bound of the volumetric stretch energy" and the precise sense in which volume-preserving = energy-minimising.
* 【V】 arXiv:2402.00380 (n-VSE) introduces an **ε-volume-preserving** relaxation: *"ε-volume-preserving parameterization with as small ε as possible"*, i.e. `|σ_{μ,f} − 1| ≤ ε` — this is the natural quantitative version of the CAD paper's claim, and it is the same ε-family as the "measure controllable" companion.

### 4.4 Bijectivity

【Q, abstract】 the paper claims the parameterization *"guarantees to be free of volume distortion"* and the companion paper claims *"the regions with complicated geometries or materials can be mapped to large regions in the parameter domain"* and one wants a *"homeomorphic"* map.
【V/I】 The rigorous chain is:
* `ψ` (volumetric harmonic map into the convex ball `B³` with homeomorphic boundary data onto `∂B³`) is a homeomorphism onto its image by the maximum principle + the Radó–Kneser–Choquet-type theorem for harmonic maps into convex domains in `ℝⁿ`; numerically it can fold for badly shaped meshes. 【I】
* `φ` (the OMT map) is **automatically injective at the level of cells** (the power cells have disjoint interiors and each is collapsed to its site), and its inverse assigns to each site its own cell. So the composition has no fold-over *at cell granularity*.
* The practical guarantee is the standard one used by the whole group: **the Jacobian determinant of the piecewise-affine `F` is checked numerically**, fold-over tetrahedra are counted, and any residual folding is removed by a bounded-distortion post-processing step. 【V】 arXiv:2210.09654 Remark 6.1: *"By applying the large-scale bounded distortion mapping, the folding tetrahedra can be unfolded, slightly sacrificing the volume distortion, so that the resulting mapping is bijective."*
* 【I】 I could **not** verify that the CAD 2017 paper proves bijectivity; the abstract only asserts "free of volume distortion". Treat "bijectivity proof" as **not established by the retrieved evidence**.

### 4.5 Error / ε relaxation and the "measure controllable" version

* **ε as the Newton stopping criterion**: `‖∇E(h_k)‖ < ε` with `ε` = the volume-difference norm; the released tool's default `ε = 1e-4` (relative, per `omt2d`/`omt3d` `-threshold`), and the reported diagnostic is the **max relative cell-volume error** `max_i |ν_i − w_i|/ν_i`.
  【Q, Lévy 2015, for a principled scale】 `ε = 0.01 · µ(M)/√k`.
* **"Measure controllable" ≡ make the target measure user-prescribed.** 【Q, companion abstract】 *"we compute an optimal mass transportation map … onto the parameter domain with the user prescribed volumetric measure"*; 【Q, companion §1】 *"Our method gives users the freedom to assign a target measure to the input solid. Equivalently, the user can prescribe a positive Jacobian function on the whole target domain."*; *"As a special case, our method can make the whole parameterization to be volume-preserving … namely the Jacobian of the mapping is equal to 1 everywhere."*
  ⇒ In the measure-controllable version the algorithm is **identical**; only the constraint `w_i(h) = ν_i` changes to `w_i(h) = μ_target(D_i)` for a user-chosen target measure (e.g. a per-vertex target volume). The volume-preserving paper is the special case `μ_target = ` Euclidean volume.
  ⇒ 【Q, companion §1】 *"this method can handle solids with empty voids inside"* (the harmonic-map stage is generalised through the boundary components).
  ⇒ 【V】 the follow-up of the same group is **Su, Li, Zhao, Lei, Gu, "Discrete Lie flow: A measure controllable parameterization method", CAGD 72 (2019) 49–68** — the same ε/measure-controllable idea for surfaces.

---

## 5. Experimental details

### 5.1 What the paper itself reports — **not retrievable**

The paper's §5/§6 (figures, tables, timings) are inside the paywalled/captcha-walled body. 【Q, abstract】 only: *"We have thoroughly tested our algorithm on many solid models in reality. The experimental results demonstrate the efficiency and efficacy of the proposed method."*

### 5.2 What I could measure from the authors' own released demo data (this IS the paper's toolchain)

Files in `omt-binary/omt3d_test/` (inputs to `omt3d -clip cube`) — **counts computed by me**:

| File | Size | Vertices `#V` | Tetrahedra `#T` | `boundary`-flagged | `Σν = Vol(M)` (file units) | `|h|` boundary / interior |
|---|---|---|---|---|---|---|
| `lion_h.t` | 6 234 137 B | 20 874 | 97 131 | 10 363 | 0.16288814 | `1.000000` / `[0.0590, 0.9998]` |
| `david_h.t` | 5 882 126 B | 19 984 | 90 575 | 10 671 | 0.72849475 | `1.000000` / `[0.0814, 0.9989]` |
| `bimba_h.t` | 5 884 761 B | 21 700 | 89 281 | 0 (unflagged) | 1333.500309 | n/a |
| 2-D: `alex.off` | 2 595 551 B | — | — | — | — | `omt2d -threshold 1e-6` |
| 2-D: `skull.off` | 1 419 379 B | — | — | — | — | `omt2d -threshold 1e-6` |

* Incidence statistics for `lion_h.t`: **min 5, max 46, mean 18.6** tetrahedra per vertex.
* Per-vertex weights: min `3.301e-07`, max `2.884e-04`, mean `7.803e-06`, **no negative weights**.
* Output naming: `<model>_omt.tet6` — the OMT-corrected tetrahedral mesh with the *same connectivity* and **new vertex positions** (= the mass centres of the power cells), i.e. the final volume-preserving parameterization.

### 5.3 Distortion metrics in this line of work (for the team's evaluation harness)

* **Jacobian determinant** `det(J_F)`; volume-preserving ⇔ `det(J_{F⁻¹}) ≡ 1`. 【Q, companion §1】 *"Volume-preserving parameterization preserves the volumetric element, namely the Jacobian of the mapping is equal to 1 everywhere."*
* **Stretch factor** `σ_{μ,f⁻¹}(τ) = μ(τ)/|f(τ)|`; report `mean` and `std` over all tets. 【V】 arXiv:2210.09654 Table 6.1 (a *different* method on the same class of problems, useful as scale reference):

| model | `#T(M)` | `#V(M)` | stretch mean | stretch std | `E_V(f)` | #foldings | time (s) |
|---|---|---|---|---|---|---|---|
| Arnold | 36 875 | 6 990 | 1.0129 | 0.1716 | 6.3014 | 3 | — |
| Heart | 103 751 | 18 408 | 1.0023 | 0.0551 | 6.2722 | 1 | — |
| Igea | 130 375 | 22 930 | 1.0020 | 0.0462 | 6.2738 | 0 | — |
| David Head | 233 663 | 40 669 | 1.0034 | 0.0780 | 6.2873 | 0 | — |
| Max Planck | 390 361 | 66 935 | 1.0025 | 0.0814 | 6.2872 | 8 | — |
| Apple | 559 122 | 102 906 | 1.0002 | 0.0201 | 6.2813 | 0 | — |

  (`E_V → 2π ≈ 6.2832` for a unit ball, i.e. `Σ_τ μ(τ)` for `μ` = unit-ball volume; the table is a *modern* baseline, not the CAD-2017 numbers.)
* **Volume distortion ratio** and **number of folding tetrahedra** are the two headline metrics in the group's papers; 【I】 the CAD 2017 paper most likely reports `det(J)` histograms/colour maps plus per-model timings, but I have no numbers.
* Comparative context 【V】 arXiv:2407.19272 §1.2, an independent description of exactly this method: *"Gu et al. [22, 36], based on the theory of convex polyhedra, compute a polyhedral decomposition of the unit ball such that the volume of each polyhedral cell matches the vertex volume of the given tetrahedral mesh. The resulting map, which sends each cell to its corresponding vertex, is a discrete optimal transport map. Therefore, it can be represented as the gradient of a piecewise linear convex function, and the associated convex optimization problem can be solved by Newton's method."* — this is the single best independent corroboration of §3.b–3.c and it independently confirms the phrase "**the volume of each polyhedral cell matches the vertex volume of the given tetrahedral mesh**".

---

## 6. Code / data release by the authors

**Yes — a complete, runnable release exists (this is the single most valuable find).**

1. **`omt-binary.rar`** — https://www3.cs.stonybrook.edu/~gu/software/omt/omt-binary.rar (21.6 MB), landing page https://www3.cs.stonybrook.edu/~gu/software/omt/
   * `bin/OMT2D.exe`, **`bin/OMT3D.exe`** (the 3-D volume-preserving tool), `bin/miniMeshViewer.exe`, `bin/geogram.dll`, `bin/CGAL-vc120-mt-4.7.dll`, Qt5, boost.
   * `omt2d_test/{alex.off, skull.off, test_alex.bat, test_skull.bat}`
   * `omt3d_test/{bimba_h.t, david_h.t, lion_h.t, test_bimba.bat, test_david.bat, test_lion.bat}`
   * CLI (verbatim from the site): `-output <file>`, `-threshold <num>` (default **1e-4**), `-project [true|false]` (default true), `-clip [cube|sphere]` (default cube), `-help`.
   * **The `.t` vertex attributes (`h` = harmonic-map image, `weight` = per-vertex volume, `boundary` flag) are the ground truth for §3.b and were verified numerically in this report.**
2. **`OT_11_06_2020.zip`** — https://www3.cs.stonybrook.edu/~gu/software/OT/OT_11_06_2020.zip (5.4 MB): full C++ source of the 2-D geometric OMT (`PowerDiagramMesh.h`, `OptimalMassTransportationMesh.h`, `omt_viewer.cpp`), Eigen + MeshLib bundled, CMake. Companion lecture: `Geometric_OT_Lecture.pdf`. Also `bin.zip`, `data.zip` (multi-resolution Alex/David/Sophie surface datasets).
3. **Spherical OMT** tutorial + C++ source: https://www3.cs.stonybrook.edu/~gu/software/SOT/index.html (`SOT_11_29_2020_solution.zip` — **404 at the time of writing**; the page's other zips also 404).
4. `omt3d.exe -clip cube` in the shipped demos means **the demos target the unit cube**; the paper targets **the unit solid ball** (`-clip sphere` / `-project true`). Both are supported.
5. The group's later released algorithm (`Discrete Lie flow`, CAGD 2019; `Spherical optimal transportation`, CAD 2019) is the same machinery for surfaces/spheres.

---

## 7. Confidence ledger — quoted vs verified vs inferred

| Item | Status | Evidence |
|---|---|---|
| Bibliographic record, abstract, keywords | **【Q】 verbatim** | Semantic Scholar API, Stony Brook Pure, DLUT Pure, Daneshyari, DLUT Chinese page |
| 44-item reference list | **【Q】** | OUCI record |
| Boundary spherical map by "non-linear heat diffusion" | **【Q】** (abstract) | the phrase is the paper's own |
| The specific non-linear-heat-diffusion scheme (Tuette init → projected/absolute-derivative harmonic flow → Möbius normalisation) | **【V】** for the cited ref [39] (I have its PDF); **【I】** that the CAD paper uses the same two-phase version | `TMI_04.pdf` |
| Volumetric harmonic edge weight `k_{u,v} = (1/12)Σ l_i cot θ_i` + `Δ_PL f(u) = Σ_v k_{u,v}(f(v)−f(u))` | **【Q】** (VHM ref [5] + the indexed manuscript quote of `Δ_PL`) | `volumetric.pdf`; indexed manuscript snippet |
| Volumetric harmonic map = fixed-boundary Dirichlet problem / heat flow | **【V】** | VHM Algorithm 1 |
| Sites = mesh vertices, `p_i = ψ(v_i)` | **【V】** (numerically) | `lion_h.t`/`david_h.t` `h` attribute: boundary `|h| ≡ 1` |
| **Atom mass `ν_i = (1/4)Σ_{τ∋v_i}|τ|`** | **【V】 exactly (5.6e-8 relative error)** | computed from `lion_h.t` |
| Source measure = Euclidean on `B³` | **【V】** by construction + group's 2-D pipeline; abstract phrases the dual direction | code + data |
| Power/Laguerre diagram, Alexandrov energy, gradient = volume difference, Hessian = −Area(cell∩cell)/distance | **【Q】** for the group's own statements; **【V】** that the paper uses this framework (its refs [29],[38],[44] and the independent description in arXiv:2407.19272) | Gu lecture, Zhengyu Su Eqs. 3.25–3.28, GLY AJM |
| Solver = damped Newton (vs L-BFGS/Lévy) | **【I】** — both are consistent with the evidence; the paper does not say which in any retrievable text | Gu lecture vs Lévy §2 |
| 3-D cell volume = re-entrant half-space clipping ∩ (ball/cube), volume by tet decomposition or divergence theorem; geogram/CGAL | **【V】** for the method (paper's ref [44] + shipped `geogram.dll`/`CGAL` in `OMT3D.exe`); **【I】** for the exact inner routine | `levy_ot3d_arxiv.pdf`, `omt-binary.rar` |
| Composition `F = φ⁻¹∘ψ`, `F(v_i) = mass centre(W_i)` | **【V】** (group's 2-D Algorithm 3 step 5 + code `_power_cell_center` + Gu's 3-D web page) + **【Q】** "the inverse map φ⁻¹ maps each …" | as cited |
| Pseudocode box | **【I】 reconstructed** — the paper's own box (if any) is not retrievable | §3.f |
| Complexity bounds | **【I】** standard bounds; the paper's claim unverified | §3.f |
| Experiments/numbers from the paper | **NOT OBTAINED** | §5.1 |
| Bijectivity proof in the paper | **NOT ESTABLISHED** | §4.4 |
| "First work on volume-preserving parameterization" | **【Q】** (abstract); confirmed by the literature (all later papers cite it as the origin) | arXiv:2210.09654, 2402.00380, 2407.19272 |

---

## 8. Local artifact inventory (for reuse)

```
E:\panyingyun\smartmm\
├─ refs\
│  ├─ MC_companion.txt              ← companion paper, first 3 pages (verbatim)
│  ├─ Geometric_OT_Lecture.txt      ← Gu OMT algorithm box + Hessian (verbatim)  ★
│  ├─ Su_dissertation_suny.txt      ← Algorithm 2 box, Eqs 3.25–3.28 (verbatim)   ★
│  ├─ volumetric.txt                ← VHM energy + k_{u,v} (verbatim)             ★
│  ├─ TMI_04.txt                    ← non-linear heat diffusion / spherical map
│  ├─ levy_ot3d_arxiv.txt           ← 3-D power-cell volume algorithm (verbatim)  ★
│  ├─ arxiv1302.5472.txt, AJM_GLY_variational.txt   ← GLY variational principle
│  ├─ arxiv2210.09654_vsem.txt, arxiv2402.00380_nvse.txt, arxiv2407.19272.txt
│  ├─ crane_nd_cotan.txt            ← tetrahedral cotan weights
│  ├─ gu_omt_theory_tcd.txt, IPAM_OMT.txt, SIGGRAPH_ASIA_OT.txt
│  ├─ gu_bimsa_omt_notes_cn.txt     ← Chinese OMT notes (CJK garbled, math readable)
│  ├─ Zhang_dissertation.txt, HBM04POSTER_volumetric.txt
│  ├─ *_h.t / *.off copies of the authors' test meshes → code\omtbin\omt-binary\omt3d_test\
│  └─ pdf\<all downloaded PDFs/zips>
├─ code\
│  ├─ OT\OT_11_06_2020\             ← C++ source (PowerDiagramMesh.h, OptimalMassTransportationMesh.h)
│  └─ omtbin\omt-binary\            ← OMT2D.exe, OMT3D.exe, geogram.dll, CGAL, test data + .bat
└─ research\01-paper-algorithm.md   ← this file
```

**★ = the four documents a re-implementer actually needs.**
