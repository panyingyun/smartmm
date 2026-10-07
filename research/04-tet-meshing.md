# 04 — 从保体积映射到四面体网格：下半段流水线调研

> 调研范围：给定保体积（共形 + 最优传输）映射 `f: M → B`（单位球）及其逆映射，如何输出一个
> **质量良好、单元尺寸可控（均匀 / 自适应）的四面体网格**。
>
> 工作目录：`E:\panyingyun\smartmm`。本文件是 `research/` 下的第 4 篇。
>
> **证据等级约定**：文中每条结论标注来源等级
> - `[原文]` = 我下载并抽取了 PDF 正文（文件在 `refs/`）
> - `[摘要]` = 只拿到官方摘要 / 元数据（Crossref / Semantic Scholar / OpenAlex / 作者主页）
> - `[实测]` = 我在本机跑通并给出可复现数字（脚本在 `refs/`、`vendor/`）
> - `[未验证]` = 推断或未能证实，请不要直接用于论文／交付

---

## 0. 结论速览（TL;DR）

1. **Su et al. CAD 2016「Measure controllable」做的事不是生成四面体网格**，而是把"体积可控"这个
   能力加到体积参数化里：先把实体经**体积调和映射**映到单位球，再在球上解一个**最优传输映射**把
   调和映射诱导的推前体积元搬成用户指定的目标测度，两者复合 = 测度可控体积参数化。它保证**精确**
   达到指定测度、并且**减少未知量个数**、收敛更快。它输出的是**映射**，不是网格 `[摘要]`。
   下游"重采样成四面体网格"这一步在论文里**没有被强调**（见 §1.4，明确标注为未验证）。

2. **真正把"参数域采样 + 逆映射"变成四面体网格的成熟范式是 Alliez et al. SIGGRAPH 2005
   「Variational tetrahedral meshing」**。它的能量 `E_ODT` 是**抛物面与其分段线性插值的 L¹ 距离**，
   `E_CVT` 是抛物面与其**下**包络（相切平面）的 L¹ 距离——换句话说 `E_CVT` 的最小化**就是**一个
   带"容量约束"（capacity-constrained）的乘方图/Laguerre 图问题，**即半离散最优传输问题**。
   Alliez 明确指出 `E_CVT` 在 3D 会大量产生 sliver，所以必须换成 `E_ODT` `[原文]`。
   这是本报告最重要的一个概念连接（§2.3）。

3. **工程上最省事的推荐路线（见 §4，Recipe A）**：
   `球 B 内均匀/自适应采样 → 用 OMT 单元结构求逆映射 ω⁻¹ → 得到 M 内的点集 → 把 M 的原始边界
   三角面片 + 这些内点喂给 fTetWild / TetGen 做约束 Delaunay 四面体化，单元尺寸用 fTetWild 的
   `--bg-mesh` 或 TetGen 的 `-m` + `.mtr` 背景网格控制`。
   **关键点：不要在球里做网格再拉回来**，而是**在球里只做采样（点集），在 M 里重新四面体化**。
   这自动解决了边界兼容性、sliver 和 sharp feature 三个难题（§4.4）。

4. **不要相信"近似保体积"能保证单元质量**：逆映射是分片光滑、可能不连续的（OT 映射的奇异集），
   直接把球内的正四面体网格拉回去会在奇异集附近把单元压扁/翻转。要么只用它做**点采样**（Recipe A），
   要么把"尺寸场"也一并传输过去让网格生成器去管质量（Recipe B），要么直接在体里做
   **r-自适应（移动节点）**（Recipe C，An–Lei–Zhao–Si–Gu IMR 2021）`[原文]`。

5. **工具许可注意（最容易踩坑）**：TetGen ≥1.5 是 **AGPL-3.0-or-later**（商业需向 WIAS 购买许可）；
   CGAL `Mesh_3` 是 **GPL-3.0-or-later OR 商业许可**；Gmsh 是 **GPL-2.0-or-later**；
   Mmg 是 **LGPL-3.0-or-later**（商用友好）；fTetWild / wildmeshing 是 **MPL-2.0**（商用友好）；
   TetWild 是 **GPL-3.0**；**PyMesh 仓库里根本没有 LICENSE 文件**，法务风险最高。
   详见 §3 与 §3.9。

---

## 1. 「Measure controllable volumetric mesh parameterization」精读

### 1.1 元数据（已用 Crossref 逐字核对）`[摘要]`

| 字段 | 值 |
|---|---|
| 标题 | Measure controllable volumetric mesh parameterization |
| 作者 | Kehua Su; Wei Chen; Na Lei; Li Cui; Jian Jiang; **Xianfeng David Gu** |
| 期刊 | Computer-Aided Design |
| 卷/页/年 | **vol. 78, pp. 188–198, 2016**（2016-09 出版；在线 2016-06-06） |
| DOI | [10.1016/j.cad.2016.04.007](https://doi.org/10.1016/j.cad.2016.04.007) |
| 专刊 | Solid and Physical Modeling (SPM) 2016 特刊（DLUT 页面记录 `卷号：78 期号：,SI`） |
| 关键词 | Parameterization; Volume; Controllable; Optimal mass transportation |
| 收录 | Scopus, CPCI-S, EI, SCIE |
| 被引 | Semantic Scholar 记 11 次（估 16.3） |

> ✅ **任务书里的元数据基本正确**，只需补全：作者共 6 人（含 Li Cui、Jian Jiang），卷 78 页 188–198。
> Crossref 记录的末位作者是 "Gu, Xianfeng David"。

**可读副本情况（合规来源排查结果）**：
- 我尝试并失败的合规途径：ScienceDirect 直链 PDF（403）、`daneshyari.com/article/preview/440675.pdf`
  （反复 504）、Unpaywall（无仓库全文，`any_repository_has_fulltext: false`）、
  CORE API（500）、Stony Brook ResearchConnect（只有元数据页，无 PDF）、
  Gu 主页 `www3.cs.stonybrook.edu/~gu/`（只有 `publications.pdf` 列表，无该文全文）。
- **因此我没有拿到该文正文**。下述算法描述来自**两处官方摘要**（DLUT 教师主页、Semantic Scholar
  论文页）交叉比对，两者文字一致，可信度高；但**图/表/质量数字我无法给出**，已在 §1.4 明确标注。
- 供你们后续取全文的**唯一合规入口**：机构订阅（Elsevier ScienceDirect）、
  Stony Brook 图书馆、或直接联系作者（Na Lei: `nalei@dlut.edu.cn`，Gu: `gu@cs.stonybrook.edu`）。
  **不要用 sci-hub/libgen 镜像**。

### 1.2 论文摘要原文（DLUT 中文主页录入，与 Semantic Scholar 一致）`[摘要]`

> "Volumetric parameterization is a fundamental problem in solid and physical modeling. In practice, it is
> highly desirable to control the volumes of the regions of interest in the parameter domain. This work
> introduces a novel volumetric parameterization method, which allows users to prescribe the target volumetric
> measure of the input solid.
>
> **Given a simply connected tetrahedral mesh with a single boundary surface, we first compute a volumetric
> harmonic map to parameterize the solid onto the unit solid ball; then we compute an optimal mass
> transportation map from the unit solid ball with the push-forward volume element induced by the harmonic map
> onto the parameter domain with the user prescribed volumetric measure. The composition of the volumetric
> harmonic map and the optimal mass transportation map gives a measure controllable volumetric
> parameterization.** Furthermore, this method can handle solids with empty voids inside.
>
> The method has solid theoretic foundation, and is based on conventional algorithms in computational geometry,
> and easy to implement. The experimental results demonstrate the efficiency and efficacy of the proposed method."

Semantic Scholar 版额外给出四条贡献（DLUT 版没有）：
> "The proposed method can handle solids with complicated topologies. **Allows users to prescribe measures on
> the target and guarantees to reach it exactly.** Comparing to other algorithms the proposed method
> **reduces the number of unknowns**. The algorithm in this paper **achieves higher convergence rate**."

### 1.3 「Measure controllable」比「保体积」多了什么？（逐步拆解）

任务书要的"到底多了什么"，我按可验证的信息严格拆开讲：

**符号**：`M` 为输入实体（单边界、单连通四面体网格），`ω: M → B` 为体积调和映射，
`μ := ω_#(dV_M)` 为 `ω` 诱导的**推前体积元**（push-forward volume element），
`ν` 为用户指定的**目标体积测度**（在参数域上，即"每个区域应该占多少体积"）。

| 步骤 | 保体积版本（Volume preserving, Su et al. CAD 2017） | Measure controllable 版本（CAD 2016） |
|---|---|---|
| ① 归约到球 | 体积调和映射 `ω: M → B` | 同左（相同起点） |
| ② 球上的传输 | 求 OT 映射把 `μ` 搬成**均匀测度**（`ν = const`） | 求 OT 映射把 `μ` 搬成**用户预设测度 `ν`** |
| ③ 复合 | `f = T ∘ ω`，`f_#(dV_M) = 均匀` | `f = T ∘ ω`，**`f_#(dV_M) = ν`（逐区域精确命中）** |
| ④ 用户自由度 | 无（只能得到等体积分布） | **每个感兴趣区域的体积可任意指定**——这是"measure controllable"的全部含义 |
| ⑤ 额外拓扑能力 | 单连通 | 摘要称可处理 **empty voids（内部空腔）** 与 **complicated topologies** |

**关键点 1：测度是怎么"被传输"的。** 摘要明确写了 OT 的**源测度**不是球上的均匀勒贝格测度，
而是 `ω_#(dV_M)`——即"体积调和映射推过去的那个不均匀体积元"。这是一个很重要的工程细节：
它意味着不需要假设调和映射近似保体积；OT 步骤会**精确地**把调和映射造成的体积畸变纠正掉。
参考该组同一时期的姊妹工作 `Discrete Lie flow: A measure controllable parameterization method`
(CAGD 72, 2019, pp. 49–68, DOI [10.1016/j.cagd.2019.05.003](https://doi.org/10.1016/j.cagd.2019.05.003))
可以看他们的技术路线 `[摘要]`。

**关键点 2：计算上用的是半离散 OT（乘方图 / Laguerre 图 + Newton）。** 摘要说"基于计算几何的常规
算法""减少未知量个数""收敛更快"，这与该组一贯的做法完全吻合：把 Brenier 势 `u` 表示为支撑平面
`π_i(x) = ⟨p_i, x⟩ + h_i` 的**上包络**，其投影是**乘方图** `{W_i(h)}`；能量
`E(h) = Σ_i ∫_0^{h_i} vol(W_i) dh_i − Σ_i ν_i h_i`，梯度 `∇E(h) = (vol(W_i(h)) − ν_i)`，
Hessian 是加权 Laplacian（`∂²E/∂h_i∂h_j = −μ(W_i ∩ W_j)/|p_i − p_j|`），用 Newton 求解。
这套算法的**完整可运行形式**见该组后来的 IMR 2021 论文（§2.4），我在下面给出逐行算法 `[原文]`。

**关键点 3：`ν` 从哪里来。** 这是"measure controllable"工程上真正要回答的问题。摘要只说
"user prescribed volumetric measure"，没说默认构造。合理实现是：用户在球内给一个标量密度
`ρ(x)`（例如按到球心的距离、按解剖区域、按曲率/厚度），离散成 `ν_i = ρ(p_i)` 再归一化。
**这一点论文里怎么做的，我无正文，属 `[未验证]`。**

### 1.4 四面体网格如何"重采样"，以及质量数字——**诚实说明**

任务书要求"how the tetrahedral mesh is then RESAMPLED, and what the resulting tetrahedral mesh
quality is. Give the algorithm and any figures/tables with quality numbers."

**我必须明确说明：我拿不到该文正文，所以无法给出它的图和表。** 我在摘要里也没有看到任何
"resample / remesh / tetrahedral quality" 的表述。基于可验证证据我能说的是：

- 该文的输出被定义为**映射（parameterization）**，标题和摘要通篇是 parameterization 而非 meshing；
- 该文**输入**本身就要求是一个四面体网格（"simply connected tetrahedral mesh"），
  所以它的目标不是"从零生成四面体网格"；
- 该组同期的**密度控制 → 网格**的证据在他们**后来**的工作里才明确出现：把参数化后的三角形网格
  通过 OT 改变密度、再由参数化拉回生成网格；四边/六面版本见 IMR 2025（§2.5）。
  **因此"CAD 2016 里具体怎么重采样、质量数字是多少"我标记为 `[未验证]`，需要拿到正文才能回答。**
- **可用的替代证据**：如果目标只是"球（参数域）内采样 → 逆映射 → 高质量四面体网格"，
  Alliez et al. 2005 给出了**可直接引用的质量表**（§2.2 表），而 IMR 2021 给出了**可运行的 OT 算法**（§2.4）。
  这两篇加起来可以完整替代 CAD 2016 在这一步的（未公开的）细节。

---

## 2. 保体积参数化如何真正用于生成四面体网格

### 2.1 范式：参数域采样 + 逆映射 = 新网格

核心恒等式（贯穿本报告）：

```
    M  ──── f (保体积映射) ────►  B (单位球)
    ▲                              │
    │                              │ ① 在 B 里定义采样/网格（规则网格、CVT、或尺寸场驱动）
    └────── f⁻¹ (逆映射) ──────────┘ ② 把样本点/单元拉回 M
```

因为 `f` 近似保体积，所以"在 `B` 里采样的**体积分布**"≈"在 `M` 里的体积分布"，
于是**在 `B` 里指定的单元尺寸场可以直接搬到 `M`**：若 `B` 里球坐标下的目标边长为 `h_B`，
则 `M` 中该处的目标边长近似为 `h_M(x) ≈ h_B(f(x)) · (dV_B/dV_M)^{1/3} ≈ h_B(f(x))`（保体积时）。
**这是"用保体积映射做尺寸场传输"的全部理论基础。**（`[未验证]`：我没有找到把这个式子在
CAD 2016 语境下写出来的论文；这是我的工程归纳，数学上是 Jacobian 行列式换元法的直接推论。）

### 2.2 Alliez, Cohen-Steiner, Yvinec, Desbrun — Variational Tetrahedral Meshing

**元数据（Crossref 核对）**`[摘要]`
- Pierre Alliez, David Cohen-Steiner, Mariette Yvinec, Mathieu Desbrun.
  "Variational tetrahedral meshing." **ACM Transactions on Graphics 24(3):617–625, July 2005**
  (SIGGRAPH 2005). DOI: [10.1145/1073204.1073238](https://doi.org/10.1145/1073204.1073238)
  （SIGGRAPH proceedings 版 DOI：10.1145/1186822.1073238；course 版 10.1145/1198555.1198669）
- **开放获取副本（合规）**：CaltechAUTHORS 机构库
  [https://authors.library.caltech.edu/records/c9exj-3ed14](https://authors.library.caltech.edu/records/c9exj-3ed14)
  → 提交版 PDF `ACSYD.pdf`（778 KB，已存 `refs/alliez2005-variational-tet.txt`）`[原文]`

**算法（原文 §3 伪代码，逐字）**`[原文]`
```
Read the input boundary mesh ∂Ω
Setup Data Structure & Preprocessing
Compute sizing field μ
Generate initial sites xi inside Ω
Do
    Construct Delaunay triangulation({xi})
    Move sites xi to their optimal positions x*_i
Until (convergence or stopping criterion)
Extract interior mesh
```
- **连通性**：直接取 Delaunay 三角化——原文说这对 `E_ODT` 是**全局最优连通性**
  （引 Chen & Xu 2004）。这是与只做局部翻边的传统方法最本质的区别。
- **顶点位置**：解析求出 `x*_i`（原文 §3.3 附近的"optimal position"），不是启发式平滑。
- **边界**：用**受约束质心 Voronoi 剖分（CCVT）**处理边界顶点（`between the 3D Voronoi cell of the
  boundary vertex and the input surface`），并对边界顶点做 **jittering** 来去掉残余 sliver。
- **尺寸场**：`μ(x) = inf_{s∈∂Ω} [K·d(s,x) + lfs(s)]`（原文式 (9)），
  即"不超过局部特征尺寸 `lfs` 的**极大 K-Lipschitz 函数**"，`K` 就是 gradation 参数
  （`K=0` 对应均匀）。`lfs` 由**中轴（medial axis / poles）**的离散近似给出：
  每个 Delaunay 顶点取最远 Voronoi 顶点作 pole，`lfs` = 到最近 pole 的距离（原文 §3.2）。
  **边界上的求积样本权重直接编码尺寸场**：面积样本 `ds/μ(x)^4`，特征线样本 `dl/μ(x)^3`
  ——这就是把尺寸场变成"mass density"的关键一步（也是与最优传输挂钩的地方）。

**质量数字（原文 Table 1 + Fig. 12 说明）**`[原文]`

| Model | #v | #tets | min **radius ratio** | **average radius ratio** | mean L² approx. error (% of bbox) |
|---|---|---|---|---|---|
| Torus | 1K | 4K | 0.42 | 0.88 | 0.17 |
| Bunny | 49K | 275K | 0.37 | 0.89 | 0.04 |
| Hand | 36K | 174K | 0.29 | 0.86 | 0.024 |
| Gargoyle | 50K | 260K | 0.23 | 0.88 | 0.053 |
| Fandisk | 3K | 14K | 0.29 | 0.87 | 0.021 |

Hand 模型额外数字（Fig. 12 caption）：**average dihedral angle = 70°**，
worst radius ratio = 0.29，`K = 1`。
**性能**：Fig. 10（≈50 次迭代）**16 秒**（Pentium IV 3 GHz）；Hand 36K 顶点 174K 单元平均
**2 秒/迭代**（含 Delaunay、边界求积、顶点更新）；"10–20 次迭代就够，我们常做到 50 次"。
**明确局限（原文 §"Limitations"）**：① 逼近边界而非严格共形于边界；
② **对结果质量没有任何理论界**；③ 比贪婪式 Delaunay refinement 慢。

**与"保体积映射"的接口**：Alliez 的 `μ` 是**几何内在**的（lfs/曲率），而本项目的尺寸场来自
**保体积映射的 Jacobian**。二者可以相乘：`h_final(x) = h_user(f(x)) · h_geometric(x)`。
`[未验证]`：这条组合式在我查到的论文中没有明确写出，请用之前先用合成算例验证。

### 2.3 ⭐ 显式连接：Variational tetrahedral meshing **就是**半离散最优传输

这是任务书特别要求点明的，我把原文证据摆出来。

**原文证据（Alliez 2005 §2.1，逐字）**`[原文]`
> "**Centroidal Voronoi Tessellations** Du and Wang [2003] propose to generate meshes that are dual to
> optimal Voronoi diagrams. These diagrams are achieved by minimizing the quadratic energy:
> `E_CVT = Σ_{i=1..N} ∫_{V_i} ||x − x_i||² dx` (1) …
> **From the analysis of E_CVT it is well known that its minimization corresponds to minimizing the volume
> between a paraboloid f(x) = ||x||² and an underlaid, circumscribing piecewise linear approximant
> f^{dual}_{PWL}, which is formed by planar patches tangent to the paraboloid:**
> `E_CVT = || f − f^{dual}_{PWL} ||_{L¹}`."

**原文证据（同节，关于 3D 失败）**`[原文]`
> "Unfortunately, and despite Du's proposal [2003] to use CVTs for tet meshing, **there exists no proof of
> such a dual property in 3D. Our own tests show that using Du's suggestion for tet meshing gives rise to
> numerous degenerate sliver tets** (see Fig. 2). We can attribute the slivers to the fact that E_CVT tends
> to optimize the compactness of the dual Voronoi cells, but not the compactness of simplices in the primal
> Delaunay triangulation: therefore, **the presence of a sliver is not penalized by this energy.**"

**原文证据（ODT）**`[原文]`
> "Recently, Chen [2004] proposed an approach 'dual' to the above … He used the following energy:
> `E_ODT = || f − f^{primal}_{PWL} ||_{L¹}` … `E_ODT = 1/(n+1) Σ_{i=1..N} ∫_{Ω_i} ||x − x_i||² dx` (2) …
> The integral is now taken over each 1-ring region Ω_i (also called the star of the vertex x_i).
> Notice that these regions overlap. … **Chen's E_ODT energy measures a quality of the simplicial mesh, not
> its dual. It is thus more prone to generate well-shaped primal elements, while E_CVT was maximizing the
> compactness of the dual Voronoi cells.**"

**为什么这"就是"半离散 OT（本项目要的显式连接）**：

1. `E_CVT` 的最优划分，是在**给定站点**下取 **Voronoi 图**；把站点加权（`w_i`）就得到
   **乘方图 / Laguerre 图（power diagram）** `{W_i(w)}`，其中
   `W_i = { x : ||x − x_i||² − w_i ≤ ||x − x_j||² − w_j, ∀j }`。
   这正是 Kantorovich 对偶里 Brenier 势的上包络投影。
2. **"每个 cell 承载指定质量"= 容量约束 Voronoi 剖分（capacity-constrained CVT, CCVT）**。
   若要求 `μ(W_i(w)) = ν_i`（`μ` 为密度测度），求合适的权重 `w`，这就是
   **L² 半离散最优传输**：源 = 连续密度 `ρ`，目标 = Dirac 质量 `{ν_i}`，
   最优传输映射恰是 `T: W_i ∋ x ↦ x_i`，由 Aurenhammer 等人的经典结果保证乘方图权重唯一决定。
   `[原文]` An et al. IMR 2021 直接把这句话写出来了：
   > "The relationship between **the capacity-constrained Voronoi tessellation** and the power diagram …"
   （见 `refs/imr2021-moving-mesh-ot.txt` 第 257 行附近，原文 §2 Related Work）
3. Alliez 的算法里**两次**用到这个结构：
   - **连通性**：Delaunay 三角化 = 乘方图（权重=0 时）的对偶。Alliez 说"remarkably, the Delaunay
     triangulation is again the optimal connectivity which minimizes `E_ODT`"——这是 Brenier 势的
     下包络直接给出的。
   - **顶点位置**：`x*_i` 的解析解由 `∫_{Ω_i} (x − x_i) dx = 0` 型条件给出，
     本质上就是"cell 的质心"条件（CVT 的 centroid 条件），但作用在**重叠 star** 而不是 disjoint cell 上。
4. **Alliez 自己在结论里承认下一站是各向异性**：
   > "In future work, we wish to explore how to extend our approach to anisotropic meshing
   > **using not just a mass density, but a tensor field.**" `[原文]`
   ——`mass density` 正是 OT 里的 `ρ`。

> **给团队的结论**：如果论文里要写"variational tetrahedral meshing is an optimal-transport-flavoured
> problem"，请写成更强、更准确的说法：
> **`E_CVT` 的最小化（在容量约束下）就是 L² 半离散最优传输问题；`E_ODT` 是同一个抛物面逼近
> 但换成上包络（primal PWL）的版本，用它可以避免 3D 中 CVT 必然产生的 sliver。**
> 这个表述有 Alliez 2005 §2.1 的逐字支持。

### 2.4 可直接落地的 OT 网格算法：An, Lei, Zhao, Si, Gu, IMR 2021 `[原文]`

**元数据**（IMR 官网论文集页 + OpenAlex 核对）：
- Dongsheng An, **Na Lei**, Tong Zhao, **Hang Si**, **Xianfeng Gu**.
  "A Moving Mesh Adaptation Method by Optimal Transport."
  **Proceedings of the 29th International Meshing Roundtable (IMR 2021)**, pp. 130–138.
  开放 PDF：[https://www.meshingroundtable.com/assets/papers/2021/10-An.pdf](https://www.meshingroundtable.com/assets/papers/2021/10-An.pdf)
  （Zenodo DOI: [10.5281/zenodo.5559175](https://doi.org/10.5281/zenodo.5559175)）
  本地文本：`refs/imr2021-moving-mesh-ot.txt`
  ⚠️ **注意**：官网列表页把它标为 `10-An.pdf`，OpenAlex 的标题是
  "A Moving Mesh Adaption Method By Optimal Transport"（Adaption 拼写），**无 DOI**（SIAM 之前的 IMR
  用 Zenodo badge）。**期刊版本我未找到。**

**这篇为什么重要**：它是"用 OT 做网格密度控制"的**完整可复现算法 + 复杂度/收敛数据**，
且作者组与 CAD 2016 高度重合（Na Lei、Gu 都在；Hang Si 是 TetGen 作者）。
它解决的是 **r-自适应（移动节点，不改变连通性）**，正是"把密度场施加到已有网格上"的标准手段。

**算法（原文 Algorithm 1，逐字转录，2D 版但结构对 3D 完全一致）**`[原文]`
```
Input:  source domain Ω with uniform distribution μ;
        density function g and its corresponding rectangular domain Ω*
Output: the map T̂ and the generated adaptive mesh

Construct the target measure ν = Σ_{i=1..n} ν_i δ(y − p_i) on Ω* with grid
    vertices {p_i}, set ν_i = g(p_i) and normalize as ν_i = μ(Ω)·ν_i / Σ_j ν_j
Initialize h_i = ½(|p_i|² − 1)
while true do
    Compute the lower convex hull of {(p_i, −h_i)}_{i=1..n}
    Compute the upper envelope of the planes {⟨p_i, x⟩ + h_i}
    Project the upper envelope to the plane to get a power diagram Ω = ∪ W_i(h)
    Compute the μ-volume of each finite cell w_i(h) = μ(W_i(h))
    Compute the gradient ∇E(h) = (ν_i − w_i(h))
    if ‖∇E(h)‖ < ε then return h
    Compute the μ-lengths of the power Voronoi edges W_i(h) ∩ W_j(h)
    Construct Hessian:
        Hess(E)_ij = − μ(W_i ∩ W_j) / |p_i − p_j|
        Hess(E)_ii =  Σ_{j~i} μ(W_i ∩ W_j) / |p_i − p_j|
    Solve Hess(E) d = ∇E(h)
    λ ← 1
    repeat
        if h + λd ∉ H then λ ← ½ ; continue
        Compute the power diagram D(h + λd)
    until no empty power cell
    Update h ← h + λd
end
Compute cell centers {m_i} of {W_i}, giving T̂(m_i) = p_i
Generate the adaptive mesh of {m_i} according to the connectivity of {p_i}
```
**要点解读（工程上最有用的几句）**：
- `H = {h : all power cells non-empty}` 是**可行域**；线搜索就是拿来保住"无空 cell"这个不变量的。
  这是整个算法的鲁棒性核心。**移植到 3D 时这一步必须保留**（否则 Newton 步会杀死 unit）。
- `∇E(h) = (ν_i − w_i(h))`：**梯度 = 目标质量 − 当前质量**。这就是"容量约束"的直接体现。
- Hessian 是权重挂在**乘方图棱长**上的加权 Laplacian，稀疏正定，可用 Cholesky。
- 最后一步 `T̂(m_i) = p_i` 是**把 cell 质心映到站点**——这是把离散 OT 映射变成**网格节点位移**的关键，
  也是"逆映射求值"在实践中的廉价实现：**不需要连续逆映射，只需要算每个 cell 的质心。**

**报告的数字（原文 Table 1 + Fig. 7/8，Franke 函数）**`[原文]`

| 顶点数 | 迭代次数 | 时间 (s) |
|---|---|---|
| 1 681 (41×41) | 8 | 0.227 |
| 10 201 (101×101) | 10 | 1.544 |
| 40 401 (201×201) | 12 | 7.938 |
| 160 801 (401×401) | 12 | 53.307 |

- 停机条件：`max_i(|w_i − ν_i|/ν_i) < 1e-8`。
- 原文："**Only after seven iterations, the algorithm converges well and nearly all of
  `log(w_i/ν_i)` concentrates on 0.**"
- 测试密度：ring `g(y)=1+10·sech²(200||y||²−0.25²)`；bell `g(y)=1+50·sech²(100||y||²)`；
  3-峰高斯混合；Franke 函数；真实灰度图（400×400 采样）。
- ⚠️ **重要局限**：**所有实验都是 2D 平面**。原文 §2 把 "three-dimensional r-adaptivity" 列为
  related work 而非本文工作。**所以"3D 版本"需要你们自己做，风险点是 3D 乘方图求交与积分**
  （见 §2.6 的 Lévy / Mérigot–Meyron–Thibert 工具）。

**姊妹论文（同一批人，同一年 IMR）**：
- Dongsheng An, Na Lei, Wei Chen, Zhongxuan Luo, Tong Zhao, Hang Si, Xianfeng Gu.
  "Efficient Approximation of Optimal Transportation Map by Pogorelov Map."
  IMR 2021, [PDF](https://www.meshingroundtable.com/assets/papers/2021/09-An.pdf),
  Zenodo DOI [10.5281/zenodo.5559179](https://doi.org/10.5281/zenodo.5559179)。`[摘要]`

### 2.5 该组的密度控制 → 网格（后来才明确出现）

- **Quad Mesh Density Control using Optimal Transport** — Yuyang Cao, Yue Wang, Yiming Zhu,
  **Hang Si**, **Na Lei**, Zhongxuan Luo, **Xianfeng Gu**.
  *Proceedings of the 2025 SIAM International Meshing Roundtable*, **pp. 45–57**.
  DOI: [10.1137/1.9781611978575.5](https://doi.org/10.1137/1.9781611978575.5)
  开放 PDF: [meshingroundtable.com/assets/papers/2025/1015-compressed.pdf](https://www.meshingroundtable.com/assets/papers/2025/1015-compressed.pdf)
  `[原文]`（`refs/imr2025-quad-density-ot.txt`）
  **这是最接近"measure controllable → 网格"这一思想的公开文献**，只是把 quad 换成你关心的 tet 就能套。
  摘要逐字："In this work, a novel method for controlling the density of quad meshes is proposed based on
  the Optimal Transport (OT) theory. … It provides the designers with the flexibility to control the meshing
  density by **prescribing arbitrary probability density functions** on the surface."
  流程（原文 §1 与 §4）：**Ricci flow 共形映射到平面矩形 → 在参数域跑半离散 OT（Algorithm 4.3，与
  §2.4 同构）→ 用 dual face 的质心作为新顶点位置 → 用 Zheng 的方法生成 quad mesh**。
  Algorithm 4.3 的关键算子（原文逐字）：
  ```
  ∇E(hⁿ) = v_i.dualface.area − v_i.area
  ∂²E/∂h_i∂h_j = − e_ij.length / e_ij.dualedge.length
  SolveHess(hⁿ)d = ∇E(hⁿ);  λ=1
  hⁿ⁺¹ ← hⁿ + λd
  if Delaunay check fails → Roll back, λ ← λ/2
  if any power cell is empty → Roll back, λ ← λ/2
  ε = Σ |v.targetarea − v.area|;  until ε ≤ threshold
  v.coordinates = v.dualface.center        // ← 质心即新顶点！
  ```
  **实现环境（原文 §5）**：generic C++ / Visual Studio 2022 / Windows；
  Intel i9-12900H 14 核 / 32 GB。
- **Optimal Surface Quadrilateral Mesh Generation** — Zhou Zhao, Siyu Fang, Na Lei, Yuanpeng Liu,
  Yiming Zhu, Chander Sadasivan, Apostolos Tassiopoulos, Shikui Chen, Xianfeng Gu.
  IMR 2024, pp. 14–27. DOI: [10.1137/1.9781611978001.2](https://doi.org/10.1137/1.9781611978001.2)
  PDF: [meshingroundtable.com/assets/papers/2024/1010.pdf](https://www.meshingroundtable.com/assets/papers/2024/1010.pdf)
- **Robust Surface Remeshing Based on Conformal Welding** — Wei Chen, Siquan Sun, Yue Wang, Na Lei,
  Chander Sadasivan, Apostolos Tassiopoulos, Shikui Chen, Hang Si.
  IMR 2024, DOI: [10.1137/1.9781611978001.3](https://doi.org/10.1137/1.9781611978001.3)
  PDF: [meshingroundtable.com/assets/papers/2024/1009.pdf](https://www.meshingroundtable.com/assets/papers/2024/1009.pdf)
- **Curvature adaptive surface remeshing by sampling normal cycle** — Kehua Su, Na Lei, Wei Chen, Li Cui,
  **Hang Si**, Shikui Chen, Xianfeng Gu. *Computer-Aided Design* 111 (2019) 22–34.
  DOI: [10.1016/j.cad.2019.01.004](https://doi.org/10.1016/j.cad.2019.01.004) `[摘要]`

### 2.6 「Higher-dimensional power diagrams for semi-discrete optimal transport」(Caplan, IMR 2021)

**元数据**（IMR 官网 + Crossref 核对）`[原文]`
- Philip Claude Caplan. "High-Dimensional Power Diagrams for Semi-Discrete Optimal Transport."
  **Proceedings of the 29th International Meshing Roundtable (IMR 2021)**, pp. 44–56.
  PDF: [https://www.meshingroundtable.com/assets/papers/2021/04-Caplan.pdf](https://www.meshingroundtable.com/assets/papers/2021/04-Caplan.pdf)
  Zenodo DOI: [10.5281/zenodo.5559221](https://doi.org/10.5281/zenodo.5559221)
  本地文本：`refs/caplan2021-power-diagrams.txt`

⚠️ **标题修正**：任务书写的是 "Higher-dimensional power diagrams for semi-discrete optimal transport"，
官网列表页写作 **"High-Dimensional Power Diagrams for Semi-Discrete Optimal Transport"**（首字母大写）。
两者同义，引用时建议用官网/PDF 标题。

**核心内容（与 Aurenhammer 的关系 + 算法）**`[原文]`
- 定位："Recent algorithms have been developed for solving the semi-discrete problem in 2d and 3d,
  however, algorithms in higher dimensions have yet to be demonstrated, which rely on the efficient
  calculation of the power diagram (Laguerre diagram) in higher dimensions. Here, we introduce an algorithm
  for computing power diagrams, which extends to any topological dimension." 关键词：
  Voronoi diagram, power diagram, Laguerre diagram, semi-discrete optimal transport, high dimensions, quantization。
- **Aurenhammer 的观察被用作全部算法基础**（原文 §1）：
  > "In 1987, Aurenhammer observed that **a power diagram can be computed via the restriction of a Voronoi
  > diagram in a higher dimensional space** [20,21]. That is, the power diagram is the intersection of a
  > higher-dimensional Voronoi diagram with our domain embedded to a higher dimensional space."
- **提升（lifting）公式，原文式 (11)——工程上可直接抄**：
  ```
  z_i = ( y_iᵀ , sqrt(max(w) − w_i) )ᵀ  ∈ R^{d+1},   i = 1..N
  ```
  然后把域也提升（坐标补 0），于是 **power cell  P_i(Y, w) = V_i(Z)**（V 是 R^{d+1} 里的 Voronoi cell）。
  这个公式来自 **Lévy**（原文引 [36]：B. Lévy, "A numerical algorithm for L² semi-discrete optimal
  transport in 3D", ESAIM: M2AN **49(6):1693–1715, 2015**,
  DOI [10.1051/m2an/2015055](https://doi.org/10.1051/m2an/2015055)）`[摘要]`。
- **算法 2（割平面）**：从单个 `d`-单纯形 `σ` 出发，按到站点 `z_i` 的**最近邻顺序**依次用 Voronoi
  平分面 `H_1, H_2, …` 去裁剪；每次裁剪后用 **security radius theorem**（引 Lévy–Bonneel）判断
  是否可以提前终止。这套"按 k-近邻顺序裁剪 + 安全半径剪枝"是 3D 幂图求交的实用做法。
- **算法 3（权重优化，原文逐字）**`[原文]`
  ```
  optimizeWeights(N, ρ(x))
  input:  number of sites N, density ρ(x)
  output: power cells with uniform mass
  1  Y ← optimizePoints(N, ρ)
  2  w ← 0
  3  ν_t ← m_t/N            // m_t is the total mass
  4  for iter = 1 : nb_iter
  5      Z ← lift Y to R^{d+1} (Eq. 11)
  6      compute Vor(Z) ∩ Ω                       (Section 3.2)
  7      compute energy and mass in Eq. 4
  8      compute gradients dE/dw_i using Eq. 6
  9      perform L-BFGS update on w with E and dE/dw
     end
  ```
  **这就是"容量约束 Voronoi 剖分 = 半离散 OT"的最干净的伪代码**，可以直接放到你们的设计文档里。
- **实测规模与结论**：维度 2–6 评测；4D 用 `N = 1000` 个站点；两个密度
  `ρ_u(x)=1` 与 `ρ_s(x)=1+100||x−μ||²`，`μ=(0.5,0.5,0.5,0.5)ᵀ`；
  **约 75 次迭代收敛**（Fig. 12：梯度模快速下降，各 cell 归一化质量收敛到 1）。
  原文结论："the cost of performing the numerical integration to compute gradients significantly outweighed
  the cost of computing the power diagram, especially for four-dimensional polytopes."
  → **给你们的直接启示**：3D 提升到 4D 之后，瓶颈是**积分（求 cell 体积）而不是求交**。
- 参考实现生态（原文 §1 明确点出，2D/3D 可用）：**geogram、Voro++、CGAL**；
  更高的维度此前无人实现。**`[实测建议]`**：你们的 3D 情形直接用 **geogram**
  (`GEO::PeriodicDelaunay`/`ConvexCell`) 或 **CGAL** 的 `Regular_triangulation_3`
  + `Power_diagram_2`，不必自己写割平面。

### 2.7 「用乘方图做网格」的其他 IMR 论文（含引 Aurenhammer 的）

| 论文 | 出处 | 与乘方图/Aurenhammer 的关系 |
|---|---|---|
| Engwirda, Liao. "**'Unified' Laguerre-Power Meshes for Coupled Earth System Modelling**" | IMR 2021, pp. 297–314, [PDF](https://www.meshingroundtable.com/assets/papers/2021/23-Engwirda-compressed.pdf), DOI [10.5281/zenodo.5558988](https://doi.org/10.5281/zenodo.5558988) | **Laguerre/power 图直接生成网格**（球面/地球系统），是"power diagram = 网格"最直白的工程案例 `[摘要]` |
| Hummel. "**Analytic Formula for the Difference of the Circumradius and Orthoradius of a Weighted Triangle**" | IMR 2021, pp. 244–254, [PDF](https://www.meshingroundtable.com/assets/papers/2021/16-Hummel.pdf), DOI [10.5281/zenodo.5559052](https://doi.org/10.5281/zenodo.5559052) | 加权（乘方/Laguerre）三角形的半径差解析公式——**乘方图网格的质量度量基础** `[摘要]` |
| Caplan. "Parallel Four-Dimensional Anisotropic Mesh Adaptation" | IMR 2022, [PDF](https://www.meshingroundtable.com/assets/papers/2022/09-Caplan-compressed.pdf), DOI [10.5281/zenodo.6562412](https://doi.org/10.5281/zenodo.6562412) | 同一作者的 4D 各向异性网格（时空网格） `[摘要]` |
| Marot, Verhetsel, Remacle. "Reviving the Search for Optimal Tetrahedralizations" | IMR 2019, [PDF](https://www.meshingroundtable.com/assets/papers/2019/22-Marot.pdf), DOI [10.5281/zenodo.3653420](https://doi.org/10.5281/zenodo.3653420) | 最优四面体化（对偶于最优 Delaunay） `[摘要]` |
| Estrems, Gargallo-Peiró, Roca. "High-Order Metric Interpolation for Curved R-Adaption by Distortion Minimization" | IMR 2022, DOI [10.5281/zenodo.6562456](https://doi.org/10.5281/zenodo.6562456), [PDF](https://www.meshingroundtable.com/assets/papers/2022/01-Estrems-compressed.pdf) | **度量场（metric field）驱动的 r-自适应**——与"把尺寸场搬到 M 上"直接相关 `[摘要]` |
| Timalsina, Knepley. "Tetrahedralization of a Hexahedral Mesh" | IMR 2023, [PDF](https://www.meshingroundtable.com/assets/papers/2023/14-Timalsina-compressed.pdf) | 六面体→四面体（如果你们上游走 hex 路线） `[摘要]` |

**未找到**："Tetrahedral mesh generation via optimal transport" 这个确切标题的论文。
任务书里这个说法在文献中**没有直接对应的标题**；最接近的就是
Alliez 2005（§2.2/§2.3）、An et al. IMR 2021（§2.4）、Caplan IMR 2021（§2.6）。
**请勿在论文中引用不存在的标题。**

### 2.8 等几何分析 / 三变量样条拟合作为下游消费者

体积参数化的经典下游消费场景（引 CAD 2016 的常见引用）：

- **T. Martin, E. Cohen, M. Kirby. "Volumetric parameterization and trivariate B-spline fitting using
  harmonic functions." SPM '08 (ACM Symposium on Solid and Physical Modeling), 2008.** `[摘要]`
  ——用离散体积调和函数参数化，再拟合**单个三变量 B 样条**，并对股骨做弹性静力学
  **等几何分析（IGA）**。这是"体积参数化 → trivariate spline → IGA"的标准范式。
  （Semantic Scholar 显示被引 268 次）
- **Hughes, Cottrell, Bazilevs. "Isogeometric analysis: CAD, finite elements, NURBS, exact geometry and
  mesh refinement." CMAME 194(39–41):4135–4195, 2005.** `[未验证]`（未逐一核对元数据，请自行 Crossref）
- 该组内部相关：**Prismatic mesh generation based on anisotropic volume harmonic field** —
  Yiming Zhu, Shengfa Wang, Xiaopeng Zheng, Na Lei, Zhongxuan Luo, Bo Chen.
  *Advances in Aerodynamics* 3, 2021. DOI: [10.1186/s42774-021-00065-y](https://doi.org/10.1186/s42774-021-00065-y)
  `[摘要]` ——**这是该组少见的"体网格生成"论文**（棱柱体网格），值得一读。
- **AI for Mesh Generation and Model Preparation: A Review** — Owen, Brown, Chrisochoides, Garimella,
  Gu, Ledoux, **Lei**, Quadros, Ray, Winovich, Zhang.
  IMR 2026, DOI [10.1137/1.9781611979138.1](https://doi.org/10.1137/1.9781611979138.1)；
  预印：*Engineering with Computers* 2026, DOI [10.1007/s00366-026-02379-1](https://doi.org/10.1007/s00366-026-02379-1)
  `[摘要]` ——最新综述，作者含 Gu 和 Lei，**推荐作为你们论文的 related-work 引子**。

---

## 3. 重采样工具清单（许可 / 平台 / Python / 尺寸场 / 安装）

### 3.0 总表

| 工具 | 许可（已核实） | Windows | pip 绑定 | 全局 max volume | **真正的逐点尺寸场 / 背景网格** |
|---|---|---|---|---|---|
| **TetGen** 1.6.x | **AGPL-3.0-or-later**（或向 WIAS 购商业许可） | ✅ MSVC/MinGW | ✅ `tetgen` (wheel) | ✅ `-a#` | ✅ **`-m` + `.mtr`（PLC 节点或背景网格）** |
| **Mmg / Mmg3d** | **LGPL-3.0-or-later** | ✅ | ✅ `pymmg` / `mmgpy` (wheel) | ✅ `-hmax` | ✅ **metric `.sol`（各向异性）** |
| **Gmsh** | **GPL-2.0-or-later** | ✅ 官方 exe + wheel | ✅ `gmsh` (wheel) | ✅ `Mesh.MeshSizeMax` | ✅ **Field 栈；PostView 需 list-based `.pos`** |
| **CGAL `Mesh_3`** | **GPL-3.0-or-later OR 商业许可** | ✅（header-only，需 Boost/GMP） | ⚠️ `pycgal` 很窄 | ✅ `cell_size` | ✅ `Mesh_domain_field_3` / 函数对象（逐点） |
| **fTetWild** | **MPL-2.0** | ✅（需 mpir.dll） | ✅ `wildmeshing` (MPL-2.0) / `pytetwild` | ✅ `--edge-length-abs/-r` | ✅ **`--bg-mesh`（.msh + 节点 "values" 场）** |
| **TetWild**（原版） | **GPL-3.0** | ✅ | 间接 | ✅ | `[未验证]` |
| **PyMesh** | ⚠️ **仓库无 LICENSE 文件** | ❌ 主要 Linux/Docker | ✅ 但需自编译 | 透传 TetGen | 透传 TetGen/MMG |
| **Netgen/NGSolve** | **LGPL-2.1-only** | ✅ | `netgen-mesher` | ✅ | ✅ meshsize 函数 |
| **SeismicMesh** | GPL-3.0-or-later (PyPI) | ✅ | ✅ | ✅ `hmin/hmax` | ✅ `SizeFunction` 类（可接背景网格） |
| **Instant Meshes** | **BSD-3-Clause**（`LICENSE.txt` 实测） | ✅ 预编译 zip | ❌（GUI/CLI） | — | — **注意：这是 hex/quad，不是 tet** |

### 3.1 TetGen ⭐（推荐主用）

**版本 / 许可（逐字核实）**`[摘要]`
- 最新 **v1.6.1, 2026-08**（我还从 `main` 分支 clone 并成功编译出 1.6.1，见 §3.1 末）。
- 许可历史（**重要，容易搞错**）：`refs/` 里我保存了 GitHub README 的原文：
  - **v1.4.3 (2011)**：MIT license **with noncommercial clause**
  - **v1.5.0 (2013) / v1.5.1 (2018) / v1.6.0 (2020) / v1.6.1 (2026)**：**AGPLv3**
  - 官方 **Licensing FAQ**（[wias-berlin.de/software/tetgen/FAQ-license.jsp](https://www.wias-berlin.de/software/tetgen/FAQ-license.jsp)）
    逐字："TetGen is free software: you can redistribute it and/or modify it under the terms of the
    **GNU Affero General Public License** … either version 3 of the License (AGPLv3), or (at your option)
    any later version."；商业许可 "for **one of the versions 1.4 – 1.6.0**"，向 `tetgen@wias-berlin.de` 询价。
  - **工程结论：闭源商业产品里不能用 AGPL 版 TetGen，必须买许可或换 LGPL/MPL 的工具。**
- 推荐稳定版：**v1.5.1 (2018)**（官方标注 "recommended as most stable version"）；
  v1.6.0 "most recent version with some rough edges"。
- 主仓库在 **Codeberg**（`https://codeberg.org/TetGen/TetGen`），GitHub 是镜像。

**平台**：C++，只用标准 C 库，MSVC/gcc/Intel 均可编译，32/64 位，Unix/Linux/Windows/macOS 全支持
（TetGen 手册 §2.1，本地文本 `refs/tetgen-manual-1.6.txt`）`[原文]`。

**Python 绑定**
- `pip install tetgen` → **`tetgen` 0.8.4**（pyvista 团队维护，**自带编译好的 `_tetgen.pyd`，无需外部 tetgen.exe**）。
  要求 Python ≥3.10；有 Windows/macOS/Linux wheel（cp310/311 + cp312-abi3）。`[实测]`
- `pip install meshpy` → **`meshpy` 2026.1.1**，含 TetGen 的 SWIG 绑定（`meshpy.tet.MeshInfo`，
  **有 `load_mtr` / `point_metric_tensors`**，可传尺寸场）。`[实测]`
- conda：`conda install -c conda-forge tetgen`（**1.6.0, AGPL-3.0-or-later**）或
  `conda install -c conda-forge python-tetgen`（0.8.4）。`[实测]`（anaconda.org API）

**尺寸场机制（本报告的核心结论之一，来自 TetGen 手册 §1.2.6 / §4.2.4）**`[原文]`
TetGen 手册逐字说明它支持四种定义 sizing function 的方式，且**可以同时使用，TetGen 自动取最小**：
> - "By default, TetGen uses the **local feature size** of the PLC."
> - "One can apply a **maximum volume bound (the `-a` switch)** on every tetrahedra … It is the same as
>   defining a global constant sizing function. For a domain consisting of multiple sub-domains (materials),
>   the volume bound can vary in each sub-domain, see the `.poly` and `.smesh` file formats."
> - "One can define a mesh sizing function **directly on the input PLC**. In this way, to each vertex of the
>   PLC, a value for the mesh sizing function is to be assigned, see the **`.mtr` file format**."
> - "One can use a **background mesh** to define a sizing function. The background mesh can be any tetrahedral
>   mesh. Its underlying space must cover the input PLC. To each mesh node of the background tetrahedral mesh,
>   a value for the mesh sizing function is assigned, which contains the desired edge length at that location
>   in the PLC."
> - "All the above ways of defining a sizing function can be used at the same time. **TetGen will automatically
>   choose the smallest mesh element size.** In the last two ways … it is possible to set the size to zero.
>   In this case, the mesh element size at this location is ignored."

`-m` 的文件命名（手册 §4.2.4）`[原文]`：
- 尺寸直接定义在**输入 PLC 的节点**上 → 文件名 `xxx.mtr`（`xxx` = 输入 `.poly`/`.smesh`/`.node` 的 basename）；
- 尺寸定义在**背景网格**的节点上 → 背景网格文件为 `xxx.b.node`、`xxx.b.ele`，尺寸文件为 `xxx.b.mtr`。
- 尺寸值 = **期望边长**（desired edge length），**当前只支持各向同性**。
- `-m` 必须和 `-q` 一起用（"used in mesh refinement, i.e., together with the `-q` switch"）。

完整 meshing 流程（手册 §1.3，12 步）第 4 步就是背景网格插值 `[原文]`：
> "4. Read the background mesh from (`.b.node`, `.b.ele`, `.b.mtr` ...) files (if it is provided) and
> interpolate the mesh element size from the background mesh to the current mesh ( **`-m`** )."

**`.mtr` 格式**（手册 §5.2.8，本地文本可查）：第一行 `节点数 1`，随后每行一个尺寸值。

**Python 侧怎么传尺寸场（实测结论）**
- `tetgen.TetGen(pts, tris).tetrahedralize(bgmesh=grid, ...)` 要求 **pyvista UnstructuredGrid
  只含线性四面体**，尺寸放在 **`point_data['target_size']`**。签名文档逐字：
  > "`bgmesh` : pyvista.UnstructuredGrid — Background mesh to be processed. Must be composed of only linear
  > tetra with the sizing contained in the `point_data` of the mesh within the `'target_size'` key."
  > "`bgmeshfilename` : str — Filename of the background mesh with the target size associated with the nodes."
  （`pip install tetgen` 后 `inspect.getsource(TetGen)` 可见）
- ⚠️ **实测踩坑（必须知道）**：我用 `pip` 版 `tetgen` 0.8.4 调用 `bgmesh=` 时**进程直接崩掉
  （Windows 上 EXIT = -1073741819 = 0xC0000005 access violation）**；
  用 `switches="pq1.4m"` + `bgmeshfilename=` 时 TetGen 本体跑完了
  （stdout 里明确出现 `Size interpolating seconds: 0.09`，证明 **`-m` 路径确实被触发**），
  但随后包装层在回传结果时同样崩掉。**结论：Python 包装器的背景网格路径在 Windows 上不可靠。**
  生产上请用**命令行 tetgen.exe**（见下面构建方法）或 **meshpy**。

**从源码构建（我在这台机器上实测成功）**`[实测]`
```powershell
git clone --depth 1 https://github.com/TetGen/TetGen.git E:\panyingyun\smartmm\vendor\tetgen
cmake -S E:\panyingyun\smartmm\vendor\tetgen -B E:\panyingyun\smartmm\vendor\tetgen\build `
      -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build E:\panyingyun\smartmm\vendor\tetgen\build -j 4
# → build\tetgen.exe (3.99 MB) + build\libtet.a ；`tetgen -h` 输出 "TetGen 1.6.1, July 2026"
```
依赖：`cmake`（本机 `C:\Program Files\CMake\bin\cmake.exe`）+ `g++`/`gcc`
（本机 `C:\TDM-GCC-64\bin\g++.exe`，GNU 10.3.0）。**两者都在，能编过。**

**⚠️ `.smesh` 输入格式的两个致命陷阱（我实测踩到）**`[实测]`
1. **小平面那一节必须写 `3 i j k`，即每行开头要有顶点个数 `3`**。写成 `i j k` 会报
   `Error: Wrong number of vertex in facet 1`（exit 10）。facet 顶点是 **0-based**。
2. **法向必须朝外，且不能有重复极点顶点**。我的第一版 UV 球把两个极点各复制了 `n_lon` 份
   → 出现退化三角形 → **TetGen 不报错，静默输出 0 个四面体**（"Mesh tetrahedra: 0"），
   极易误判为"算法没跑"。修好极点 + 统一朝外后正常出网格。
   诊断方法：算 `Σ_facet (n · centroid)`，必须显著为正。

### 3.2 Mmg / Mmg3d ⭐（推荐的商用友好 remesher）

**元数据 / 许可**`[原文]`
- 官网 [https://www.mmgtools.org](https://www.mmgtools.org)；主仓库
  [https://github.com/MmgTools/mmg](https://github.com/MmgTools/mmg)。
- `LICENSE` 逐字（我抓了 raw 文本）：
  > "mmg is free software: you can redistribute it and/or modify it under the terms of the
  > **GNU Lesser General Public License** as published by the Free Software Foundation, either version 3
  > of the License, or (at your option) any later version."
  → **LGPL-3.0-or-later**。Copyright © Bx INP/Inria/UBordeaux/UPMC, 2004–。
  **商用友好（动态链接即可）。**
- 最新活动：仓库 2026-09-25 仍有 push。

**平台**：Windows / Linux / macOS，CMake 构建，有 CI。

**Python 绑定（实测 PyPI 上有两个，选哪个要谨慎）**
| 包 | 版本 | 许可 | wheel 情况 | 说明 |
|---|---|---|---|---|
| `pymmg` | 1.0.0 | **LGPL-3.0-or-later** | `py3-none-win_amd64` / manylinux / macOS（**但 sdist 存在，说明需编译**） | gnikit 维护，[github.com/gnikit/pymmg](https://github.com/gnikit/pymmg) |
| `mmgpy` | 0.17.0 | **MIT**（绑定层） | win_amd64 / manylinux / macOS，cp310–cp314 | [github.com/kmarchais/mmgpy](https://github.com/kmarchais/mmgpy) |
| `mmg` (PyPI) | 2.2.0 | — | ❌ | **这是 i18n 文档生成库，不是 Mmg！常见的名字陷阱** |

conda-forge 上的包名是 **`mmgsuite`**（5.8.0）而不是 `mmg`；另有 `conda-forge/mmgpy` 0.17.0。`[实测]`

```bash
pip install mmgpy            # 推荐：有 win wheel + MIT
pip install pymmg            # LGPL，可能需要本地编译器
conda install -c conda-forge mmgsuite mmgpy
```

**尺寸场机制**：Mmg 用 **(n+1)²/2 分量的各向异性 metric 场**（3D 时 6 个分量/点），
通过 `-sol`/metric 文件或 API `MMG3D_Set_solSize` + `MMG3D_Set_solsAtVertices` 传入；
也支持 `-hmin/-hmax/-hsiz/-hausd` 等标量控制。
`[未验证]`：我**没能**从 Mmg 官方文档逐字核实 `-ls` 与 `.sol` 文件格式的确切字段
（`mmgtools.org` 与 `gitlab.onelab.info` 都被反爬/Anubis 挡住，`www.mmgtools.org/page/11` 未取到正文）。
**上线前请用 `mmg3d -h` 和官方 `mmg3d` 文档核对一遍。**

### 3.3 Gmsh ⭐（最容易做出"背景网格尺寸场"）

**元数据 / 许可**`[实测]`
- 版本 **4.15.2**；官方 wheel 覆盖 `win_amd64` / `manylinux_2_24_x86_64` / `macosx`；
  `requires_python = None`。
- 许可：PyPI 元数据 `GPLv2+`；conda-forge `GPL-2.0-or-later`。→ **GPL-2.0-or-later**。

```bash
pip install gmsh meshio
conda install -c conda-forge gmsh
```

**尺寸场机制（我实测跑通两条路）**`[实测]`

**路线 A：Field 栈（不需要背景网格）——4 行就能用，且已验证生效**
```python
import gmsh
gmsh.initialize(); gmsh.model.add("ball")
gmsh.model.occ.addSphere(0,0,0,1.0, tag=1); gmsh.model.occ.synchronize()

gmsh.model.mesh.field.add("Box", 1)          # 一个立方体区域内的目标尺寸
for k,v in dict(VIn=0.05, VOut=0.30,
                XMin=-0.4,XMax=0.4, YMin=-0.4,YMax=0.4, ZMin=-0.4,ZMax=0.4).items():
    gmsh.model.mesh.field.setNumber(1,k,v)

gmsh.model.mesh.field.add("Threshold", 2)    # 平滑过渡
gmsh.model.mesh.field.setNumber(2,"InField",1)
for k,v in dict(SizeMin=0.05, SizeMax=0.30, DistMin=0.1, DistMax=0.7).items():
    gmsh.model.mesh.field.setNumber(2,k,v)

gmsh.model.mesh.field.setAsBackgroundMesh(2)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)  # ← 这三个必须关掉，否则曲率会覆盖你的场
gmsh.model.mesh.generate(3)
```
**实测结果**（球半径 1，脚本 `refs/gmsh_bgmesh_demo.py`，单元"尺寸"用体积等价边长
`h_eff = (6√2·V_tet)^{1/3}` 度量）：

| 区域 | 单元数 | h_eff 中位数 | p10 | p90 |
|---|---|---|---|---|
| core `r<0.25` | 2 247 | **0.0623** | 0.0552 | 0.0699 |
| mid `0.25–0.5` | 13 273 | 0.0627 | 0.0553 | 0.0739 |
| mid `0.5–0.75` | 5 640 | 0.1018 | 0.0611 | 0.1673 |
| shell `r>0.75` | 5 255 | **0.1542** | 0.1355 | 0.1775 |

→ 网格数 26 415，体积 4.16191（精确 `4π/3 = 4.18879`，误差 0.64%）。
**尺寸场确实被遵守**（从内到外单调递增），且**体积守恒到 0.6%**。

**路线 B：真·背景网格**。Gmsh **拒绝**把"基于当前网格的 view"当背景网格，报错逐字：
> `Cannot use view based on current mesh for background mesh: you might want to use a list-based view (.pos file) instead`

所以必须写一个 **list-based `.pos`**（`SP(x,y,z){s};` 形式）再 `gmsh.merge(...)` 回来，
用 `gmsh.view.getTags()[-1]` 拿 tag，再 `field.add("PostView", ...)` + `setAsBackgroundMesh`。
**实测通过**，脚本 `refs/gmsh_bgmesh_demo.py`。结果（半径线性 `h: 0.05→0.30`）：

| 区域 | 单元数 | h_eff 中位数 | 处方值 |
|---|---|---|---|
| core `r<0.25` | 501 | 0.0922 | ~0.05–0.11 |
| mid `0.25–0.5` | 676 | 0.1743 | ~0.11–0.175 |
| mid `0.5–0.75` | 654 | 0.2342 | ~0.175–0.24 |
| shell `r>0.75` | 908 | 0.2836 | ~0.24–0.30 |

→ 与处方场**逐段吻合**，证明"背景网格 + 逐点尺寸"这条路可用。
（这条路正是"保体积映射 → 在 B 里造背景网格 → 搬到 M"所需要的。）

**可用于尺寸场的 Field 类型**（Gmsh 文档口径，`[未验证]`——官方 texinfo 文档站被 Anubis 挡住，
我只从 Python API 与实测确认了 `Box`/`Threshold`/`PostView`）：
`Box`, `Distance`, `Threshold`, `MathEval`（表达式！最灵活）, `MathEvalField`, `Frustum`,
`Attractor`, `Curvature`, `PostView`, `Min`, `Max`, `Restrict`, `Param`, `Constant`, `External`.
**`MathEval` 是你们真正需要的**：可以直接把"保体积映射的 Jacobian"写成坐标的解析表达式。
上线前请用 `gmsh.model.mesh.field.list()` 和 `gmsh -help` 核对。

### 3.4 CGAL `Mesh_3`

**许可（从源码文件头逐字核实）**`[实测]`
- `Mesh_3/package_info/Mesh_3/license.txt` 内容就是一行：**`GPL (v3 or later)`**
- `Mesh_3/include/CGAL/Mesh_criteria_3.h` 文件头 SPDX：
  **`SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial`**
- conda-forge `cgal` 的 license 字段：`GPL-3.0-or-later`（最新 6.0.1）。
- `Mesh_3/package_info/Mesh_3/description.txt`："The Mesh_3 package is a mesh generator to build
  **isotropic simplicial meshes** of 3D domains."
- **结论：`Mesh_3` 是 GPL，不是 LGPL。** 闭源商用必须买 CGAL 商业许可，或改用净室实现。

**平台 / 绑定**：header-only，需要 Boost、GMP/MPFR；MSVC 可编。
Python：**没有官方绑定**。`pycgal`（0.4.3，[gitlab.com/brgm/geomodelling/geometry/pycgal](https://gitlab.com/brgm/geomodelling/geometry/pycgal)）
有 win_amd64 wheel，但只包了很小一部分 CGAL 功能（几何基元、triangulation），
**`Mesh_3` 的 criteria API 我未在其中确认** `[未验证]`。

**尺寸场机制**：`CGAL::Mesh_criteria_3` 接受
`edge_size` / `facet_size` / `facet_distance` / `cell_size` / `cell_radius_size` / `cell_radius_edge_ratio`，
每个都可以是**函数对象**（`FT (const Point&)`）或 `CGAL::Mesh_constant_domain_field_3<Gt,Index>`，
于是可以做**逐点尺寸场**；域用 `Mesh_domain_with_polyline_features_3` 支持 sharp features。
`[未验证]`：我没能抓到 `Mesh_criteria_3.h` 的模板参数注释全文（只抓到文件头 + 类声明），
**API 细节请查 `doc.cgal.org/latest/Mesh_3/`。**

**结论：不推荐作为首选**——GPL 许可 + 无 Python + 编译重，对"快速迭代出网格"性价比低。

### 3.5 fTetWild / tetwild 家族 ⭐（推荐的"最鲁棒"方案）

**元数据 / 许可（GitHub API 精确核对）**`[实测]`
| 仓库 | SPDX（GitHub API） | 说明 |
|---|---|---|
| `wildmeshing/fTetWild` | **MPL-2.0** | 2026-09-23 仍有 push；607 stars |
| `Yixin-Hu/TetWild`（原版） | **GPL-3.0** | 723 stars，2023 后停更 |
| `wildmeshing/wildmeshing-python` | NOASSERTION → `LICENSE` 内容为 **MPL-2.0** | 57 stars |

> ⚠️ 任务书猜测 fTetWild 是 MIT——**不对，是 MPL-2.0**（GitHub API 的 license 字段就是 `MPL-2.0`；
> 我同时抓了 `wildmeshing-python` 的 `LICENSE`，正文是 Mozilla Public License Version 2.0，
> 并列出 bundle 的依赖许可：triwild MPL2、tetwild MPL2、pybind11 BSD、aabbcc zlib、cli11 BSD-3、
> geogram BSD-3、libigl MPL2、nlopt MIT、spdlog MIT、tbb Apache-2.0、LM6 LGPL、rply MIT、zlib MIT）。
> **MPL-2.0 是商用友好的（file-level copyleft），比 TetGen 的 AGPL 好得多。**

**论文**
- Yixin Hu, Teseo Schneider, Bolun Wang, Denis Zorin, Daniele Panozzo.
  "**Fast tetrahedral meshing in the wild**." *ACM Trans. Graph.* **39(4)**, Article 117, 2020, 18 pp.
  DOI: [10.1145/3386569.3392385](https://doi.org/10.1145/3386569.3392385)
- Yixin Hu, Qingnan Zhou, Xifeng Gao, Alec Jacobson, Denis Zorin, Daniele Panozzo.
  "**Tetrahedral meshing in the wild**." *ACM Trans. Graph.* **37(4)**:1–14, 2018.
  DOI: [10.1145/3197517.3201353](https://doi.org/10.1145/3197517.3201353)

**Python 绑定**
```bash
pip install wildmeshing     # 0.4.1，MPL-2.0；**只发了 Linux/macOS wheel，没有 win_amd64！**
pip install pytetwild       # 0.4.2，**有 win_amd64 wheel**（cp310–cp312）
```
`wildmeshing` 0.4.1 的 wheel 列表（我逐条查过 PyPI）：`macosx_13_0_x86_64`、`macosx_14_0_arm64`、
`manylinux_2_17_x86_64`——**没有 Windows wheel**。`environment.yml` 里写的是
`conda install -c conda-forge wildmeshing`（那里应该有 Windows 构建）`[实测]`。

**API（实测，从源码 `src/binding.cpp` / `tetrahedralize.cpp` / `test/tet_test.py`）**
```python
import wildmeshing as wm, numpy as np
# 路线 1：文件
wm.tetrahedralize("in.stl", "out.msh", mute_log=True, stop_quality=1000)
# 路线 2：数组 + 类
tetra = wm.Tetrahedralizer(stop_quality=1000)
tetra.set_mesh(V, F)          # V:(n,3) float, F:(m,3) int
tetra.tetrahedralize()
VT, TT, _ = tetra.get_tet_mesh()
# 还有 wm.boolean_operation(json_csg_tree, out, stop_quality=...)（Windows 上测试被跳过）
```
`Tetrahedralizer.__init__` 的真实签名（源码逐字）：
```cpp
Tetrahedralizer(double stop_quality, int max_its, int stage, int stop_p,
                double epsilon, double edge_length_r,
                bool skip_simplify, bool coarsen)
```
内部映射到 fTetWild 的 `params.stop_energy / max_its / stage / stop_p / eps_rel / ideal_edge_length_rel`。

**尺寸场机制（关键，任务书猜测"不支持一般密度场"是错的）**`[原文]`+`[实测]`
fTetWild `src/main.cpp` 里有：
```
std::string background_mesh = "";
command_line.add_option("--bg-mesh", background_mesh,
                        "Background mesh for sizing field (.msh file).")->check(CLI::ExistingFile);
...
if (!background_mesh.empty()) {
    PyMesh::MshLoader mshLoader(background_mesh);
    V_in   = mshLoader.get_nodes();
    T_in   = mshLoader.get_elements();
    values = mshLoader.get_node_field("values");     // ← 节点标量场名固定为 "values"
}
if (V_in.rows()!=0 && T_in.rows()!=0 && values.rows()!=0) {
    params.apply_sizing_field      = true;
    params.V_sizing_field          = V_in;
    params.T_sizing_field          = T_in;
    params.values_sizing_field     = values;
}
```
实现（`src/MeshImprovement.cpp`）：用背景四面体的**重心坐标插值**得到点 `p` 处的值，然后
```
p.sizing_scalar = value / mesh.params.ideal_edge_length;   // 第1346行
```
后续 `sizing_scalar` 通过 `int n = -int(log2(v.sizing_scalar) - 0.5)` 转成**2 的幂次细分层数**
（第659行）——**这意味着 fTetWild 的尺寸场控制是"每点一个倍率，量化到 2 的幂"**，
比 Gmsh/TetGen 的连续尺寸场粗糙，但鲁棒性极强。
- 其它相关开关：`--edge-length-abs`（绝对值）、`-l/--ideal-edge-length-r`（默认 `1/20`，相对 bbox 对角线）、
  `-e/--epsr`（默认 `1e-3`）、`--use-floodfill`、`--use-general-wn`、`--use-input-for-wn`、
  `--disable-filtering`、`--epsr-tags`（逐面 envelope 大小）、`--tetgen`（同时跑 tetgen）。
- `.msh` 里字段名必须是 **`values`**（`PyMesh::MshLoader::get_node_field("values")`）。

**输入要求**：`.off/.obj/.stl/.ply` 三角面片，**朝向必须一致且有正确法向**
（原文："The orientation of input faces is as important as the position of their vertices"）。
输出默认 `.msh`（Gmsh 可看），并把最小二面角写成元素标量场。

**构建（Windows 注意）**：需 `gmp`；Windows 上用 `conda install -c conda-forge mpir`，
并把 `mpir.dll`（例如 `<conda_dir>\Library\bin`）拷到 `FloatTetwild_bin.exe` 同目录。
若 cmake 找不到 gmp，设 `GMP_INC` / `GMP_LIB` 环境变量。

### 3.6 PyMesh ⚠️（不推荐，许可不清）

- 仓库 [github.com/PyMesh/PyMesh](https://github.com/PyMesh/PyMesh)，Qingnan Zhou (NYU)。
- ⚠️ **GitHub API 的 `license` 字段是 `null`，仓库根目录没有 `LICENSE`/`COPYING`/`LICENSE.TXT`
  （我逐个 404 试过），README 里也没有许可声明。** 这在法务上是**明确的未授权**状态。
  **如果你的团队要做商业交付，请不要用 PyMesh。**（社区长期以为它是 MPL-2.0，但仓库里没有任何证据。）
- 平台：主要 Linux/macOS；**推荐用 Docker `pymesh/pymesh`**；Windows 原生构建没有官方支持。
- 功能定位：几何处理平台（读写、布尔、自交检测、四面体化封装 TetGen、外包 Quartet 的
  "isosurface stuffing"）。**它自己没有独立的尺寸场机制**，透传 TetGen 的 `-a`/`-m`。
- 安装（若要试）：`docker run -it pymesh/pymesh`；或源码 `git clone` + `git submodule update --init` + CMake 构建。

### 3.7 其它：Netgen、SeismicMesh、Instant Meshes

- **Netgen / NGSolve**：PyPI `netgen-mesher` 6.2.2607；conda-forge `netgen` **LGPL-2.1-only**
  （商用友好）。支持逐点 meshsize 函数（`meshsize = lambda p: ...`）。适合做"球内自适应网格"。
- **SeismicMesh**：PyPI `SeismicMesh` **3.6.2，GPL-3.0-or-later**（注意是 GPL）。
  有 `SizeFunction` 类，支持 `hmin/hmax/grade/apply`，**可接背景网格**；conda-forge 上我**没查到**
  （`anaconda.org/search?name=seismicmesh` 返回空）。`pip install SeismicMesh`。
- **Instant Meshes**（作者 wjakob）：**这是 quad/hex 重网格，不是 tet**。
  许可 **BSD-3-Clause**（我抓了 `LICENSE.txt` 全文：三条款 BSD，带 "Enhancements" 授权条款）；
  有官方 **Windows 预编译 zip**（`instant-meshes-windows.zip`），需 CMake + VS2015 自编；
  **无 Python 绑定**（GUI + CLI）。论文：Jakob, Tarini, Panozzo, Sorkine-Hornung,
  "Instant Field-Aligned Meshes", ACM TOG (SIGGRAPH Asia) 2015。
  **如果最终目标是四面体网格，不要走这条路**；除非上游要 hex。
- **optimesh**（PyPI 0.12.6）：网格优化/平滑，可用于 tet 质量后处理 `[未验证许可证]`。

### 3.8 conda 环境一行装（推荐组合）

```bash
# 商用友好组合（无 AGPL/GPL）
conda create -n tet python=3.12 -c conda-forge
conda activate tet
conda install -c conda-forge mmgsuite mmgpy          # LGPL remesher
pip install pytetwild                                 # MPL-2.0，Windows 有 wheel
pip install gmsh meshio numpy scipy pyvista           # GPL-2（工具链可接受则用）

# 若团队能接受 AGPL（研究/开源）
conda install -c conda-forge tetgen python-tetgen
pip install tetgen meshpy

# 可选：CGAL 路线（GPL）
conda install -c conda-forge cgal boost-cpp
```

### 3.9 许可总结（给法务看的一页）

| 工具 | 许可 | 闭源商用免费？ |
|---|---|---|
| **Mmg / mmgsuite / pymmg** | LGPL-3.0-or-later | ✅ 是（动态链接） |
| **fTetWild / wildmeshing** | MPL-2.0 | ✅ 是 |
| **Instant Meshes** | BSD-3-Clause | ✅ 是（但非 tet） |
| **Netgen** | LGPL-2.1-only | ✅ 是 |
| **mmgpy**（绑定层） | MIT | ✅ 是 |
| **TetGen ≥1.5** | **AGPL-3.0-or-later** | ❌ 否（需向 WIAS 购买） |
| **TetGen 1.4.3** | MIT + noncommercial | ❌ 否 |
| **CGAL `Mesh_3`** | **GPL-3.0-or-later OR 商业** | ❌ 否 |
| **Gmsh** | GPL-2.0-or-later | ❌ 否 |
| **TetWild（原版）** | GPL-3.0 | ❌ 否 |
| **SeismicMesh** | GPL-3.0-or-later | ❌ 否 |
| **PyMesh** | **无许可文件** | ❌ 否（未授权） |

---

## 4. 端到端方案设计（3 条可落地方案）

**共同前提（`P0`）**：已有
- `f = T ∘ ω : M → B`（保体积映射到单位球）及其**逆** `f⁻¹`；
- 若用半离散 OT 实现，则手上还有 **OT 单元结构**：
  - 站点 `{p_i}`（球内 / 参数域上的目标点）
  - 权重 `{h_i}`（Brenier 势的支撑平面高度）
  - 乘方图 cell `{W_i(h)}`（Laguerre 单元）
  - **逆映射的廉价实现**：`f⁻¹(p_i) ≈ centroid(W_i)`（这是 IMR 2021 Algorithm 1 的最后一行，注意方向）

### 4.0 三个必须事先想清楚的失败模式

| 失败模式 | 机理 | 缓解 |
|---|---|---|
| **F1 映射只是"近似"保体积** | 数值 OT 的 `ε`（`max_i|w_i−ν_i|/ν_i`）不为 0；体积畸变直接变成单元尺寸畸变 `h ∝ J^{1/3}` | 把"目标尺寸"和"实际尺寸"**解耦**：只用映射做**采样/尺寸场**，让下游 mesher 自己保证质量（Recipe A/B）。或者用 `-a`/`maxvolume` 兜底。 |
| **F2 逆映射求值贵** | 连续逆映射（Newton 迭代求 Brenier 势梯度）每次 O(迭代数 × 单元定位) | **不要逐点求连续逆映射**。用 cell 质心（`f⁻¹(p_i)=centroid(W_i)`，O(1) 每 cell）；或建 KD-tree / BVH 缓存在 `{p_i}→{centroid}` 上做插值。 |
| **F3 边界不兼容** | 球内采样的边界点拉回来**不落在** `∂M` 上（数值误差 + 离散化），产生"壳层"或缝隙 | **永远把 `∂M` 的原始三角面片作为硬约束交给 mesher**（fTetWild / TetGen `-p` / Gmsh / CGAL 都支持）。**不要**试图用拉回来的点当边界。 |
| **F4 sliver** | 3D 中 `E_CVT` 不惩罚 sliver（Alliez 原文）；拉回映射会放大畸变 | 用 `E_ODT` 类方法（Alliez）或 Delaunay refinement（TetGen `-q`）；fTetWild 内建质量优化。 |
| **F5 sharp feature 丢失** | 保体积映射把尖角磨圆；采样不到位 | 在 `∂M` 上标出 sharp edge/corner（二面角阈值），作为 **polyline feature** 传给 mesher；Alliez 原文用 "boundary vertex jittering" 去 sliver。 |
| **F6 奇异集** | OT 映射在不连续处（singular set）会产生"撕裂"，把单元压扁 | 用 `robust and accurate optimal transportation map by self-adaptive sampling` 的思路（Wang et al., FITEE 2021）加密奇异集附近；或只采样不拉网格（Recipe A）。 |

### 4.1 Recipe A ⭐ 推荐：**球内采样 + 在 M 里重新四面体化**

**核心思想**：保体积映射只用来**生成 M 内部的好点集**，四面体化交给成熟的鲁棒 mesher。
这样 F1/F2/F3/F4/F5 全部变成 mesher 的责任，而 mesher 已经解决了这些问题。

**需要的输入数据**
- `∂M` 的三角面片（原始，保留 sharp feature 标注）
- 球内采样点集 `{q_j} ⊂ B`（均匀：规则网格；自适应：CVT 或 Poisson-disk，密度 ∝ 用户 `ρ`）
- 逆映射求值器 `f⁻¹`（用 OT cell 质心，或缓存的 `{p_i → centroid(W_i)}` + KD-tree 插值）

**算法**
```
1  在 B 内按用户密度 ρ 生成采样点 {q_j}：
   均匀 → 球坐标规则网格（r 方向按 r^{1/3} 分层，保证体积均匀！）
   自适应 → 在 B 上跑带容量约束的 CVT（即半离散 OT，§2.4 Algorithm 1）
           或 Poisson-disk sampling，半径 ∝ ρ(q)^{-1/3}
2  对每个 q_j 求 x_j = f⁻¹(q_j)
   —— 廉价实现：用 {p_i} 的 KD-tree 找最近的 OT cell，取 centroid(W_i)
   —— 精确实现：f⁻¹(q) = ∇u*(q)，Newton 解 ∇u(x)=q（只在需要时做）
3  去除落在 ∂M 外侧/过近的点（|signed_dist(x_j, ∂M)| < 0.5·h_local 的丢掉）
4  把 X = {x_j} 作为 Steiner 点，与 ∂M 一起喂给 mesher：
   (a) fTetWild（推荐）：
       FloatTetwild_bin -i M_surface.obj -o out.msh \
                        --bg-mesh size_field.msh -e 1e-3 --use-floodfill
   (b) TetGen：把 X 写进 .a.node，用 -i 开关插入
       tetgen -pq1.4a<maxvol>i M.smesh
   (c) Gmsh：把 X 写成 .geo 里的 Point + 用 Field 控制尺寸
5  质量后处理：TetGen -O / optimesh / Mmg 的 -optim
```

**关键工程细节：球内均匀采样必须按 `r^{1/3}` 分层！**
球内**体积均匀**要求半径分位数 `r_k = R·(k/N)^{1/3}`；如果按 `r` 等距分层，
得到的点在球心附近会**过密**（体积偏差 ∝ r²）。这一步错了，整个"保体积"的收益就没了。

**数据要存**
- `∂M` 三角面片 + feature 标记：O(n_surf)
- `{p_i}`、`{h_i}`、`{centroid(W_i)}`：O(N_cell)，N_cell = 目标单元数
- KD-tree over `{p_i}`：O(N_cell log N_cell)
- **不需要存** `f` 的解析形式，也不需要存球内的网格

**复杂度**
- 采样：O(N_cell)
- 逆映射（cell 质心路径）：O(N_cell log N_cell)（建 KD-tree）+ O(N_sample log N_cell)（查询）
- 四面体化：fTetWild 约 O(N log N)～O(N^1.2)；TetGen Delaunay refinement 类似。
  实测：球面 1302 点 → `-pq1.4a0.002` 出 10321 点 / 43955 四面体，**总运行 0.81 秒**（本机 1.6.1）`[实测]`
- **总计：分钟级（10^5–10^6 单元）**

**优点**：最鲁棒；边界天然共形（因为 `∂M` 是输入）；sharp feature 可保留；避开 F1–F5。
**缺点**：Steiner 点可能被 mesher 拒掉/移动；得到的网格与 `f` 的"体积保真"不再严格对应
（如果你的应用依赖体积精确，需要再加一步 §4.3 的校验）。

### 4.2 Recipe B：**在 B 里造背景网格尺寸场，搬到 M，交给 mesher**

**核心思想**：不改采样策略，而是把"用户密度 `ρ`"翻译成 M 上的**逐点目标边长场**
`h_M(x)`，让 mesher 用它做自适应。

**尺寸场的推导（分两种情况）**
- **保体积时**（`J_f ≡ 1`）：`h_M(x) = h_B(f(x))`，其中 `h_B` 是球内目标边长。
  若用户给的是密度 `ρ_B(q)`（单位体积的单元数比例），则
  `h_B(q) ∝ ρ_B(q)^{-1/3}`（因为单元数 ∝ 1/h³）。
- **近似保体积时**（Jacobi `J_f(x) = det(∂f/∂x)`）：
  `h_M(x) = h_B(f(x)) · J_f(x)^{1/3}`。
  ⚠️ 这是我用换元法推的工程公式，`[未验证]`——**请用合成算例验证**：
  取已知解析映射（例如 `M = B`、`f = Id`，或一个已知 Jacobian 的椭球→球映射），
  比较处方 `h` 与实际 `h_eff` 的中位数比值，应该接近 1。

**算法**
```
1  在 B 里构造"尺寸场载体"：
   (a) 一个四面体背景网格 T_B（可以用 Gmsh 在球里生成，尺寸场就是 h_B）
       每个节点带一个标量 = 目标边长 h_B
   (b) 或一张散点 .pos（Gmsh） / 一组 (point, h) 元组
2  把载体从 B 搬到 M：对每个节点 q_k ∈ T_B，令 x_k = f⁻¹(q_k)，
   得到 M 上的背景网格 T_M（连通性照抄 T_B！保体积映射是微分同胚，拓扑不变）
   节点值照抄（尺寸是"边长"量，保体积时不变；近似保体积时乘 J^{1/3}）
3  交给 mesher：
   - TetGen：写 xxx.b.node / xxx.b.ele / xxx.b.mtr，跑 tetgen -pq1.4m
   - Gmsh：写成 list-based .pos（SP 格式）+ PostView field（§3.3 路线 B，实测通过）
   - fTetWild：写成 .msh（节点场名必须是 "values"）+ --bg-mesh
   - Mmg：写成 metric .sol，跑 mmg3d
4  质量后处理
```

**为什么连通性可以直接照抄？** 因为 `f` 是微分同胚（保体积 + 共形 ⇒ 双射且光滑），
球里一个合法的四面体网格拉回去**拓扑上仍然合法**（可能退化但不会翻转，除非数值失败）。
但**几何质量会变差**，所以第 4 步的后处理不能省。

**数据要存**
- `T_B` 的连通性 + 逐节点 `h_B`：O(N_bg)
- `{p_i}`、`{h_i}` 的 OT 单元结构（用于求 `f⁻¹`）：O(N_cell)
- `J_f` 的估计（若要用 `J^{1/3}` 修正）：每个背景节点一个值，可用 `f` 在四面体上的
  **体积比** `vol(f(τ))/vol(τ)` 直接算，不需要解析微分

**复杂度**
- 搬背景网格：O(N_bg log N_cell)（逆映射查询）
- mesher：同 Recipe A
- **总计：分钟级**

**优点**：与现有 mesher 工具链零摩擦（每个 mesher 都吃尺寸场）；可实现真正的自适应。
**缺点**：`J^{1/3}` 修正是近似的；背景网格本身的分辨率限制了尺寸场的细节；
fTetWild 的尺寸场量化到 2 的幂，会丢细节。

### 4.3 Recipe C：**体网格 r-自适应（移动节点）**——最贴近 An et al. IMR 2021

**核心思想**：先得到一个**任意**的 M 上四面体网格（甚至就用原始输入网格），
然后在**体里**解一次半离散 OT 把节点搬到与用户密度匹配的位置。**连通性完全不变。**

**为什么这可能比 Recipe A/B 更好**：如果你已经有一个高质量网格（例如 CAD/CAM 来的），
r-自适应是**代价最低**的密度控制手段，且不动拓扑、不碰边界（边界节点可以锁住）。

**算法（3D 版，由 §2.4 的 2D Algorithm 1 推广）**
```
输入：M 的四面体网格 T（n 个节点 {x_i}），用户密度函数 g: M → R+
1  ν_i ← g(x_i)，归一化使 Σν_i = vol(M)
2  h_i ← ½(|x_i|² − 1)   （用 M 的包围盒中心/尺度归一化后再算）
3  repeat
4      构造 {(x_i, −h_i)} 的 LOWER CONVEX HULL        ← 3D 提升到 4D！
5      投影得乘方图 {W_i(h)}（= 加权 Delaunay 的对偶）
6      w_i(h) ← vol(W_i(h) ∩ M)               ← 3D 乘方图与四面体网格求交 + 积分
7      ∇E(h) = (ν_i − w_i(h))
8      if max_i |w_i − ν_i|/ν_i < 1e-8: break
9      组装 Hessian（稀疏 SPD）
10     解 Hess(E) d = ∇E(h)
11     λ←1; while h+λd ∉ H (有空 cell): λ←λ/2
12     h ← h + λd
13 直到收敛
14 新节点位置：x_i ← centroid(W_i(h))      ← 关键：质心就是新位置
15 （可选）局部翻转/平滑修复退化单元
输出：移动后的网格（连通性不变）
```
**3D 化的风险点**：
- 第 4 步是 **4D 凸包**（n 个点 → 4D Delaunay 的对偶）。CGAL 的
  `Regular_triangulation_3`（带权 Delaunay）+ `Power_diagram_2`，或 **geogram** 的
  `GEO::ConvexCell`（Lévy 本人的实现，专为此设计）是最稳的选择。
- 第 6 步（幂图 cell 与四面体网格求交后积分）是**性能瓶颈**。Caplan 原文明确说
  "the cost of performing the numerical integration to compute gradients **significantly outweighed**
  the cost of computing the power diagram"（§2.6）。可行做法：
  **每步只重算"质量变化超过阈值"的 cell**（增量更新），或降低求积精度。
- 边界节点：锁住（`h` 固定）或允许切向滑动，否则网格会脱离 `∂M`。

**数据要存**
- 网格 `(X, T)`：O(n + n_tet)
- `ν`、`h`、`w`：O(n)
- 乘方图 cell 的**多面体表示**（每步重建）：O(n × 平均 cell 面数)，内存是主要开销
- Hessian 的稀疏结构：O(n × 平均邻居数)

**复杂度**
- **每次迭代**：4D 凸包 O(n log n)～O(n²)（最坏）；幂图 cell 求交 + 积分 O(n · k)，`k` = cell 复杂度；
  Newton 解 O(n^1.5)～O(n²)（稀疏 Cholesky）
- **迭代次数**：IMR 2021 实测 **8–12 次**（2D，Franke）
- **实测时间（2D）**：`[原文]`

  | 顶点 | 迭代 | 时间 |
  |---|---|---|
  | 1 681 | 8 | 0.227 s |
  | 10 201 | 10 | 1.544 s |
  | 40 401 | 12 | 7.938 s |
  | 160 801 | 12 | 53.307 s |

  → 3D 预计**贵 1–2 个数量级**（4D 凸包 + 多面体体积）。**10^5 节点 3D 估计分钟～小时级**，
  需要实测。

**优点**：连通性不变（对已有网格友好）；理论最干净（真正解 OT）；
尺寸场是"精确"满足的（`w_i = ν_i` 到 1e-8）。
**缺点**：**算法只在 2D 验证过**（`[原文]` 明确）；3D 实现难度高；边界处理需额外设计；
r-自适应**不能改变单元数量**，只改分布——如果用户要的是"从粗到细"，得先 refine。

### 4.4 方案对比与**明确推荐**

| | Recipe A 采样+重剖 | Recipe B 背景网格尺寸场 | Recipe C 体 r-自适应 |
|---|---|---|---|
| 边界共形 | ✅ 天然 | ✅（`∂M` 仍是输入） | ⚠️ 需锁边界 |
| Sharp feature | ✅ 可保留 | ✅ 可保留 | ⚠️ 可能被移动 |
| 实现难度 | **低** | 中 | **高** |
| 3D 已验证 | ✅（fTetWild/TetGen/Gmsh 都是 3D） | ✅（§3.3 我实测通过） | ❌ 论文只有 2D |
| 单元质量 | **最好**（mesher 保证） | 好 | 取决于初始网格 |
| 精确度量控制 | ❌ 近似 | ❌ 近似 | ✅（1e-8） |
| 能否改单元数 | ✅ | ✅ | ❌ |
| 计算量 | 分钟 | 分钟 | 分钟～小时 |

**推荐（分阶段）**：

> **阶段 1（立刻可做，1–2 周）：Recipe A + fTetWild。**
> 用 §2.4 的 OT 求解器（或你们已有的 OMT 代码）生成球内**按 `r^{1/3}` 分层**的均匀采样点，
> 用 cell 质心求逆映射，把点写到 `.a.node` 或直接作为 fTetWild 的 Steiner 点，
> 边界用原始 `∂M` 面片。**验收指标**：体积误差 < 1%（我实测 TetGen 是 0.6%，Gmsh 是 0.64%）、
> 最小二面角 ≥ 10°、无负体积单元、radius ratio 中位数 ≥ 0.7。

> **阶段 2（1–2 个月）：加 Recipe B 做自适应。**
> 用 Gmsh 的 list-based `.pos` + `PostView` Field（我在 §3.3 已实测通过）或 TetGen 的
> `-m` + `.b.mtr` 实现真正的逐点尺寸场。同时验证 `J_f^{1/3}` 修正公式（§4.2）。

> **阶段 3（研究性）：Recipe C 做 r-自适应。**
> 只在阶段 2 的密度场效果不够、且已有高质量初始网格时才投入。
> 直接基于 **Lévy 的 geogram** + Caplan IMR 2021 的 Algorithm 3 实现，不要从零写 4D 凸包。

> **法务线**：如果最终要闭源商用，**阶段 1/2 用 fTetWild（MPL-2.0）+ Mmg（LGPL-3.0）**，
> 把 TetGen 只当**研究期验证工具**。Gmsh（GPL-2）如果被链接进产品也有风险，
> 只用作离线网格生成器则通常没问题（但请法务确认 GPL 的"工具 vs 库"边界）。

---

## 5. Na Lei / Xianfeng Gu 的「共形几何 + 最优传输 + 网格生成」论文清单

**方法与局限说明**：DBLP 与 gmsh 官方文档站都被 Anubis 反爬挡住（我确认了 HTTP 200 但返回
"Making sure you're not a bot!"），DLUT 教师主页的分页是 JS 驱动的（`javascript:void(false)`），
拿不到第 2–17 页。**我改用 OpenAlex 作者 ID 全量抓取**（`A5086990531` = Na Lei @ DUT，190 篇；
`A5100607707` = Xianfeng David Gu，580 篇），再用正则筛标题，配 Crossref 逐条核对 DOI。
**这是我能做到的最可靠的清单。** 每条都标了来源等级。

### 5.1 核心：测度/密度可控参数化 + 网格生成

| # | 论文 | 出处 | 贡献 |
|---|---|---|---|
| 1 | Su, Chen, **Lei**, Cui, Jiang, **Gu**. "**Measure controllable volumetric mesh parameterization**" | *Computer-Aided Design* **78**:188–198, 2016. DOI [10.1016/j.cad.2016.04.007](https://doi.org/10.1016/j.cad.2016.04.007) | 体积调和映射 → 球 → OT 到**用户指定测度**。精确命中目标测度、减少未知量、收敛更快。**本报告 §1。** `[摘要]` |
| 2 | Su, Chen, **Lei**, Zhang, Qian, **Gu**. "**Volume preserving mesh parameterization based on optimal mass transportation**" | *Computer-Aided Design* **82**:42–56, 2017. DOI [10.1016/j.cad.2016.05.020](https://doi.org/10.1016/j.cad.2016.05.020) | #1 的**保体积特例/前身**：目标测度取均匀。**注意 DOI 是 2016.05.020，不是 2016.07.001**（后者是 Sanfilippo & Borgo 的 feature ontology 综述）`[实测 Crossref]` |
| 3 | Su, Cui, Qian, **Lei**, Zhang, Zhang, **Gu**. "**Area-preserving mesh parameterization for poly-annulus surfaces based on optimal mass transportation**" | *Computer Aided Geometric Design* **46**:76–91, 2016. DOI [10.1016/j.cagd.2016.05.005](https://doi.org/10.1016/j.cagd.2016.05.005) | 曲面版（多环域）保面积参数化 `[摘要]` |
| 4 | Su, Li, Zhao, **Lei**, **Gu**. "**Discrete Lie flow: A measure controllable parameterization method**" | *Computer Aided Geometric Design* **72**:49–68, 2019. DOI [10.1016/j.cagd.2019.05.003](https://doi.org/10.1016/j.cagd.2019.05.003) | **"measure controllable" 方法论的姊妹篇**，Lie flow 视角。**强烈建议找来读，很可能补上 CAD 2016 缺失的算法细节** `[摘要]` |
| 5 | Cao, Wang, Zhu, **Si**, **Lei**, Luo, **Gu**. "**Quad Mesh Density Control using Optimal Transport**" | *Proc. 2025 SIAM International Meshing Roundtable (IMR 2025)*, pp. 45–57. DOI [10.1137/1.9781611978575.5](https://doi.org/10.1137/1.9781611978575.5) | ⭐ **最直接的"OT 密度控制 → 网格"公开文献**。用户给任意概率密度函数；算法 = 半离散 OT（乘方图 + Newton + 线搜索保非空 cell），顶点 = dual face 质心。**我把它的 Algorithm 4.3 逐行抄在 §2.5** `[原文]` |

### 5.2 曲面叶状结构（foliation）→ 四边/六面网格系列（"密度控制"的另一条线）

| # | 论文 | 出处 | 贡献 |
|---|---|---|---|
| 6 | **Lei**, Zheng, Jiang, Lin, **Gu**. "**Quadrilateral and hexahedral mesh generation based on surface foliation theory**" | *CMAME* **316**:758–**781**, 2017. DOI [10.1016/j.cma.2016.09.044](https://doi.org/10.1016/j.cma.2016.09.044) | ✅ **DOI 与任务书一致**；⚠️ **页码是 758–781，不是 758–780**；作者是 Lei, Zheng, Jiang, Lin, Gu `[实测 Crossref]`。用曲面叶状结构生成四边/六面网格 |
| 7 | **Lei**, Zheng, Luo, **Gu**. "**Quadrilateral and hexahedral mesh generation based on surface foliation theory II**" | *CMAME* **321**:406–426, 2017. DOI [10.1016/j.cma.2017.04.012](https://doi.org/10.1016/j.cma.2017.04.012) | 续篇 `[摘要]` |
| 8 | **Lei**, Zheng, Si, Luo, **Gu**. "**Generalized Regular Quadrilateral Mesh Generation based on Surface Foliation**" | *Procedia Engineering* **203**:336–348, 2017. DOI [10.1016/j.proeng.2017.09.818](https://doi.org/10.1016/j.proeng.2017.09.818) | 广义正则四边网格 `[摘要]` |
| 9 | **Chen W.**, Zheng, Ke, **Lei**, Luo, **Gu**. "**Quadrilateral mesh generation I: Metric based method**" | *CMAME* **356**:39–67, 2019. DOI [10.1016/j.cma.2019.07.023](https://doi.org/10.1016/j.cma.2019.07.023) | 度量法四边网格 `[摘要]` |
| 10 | **Lei**, Zheng, Luo, Luo, **Gu**. "**Quadrilateral mesh generation II: Meromorphic quartic differentials and Abel–Jacobi condition**" | *CMAME* **366**:112980, 2020. DOI [10.1016/j.cma.2020.112980](https://doi.org/10.1016/j.cma.2020.112980) | 亚纯四次微分 + Abel–Jacobi 条件 `[摘要]` |
| 11 | Zheng, Zhu, Chen, **Lei**, Luo, **Gu**. "**Quadrilateral mesh generation III: Optimizing singularity configuration based on Abel–Jacobi theory**" | *CMAME* **387**:114146, 2021. DOI [10.1016/j.cma.2021.114146](https://doi.org/10.1016/j.cma.2021.114146) | 奇异点配置优化 `[摘要]` |
| 12 | **Lei**, Zhu, Zheng, **Si**, Luo, **Gu**. "**Why cross fields are not equivalent to quadrilateral meshes**" | *CMAME* **417**:116442, 2023-12. DOI [10.1016/j.cma.2023.116442](https://doi.org/10.1016/j.cma.2023.116442) | 理论：cross field ≠ quad mesh（对做参数化的团队很重要） `[实测 Crossref]` |
| 13 | Zheng, Wang, **Lei**, Luo. "**Quadrilateral mesh generation based on foliation and meromorphic quadratic differential**" | *Computer-Aided Design* **190**:103953, 2026. DOI [10.1016/j.cad.2025.103953](https://doi.org/10.1016/j.cad.2025.103953) | 最新（2026 卷） `[摘要]` |
| 14 | Zheng, Qi, **Lei**, Luo, Si, **Gu**. "**Hexahedral Mesh Generation Based on Meromorphic Quadratic Differentials on Surfaces**" | *Proc. 2026 SIAM IMR*, DOI [10.1137/1.9781611979138.17](https://doi.org/10.1137/1.9781611979138.17) | 六面体 `[摘要]` |
| 15 | Zhu, **Lei**, Zheng, Luo, Si, **Gu**. "**Quadrilateral Mesh Generation for Open Surfaces with Negative Euler Characteristics Based on Symmetric Abel Differentials**" | IMR 2025, pp. 58–71. DOI [10.1137/1.9781611978575.6](https://doi.org/10.1137/1.9781611978575.6) | `[摘要]` |

### 5.3 最优传输 → 球面/生成模型/自适应（与 meshing 相邻）

| # | 论文 | 出处 | 贡献 |
|---|---|---|---|
| 16 | Cui, Qi, Wen, **Lei**, Li, Zhang, **Gu**. "**Spherical optimal transportation**" | *Computer-Aided Design* **115**:181–**193**, 2019. DOI [10.1016/j.cad.2019.05.024](https://doi.org/10.1016/j.cad.2019.05.024) | ✅ **DOI 与任务书给的不同**（任务书写 10.1016/j.cad.2018.04.021，那是 Xin et al. 的测地距离论文）。⚠️ **页码 181–193**（Crossref 口径；Semantic Scholar 也一致）。第一作者是 **Li Cui**，不是 Na Lei。球面 OT 的几何变分算法，可用于 Minkowski 问题、光学设计，**以及球域参数化下的采样** `[实测 Crossref]` |
| 17 | **An, Lei, Zhao, Si, Gu**. "**A Moving Mesh Adaptation Method by Optimal Transport**" | **IMR 2021**, pp. 130–138. [PDF](https://www.meshingroundtable.com/assets/papers/2021/10-An.pdf), DOI [10.5281/zenodo.5559175](https://doi.org/10.5281/zenodo.5559175) | ⭐ **任务书没列但最该看的一篇**：2D r-自适应完整算法 + 收敛表。**§2.4 有逐行转录** `[原文]` |
| 18 | **An, Lei, Chen, Luo, Zhao, Si, Gu**. "**Efficient Approximation of Optimal Transportation Map by Pogorelov Map**" | IMR 2021, [PDF](https://www.meshingroundtable.com/assets/papers/2021/09-An.pdf), DOI [10.5281/zenodo.5559179](https://doi.org/10.5281/zenodo.5559179) | Pogorelov 映射加速 OT `[摘要]` |
| 19 | Wang, Zheng, Chen, Qi, Ren, **Lei**, **Gu**. "**Robust and accurate optimal transportation map by self-adaptive sampling**" | *Frontiers of Information Technology & Electronic Engineering* **22**:1218–1231, 2021. DOI [10.1631/fitee.2000250](https://doi.org/10.1631/fitee.2000250) | 自适应采样提高 OT 精度——**直接对应 §4.0 的 F6 奇异集问题** `[摘要]` |
| 20 | **Lei**, Chen, Luo, **Si**, **Gu**. "**Secondary Polytope and Secondary Power Diagram**" | *Computational Mathematics and Mathematical Physics* **59**(12):2024–2037, 2019. DOI [10.1134/s0965542519120121](https://doi.org/10.1134/s0965542519120121) | **二阶乘方图**——乘方图理论的推广，对高阶网格/更高维 OT 有潜在价值 `[摘要]` |
| 21 | **Lei**, Su, Cui, Yau, **Gu**. "**A geometric view of optimal transportation and generative model**" | *CAGD* **68**:1–21, 2019. DOI [10.1016/j.cagd.2018.10.005](https://doi.org/10.1016/j.cagd.2018.10.005) | 综述性：OT 的几何视角 `[摘要]` |
| 22 | Zhao, **Lei**, Cui, Su, Xu, Luo, **Gu**, Yau. "**A geometric variational framework for computing spherical optimal transportation maps II**" | *Mathematics, Computation and Geometry of Data* **2**(2), 2022. DOI [10.4310/mcgd.2022.v2.n2.a2](https://doi.org/10.4310/mcgd.2022.v2.n2.a2) | 球面 OT 变分框架 II `[摘要]` |
| 23 | **Lei**, Li, Xu, Li, **Gu**. "**What's the Situation With Intelligent Mesh Generation: A Survey and Perspectives**" | *IEEE TVCG* **30**:4997–5017, 2024-08. DOI [10.1109/tvcg.2023.3281781](https://doi.org/10.1109/tvcg.2023.3281781) | 智能网格生成综述（Lei 一作） `[实测 Crossref]` |
| 24 | Owen, Brown, Chrisochoides, Garimella, **Gu**, Ledoux, **Lei**, Quadros, ... "**AI for Mesh Generation and Model Preparation: A Review**" | IMR 2026, DOI [10.1137/1.9781611979138.1](https://doi.org/10.1137/1.9781611979138.1)；期刊版 *Engineering with Computers* 2026, DOI [10.1007/s00366-026-02379-1](https://doi.org/10.1007/s00366-026-02379-1) | 最新综述 `[摘要]` |
| 25 | Zhu, Wang, Zheng, **Lei**, Luo, Chen. "**Prismatic mesh generation based on anisotropic volume harmonic field**" | *Advances in Aerodynamics* **3**:16, 2021. DOI [10.1186/s42774-021-00065-y](https://doi.org/10.1186/s42774-021-00065-y) | **该组少见的体网格生成论文**（棱柱），用各向异性体积调和场 `[摘要]` |
| 26 | Ren, **Lei**, **Si**, **Gu**. "**A Construction of Anisotropic Meshes Based on Quasi-Conformal Mapping**" | *Lecture Notes in Computational Science and Engineering* (2019), DOI [10.1007/978-3-030-13992-6_14](https://doi.org/10.1007/978-3-030-13992-6_14) | 准共形映射构造各向异性网格——**和"共形 + 网格"最直接相关的一篇** `[摘要]` |

### 5.4 任务书里三条引用的**修正建议**

| 任务书写法 | 实际情况 | 处理 |
|---|---|---|
| CMAME 2017, DOI 10.1016/j.cma.2016.09.044 | ✅ DOI 正确；**页码应为 758–781**（不是 758–780）；作者 Lei, Zheng, Jiang, Lin, Gu | 改页码 |
| "Quad mesh density control using optimal transport" (Na Lei et al.) | ✅ 存在，但是 **IMR 2025 会议论文**（Cao, Wang, Zhu, Si, Lei, Luo, Gu），pp. 45–57，DOI [10.1137/1.9781611978575.5](https://doi.org/10.1137/1.9781611978575.5)；**不是 Na Lei 一作** | 补全作者与出处 |
| "Spherical optimal transportation" CAD 2019, DOI 10.1016/j.cad.2018.04.021 | ⚠️ **DOI 错误**。正确：DOI [10.1016/j.cad.2019.05.024](https://doi.org/10.1016/j.cad.2019.05.024)，*CAD* **115**:181–193, 2019；作者 **Cui**, Qi, Wen, Lei, Li, Zhang, Gu。`10.1016/j.cad.2018.04.021` 是 Xin, Wang, He, Zhou, Chen, Tu, Shu 的 "Lightweight preprocessing and fast query of geodesic distance via proximity graph" | **必须改 DOI** |
| "Quadrilateral and hexahedral mesh generation based on surface foliation theory" (CMAME 2017) | ✅ 正确 | 无需改 |

**未找到**：任务书提到的 "Tetrahedral mesh generation via optimal transport" 这一确切标题；
Na Lei / Gu 的著作里**没有**一篇标题同时含 "optimal transport" 与 "tetrahedral mesh" 的论文。
最接近的就是 #1（CAD 2016，体积测度可控参数化）、#17（IMR 2021，OT 移动网格，2D）、
#5（IMR 2025，OT 四边网格密度控制）。
**请不要引用不存在的标题。**

---

## 6. 复现资产（本仓库内）

| 文件 | 内容 | 状态 |
|---|---|---|
| `refs/tetgen-manual-1.6.txt` | TetGen 1.6 官方手册全文抽取（1931 行） | ✅ |
| `refs/caplan2021-power-diagrams.txt` | Caplan IMR 2021 全文抽取 | ✅ |
| `refs/imr2021-moving-mesh-ot.txt` | An–Lei–Zhao–Si–Gu IMR 2021 全文抽取（含 Algorithm 1 与 Table 1） | ✅ |
| `refs/imr2021-pogorelov-ot.txt` | IMR 2021 Pogorelov 论文 | ✅ |
| `refs/imr2025-quad-density-ot.txt` | IMR 2025 Quad Mesh Density Control 全文（含 Algorithm 4.1–4.3） | ✅ |
| `refs/alliez2005-variational-tet.txt` | Alliez 2005 SIGGRAPH 提交版全文（含 Table 1 质量数字） | ✅ |
| `refs/imr2024-optimal-quad.txt`, `refs/imr2024-conformal-welding.txt` | IMR 2024 两篇 | ✅ |
| `refs/imr2021-engwirda-laguerre.txt`, `refs/imr2021-hummel-power.txt` | 乘方图/Laguerre 网格相关 | ✅ |
| `refs/imr2019-marot-optimal-tet.txt`, `refs/imr2019-remacle-gmsh.txt` | 最优四面体化 / Gmsh 鲁棒网格 | ✅ |
| `refs/imr2022-estrems-metric.txt`, `refs/imr2022-caplan-4d.txt`, `refs/imr2024-odeco.txt` | 度量场 r-自适应 / 4D 各向异性 | ✅ |
| `refs/tools-raw-notes.txt` | 各工具 README/LICENSE 原文片段 | ✅ |
| `refs/gmsh_sizing_demo.py` | Gmsh `Box`+`Threshold` 尺寸场实验 | ✅ 跑通 |
| `refs/gmsh_bgmesh_demo.py` | Gmsh Field 栈 **vs** list-based `.pos` 背景网格对比实验 | ✅ 跑通（两条路线都验证） |
| `refs/tetgen_bgmesh_demo.py` | `tetgen` Python 包裹的 `bgmesh=` 路径 | ❌ **崩溃（0xC0000005）**，已记录 |
| `refs/tetgen_mtr_demo.py` | 同上，`bgmeshfilename=` 路径 | ⚠️ TetGen 本体跑通（stdout 有 `Size interpolating seconds`），包装层回传时崩溃 |
| `vendor/tetgen/` | TetGen 1.6.1 源码 + MinGW 构建（`build/tetgen.exe`，3.99 MB） | ✅ 编译成功 |
| `vendor/tetgen/run_sizing_test.py` | 用 TetGen CLI 跑 `-pq1.4a#` / `-m` / `-mR` 的尺寸场对照实验 | ⚠️ 脚本可用，但 `tetgen.exe` 在不同调用方式下给出**不一致的单元数**（13553 vs 43955），**根因未定位**，见 §7 遗留问题 |
| `tools/pdf2txt.py` | PDF→文本工具，已加固（连字/引号清洗 + 控制字符剥离） | ✅ |

---

## 7. 明确标注的未验证项与遗留问题

1. **CAD 2016 正文未获取** → 它的"重采样"步骤、图、表、质量数字**全部未验证**。
   合规获取途径见 §1.1。**建议**：让 Na Lei 组或 Gu 组直接发一份作者版。
2. **CAD 2016 的目标测度 `ν` 如何由用户给定**（标量场？区域常数？）未验证。建议读
   §5.1 的第 4 条（Discrete Lie flow, CAGD 2019）来补。
3. **`h_M(x) = h_B(f(x))·J_f(x)^{1/3}` 这个尺寸场换元公式**是我推导的，未见论文写出。需合成算例验证。
4. **Mmg 的 metric `.sol` 文件格式与 `-ls`/`-hmax` 的确切语义未逐字核实**（官网被反爬）。
   上线前用 `mmg3d -h` 核对。
5. **Gmsh 的完整 Field 列表**（尤其 `MathEval`）未从官方文档逐字核实（文档站被 Anubis 挡）。
   `Box`/`Threshold`/`PostView` 已实测通过。
6. **`tetgen` Python 包裹的 `bgmesh=` 在 Windows 上崩溃**（0xC0000005）——根因未定位，
   可能是 pyvista/tetgen 版本组合问题（pyvista 0.49.0 + tetgen 0.8.4）。**规避：用 CLI 或 meshpy。**
7. **TetGen CLI 在不同调用方式下给出不同单元数**（Python `subprocess` 13553 vs PowerShell 43955，
   同一输入、同一开关 `-pq1.4a0.002`）。我排除了输出文件残留和输入格式问题，
   **根因未定位**（怀疑与 `-a` 后接数字的解析/终端字符编码有关）。
   **建议**：生产脚本一律用 `subprocess.run([...])` 并**显式拆开开关**（`-p`,`-q1.4`,`-a0.002`），
   并对输出做单元数/体积的自动断言。
8. **CGAL `Mesh_3` 的 criteria 模板参数细节**未逐字核实（只拿到文件头 SPDX + 类声明）。
9. **`pycgal` 是否暴露 `Mesh_3`** 未确认（很可能没有）。
10. **IMR 2023/2022 中是否有引 Aurenhammer 的 tet meshing 论文**——我逐一扫了 IMR 2023/2022 的
    论文列表，**没有发现**以"乘方图/power diagram"为标题的四面体网格论文；
    与乘方图相关的是 IMR 2021 的 Engwirda & Liao（Laguerre-Power，地球系统）
    与 Hummel（加权三角形半径公式）。**若需要更彻底的检索**，建议用 Semantic Scholar 的
    citation graph 反查引用了 Aurenhammer 1987（DOI 10.1137/0216006）的 meshing 论文。

---

## 8. 参考文献（本报告引用到的全部条目，含核实状态）

**核心（本报告主体）**
1. Su, K.; Chen, W.; Lei, N.; Cui, L.; Jiang, J.; Gu, X.D. "Measure controllable volumetric mesh parameterization." *Computer-Aided Design* **78** (2016) 188–198. DOI: [10.1016/j.cad.2016.04.007](https://doi.org/10.1016/j.cad.2016.04.007). `[摘要]` 元数据已核实
2. Alliez, P.; Cohen-Steiner, D.; Yvinec, M.; Desbrun, M. "Variational tetrahedral meshing." *ACM Transactions on Graphics* **24**(3) (2005) 617–625. DOI: [10.1145/1073204.1073238](https://doi.org/10.1145/1073204.1073238). OA: [CaltechAUTHORS](https://authors.library.caltech.edu/records/c9exj-3ed14). `[原文]`
3. Du, Q.; Faber, V.; Gunzburger, M. "Centroidal Voronoi Tessellations: Applications and Algorithms." *SIAM Review* **41**(4) (1999) 637–676. DOI: [10.1137/s0036144599352836](https://doi.org/10.1137/s0036144599352836). `[实测 Crossref]`（注意：Alliez 正文里写的是 "Du and Wang [2003]"，实为 Du–Faber–Gunzburger 1999）
4. Du, Q.; Wang, D. "Anisotropic Centroidal Voronoi Tessellations and Their Applications." *SIAM J. Sci. Comput.* **26**(3) (2005) 737–761. DOI: [10.1137/s1064827503428527](https://doi.org/10.1137/s1064827503428527). `[实测 Crossref]`
5. Chen, L.; Xu, J. "Optimal Delaunay Triangulations." *Journal of Computational Mathematics* **22**(2) (2004) 299–308. DOI: [10.4208/jcm.v22.n2.p299](https://doi.org/10.4208/jcm.v22.n2.p299). `[实测 Crossref]`
6. An, D.; Lei, N.; Zhao, T.; Si, H.; Gu, X. "A Moving Mesh Adaptation Method by Optimal Transport." *Proc. IMR 2021*, 130–138. [PDF](https://www.meshingroundtable.com/assets/papers/2021/10-An.pdf), DOI [10.5281/zenodo.5559175](https://doi.org/10.5281/zenodo.5559175). `[原文]`
7. Caplan, P.C. "High-Dimensional Power Diagrams for Semi-Discrete Optimal Transport." *Proc. IMR 2021*, 44–56. [PDF](https://www.meshingroundtable.com/assets/papers/2021/04-Caplan.pdf), DOI [10.5281/zenodo.5559221](https://doi.org/10.5281/zenodo.5559221). `[原文]`
8. Cao, Y.; Wang, Y.; Zhu, Y.; Si, H.; Lei, N.; Luo, Z.; Gu, X. "Quad Mesh Density Control using Optimal Transport." *Proc. 2025 SIAM IMR*, 45–57. DOI: [10.1137/1.9781611978575.5](https://doi.org/10.1137/1.9781611978575.5). [PDF](https://www.meshingroundtable.com/assets/papers/2025/1015-compressed.pdf). `[原文]`
9. Aurenhammer, F. "Power Diagrams: Properties, Algorithms and Applications." *SIAM Journal on Computing* **16**(1) (1987) 78–96. DOI: [10.1137/0216006](https://doi.org/10.1137/0216006). `[实测 Crossref]`（⚠️ Caplan 文中引作 [20]，DOI 10.1137/0215043 是错的，那是 1986 年矩阵乘法论文）
10. Lévy, B. "A Numerical Algorithm for L² Semi-Discrete Optimal Transport in 3D." *ESAIM: M2AN* **49**(6) (2015) 1693–1715. DOI: [10.1051/m2an/2015055](https://doi.org/10.1051/m2an/2015055). `[实测 Crossref]`
11. Mérigot, Q. "A Multiscale Approach to Optimal Transport." *Computer Graphics Forum* **30**(5) (2011) 1583–1592. DOI: [10.1111/j.1467-8659.2011.02032.x](https://doi.org/10.1111/j.1467-8659.2011.02032.x). `[实测 Crossref]`
12. Mérigot, Q.; Meyron, J.; Thibert, B. "An Algorithm for Optimal Transport between a Simplex Soup and a Point Cloud." *SIAM J. Imaging Sci.* **11**(2) (2018) 1363–1389. DOI: [10.1137/17m1137486](https://doi.org/10.1137/17m1137486). `[实测 Crossref]`
13. Kitagawa, J.; Mérigot, Q.; Thibert, B. "Convergence of a Newton algorithm for semi-discrete optimal transport." *J. Eur. Math. Soc.* **21**(9) (2019) 2603–2651. DOI: [10.4171/jems/889](https://doi.org/10.4171/jems/889). `[实测 Crossref]`

**Alliez 相关**
14. Alliez, P.; Cohen-Steiner, D.; Devillers, O.; Lévy, B.; Desbrun, M. "Anisotropic polygonal remeshing." *ACM TOG* **22**(3) (2003) 485–493. DOI: [10.1145/882262.882296](https://doi.org/10.1145/882262.882296). `[实测 Crossref]`
15. Cohen-Steiner, D.; Alliez, P.; Desbrun, M. "Variational shape approximation." *ACM TOG* **23**(3) (2004) 905–914. DOI: [10.1145/1015706.1015817](https://doi.org/10.1145/1015706.1015817). `[实测 Crossref]`

**该组（Lei / Gu）其他**
16. Su, K.; Chen, W.; Lei, N.; Zhang, J.; Qian, K.; Gu, X. "Volume preserving mesh parameterization based on optimal mass transportation." *Computer-Aided Design* **82** (2017) 42–56. DOI: [10.1016/j.cad.2016.05.020](https://doi.org/10.1016/j.cad.2016.05.020). `[实测 Crossref]`
17. Cui, L.; Qi, X.; Wen, C.; Lei, N.; Li, X.; Zhang, M.; Gu, X. "Spherical optimal transportation." *Computer-Aided Design* **115** (2019) 181–193. DOI: [10.1016/j.cad.2019.05.024](https://doi.org/10.1016/j.cad.2019.05.024). `[实测 Crossref]`
18. Lei, N.; Zheng, X.; Jiang, J.; Lin, Y.-Y.; Gu, D.X. "Quadrilateral and hexahedral mesh generation based on surface foliation theory." *CMAME* **316** (2017) 758–781. DOI: [10.1016/j.cma.2016.09.044](https://doi.org/10.1016/j.cma.2016.09.044). `[实测 Crossref]`
19. Lei, N.; Zheng, X.; Luo, Z.; Gu, D.X. "...theory II." *CMAME* **321** (2017) 406–426. DOI: [10.1016/j.cma.2017.04.012](https://doi.org/10.1016/j.cma.2017.04.012). `[实测 Crossref]`
20. Su, K.; Li, C.; Zhao, S.; Lei, N.; Gu, X. "Discrete Lie flow: A measure controllable parameterization method." *CAGD* **72** (2019) 49–68. DOI: [10.1016/j.cagd.2019.05.003](https://doi.org/10.1016/j.cagd.2019.05.003). `[摘要]`
21. Su, K.; Cui, L.; Qian, K.; Lei, N.; Zhang, J.; Zhang, M.; Gu, X. "Area-preserving mesh parameterization for poly-annulus surfaces based on optimal mass transportation." *CAGD* **46** (2016) 76–91. DOI: [10.1016/j.cagd.2016.05.005](https://doi.org/10.1016/j.cagd.2016.05.005). `[摘要]`
22. Ren, Y.; Lei, N.; Si, H.; Gu, X. "A Construction of Anisotropic Meshes Based on Quasi-Conformal Mapping." *LNCSE* (2019). DOI: [10.1007/978-3-030-13992-6_14](https://doi.org/10.1007/978-3-030-13992-6_14). `[摘要]`
23. Wang, Y.; Zheng, X.; Chen, W.; Qi, X.; Ren, Y.; Lei, N.; Gu, X. "Robust and accurate optimal transportation map by self-adaptive sampling." *FITEE* **22** (2021) 1218–1231. DOI: [10.1631/fitee.2000250](https://doi.org/10.1631/fitee.2000250). `[摘要]`
24. Lei, N.; Li, Z.; Xu, Z.; Li, Y.; Gu, X.D. "What's the Situation With Intelligent Mesh Generation: A Survey and Perspectives." *IEEE TVCG* **30**(8) (2024). DOI: [10.1109/tvcg.2023.3281781](https://doi.org/10.1109/tvcg.2023.3281781). `[摘要]`
25. Zhu, Y.; Wang, S.; Zheng, X.; Lei, N.; Luo, Z.; Chen, B. "Prismatic mesh generation based on anisotropic volume harmonic field." *Advances in Aerodynamics* **3** (2021) 16. DOI: [10.1186/s42774-021-00065-y](https://doi.org/10.1186/s42774-021-00065-y). `[摘要]`

**工具论文**
26. Si, H. "TetGen, a Delaunay-Based Quality Tetrahedral Mesh Generator." *ACM TOMS* **41**(2) (2015) 1–36. DOI: [10.1145/2629697](https://doi.org/10.1145/2629697). `[实测 Crossref]`
27. Hu, Y.; Schneider, T.; Wang, B.; Zorin, D.; Panozzo, D. "Fast tetrahedral meshing in the wild." *ACM TOG* **39**(4) (2020), Article 117. DOI: [10.1145/3386569.3392385](https://doi.org/10.1145/3386569.3392385). `[实测 Crossref]`
28. Hu, Y.; Zhou, Q.; Gao, X.; Jacobson, A.; Zorin, D.; Panozzo, D. "Tetrahedral meshing in the wild." *ACM TOG* **37**(4) (2018) 1–14. DOI: [10.1145/3197517.3201353](https://doi.org/10.1145/3197517.3201353). `[实测 Crossref]`
29. Martin, T.; Cohen, E.; Kirby, M. "Volumetric parameterization and trivariate B-spline fitting using harmonic functions." *SPM '08*. `[摘要]`
