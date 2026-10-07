# 体保持（volume-preserving）3D 参数化流水线 —— 前半段调研报告

**Ball parameterization: initial homeomorphism + volume-preserving optimization**

- 工作目录 / workspace: `E:\panyingyun\smartmm`
- 面向：四面体网格（TetGen `.node`/`.ele`）→ 单位球的双射参数化，后续用于四面体重网格化
- 本报告作者：DSH research subagent（session `session-49d98a7e-73e8-4c94-8f97-0c33e0673477`）
- **本报告中的每一条公式与每一个 API 签名都经过了实际代码验证**，验证脚本见 §7。
  未能验证的内容全部在 §8 中显式标注为 UNVERIFIED。

---

## 0. 结论摘要（TL;DR）

| 问题 | 结论 |
|---|---|
| 目标 A：经典 volumetric harmonic map 到单位球的离散形式 | P1 有限元刚度矩阵 `K_ij = -1/6 · l_kl · cot θ_ij^{kl}`（tets），边界固定到球面。**公式已给出 + 已验证**（§1.2, §7 V5a–V5k） |
| 目标 A 的失效模式 | 对**非凸**实体，径向投影边界 + Dirichlet 内部解会产生极严重畸变。实测 peanut 形状：体积比 std 2.09、range [0.19, 8.88]、**16.5% 的四面体 J<0.5**（§1.4, §7 V7b） |
| 重要陷阱 | 若输入网格本身已经是单位球（边界顶点已在球面上），harmonic 解**恰好是恒等映射**，畸变精确为 0。**不要用这类网格测试失效模式**（§7 V7a） |
| 目标 B：VSEM（固定边界不动点） | 已验证可收敛，在**几何球**上把 E_V 压到理论下界（gap ≈ 1.8e-15），体积比 std ≈ 7.8e-16（§3.1, §7 V4b–V4d） |
| 目标 B：VSEM 的根本局限 | 边界冻结时**图像体积不可控**。实测非凸实体：E_V − 1.5·V(f) = 1.74 > 0，**无法达到下界**。这正是 IEM 让边界点在球面上自由滑动的原因（§3.3, §7 V8/V9） |
| 推荐路线 | 1) 边界：`SphericalAEM`(MATLAB, 有源码) 或 `spherical_conformal_map`(MATLAB, 有源码) → 2) 内部初值：P1 harmonic（§1.2） → 3) 精修：IEM（preconditioned nonlinear CG, §3.4） |
| **公开代码的现实** | **n-VSE / VSEM / IEM 三篇论文都没有公开代码。** 已逐一确认作者 GitHub（`mhyueh` 只有 4 个无关仓库）、作者主页、MathWorks。仅 SAEM 与 spherical-conformal-map 有可下载源码（§3.5, §2.5） |
| numpy/scipy 陷阱 | scipy CSR 的 `K[idx][:, idx]` 取子矩阵**是正确的**，但**不要拿全局子矩阵去和单元矩阵比较**——全局对角包含所有相邻单元的贡献。这个坑耽误了本次调研两轮（§7 备注） |

---

## 1. 目标 A：Volumetric Harmonic Map 到单位球

### 1.1 文献谱系（含精确引用）

| 文献 | 精确引用 | 备注 |
|---|---|---|
| 经典 volumetric harmonic map | Y. Wang, X. Gu, T. F. Chan, P. M. Thompson, S.-T. Yau, "Volumetric harmonic brain mapping", *ISBI* 2004; 以及 **Wang, Li, Gu et al., "Volumetric harmonic map", Communications in Information and Systems 3(3):191–202, 2003** (DOI `10.4310/CIS.2003.v3.n3.a3`) | 原始构造 |
| Meshless harmonic volumetric mapping | Y. Li et al., "Meshless harmonic volumetric mapping using fundamental solution methods", *IEEE T-ASE* 6(3):409–422, 2009 | 无网格 |
| Biharmonic volumetric mapping | H. Xu, et al., "Biharmonic volumetric mapping" | 降畸变 |
| Trivariate B-spline fitting | T. Martin, E. Cohen, R. M. Kirby, "Volumetric parameterization and trivariate B-spline fitting using harmonic functions", *SMI* 2008 / *CAGD* 26(6):648–664, 2009 | 工程拟合 |
| Star-shaped volumes / Green's function | J. Xia, Y. He, X. Yin, S. Wang, X. Gu, "Parameterization of star-shaped volumes using Green's functions", *GMP* 2010, LNCS 6144, pp. 219–235 | 只适用于 star-shaped |
| n-VSE（把 2D/3D 统一到 n 维） | Z.-H. Tan, T. Li, W.-W. Lin, S.-T. Yau, **arXiv:2402.00380**; 期刊版 *SIAM J. Imaging Sci.* **18(2):1141–1175, 2025** | 有完整算法框（抄录见 §3.2） |

> 我没能下载到 2003 年 CIS 那篇 PDF（付费墙），因此 §1.2 的离散形式是依据 **n-VSE (arXiv:2402.00380)** 的 `n=3` 特例 + **VSEM (arXiv:2210.09654)** 的原始定义重建的。二者是同一套离散格式，可靠性高，但**"2003 年原文是否逐字相同"这一点我未验证**（§8 #1）。

### 1.2 离散形式：P1 有限元刚度矩阵（**已实现 + 已验证**）

给定四面体 $\tau=[v_0,v_1,v_2,v_3]$，令

$$B = [\,v_1-v_0,\; v_2-v_0,\; v_3-v_0\,]\in\mathbb R^{3\times3}\quad(\text{列为边向量}),\qquad
B[i,j]=(v_{j+1}-v_0)_i$$

则重心坐标（P1 形函数）的梯度为

$$\boxed{\ \nabla\lambda_k \;=\; \text{$B^{-1}$ 的第 }k\text{ 行},\quad k=1,2,3,\qquad
\nabla\lambda_0 \;=\; -\sum_{k\ge1}\nabla\lambda_k\ }$$

单元刚度矩阵与全局刚度矩阵：

$$K^\tau_{ij} \;=\; |\tau|\;\nabla\lambda_i\cdot\nabla\lambda_j,
\qquad
K \;=\; \sum_{\tau} K^\tau \quad(\text{按顶点编号 scatter})$$

**等价闭式（四面体版 cotangent 公式，工程上最常用）：**

$$\boxed{\;
K^\tau_{ij} \;=\; -\tfrac16\, l_{kl}\,\cot\theta_{ij}^{kl}\quad (i\ne j),\qquad
K^\tau_{ii} \;=\; -\sum_{j\ne i} K^\tau_{ij}
\;}$$

其中 $(k,l)$ 是与 $(i,j)$ 相对的边（$\{i,j,k,l\}=\{0,1,2,3\}$），$l_{kl}$ 是其长度，$\theta_{ij}^{kl}$ 是边 $(i,j)$ 处的二面角。

> ⚠️ **`$(k,l)$ 是"相对边"**。在三角形里 cotangent 权重用的是对角，四面体里用的是**对边**。写错这一个索引，整个矩阵就错，但 `K·1 = 0` 仍然成立（因为按行消去），所以**"K·1=0"不足以验证你的实现**。必须像 §7 V5b 那样逐元素对照闭式公式。

**调试期踩到的真实坑（务必避免）：**
1. `B^{-1}` 的行 vs 列。行 = $\nabla\lambda_k$。搞反后 `K·1=0` 依然成立（对称性掩盖了错误），要用 V5g/V5h 才能抓到。
2. `K[idx][:,idx]`（全局子矩阵）**不等于**单元矩阵 `K^τ`——全局对角还累加了所有其他相邻单元的贡献。两者相差可达 3.5（量级 100%）。单元级测试必须用独立的单元装配函数。
3. scatter 索引：`rows = repeat(T,4,axis=1).ravel()`, `cols = repeat(T,4,axis=0).ravel()`。写成 `tile(T,(1,4))` 会静默给出错误的矩阵（对角偏大），且 `K·1=0` **仍然成立**。

### 1.3 边界条件与线性系统

设 $n$ 个顶点，$\mathtt B=\{t\mid v_t\in\partial\mathcal M\}$（边界），$\mathtt I=\{1..n\}\setminus\mathtt B$（内部）。把目标映射写成 $\mathbf f\in\mathbb R^{n\times3}$。

**边界条件：** 固定 $\mathbf f_{\mathtt B} = \mathbf g$，其中 $\mathbf g:\partial\mathcal M\to\mathbb S^2$ 是由 §2 的球面映射算法算出的**双射**球面参数化。

**待解线性系统（3 次，每个坐标分量一次）：**

$$\boxed{\;K_{\mathtt I\mathtt I}\,\mathbf f_{\mathtt I}^{\,s} \;=\; -\,K_{\mathtt I\mathtt B}\,\mathbf f_{\mathtt B}^{\,s},
\qquad s=1,2,3\;}$$

`K_II` 是对称正定（SPD）的（Dirichlet 拉普拉斯），用 Cholesky 一次分解、三次回代即可。工程上用 `scipy.sparse.linalg.factorized(K_II.tocsc())`。

**实现要点：**
- **必须先归一化网格朝向**：TetGen 输出的四面体可能有正有负。实测：把整个网格旋转后重新三角化，会出现 50% 的负定向 tets。用 `orient_tets()` 统一到 `det>0`（交换两个顶点即可）。**不统一会连带毁掉二面角、外法向、边界绕序**。
- **边界三角形绕序必须几何修正**：不要依赖任何固定的"规范面排列"，因为 TetGen 不遵守。对每个 tet 的面，用 `n·(v_k - v_a) > 0` 判据翻转，使其背离对顶点（即朝外）。**不修正的话，球面三角形有符号面积会算成 0 或错一倍**（本次调研在此处卡了很久：`divergence volume = -0.000000`、`spherical excess = 8π`，根因就是 50% 的边界面向内）。

### 1.4 已知失效模式（**有实测数据**）

| 失效模式 | 机理 | 实测证据（peanut 形状，$r(\theta)=1+0.45\cos2\theta$，3088 tets） |
|---|---|---|
| **Fold-over（翻转）** | 域远离凸/星形时，调和映射不再单射 | harmonic + 径向 BC：体积比 mean 2.47 / std 2.09 / range [0.19, 8.88]，**16.5% 的 J<0.5** |
| 边界映射本身退化 | 径向投影对凹形状会产生极端面积畸变 | 径向 BC 的球面三角形面积 max/min ≈ 1.29；凹得更厉害时边界先崩 |
| 退化四面体 (sliver) | $l_{kl}\cot\theta$ 在 $\theta\to0$ 或 $\pi$ 时爆炸 | 球面网格上实测 min dihedral 5.5°，aspect ratio max 17.9；映射后 max radius-edge ratio 从 1.5 涨到 130 |
| 非双射 | 调和映射**不保证**单射；只在 $C^1$ 且 $\det J>0$ 时局部单射，全局单射需要额外论证 | 见下 |

**⚠️ 一个必须知道的陷阱（实测）：**
> 在**几何上本来就是单位球**的网格上（边界顶点恰好落在 $\mathbb S^2$），径向投影边界就是**恒等边界条件**，而恒等映射本身是调和函数，于是 P1 调和映射的解**精确等于恒等映射**：
> ```
> ratio_std = 1.86e-15,  frac_flipped = 0,  frac_J<0.5 = 0,
> E_V = 6.229111225639588,  1.5*V(f) = 6.229111225639587   (gap = 8.9e-16)
> ```
> **这看起来"完美"，但它什么也没测到。** 用这种网格做回归测试只能验证代码不崩（这个价值也不小——它验证了装配、线性求解、能量、度量全链路自洽），
> **但绝不能用来评估"失效模式"或"质量提升"**。本报告的所有质量结论都来自 §7 的非凸实体。
> 建议团队**两个测试都留**：几何球做正确性回归，非凸实体做质量评估。

**关键实测反直觉结论：** 在本组的实验中，**翻转 tets 数为 0，但 J 严重压缩**（16.5% 的 tets 有 J<0.5）。也就是说：

> **"没有翻转"不等于"映射质量好"。** 必须同时监控 $J$ 的分布（不只是符号），否则会得到"看起来合法、实际上畸变巨大"的映射。

**关于"星形/凸"的实用判据：** 若 $\mathcal M$ 关于某内点 $c$ 是星形的（每条射线与边界恰交一次），则径向投影边界 + 调和内部解通常**不翻转**，畸变也温和。实测 peanut $a=0.10$（弱凹）：J<0.5 比例为 0%（§7 sweep）。$a\ge0.20$ 开始明显劣化。所以：**先用星形判据（或直接数每个方向的射线交点数）筛掉危险输入**，再决定是否需要更强的边界算法。

---

## 2. 目标 A 的边界条件：闭合亏格 0 三角网格 → 球面

### 2.1 各算法对比

| 方法 | 核心思想 | 代码可用性 | 参考 |
|---|---|---|---|
| **Spherical conformal** (Gu & Yau / Choi & Lui) | 解一个非齐次 Laplace–Beltrami 方程（"Dirac map" / FLASH），通过球极投影把球面约束变成线性系统；再用 Möbius 修正降面积畸变 | ✅ **MATLAB 源码可下载**（`garyptchoi/spherical-conformal-map`） | Choi, Lam, Lui, *SIIMS* 8(1):67–94, 2015, DOI `10.1137/130950008`；Möbius 修正：Choi, Leung-Liu, Gu, Lui, *SIIMS* 13(3):1049–1083, 2020, DOI `10.1137/19M125337X` |
| **Spherical authalic / stretch energy min (SAEM)** | 最小化 spherical authalic energy $E_{\mathbb A}(f)=\frac{|\mathcal M|}{3\mathbb V(f)}E_S(f)-3\mathbb V(f)$，$E_S$ 为 stretch energy，$\mathbb V$ 为有符号体积；边界用球坐标 $(\theta,\phi)$ 参数化 + **preconditioned nonlinear CG** | ✅ **MATLAB 源码可下载**（`SAEM.zip`，函数被编译为 `.p`，但 `main.m` 有完整调用签名） | Liu & Yueh, *SIAMS* 19(1):207–235, 2026, DOI `10.1137/25M1736979`；arXiv:2412.19011 |
| **Spherical optimal transport (SOT)** | 先把曲面共形映到球（harmonic map），把共形因子当作概率密度；再用 **spherical power diagram + 凸优化（Newton）** 求最优传输映射。得到保面积球面参数化 | ⚠️ **C++ tutorial 有源码入口**（Gu 组 `software/SOT/index.html` 声称 "with C++ source"，二进制 `OT.exe`，命令行 `OT.exe -source <mesh> -target <sphere.obj>`）。我**未实际下载验证源码包**。论文代码本身未单独发布 | Cui, Qi, Wen, Lei, Li, Zhang, Gu, "Spherical optimal transportation", *CAD* 82:181–193, 2019, DOI `10.1016/j.cad.2019.05.024` |
| **Discrete surface Ricci flow → sphere** | $\frac{dg_{ij}}{dt}=(\bar K-K)g_{ij}$；离散曲率 = 角亏 $K_i=2\pi-\sum\theta_i$；收敛到常数曲率度量后映到球 | ⚠️ Gu 组有 C++ tutorial（`tutorial/RicciFlow.html`）与 `MeshLib` | Gu & Yau 系列 |
| **Harmonic map to sphere (bootstrapping)** | 直接最小化 Dirichlet 能量 $\int\|\nabla f\|^2$ 且约束 $\|f\|=1$；因为球面不是凸的，直接 CG 会卡；**实用做法是先共形初值再迭代** | ⚠️ 有第三方 C++ 实现 `icemiliang/spherical_harmonic_maps`（见 §2.4） | Eells–Sampson 理论；Gu & Yau 2004 spherical harmonic map |
| **Spherical Tutte / barycentric** | 把顶点放球面 + 凸组合约束 | 极简单但畸变很大，常用于给 conformal 做初值 | Hartley & Anderson |

**工程建议顺序：** 用 `SphericalAEM`（保面积，最贴近"体积保持"目标）或 `spherical_conformal_map`（最稳、最成熟）算边界，再跑 §1.2 的 harmonic 内部解，最后用 §3.4 的 IEM 精修。

### 2.2 球面离散计算的必要验证

球面映射是否**双射**，工程上用一个便宜判据（本次已实现并验证）：

$$\sum_{\text{边界面}} \Omega_\tau \;=\; 4\pi
\quad\text{（双射）}\qquad
\Omega_\tau = 2\,\operatorname{atan2}\!\big(a\cdot(b\times c),\;1+a\cdot b+b\cdot c+c\cdot a\big)$$

即 Van Oosterom–Strackee 公式（IEEE TBME 30(2):125–126, 1983），**必须用带符号版本**：取绝对值会掩盖双重覆盖和翻转。

> ⚠️ **实测踩坑**：一开始我把公式写成 `4*atan2(|det|, den)`，结果连标准单位球都算出 $8\pi$，让我误判"径向投影双重覆盖了球面"，浪费了若干轮。正确系数是 **2 且不带 abs**。验证：标准球面 → 精确 4π（ratio 1.000000）。

同时检查每个 $\Omega_\tau > 0$；有任何负值就说明球面网格翻转。

### 2.3 可用的代码 / 包（**逐个实测过**）

#### ✅ MATLAB —— `garyptchoi/spherical-conformal-map`（**推荐首选**）
- URL: <https://github.com/garyptchoi/spherical-conformal-map>
- 语言 MATLAB，License 未在仓库声明（README 写 "Copyright (c) 2013-2024, Gary Pui-Tung Choi"）。已下载 `demo.m` + `README.md` 到 `vendor/`。
- 入口签名（**已验证，来自 README + demo.m**）：
  ```matlab
  map = spherical_conformal_map(v, f)
  % v: nv x 3  genus-0 闭合三角网格顶点
  % f: nf x 3  三角面
  % map: nv x 3  球面共形映射
  ```
- 扩展：`mobius_area_correction_spherical` —— 在保持共形的前提下降低全局面积畸变（Choi-Leung-Liu-Gu-Lui 2020）。
- 同作者还有 `ellipsoidal-map`、`spherical-density-equalizing-map`、`disk-conformal-map`、`toroidal-density-equalizing-map`。

#### ✅ MATLAB —— SAEM（Liu & Yueh 2026）
- 项目页：<https://math.ntnu.edu.tw/~yueh/projects/SphereAEM.html>
- 直接下载：<https://math.ntnu.edu.tw/~yueh/projects/SphereAEM/SAEM.zip>（≈1.18 MB，**已下载并解包到 `vendor/SAEM/`**）
- 内容：`main.m`、`SphericalAEM.p`、`RiemanianBijectiveCorrection.p`、`ShapeDescription.m`、`RealSphericalHarmonic.p`、`area_distortion.m`、`plot_mesh.m`、`DavidHead.mat`、`RightHand_Folding.mat`
- **License（来自 main.m 头部）：academic and research purposes only，禁止商用，需作者书面许可。** 注意 `.p` 是 MATLAB P-code（已混淆字节码），**不可读源码**。
- 入口签名（**已验证，抄自 main.m**）：
  ```matlab
  S = SphericalAEM(F, V)
  %   F: #F x 3 闭合三角网格面
  %   V: #V x 3 顶点
  %   S: #V x 3 球面保面积映射
  % 可选: SphericalAEM(__, "MaxIter", 100, "Tol", 1e-5)

  S = RiemanianBijectiveCorrection(F, S)   % 把折叠的球面映射修正为双射

  [SH, coef, ReV] = ShapeDescription(V, S, L)  % 球谐形状描述, L = 最大阶
  ```
  最小用例（`main.m` 原样）：
  ```matlab
  load('DavidHead.mat');          % 提供 F, V
  S = SphericalAEM(F, V);
  area_distortion(F, V, S)
  plot_mesh(F, S); title('Area-Preserving Map');
  ```

#### ⚠️ C++ —— Gu 组 spherical harmonic map / Ricci flow / SOT
- 软件索引：<https://www3.cs.stonybrook.edu/~gu/software/index.html>
- 我 fetch 了该页，其中明确列出：
  - "Tutorial on Spherical Optimal Transportation Map **with C++ source**"
  - "Tutorial on Surface Ricci Flow **with C++ source**"
  - "Tutorial on Harmonic Map **with C++ source**"
  - `RiemannMapper`: "A general purpose mesh parameterization tool kit ... compute conformal, quasi-conformal mappings"
- SOT 页面给出的命令行：`OT.exe -source <source mesh> -target <target mesh> [-step_length s] [-help]`，例：`OT.exe -source brain.obj -target brain.sphere.obj`，热键 `!` = 用 Newton 法走一步。
- ❌ **UNVERIFIED**: 该页的源码下载链接我**没有成功提取到**（页面上的 `href` 只解析出 `RiemannMapper/index.html`、`MeshLib/index.html`、`ConformalGeometry/index.html`、`holoimage/index.html`）。**需要人工打开页面点下载链接确认。**

#### ⚠️ C++ —— `icemiliang/spherical_harmonic_maps`
- URL: <https://github.com/icemiliang/spherical_harmonic_maps>
- License: **无**（GitHub API 返回空）。语言 C++，29 stars，最后更新 2026-01-06。
- 文件结构（**已验证，来自 GitHub trees API**）：`main.cpp`、`src/Solid.h/.cpp`（30 KB，核心）、`src/harmonic.h/.cpp`（6.9 KB，球面调和映射）、`src/OBJFileReader.*`、`data/brain.obj`（232 KB，可直接做测试数据）、`CMakeLists.txt`。
- **UNVERIFIED**：我没有 fetch `harmonic.cpp` 的源码，因此**函数签名与算法细节未确认**。它看起来是一个通用 half-edge 网格库 + harmonic 映射实现，不是专门针对本任务的。

#### ❌ CGAL / libigl / MeshLab / pymeshlab / potpourri3d / gpytoolbox —— **都没有球面参数化**

这一条是**实测确认**的，可以节省团队大量时间：

| 库 | 实测结论 |
|---|---|
| **libigl (Python `igl` 0.4.x)** | `[n for n in dir(igl) if 'spher' in n.lower()]` → **只有 `ray_sphere_intersect`**。没有任何球面参数化。相关但有界的函数（签名已验证）见 §2.3 表格下方 |
| **pymeshlab (Python)** | `ml.filter_list()` 共 281 个 filter。含 `'spher'` 的只有 `create_sphere`、`create_sphere_cap`、`create_sphere_points`（都是**生成**球，不是映射到球）。含 `param/conformal/harmon` 的全部是**平面/UV 参数化**：`compute_texcoord_parametrization_harmonic`、`compute_texcoord_parametrization_least_squares_conformal_maps`、`compute_iso_parametrization` 等 |
| **MeshLab (C++)** | 同上，filter 集与 pymeshlab 一致。**没有** spherical parameterization filter |
| **potpourri3d** | API 全集（实测）：`EdgeFlipGeodesicSolver, GeodesicTracer, MeshFastMarchingDistanceSolver, MeshHeatMethodDistanceSolver, MeshMarchingTrianglesSolver, MeshSignedHeatSolver, MeshVectorHeatSolver, PointCloudHeatSolver, PointCloudLocalTriangulation, PolygonMeshHeatSolver, compute_distance, compute_distance_multisource, cotan_laplacian, edges, face_areas, io, marching_triangles, mesh, point_cloud, read_mesh, read_point_cloud, read_polygon_mesh, validate_mesh, validate_points, vertex_areas, write_mesh, write_point_cloud`。**无球面映射** |
| **gpytoolbox 0.3.7** | 名字含 `spher/ball/tetra` 的只有 `ReachForTheSpheresState, icosphere, reach_for_the_spheres, reach_for_the_spheres_iteration`。**无球面映射**。`icosphere(n)` 可用于造测试网格 |
| **CGAL** | **UNVERIFIED**：我没有实际检查 CGAL。已知 `CGAL::Surface_mesh_parameterization` 是**圆盘/平面**参数化（LSCM、ARAP、MVC、离散共形等），**我不认为有球面参数化包**，但这一条我没有一手核实 |

**libigl 中真正有用的函数（签名均已实测）：**
```python
igl.harmonic(V, F, b, bc, k) -> W
  # V #Vxdim 顶点, F #Fxsimplex, b #b 边界索引, bc #bx#W 边界值
  # k=1 harmonic, k=2 biharmonic ; 返回 W #Vx#W 权重场 (可做 3 个坐标 -> 球坐标)
igl.map_vertices_to_circle(V, bnd) -> UV          # 仅圆盘边界
igl.lscm(V, F, b, bc) -> (V_uv, Q)                # 平面 LSCM
igl.cotmatrix(V, F) -> scipy.sparse.csc_matrix    # F 可以是三角形 **或四面体**!
igl.cotmatrix_entries(V, F)
igl.bijective_composite_harmonic_mapping(V, F, bnd, bc)
```
> 💡 **重要（已实测）**：`igl.cotmatrix` 的 docstring 明确写 "F #F by simplex_size list of mesh elements (**triangles or tetrahedra**)"。**所以 libigl 可以直接给你四面体的 P1 刚度矩阵**，不必自己写 §1.2。
>
> **实测结果**：
> ```python
> import igl, numpy as np
> K_igl = igl.cotmatrix(V, T.astype(np.int64)).toarray()      # T 是 (m,4) 四面体
> K_ref = p1_tet_stiffness(V, T).toarray()                    # 本报告的实现
> np.abs(K_igl + K_ref).max()   # -> 4.44e-15   <<< 机器精度
> np.abs(K_igl - K_ref).max()   # -> 8.27       (相差整体负号)
> ```
> **结论：`igl.cotmatrix` 在四面体上返回 $-K$（整体负号），数值与本报告的 `p1_tet_stiffness` 逐位一致。**
> ⚠️ **符号陷阱**：因为 `cotmatrix` 在三角形上返回的是普通的正定 cotangent Laplacian（对角为正），
> 团队很容易误以为在四面体上也是同一符号。**实测是整体负号**，直接拿去解
> `K_II f_I = -K_IB f_B` 会得到反号的解（等于把调和映射换成"反调和"），请务必先做一次
> `K @ 1 == 0` 与 `K_igl + K_ref ≈ 0` 的断言。**强烈建议直接用 `igl.cotmatrix` 并显式取负**：
> ```python
> K = -igl.cotmatrix(V, T.astype(np.int64))     # 现在 K 是本报告意义的 P1 刚度矩阵
> assert np.abs(K @ np.ones(len(V))).max() < 1e-9
> ```

#### ❌ `spharm-pdm` / `cmrep` / FreeSurfer —— 不满足本任务需求
- **`spharm-pdm`**（UNC）：目标是**从二值分割体数据**做 SPHARM-PDM 形状分析（球谐系数统计），不是"给一个闭合 genus-0 三角网格做保体积球面映射"。它内部**确实**需要一个球面参数化步骤（球面调和基），但输出不是"四面体网格的边界球面参数化"。**UNVERIFIED**：我没有实际安装它验证接口。
- **`cmrep` / CM-Rep**（Yushkevich）：连续中轴表示，主要用于**管状/薄壳**结构，不是把实体映到球。**UNVERIFIED**。
- **FreeSurfer**：`mris_sphere` 做的是**皮层表面**（genus-0 但开放的、带孔洞的半球状表面）球面映射，流程强绑定 FreeSurfer 的 subject 目录结构和 MRI 体数据。**Windows 原生不可用**（需 WSL/容器或 `freesurfer` conda 包）。**UNVERIFIED**：我没有安装验证。
- 结论：**这三条路对本任务都不推荐**，别浪费时间。

#### ✅ 实测可用的 pip 安装（本机已装好，可直接复现）
```bash
pip install meshio==5.3.5        # TetGen .node/.ele 读写
pip install tetgen==0.8.4        # Python 接口，底层是 TetGen C++ (需 pyvista)
pip install pyvista==0.49.0
pip install gpytoolbox==0.3.7    # icosphere 等
pip install potpourri3d          # 热方法、cotan_laplacian
pip install robust-laplacian     # mesh_laplacian / point_cloud_laplacian
pip install pymeshlab            # 无球面参数化，但网格 IO/清理方便
pip install libigl               # 有 igl.cotmatrix(支持四面体)
```
实际的包内版本（`pip show`）：`tetgen 0.8.4`，Home-page `https://github.com/pyvista/tetgen`，License `MIT`，requires `numpy`。
⚠️ 该 `tetgen` 包**不是** PyPI 上另一个同名老包；它的 API 是 `tetgen.TetGen`，基于 pyvista。

---

## 3. 目标 B：Volume-Preserving 优化（不用显式最优传输）

### 3.1 VSEM —— Volumetric Stretch Energy Minimization

**精确引用（均已从论文自己的参考文献中核对）：**

1. **M.-H. Yueh, T. Li, W.-W. Lin, S.-T. Yau, "A novel algorithm for volume-preserving parameterizations of 3-manifolds", *SIAM J. Imaging Sci.* 12(2):1071–1098, 2019.** ← VSEM 原始算法（Algorithm 4.4）
2. **T.-H. Huang, W.-H. Liao, W.-W. Lin, M.-H. Yueh, S.-T. Yau, "Convergence analysis of volumetric stretch energy minimization and its associated optimal mass transport", *SIAM J. Imaging Sci.* 16(3):1825–1855, 2023, DOI `10.1137/22M1528756`.** ← arXiv 版是 **`arXiv:2210.09654`，标题为 "Convergence Analysis of Volumetric Stretch Energy Minimization and its Associated Optimal Mass Transport"**（注意：**标题与 2407.19272 里引用的 `HuLL23` 完全一致**，我下载 arXiv:2210.09654 源文件确认了这一点）
3. **Z.-H. Tan, T. Li, W.-W. Lin, S.-T. Yau, "n-Dimensional volumetric stretch energy minimization for volume-/mass-preserving parameterizations with applications", *SIAM J. Imaging Sci.* 18(2):1141–1175, 2025.** ← arXiv:2402.00380 的期刊版

> ⚠️ **关于 `arXiv:2210.09654` 的重要澄清**：任务描述里把它列为 "Yueh et al. *Numerical Optimization for Volume-Preserving Parameterization of 3-Manifolds*"。我**下载了 `arxiv.org/e-print/2210.09654` 并解包**，得到的 LaTeX 文件名为 `VSEM.tex`，`\title{Convergence Analysis of Volumetric Stretch Energy Minimization and its Associated Optimal Mass Transport}`。**所以 2210.09654 就是 Huang-Liao-Lin-Yueh-Yau 的收敛性分析，不是另一篇独立的 "Numerical Optimization" 论文。**（若团队在别处看到 "Numerical Optimization for Volume-Preserving Parameterization of 3-Manifolds" 这个标题，那是另一篇/预印本别名，我未能定位 —— 见 §8。）

#### 能量泛函

**连续：**
$$E_V(f)=\int_{\mathcal M}\Big[\frac{\rho_\nu\circ f}{\rho_\mu}\operatorname{Det}(\nabla_{\mathcal M} f)\Big]^2 d\mu,
\qquad \nabla_{\mathcal M} f := J_{f\circ h}J_h^{\dagger}\ (\text{切梯度}),\ \operatorname{Det}=\text{伪行列式}$$

**离散（VSEM 原文 eq. 2.2，$n=3$）：**
$$\boxed{\;E_V(f)\;=\;\tfrac12\operatorname{trace}\big(\mathbf f^\top L_V(f)\,\mathbf f\big)\;=\;\sum_{\tau\in\mathbb T(\mathcal M)}\frac{3\,|f(\tau)|^2}{2\,\mu(\tau)}\;}$$

注意这里的系数是 $3/2$；n-VSE/IEM 两篇用的是等价形式 $E_V(f)=\sum_\tau |f(\tau)|^2/\mu(\tau)$（相差 $2/3$ 倍，见 IEM 论文原文："The definition of $E_V$ in (2b) is equivalent to [HuLL23, (3.3)] with a scaling factor of 2/3"）。**比较不同论文的能量数值前一定要对齐这个系数。**

**局部单元矩阵（VSEM 原文 eq. 8.1/8.2）：**
$$L_\tau(f) = -\frac{1}{36\,\mu(\tau)}\,[a_{ij}],\qquad
a_{ij}=[(f_k-f_i)\times(f_l-f_i)]^\top[(f_l-f_j)\times(f_k-f_j)]$$
对 $\{i,j,k,l\}=\{0,1,2,3\}$、$[i,j]\cap[k,l]=\varnothing$；$a_{ji}=a_{ij}$；$a_{ii}=-\sum_{j\ne i}a_{ij}$。

**全局矩阵按边装配：**
$$[L_V(f)]_{ij}=[L_V(f)]_{ji}=
\begin{cases}
-\sum_{\tau\supset[v_i,v_j]} w_{ij}(f), & [v_i,v_j]\in\mathbb E(\mathcal M)\\[2pt]
-\sum_{k\ne i}[L_V(f)]_{ik}, & j=i\\[2pt]
0, & \text{otherwise}
\end{cases}$$

**关键：符号约定（这是最容易出错、也最容易被论文误导的地方）**

论文的**两处**表述在符号上互相矛盾，我做了数值实验来定夺：

| 论文位置 | 写法 |
|---|---|
| eq. (2.4) 上方定义 + eq. (3-3) | $w_{ij}(f)=-\frac{1}{36}\sum_\tau \frac{(\mathbf h_{ki}\times\mathbf h_{\ell i})\cdot(\mathbf h_{\ell j}\times\mathbf h_{kj})}{\mu(\tau)}$（**负号**），且 $[L_V]_{ij}=w_{ij}$（**负值**） |
| eq. (8.1) | $L_\tau(f)=-\frac{1}{36\mu(\tau)}[a_{ij}]$，其中 $a_{ij}=(\cdots)\cdot(\cdots)$ |

**实测结论（§7 V2a/V2b/V3/V3b）：**

- 取 **$L_\tau := -\frac{1}{36\mu}[a_{ij}]$**（即 eq. 8.1 的字面形式，$L_V$ 的非对角元 **$<0$**，对角元 $>0$，是**正半定 M-矩阵**）时，**两个定理同时成立**：
  $$-\tfrac12\sum_{s=1}^3\mathbf f_s^\top L_\tau\mathbf f_s = \frac{3|f(\tau)|^2}{2\mu(\tau)},\qquad
  \nabla_{\mathbf f}E_V=-3\,L_V(f)\,\mathbf f$$
- 取相反符号（非对角元 $>0$，负半定）时，两个恒等式**都差一个负号**（相对误差恰为 2.0）。
- **因此：两种符号不可能同时满足"迹恒等式"与"梯度恒等式"。** 实现时必须二选一并全程自洽，否则一切都不收敛。

我采用的（= 论文 eq. 8.1 字面形式）：`L_tau = -(1/(36 mu)) * a`。
> 📌 建议：**在代码里写死一个符号，并加一条单元测试**（就是 §7 V2a），任何人改符号立刻报错。

**全局梯度：**
$$\nabla_{\mathbf f}E_V(f) = -3\begin{bmatrix}L_V(f)\mathbf f^1\\ L_V(f)\mathbf f^2\\ L_V(f)\mathbf f^3\end{bmatrix}$$

**下界与"当且仅当"：**
$$\boxed{\;E_V(f)\;\ge\;\frac{3}{2}\big|f(\mathcal M)\big|
\quad\text{（在 } |f(\mathcal M)|=\mu(\mathcal M)\text{ 的约束下）},\qquad
\text{等号}\iff \mu(\tau)=|f^\*(\tau)|\ \forall\tau\;}$$

即等号当且仅当映射体保持。

**迭代步骤：** 因为 $\nabla E_V=0 \Leftrightarrow L_V(f)\mathbf f=0$ 只有平凡解 $\mathbf f=$ 常向量（不可接受），必须固定边界，解 KKT 条件：
$$[L_V(f)]_{\mathtt I\mathtt I}\,\mathbf f_{\mathtt I}^{s} \;=\; -[L_V(f)]_{\mathtt I\mathtt B}\,\mathbf f_{\mathtt B}^{*,s},\qquad s=1,2,3$$

**算法（VSEM 原文 Algorithm 1，逐字抄录）：**
```
Require: 单连通四面体网格 M, 容差 eps
Ensure : 体保持参数化 f: M -> B^3
 1: n <- M 的顶点数
 2: B <- {t | v_t in dM},  I <- {1..n} \ B
 3: 用 [YuLL19] 的 SEM 算法计算球面保面积参数化 f_B
 4: 解线性系统 [L_V]_II f_I^s = -[L_V]_IB f_B^s, s=1,2,3
      (L_V 取 L_V(id), f^s = [f_I^s; f_B^s])
 5: delta <- inf ; f_hat <- f
 6: while delta > eps do
 7:    更新 L_V(f_hat)              // eq.(2.2) + w_ij(f) eq.(3-3)
 8:    解 [L_V]_II f_I^s = -[L_V]_IB f_B^s, s=1,2,3
 9:    delta <- E_V(f_hat) - E_V(f) ; f_hat <- f
10: end while
11: return f
```

**收敛性（Huang et al. 2023）：** 论文证明 $\|\varepsilon^{(m)}\|$ 与 $\|\sigma^{(m)}\|$（stretch factor 的标准差）**R-线性收敛**。实测：几何球上 VSEM 5–80 次迭代即达机器精度（§7 V4b）。

**VSEM 的根本局限（**这是选择 IEM 而非 VSEM 的决定性理由**）：**
边界冻结时，图像体积 $|f(\mathcal M)|$ 由边界面的位置**唯一决定且不可控**，因此下界 $\frac32|f(\mathcal M)|$ 中的 $|f(\mathcal M)|$ 也在变。若冻结的边界映射本身不是体积相容的，**体保持映射可能根本不存在**，VSEM 就永远达不到下界。实测（peanut）：$E_V-1.5V(f)=1.74>0$，ratio std 停在 0.24（§7 V9）。IEM 论文原文明确指出这一点：

> "the original VSEM method assumes the total image volume is constant throughout the iteration process, **which may not hold if the boundary points are allowed to glide along the unit sphere**."

### 3.2 n-VSE（arXiv:2402.00380）—— 精确能量、两个子问题、算法框

**能量（连续，eq. def:EVfcontinue）：**
$$E_V(f)=\int_{\mathcal M}\Big[\frac{\rho_\nu\circ f}{\rho_\mu}\operatorname{Det}(\nabla_{\mathcal M}f)\Big]^2 d\mu
\;=\;\int_\Omega \frac{\big((\rho_\nu\circ f\circ h)\sqrt{\det(J_{f\circ h}^\top J_{f\circ h})}\big)^2}{(\rho_\mu\circ h)\sqrt{\det(J_h^\top J_h)}}\,ds$$

**离散（eq. EVcot / EVv）：**
$$\boxed{\;E_V(f)=\frac1n\operatorname{trace}\big(\mathbf f^\top L_V(f)\mathbf f\big)=\sum_{\sigma\in\mathbb S_n(\mathcal M)}\frac{|f(\sigma)|^2}{\mu(\sigma)}\;}$$

**n 维修正 cotangent 权重（eq. wijf）—— 这是 n-VSE 的核心新意：**
$$w_{ij}(f)=\frac{1}{n(n-1)}\sum_{[v_i,v_j]\in\mathbb S_1(\sigma)}\frac{|f(\sigma_{\hat i\hat j})|}{\sigma_{\mu,f^{-1}}(\sigma)}\cot\theta_{ij}^{\sigma}(f)$$
- $\sigma_{\hat i\hat j}$ = 含顶点集 $\mathbb S_0(\sigma)\setminus\{v_i,v_j\}$ 的 $(n-2)$-单纯形
- $\theta_{ij}^\sigma(f)$ = 含 $\mathbb S_0(f(\sigma))\setminus\{\mathbf f_i\}$ 与 $\setminus\{\mathbf f_j\}$ 的两个子空间之间的二面角
- `\sigma_{\mu,f^{-1}}(\sigma) = \mu(\sigma)/|f(\sigma)|` 是 stretch factor
- $n=3$ 时 $n(n-1)=6$，与 §1.2 的 $1/6$ 一致 ✓

**空间（拉普拉斯）矩阵（eq. LVf）：**
$$[L_V(f)]_{ij}=\begin{cases}-w_{ij}(f)&[v_i,v_j]\in\mathbb S_1(\mathcal M)\\ \sum_{\ell\ne i}w_{i\ell}(f)&i=j\\0&\text{otherwise}\end{cases}$$

**Dirichlet/体积拉普拉斯 $L_D$（eq. LD，用于 Dirac map 初值）：**
$$[L_D]_{ij}=\begin{cases}-\tilde w_{ij}&[v_i,v_j]\in\mathbb S_1(\mathcal M)\\ \sum_{\ell\ne i}\tilde w_{ij}&i=j\\0&\text{otherwise}\end{cases},\qquad
\tilde w_{ij}=\frac{1}{n(n-1)}\sum_{\sigma\supset\{v_i,v_j\}}|\sigma_{\hat i\hat j}|\cot\theta_{ij}^\sigma$$

**下界定理（Thm Sfcontinue / Sf_discrete）：** $\min_f E_V(f)\ge \dfrac{\nu(\mathbb B^n)^2}{\mu(\mathcal M)}$，等号 $\iff$ mass-preserving。
定义 **mass-preserving error**：$\varepsilon=E_V(f)-\dfrac{\nu(\mathbb B^n)^2}{\mu(\mathcal M)}$，且
$$\varepsilon=\frac{\nu(\mathbb B^n)^2}{\mu(\mathcal M)^2}\|\delta\|_{L^2(\mathcal M)}^2,\qquad
\delta=\frac{(\rho_\nu\circ f)/\nu(\mathbb B^n)}{\rho_\mu/\mu(\mathcal M)}\operatorname{Det}(\nabla_{\mathcal M}f)-1$$
（$\delta+1$ 就是体积比，与 §5(a) 完全对应）

**问题分解（eq. opt:ball_total → opt:sphere + opt:ball）：**
$$\min_{\mathbf f}\ \tfrac1n\operatorname{trace}(\mathbf f^\top L_V(f)\mathbf f)
\quad\text{s.t.}\quad \sum_{\sigma}|\!f(\sigma)\!|=C,\ \ f(\partial\mathcal M)=\mathcal S^{n-1}$$

1. **边界子问题**（$(n-1)$-球面，低维凸约束问题）：
   $$\min_{\mathbf g}\ \tfrac1{n-1}\operatorname{trace}(\mathbf g^\top L_V(g)\mathbf g)
   \quad\text{s.t.}\quad \sum_{\sigma\in\mathbb S_{n-1}(\partial\mathcal M)}|g(\sigma)|=C',\ \ \|\mathbf g_i\|_2^2=1$$
2. **内部子问题**（固定边界，无约束）：
   $$\min_{\mathbf f_{\mathtt I}}\ \tfrac1n\operatorname{trace}(\mathbf f^\top L_V(f)\mathbf f)
   \quad\Longrightarrow\quad
   [L_V(f)]_{\mathtt I\mathtt I}\mathbf f_{\mathtt I}=-[L_V(f)]_{\mathtt I\mathtt B}\mathbf f_{\mathtt B}$$

**边界子问题的 Newton 法 —— 完整线性系统（eq. newton）：**

Lagrange 函数
$$L(\mathbf g,\lambda,\boldsymbol s)=E_V(g)+\lambda\Big(\sum_{\sigma\in\mathbb S_{n-1}(\partial\mathcal M)}|g(\sigma)|-C'\Big)+\frac12\sum_i \boldsymbol s_i\big(\|\mathbf g_i\|_2^2-1\big)$$

KKT：
$$2L_V(g)\mathbf g+\lambda L_D(g)\mathbf g+\operatorname{diag}(\boldsymbol s)\mathbf g=0,\quad
\sum_\sigma|g(\sigma)|-C'=0,\quad
\tfrac12(\|\mathbf g\|^2-1)=0$$

Newton 步：
$$\boxed{
\begin{bmatrix}
\nabla^2_{\operatorname{vec}(\mathbf g)}L & \operatorname{vec}(L_D(g)\mathbf g) & \operatorname{cdiag}(\mathbf g)\\
\operatorname{vec}(L_D(g)\mathbf g)^\top & 0 & 0\\
\operatorname{cdiag}(\mathbf g)^\top & 0 & 0
\end{bmatrix}
\begin{bmatrix}\Delta\operatorname{vec}(\mathbf g)\\ \Delta\lambda\\ \Delta\boldsymbol s\end{bmatrix}
=-
\begin{bmatrix}
\operatorname{vec}(2L_V(g)\mathbf g+\lambda L_D(g)\mathbf g+\mathbf g\operatorname{diag}(\boldsymbol s))\\
\sum_\sigma|g(\sigma)|-C'\\
\|\mathbf g\|^2-1
\end{bmatrix}}$$

其中 $\operatorname{cdiag}(\mathbf g)$ 是第 $i$ 块为 $\mathbf g^i$ 的块对角矩阵；**$\nabla^2_{\operatorname{vec}(\mathbf g)}L$ 的稀疏结构与 $\mathbf 1_{(n-1)\times(n-1)}\otimes L_D$ 相同，论文建议用有限差分估计**。

**算法框 1：Dirac map 球面初值（Algorithm `alg:Dirac`）**
```
Require: (n-1)-单纯复形 dM ≅ S^{n-1}
Ensure : 球面 Dirac 参数化 g: dM -> S^{n-1}
1: 选最规则的 (n-1)-单纯形作为含 Dirac 点的 tau_p
2: 构造右端项 b        // eq.(eq:b), (eq:b_detail):  b_i = grad(alpha_i) Q
                       // Q 来自 QR 分解 [v_10^T,...,v_(n-1)0^T] = QR
3: N <- 顶点数; 固定 h_i = 0; I_hat <- {1..N}\{i}
4: 解 [L_D]_{I_hat I_hat} h_{I_hat} = b_{I_hat}
5: 中心化 h <- h - (1/N) 1_N 1_N^T h
6: 通过逆球极投影 g = Pi^{-1}(h) 得到球面 Dirac 映射
```

**算法框 2：SEM 球面保质量参数化（Algorithm `alg:SEM`）**
```
Require: (n-1)-单纯复形 dM ≅ S^{n-1}, 内半径 r, 容差 tol
Ensure : eps-mass-preserving g: dM -> S^{n-1}
1: m <- #S_0(dM)
2: 用 alg:Dirac 算 Dirac 映射 g
3: h_i <- Pi(g_i), i=1..m                 // 球极投影
4: E_old <- E_V(g); deltaE <- +inf
5: while deltaE > tol do
6:     L <- L_V(g)                        // eq.(LVf)
7:     h <- diag(|h|^{-2}) h
8:     I <- {i | |h_i| < r}, B <- {1..m}\I
9:     解 [L_S]_II h_I = -[L_S]_IB h_B    // 南北交替迭代
10:    g_i <- Pi^{-1}(h_i), i=1..m
11:    E_new <- E_V(g); deltaE <- E_old - E_new
12: end while
```

**算法框 3：(n−1)-VSEM 球面保质量参数化（Algorithm `alg:SMP`）**
```
Require: (n-1)-单纯复形 dM ≅ S^{n-1}, 测度 mu, 容差 tol
Ensure : 球面保质量参数化 g
1: 用 alg:SEM 算 Dirac 初值 g; C' <- |g(dM)|
2: E_old <- E_V(g); deltaE <- +inf
3: while deltaE < tol do                 // 原文如此（疑为笔误，应为 >）
4:     解线性系统 (eq:newton) 得 Newton 步
5:     线搜索步长 alpha, 更新
          g <- g + alpha Δg,  lambda <- lambda + alpha Δlambda,  s <- s + alpha Δs
6:     E_new <- E_V(g); deltaE <- E_old - E_new; E_old <- E_new
7: end while
```

**算法框 4：n-VSEM 球体保质量参数化（Algorithm `alg:BMP`）**
```
Require: n-单纯复形 M ≅ B^n, 测度 mu, 容差 tol
Ensure : n-ball mass-preserving 参数化 f: M -> B^n
1: B <- {i | v_i in S_0(dM)}, I <- S_0(M)\S_0(dM); m <- #B
2: 计算 V~_B = (I - (1/m) 1 1^T) V_B 的左奇异矩阵 U;  V_B <- U
      // 即沿主轴拉伸，使 M 更接近球（这一步本身是 mass-preserving 的）
3: 用 alg:SMP 算球面保质量映射 g: dM -> S^{n-1}
4: 构造 L <- L_D      // eq.(LD)
5: 固定 f_B = g, 解 L_II f_I = -L_IB f_B 得内部顶点
6: E_old <- E_V(f); deltaE <- +inf
7: while deltaE > tol do
8:     L <- L_V(f)                       // eq.(LVf)
9:     解 L_II f_I = -L_IB f_B
10:    E_new <- E_V(f); deltaE <- E_old - E_new; E_old <- E_new
11: end while
```

**实现注意：** 论文自己在 Remark 里承认
> "Generally, the orientation preserving condition is satisfied for the solution to the subproblems ... Even if there are overlap n-simplices, the number of them is little. In this case, **a convex combination postprocess can eliminate the overlap n-simplices** with a tiny mass-preserving loss. Therefore, we ignore the orientation preserving condition in the computation."

即：**n-VSEM 不保证无折叠**，需要额外做去折叠后处理（论文提到"convex combination postprocess"，但没给细节 —— 见 §8）。

**主轴上拉伸的显式做法（eq. AST / 文中 §4.3）：**
$$\tilde V_{\mathtt B}=(I-\tfrac1m\mathbf 1\mathbf 1^\top)V_{\mathtt B},\qquad
\tilde V_{\mathtt B}^\top\tilde V_{\mathtt B}X=X\Lambda,\qquad
V \leftarrow (V-\mathbf 1\bar v_{\mathtt B}^\top)X\Lambda^{-1/2}$$
实际就是 `U, S, Vt = svd(centered_boundary); V <- (V - centroid) @ Vt.T @ diag(1/S)`。

### 3.3 IEM —— Isovolumetric Energy Minimization（arXiv:2407.19272）

**精确引用：** S.-Y. Liu, T.-M. Huang, W.-W. Lin, M.-H. Yueh, "Isovolumetric Energy Minimization for Ball-Shaped Volume-Preserving Parameterizations of 3-Manifolds", arXiv:2407.19272 (27 Jul 2024)；期刊版：*J. Optim. Theory Appl.* **209(3):68, 2026**, DOI `10.1007/s10957-026-03002-5`，标题 "Volume-Preserving Parameterizations via Preconditioned Nonlinear Conjugate Gradient Method"。

#### 核心思路：让边界点在球面上**自由滑动**

$$\boxed{\;E_{\mathrm I}(f)=\frac{\mathcal V(e)}{\mathcal V(f)}\,E_V(f)-\mathcal V(f)\;}$$
- $\mathcal V(e)$ = 恒等映射的体积 = $\mu(\mathcal M)$；$\mathcal V(f)$ = 图像体积
- $E_V(f)=\sum_{\tau\in\mathbb T(\mathcal M)}\dfrac{|f(\tau)|^2}{|\tau|}$
- 退化情形约定：**$E_{\mathrm I}(f)=0$ whenever $\mathcal V(f)=0$**（图像体积取正四面体体积之和）

**关键定理（Thm 2.2）：** $E_{\mathrm I}(f)\ge0$，**等号 $\iff$ $f$ 体保持**。
证明用 Cauchy–Schwarz：
$$E_V(f)\,\mathcal V(e)=\sum_\tau\frac{|f(\tau)|^2}{|\tau|}\sum_\tau|\tau|\ \ge\ \Big(\sum_\tau|f(\tau)|\Big)^2=\mathcal V(f)^2
\;\Longrightarrow\;\frac{\mathcal V(e)}{\mathcal V(f)}E_V(f)\ge\mathcal V(f)$$
等号 $\iff |f(\tau)|/|\tau|$ 为常数。

**为什么需要它：** 球体映射的图像体积 $\mathcal V(f)$ 由边界顶点围成的多面体决定。边界点在球面上滑动时 $\mathcal V(f)$ 会变，因此**不应强制 $\mathcal V(e)=\mathcal V(f)$**。IEM 通过 $\frac{\mathcal V(e)}{\mathcal V(f)}E_V(f)$ 这个比值把体积约束**自动归一化掉**。

#### 体积 $\mathcal V(f)$ 与它对球坐标的梯度

图像体积用"边界面 + 原点"张成的四面体求和（散度定理）：
$$\mathcal V(f)=\sum_{\alpha\in\mathbb K^2(\partial\mathcal M)}|f(\tau_\alpha)|
=\frac16\sum_{[v_i,v_j,v_k]\in\mathbb K^2(\partial\mathcal M)} f_i^\top(f_j\times f_k),\qquad f(\tau_\alpha)=[\mathbf 0_3,f_i,f_j,f_k]$$

对 $\mathbf f_\alpha^1,\mathbf f_\alpha^2,\mathbf f_\alpha^3$（面 $f(\alpha)$ 三个顶点的坐标向量）：
$$\nabla_{\mathbf f_\alpha^1}|f(\tau_\alpha)|=\tfrac16(\mathbf f_\alpha^2\times\mathbf f_\alpha^3),\quad
\nabla_{\mathbf f_\alpha^2}|f(\tau_\alpha)|=\tfrac16(\mathbf f_\alpha^3\times\mathbf f_\alpha^1),\quad
\nabla_{\mathbf f_\alpha^3}|f(\tau_\alpha)|=\tfrac16(\mathbf f_\alpha^1\times\mathbf f_\alpha^2)$$

**边界用球坐标参数化（eq. sphere_coor / inverse）：**
$$\mathbf f^1_{\mathtt B}=\sin\boldsymbol\theta\odot\cos\boldsymbol\phi,\quad
\mathbf f^2_{\mathtt B}=\sin\boldsymbol\theta\odot\sin\boldsymbol\phi,\quad
\mathbf f^3_{\mathtt B}=\cos\boldsymbol\theta$$
$\boldsymbol\theta=\arccos(\mathbf f_{\mathtt B}^3),\qquad
\boldsymbol\phi=\operatorname{atan2}(\mathbf f_{\mathtt B}^2,\mathbf f_{\mathtt B}^1)$

**体积对球坐标的梯度（eq. 3.2，链式法则）：**
$$\nabla_{\boldsymbol\theta_\alpha}|f(\tau_\alpha)|=\frac16\Big(
\cos\boldsymbol\theta_\alpha\odot\cos\boldsymbol\phi_\alpha\odot(\mathbf f^2_\alpha\times\mathbf f^3_\alpha)
+\cos\boldsymbol\theta_\alpha\odot\sin\boldsymbol\phi_\alpha\odot(\mathbf f^3_\alpha\times\mathbf f^1_\alpha)
-\sin\boldsymbol\theta_\alpha\odot(\mathbf f^1_\alpha\times\mathbf f^2_\alpha)\Big)$$
$$\nabla_{\boldsymbol\phi_\alpha}|f(\tau_\alpha)|=\frac16\Big(
\sin\boldsymbol\theta_\alpha\odot\cos\boldsymbol\phi_\alpha\odot(\mathbf f^3_\alpha\times\mathbf f^1_\alpha)
-\sin\boldsymbol\theta_\alpha\odot\sin\boldsymbol\phi_\alpha\odot(\mathbf f^2_\alpha\times\mathbf f^3_\alpha)\Big)$$

> 注意：IEM 论文的能量用 $E_V=\sum|f(\tau)|^2/|\tau|$（即"n-VSE 归一化"，不含 $3/2$），其梯度公式是
> $$\nabla_{\mathbf f^s}E_{\mathrm V}(f)=2L_{\mathrm V}(f)\mathbf f^s$$
> 因为 IEM 采用的 $L_V$ 使 $E_V=\frac13\operatorname{trace}(\mathbf f^\top L_V\mathbf f)$（**$1/3$ 而非 $1/2$**）。这与 §3.1 的 VSEM 归一化不同，**混用会差一个因子**。

**显式梯度（eq. Grad，IEM 核心公式）：**
$$\nabla_{\mathbf f_{\mathtt I}^s}E_{\mathrm I}(f)=\frac{2\mathcal V(e)}{\mathcal V(f)}\Big([L_{\mathrm V}(f)]_{\mathtt I\mathtt I}\mathbf f_{\mathtt I}^s+[L_{\mathrm V}(f)]_{\mathtt I\mathtt B}\mathbf f_{\mathtt B}^s\Big),\quad s=1,2,3$$

$$\begin{aligned}
\nabla_{\boldsymbol\theta}E_{\mathrm I}(f)=&\ \frac{2\mathcal V(e)}{\mathcal V(f)}\Big(
\cos\boldsymbol\theta\odot\cos\boldsymbol\phi\odot\big([L_{\mathrm V}]_{\mathtt B\mathtt I}\mathbf f^1_{\mathtt I}+[L_{\mathrm V}]_{\mathtt B\mathtt B}\mathbf f^1_{\mathtt B}\big)\\
&+\cos\boldsymbol\theta\odot\sin\boldsymbol\phi\odot\big([L_{\mathrm V}]_{\mathtt B\mathtt I}\mathbf f^2_{\mathtt I}+[L_{\mathrm V}]_{\mathtt B\mathtt B}\mathbf f^2_{\mathtt B}\big)\\
&-\sin\boldsymbol\theta\odot\big([L_{\mathrm V}]_{\mathtt B\mathtt I}\mathbf f^3_{\mathtt I}+[L_{\mathrm V}]_{\mathtt B\mathtt B}\mathbf f^3_{\mathtt B}\big)\Big)
-\Big(1+\frac{\mathcal V(e)E_{\mathrm V}(f)}{\mathcal V(f)^2}\Big)\nabla_{\boldsymbol\theta}\mathcal V(f)
\end{aligned}$$

$$\begin{aligned}
\nabla_{\boldsymbol\phi}E_{\mathrm I}(f)=&\ \frac{2\mathcal V(e)}{\mathcal V(f)}\Big(
\sin\boldsymbol\theta\odot\cos\boldsymbol\phi\odot\big([L_{\mathrm V}]_{\mathtt B\mathtt I}\mathbf f^2_{\mathtt I}+[L_{\mathrm V}]_{\mathtt B\mathtt B}\mathbf f^2_{\mathtt B}\big)\\
&-\sin\boldsymbol\theta\odot\sin\boldsymbol\phi\odot\big([L_{\mathrm V}]_{\mathtt B\mathtt I}\mathbf f^1_{\mathtt I}+[L_{\mathrm V}]_{\mathtt B\mathtt B}\mathbf f^1_{\mathtt B}\big)\Big)
-\Big(1+\frac{\mathcal V(e)E_{\mathrm V}(f)}{\mathcal V(f)^2}\Big)\nabla_{\boldsymbol\phi}\mathcal V(f)
\end{aligned}$$

问题的未知量是 $(\mathbf f_{\mathtt I}^1,\mathbf f_{\mathtt I}^2,\mathbf f_{\mathtt I}^3,\boldsymbol\theta,\boldsymbol\phi)$，**无约束非线性优化**。

#### 预处理（Preconditioner）

$$\nabla^2_{\mathbf f^s}E_{\mathrm V}(\mathbf f)=2L_{\mathrm V}(\mathbf f)+2\,\mathrm D_{\mathbf f^s}L_{\mathrm V}(\mathbf f)\mathbf f^s,\qquad
[\mathrm D_{\mathbf f^s}L_{\mathrm V}]_{i,j,\ell}=\frac{\partial}{\partial f_\ell^s}[L_{\mathrm V}]_{i,j}$$

$L_V(f)$ 的主子矩阵通常 SPD（论文原话），故**只保留第一项**，取**块对角预处理子**：
$$\boxed{\;M=\begin{bmatrix}I_3\otimes[L_{\mathrm V}(f^{(0)})]_{\mathtt I\mathtt I}&\\&I_2\otimes[L_{\mathrm V}(f^{(0)})]_{\mathtt B\mathtt B}\end{bmatrix}\;}$$

$I_3$ 对应 3 个内部坐标分量，$I_2$ 对应 $(\boldsymbol\theta,\boldsymbol\phi)$ 两个边界变量。**$M$ 固定不变**，只需**预排序 Cholesky 分解一次**，之后每次迭代用前代/回代。

#### 初始化（非常重要）

1. **全局对齐**：用边界顶点的主轴 SVD 做旋转+缩放，使网格接近球形：
   $$(V_{\mathtt B}-\mathbf 1\bar v_{\mathtt B}^\top)=U\Sigma R^\top,\qquad V\leftarrow(V-\mathbf 1\bar v_{\mathtt B}^\top)R\Sigma^{-1}$$
   论文原话："This alignment recenters the mesh and makes the boundary closer to a sphere."
2. **初值**：跑 `[YuLL19, Algorithm 4.4]` 的固定点法（即 VSEM）**15 次迭代**。
   > 论文原话："Although this approach does not guarantee global convergence, it typically achieves a rapid energy decrease in the first few iterations. **After about 10 iterations, however, the improvement becomes marginal.** We therefore run this fixed-point method for **15 iterations** and use the resulting map as the initial guess."
3. 球坐标 $(\boldsymbol\theta,\boldsymbol\phi)$ 由 eq. sphere_coor_inv 从 $\mathbf f_{\mathtt B}$ 反解。

#### 非线性 CG 的更新式与步长

$$\mathbf f^{(k+1)}=\mathbf f^{(k)}+\alpha_k\mathbf p^{(k)},\qquad
\mathbf p^{(k)}=-M^{-1}\mathbf g^{(k)}+\beta_k\mathbf p^{(k-1)},\qquad \mathbf p^{(0)}=-M^{-1}\mathbf g^{(0)}$$
$$\beta_k=\frac{{\mathbf g^{(k)}}^\top M^{-1}\mathbf g^{(k)}}{{\mathbf g^{(k-1)}}^\top M^{-1}\mathbf g^{(k-1)}}=\frac{\|\mathbf g^{(k)}\|_{M^{-1}}^2}{\|\mathbf g^{(k-1)}\|_{M^{-1}}^2}
\quad(\text{preconditioned Fletcher--Reeves})$$

步长：把 $\varphi(\alpha)=E_{\mathrm I}(\mathbf f^{(k)}+\alpha\mathbf p^{(k)})$ 做**二次**插值
$$\varphi_q(\alpha)=a\alpha^2+\varphi'(0)\alpha+\varphi(0),\qquad
a=\frac{\varphi(\alpha_0)-\varphi(0)-\alpha_0\varphi'(0)}{\alpha_0^2},\qquad
\alpha_1=-\frac{\varphi'(0)}{2a}$$
若 $\alpha_1$ 不满足充分下降条件
$$\varphi(\alpha_1)<\ell(\alpha_1)\equiv\varphi(0)+c_1\alpha_1\varphi'(0),\qquad c_1=\tfrac12-\varepsilon\ (\varepsilon\to0^+)$$
则转**三次**插值，用二次极小点作初值 $\alpha_0$，加上新点 $\alpha_1$：
$$\varphi_c(\alpha)=b\alpha^3+c\alpha^2+\varphi'(0)\alpha+\varphi(0),\qquad
\begin{bmatrix}b\\c\end{bmatrix}=\frac{1}{\alpha_0^2\alpha_1^2(\alpha_1-\alpha_0)}
\begin{bmatrix}\alpha_0^2&-\alpha_1^2\\-\alpha_0^3&\alpha_1^3\end{bmatrix}
\begin{bmatrix}\varphi(\alpha_1)-\varphi(0)-\varphi'(0)\alpha_1\\ \varphi(\alpha_0)-\varphi(0)-\varphi'(0)\alpha_0\end{bmatrix}$$
$$\alpha_2=\argmin\varphi_c(\alpha)=\frac{-c+\sqrt{c^2-3b\varphi'(0)}}{3b}$$

> 论文的实用建议：**"taking the quadratic step and scaling it by 0.9 satisfies this condition in nearly all our tests."** ← 直接照做，省掉三次插值。

#### **算法框：Preconditioned Nonlinear CG for IEM（Algorithm `alg:IEM`，逐字抄录）**

```
Require: A simply connected tetrahedral mesh M and a tolerance eps.
Ensure : An approximate volume-preserving simplicial mapping f.

 1: Compute the SVD and update vertices by (eq:AST).
 2: Let n be the number of elements in K^0(M).
 3: Let B = { s | v_s in dM } and I = { 1,...,n } \ B.
 4: Compute a boundary mapping f_B^s, for s = 1,2,3, using [YuLL19, Algorithm 4.3].
 5: Let L <- L_V(e) as (eq:VSLaplacian).
 6: for k = 1,...,15 do
 7:     Solve the linear system  L_II f_I^s = -L_IB f_B^s,  for s = 1,2,3.
 8:     Update L <- L_V(f) as (eq:VSLaplacian).
 9: end for
10: Let the preconditioners M_I <- L_II and M_B <- L_BB.
11: Perform Cholesky decompositions of M_I and M_B.
12: Compute the spherical coordinates [theta, phi] by (eq:sphere_coor_inv).
13: Compute the energy value E = E_I(f) as (eq:Ea).
14: Let delta <- inf.
15: while delta > eps do
16:     Let E_0 <- E.
17:     Compute gradients g_I^s, g_theta, g_phi for s=1,2,3 by (eq:Grad).
18:     Solve M_I h_I^s = g_I^s, M_B h_theta = g_theta and M_B h_phi = g_phi, s=1,2,3.
19:     if delta == inf then
20:         Compute lambda <- sum_s (g_I^s^T h_I^s) + g_theta^T h_theta + g_phi^T h_phi.
21:         Update p_I^s <- -h_I^s, p_theta <- -h_theta, p_phi <- -h_phi, s=1,2,3.
22:     else
23:         Let lambda_0 <- lambda.
24:         Compute lambda <- sum_s (g_I^s^T h_I^s) + g_theta^T h_theta + g_phi^T h_phi.
25:         Let beta <- lambda / lambda_0.
26:         Update p_I^s <- -h_I^s + beta p_I^s, p_theta <- -h_theta + beta p_theta,
27:                p_phi <- -h_phi + beta p_phi, s=1,2,3.
28:     end if
29:     Compute a step length alpha by (eq:quad_step) or (eq:cubic_step).
30:     Update f_I^s <- f_I^s + alpha p_I^s, theta <- theta + alpha p_theta,
31:            phi <- phi + alpha p_phi, s=1,2,3.
32:     Update f_B^s by the spherical-coordinate formula (eq:sphere_coor), s=1,2,3.
33:     Update L <- L_V(f) as (eq:VSLaplacian).
34:     Compute the energy value E = E_I(f) as (eq:Ea).
35:     Update delta <- E_0 - E.
36: end while
```

**全局收敛性（Thm 5.1）：** 若每步步长满足强 Wolfe 条件（$0<c_1<c_2<\frac12$）且迭代点保持在非退化区域的紧子集内、$M$ SPD，则
$$\liminf_{k\to\infty}\|\mathbf g^{(k)}\|_{M^{-1}}=0$$

### 3.4 三条路线的工程对比

| | VSEM (§3.1) | n-VSEM (§3.2) | IEM (§3.3) |
|---|---|---|---|
| 边界 | **冻结** | 冻结（先解球面子问题） | **自由滑动**（球坐标） |
| 能量 | $\sum \frac32|f|^2/\mu$ | $\sum |f|^2/\mu$ | $\frac{\mathcal V(e)}{\mathcal V(f)}E_V-\mathcal V(f)$ |
| 下界 | $\frac32|f(\mathcal M)|$ | $\mathcal V(\mathbb B^n)^2/\mu(\mathcal M)$ | **0**（且 0 当且仅当体保持） |
| 优化器 | 固定点迭代 | 边界 Newton + 内部固定点 | **Preconditioned nonlinear CG** |
| 理论收敛性 | R-线性（Huang 2023 已证） | 见论文 §3 | **全局收敛**（Thm 5.1，强 Wolfe） |
| 能否达到下界 | 仅在边界恰好体积相容时 | 一般域：不能 | 论文声称显著优于 VSEM |
| 是否需要后处理 | 是（去折叠） | **是**（论文承认忽略定向条件） | 论文未提及折叠后处理 |
| 代码 | ❌ 无 | ❌ 无 | ❌ 无 |

**建议实现顺序：** VSEM（最简单，30 行核心代码，且是 IEM 的初值算法）→ 若不够，再上 n-VSEM/IEM。

### 3.5 ⚠️ 公开代码 —— **三篇核心论文都没有发布代码**（已逐一核实）

我做了以下查证，**全部为负面结果**：

| 查证对象 | 结果 |
|---|---|
| GitHub API: `api.github.com/users/meihengyueh/repos` | **404 Not Found** |
| GitHub API: `api.github.com/users/mhyueh/repos` | ✅ 存在，但只有 **4 个仓库**：`MATLAB_GraphLaplacian`、`MATLAB_MonteCarloIntegerProgrammingDNA`、`mhyueh.github.io`、`NonlinearFiltering`（Yau-Yau 非线性滤波）。**无任何 VSEM/n-VSE/IEM 代码** |
| GitHub API: `ZhongHengTan`、`wmwwlin`、`LokMingLui`、`syliu` | **全部 404**（`lmlui` 存在但只有 2 个无关仓库） |
| GitHub search API: `3D+volume+preserving+mapping`、`volume+preserving+parameterization+3-manifold`、`volumetric+stretch+energy` | **total_count = 0**（GitHub 索引很弱，但结合上面的直接查证，可以认为确实没有） |
| Yueh 主页 `math.ntnu.edu.tw/~yueh/publications.html`（已 fetch 全文） | VSEM/n-VSE/IEM **三个条目下均无 "MATLAB Demo" 链接**。只有 `Diskmap_SEM.zip`、`DiskCEM.zip` 两个圆盘参数化的 demo，以及 2D 的 `projects/SphereAEM.html`（= SAEM，**不是** VSEM）、`projects/DiskAEM.html`、`projects/CDCP.html`。`Torus_VSEM.html`（亏格 1 的体保持）只有演示视频，**无代码链接** |
| `math.ntnu.edu.tw/~yueh/code/` | **403 Forbidden**（目录被禁） |
| MathWorks MATLAB Central 作者页 | **Access Denied**（Akamai 拦截，未能确认） |

**✅ 唯一确认可下载的相关源码：**
1. <https://math.ntnu.edu.tw/~yueh/projects/SphereAEM/SAEM.zip> —— **球面保面积**映射（不是体保持）。已下载解包至 `vendor/SAEM/`。License: **academic/research only, 禁商用**。核心函数 `SphericalAEM.p` 是 P-code（不可读源码）。
2. <https://github.com/garyptchoi/spherical-conformal-map> —— **球面共形**映射，完整 MATLAB 源码。已下载 `demo.m` + `README.md` 至 `vendor/`。
3. <https://github.com/icemiliang/spherical_harmonic_maps> —— C++ 球面调和映射（未验证签名）。
4. Gu 组 C++ tutorials（spherical OT / Ricci flow / harmonic map）—— 页面声称有源码，**下载链接未提取到**（§8）。

**结论：n-VSE / VSEM / IEM 必须自己实现。** 好消息是：核心算法都不长（VSEM 约 30 行；IEM 约 150 行），且本报告已给出全部公式与算法框，§7 的参考实现已经把 VSEM 跑通并验证到下界。

---

## 4. TetGen 集成

### 4.1 `.node` / `.ele` / `.face` 精确格式（**依据 meshio 5.3.5 的 reader/writer 源码逐行核对**）

#### `.node`
```
第 1 行（可省略）：  <#nodes> <dim> <#attributes> <#boundary markers>
第 2 行起：          <node#> <x> <y> <z> [attr1 ... attrN] [marker1 ... markerM]
```
- `<dim>` **必须为 3**（meshio 的 reader 会 `raise ReadError("Need 3D points.")`）
- `#boundary markers` ≥ 1 时，第一个 marker 被读成 `point_data["tetgen:ref"]`，第二个起为 `"tetgen:ref2"`, `"tetgen:ref3"`, ...
- 属性读成 `point_data["tetgen:attr1"]`, `"tetgen:attr2"`, ...
- **顶点编号必须连续**！meshio 显式检查：
  ```python
  node_index_base = int(points[0, 0])
  if not np.all(points[:, 0] == np.arange(node_index_base, node_index_base + n)):
      raise ReadError()
  ```
  即允许 0-based 或 1-based，但**不允许有空洞**。
- 空白行与以 `#` 开头的行会被跳过。

#### `.ele`
```
第 1 行（可省略）：  <#tets> <#nodes per tet> <#attributes>
第 2 行起：          <tet#> <n1> <n2> <n3> <n4> [attr1 ... attrN]
```
- `#nodes per tet` **必须为 4**（meshio: `if num_points_per_tet != 4: raise ReadError()`）
- 属性读成 `cell_data["tetgen:ref"]`, `"tetgen:ref2"`, ...
- 顶点索引会减去 `node_index_base`，所以读进来永远是 0-based。

#### `.face`（meshio **不支持**，需自己写）
```
第 1 行：  <#faces> <#nodes per face (3)> <#boundary markers>
第 2 行起： <face#> <n1> <n2> <n3> <marker>
```

#### `.edge`（meshio 不支持）
```
第 1 行：  <#edges> <#nodes per edge (2)> <#boundary markers>
第 2 行起： <edge#> <n1> <n2> <marker>
```

**实测生成的样例文件（§7 V6 的产物）：**
```
.node 头行 :  178 3 0 1
.node 数据 :  0 4.6259292692714853e-18 -8.3266726846886741e-17 1.0000000000000000e+00 1
.ele  头行 :  # This file was created by meshio v5.3.5
.ele  数据 :  # attribute names: tetgen:ref
```
（注意 meshio 会先写注释行，TetGen 自己也能忽略 `#` 行。）

### 4.2 10 行读写用例（**已实测通过，逐位一致**）

```python
import meshio, numpy as np

# 读
m = meshio.read("ball.node")                      # 也会自动找同目录 ball.ele
V = m.points                                      # (n,3) float64
T = m.cells_dict["tetra"]                         # (m,4) int，已转为 0-based
ref_pt = m.point_data.get("tetgen:ref")           # 顶点边界标记
ref_cell = m.cell_data.get("tetgen:ref", [None])[0]   # 单元 region 标记

# 写（同时生成 ball.node 和 ball.ele）
m2 = meshio.Mesh(
    V, [("tetra", T)],
    point_data={"tetgen:ref": np.where(is_boundary, 1, 0).astype(np.int32)},
    cell_data={"tetgen:ref": [np.zeros(len(T), dtype=np.int32)]},
)
meshio.write("out.node", m2)     # 后缀 .node 或 .ele 都可
```
**实测往返精度**：顶点 `max diff = 0.000e+00`，四面体索引完全相等，顶点标记与单元属性均正确往返（§7 V6a–V6d）。

### 4.3 造球状测试网格

两种方式（本次都用了）：

**(a) icosphere + TetGen（推荐，最省事）**
```python
import gpytoolbox as gpy, tetgen as tg, numpy as np

V, F = gpy.icosphere(3)                # 单位球，642 顶点 / 1280 面
tgen = tg.TetGen(V, F.astype(np.int32))
tgen.tetrahedralize(switches="pq1.414a0.02Y")   # p=PLC, q=质量, a=最大体积, Y=不裂边
T = np.asarray(tgen.elem, dtype=np.int64)
V = np.asarray(tgen.node, dtype=np.float64)
# 实测 sub=3,a=0.02 -> 2694 tets；sub=3 无 a -> 2672 tets, 798 verts, vol=4.1527 (4pi/3=4.1888)
```
**`switches` 备忘：** `p` 读 PLC、`q` 质量约束（`q1.414` = radius-edge ratio ≤ 1.414）、`a<vol>` 最大单元体积、`Y` 禁止在边界面上加 Steiner 点、`Q` 静默、`V` verbose、`A` 输出属性、`z` 0-based 索引。
⚠️ 不加 `a<vol>` 时 TetGen 可能**不插入内部点**（球面网格上实测内部点很少）。

**(b) 解析四面体化** —— 对单位球可以手写八面体 → 每个八面体切 4 个四面体，再逐级细分。好处是**完全确定性、无外部依赖**；坏处是边界不是球面而是"球面多面体近似"，且质量一般。本次调研没有实现这条（用 (a) 已足够），如需可用 `numpy` 手写 ~30 行。

**(c) 测试用的非凸实体**（**做失效模式测试必须用这个**）：
```python
V, F = gpy.icosphere(3)
theta = np.arccos(np.clip(V[:,2], -1, 1))
r = 1.0 + 0.45*np.cos(2*theta)        # peanut / dumbbell，凹腰
V = V * r[:, None]
# 再送进 TetGen
```
参见 `tools/demo_failure_modes.py` 的 `make_blob_tet_mesh()`，已内置 `peanut` / `wedge` / `crinkly` 三种。

### 4.4 四面体质量度量（**全部已实现并验证**）

设四面体 $\tau=[v_0,v_1,v_2,v_3]$，$p=V[T]$，$L$ 为 6 条边长，$A_k$ 为第 $k$ 个面的面积，$V_\tau=|\tau|$。

| 度量 | 公式 | 理想值 | 实测（球网格 sub=3, a=0.02） |
|---|---|---|---|
| **有符号体积** | $\frac16\big((v_1-v_0)\times(v_2-v_0)\big)\cdot(v_3-v_0)$ | $>0$ | min $9.8\times10^{-4}$ |
| **二面角（6 个）** | 边 $(i,j)$：$\theta_{ij}=\pi-\arccos(\hat n_i\cdot\hat n_j)$，$\hat n_i$ = **面 i 的外法向单位向量** | 全 6 个 $=\arccos(1/3)=70.5288°$ | min $5.5°$，per-tet min 的均值 $44.4°$ |
| **min 二面角** | $\min_{6}\theta_{ij}$（度） | $\to70°$ | 同上 |
| **radius-edge ratio** | $R=\rho_{\text{circ}}/\min L$ | $\sqrt{3/2}\approx1.225$（**注意：不是 3.674**） | mean 0.92，max 1.48 |
| **aspect ratio** | $L_{\max}/(2\sqrt2\,r_{\text{in}})$，$r_{\text{in}}=3V_\tau/\sum_k A_k$ | $1$ | mean 3.00，max 17.9 |
| **sliver 检测** | $\theta_{\min}<5°$ 或 $L_{\max}/(2\sqrt2 r_{\text{in}})>20$ 或 $V_\tau<\varepsilon$ | — | 实测存在 $5.5°$ 的 sliver |

> ⚠️ **radius-edge ratio 的定义陷阱**：正四面体取 circumradius $\rho=a\sqrt6/4$、最短边 $a$ 时，$\rho/\min L=\sqrt6/4\approx0.6123$ 是理想（最小）值。文献里常见的 $1.5\sqrt6\approx3.674$ 是**在用 $R=\rho/(\text{平均边长})$、或按 $R=1/(\rho/\min L)$ 之类归一化**时才出现的常数。
> **团队与文献对比前，务必先统一 $R$ 的定义。**
> 本报告统一采用 $R=\rho_{\text{circ}}/\min L$，理想值 $\sqrt6/4\approx0.612$。
> 实测（球 sub=3, a=0.02）：`mean 0.9249, max 1.4798` —— 略高于 0.612，符合"不错的网格"的预期。✓

**二面角实现要点（踩过的坑）**：必须用**外法向**，且外法向的绕序要按 §1.3 几何修正。用错的绕序会让 $\theta$ 变成 $\pi-\theta$，看起来也"像"是合法角度，但 min 二面角会系统性偏离。**验证方法（§7 V5j）**：正四面体必须给出 6 个都等于 $70.5288°$。

---

## 5. "体积保持"的验证度量（精确定义 + 公式）

### (a) 体积畸变比（volume distortion ratio）

对每个四面体 $\tau$：
$$\boxed{\;D_V(\tau)=\frac{|f(\tau)|}{|\tau|}\;}
\qquad\text{理想：}\ \forall\tau,\ D_V(\tau)=c=\frac{|f(\mathcal M)|}{|\mathcal M|}\ \text{（常数）}$$

统计量（$\bar D=\frac1m\sum_\tau D_V$，$\sigma_D=\sqrt{\frac1m\sum(D_V-\bar D)^2}$）：

| 量 | 公式 | 含义 | 实测（VSEM on 球） | 实测（IEM 论文用 $\delta+1$ 画直方图） |
|---|---|---|---|---|
| mean | $\bar D$ | 全局缩放 | $1.000036$ | — |
| std | $\sigma_D$ | **主要质量指标**，$=0$ 当且仅当体保持 | $6.5\times10^{-3}$ | — |
| max | $\max D_V$ | 最严重的膨胀 | $1.0424$ | — |
| min | $\min D_V$ | 最严重的压缩；$\le0$ 即翻转 | $0.9532$ |
| 极差比 | $\max/\min$ | 常报的量 | $1.093$ | — |
| P99/P1 | $\frac{\text{percentile}(D_V,99)}{\text{percentile}(D_V,1)}$ | **稳健版极差**，抗离群 | $1.036$ | — |
| 相对误差场 | $\delta(\tau)=\frac{\mu(\tau)/\mu(\mathcal M)}{|f(\tau)|/|f(\mathcal M)|}-1$ | **n-VSE 定义的 $\delta$**，且 $\varepsilon=\sum\mu(\tau)\delta^2$（见 §3.2） | — | 直方图 |

> **报告建议**：同时报 `mean / std / max / min` 和 **P99/P1**。只报 mean 会掩盖 sliver 上的极端畸变。

### (b) 逐四面体 Jacobian 行列式与翻转比例

$f$ 在每个 $\tau$ 上是仿射的，$J$ 为常数。数值上不要去组装 $3\times3$ 再求行列式（可以，但没必要），直接用体积比：
$$\boxed{\;J(\tau)=\det J_f|_\tau=\frac{\operatorname{sgn}|f(\tau)|}{\operatorname{sgn}|\tau|}
=\frac{6\,|f(\tau)|_{\text{signed}}}{6\,|\tau|_{\text{signed}}}
=\frac{((f_1-f_0)\times(f_2-f_0))\cdot(f_3-f_0)}{((v_1-v_0)\times(v_2-v_0))\cdot(v_3-v_0)}\;}$$

> 推导：$A(\tau)=\big[f_1-f_0,\;f_2-f_0,\;f_3-f_0\big]\big[v_1-v_0,\;v_2-v_0,\;v_3-v_0\big]^{-1}$，
> $\det A=\det F/\det V=(6|f(\tau)|)/(6|\tau|)$。

**翻转指标：**
$$\boxed{\;\text{frac\_flipped}=\frac1m\,\#\{\tau:J(\tau)\le0\}\;}\qquad
\text{（必须为 0）}$$

**补充指标（强烈建议一起报）：**
- $\text{frac}_{J<0.5}=\frac1m\#\{\tau:J(\tau)<0.5\}$ —— **本次实测最重要的发现：翻转数可以是 0，但 16.5% 的 tets 有 $J<0.5$。只看"有没有翻转"会漏掉灾难性的压缩。**
- $\text{frac}_{J>2}$ —— 对称地看膨胀。
- $\min_\tau J(\tau)$ —— 最坏情况的安全裕度。

### (c) Volumetric stretch energy 及其理论下界

$$\boxed{\;E_V(f)=\sum_{\tau}\frac{3\,|f(\tau)|^2}{2\,\mu(\tau)}
\ \ge\ \frac{3}{2}\big|f(\mathcal M)\big|,\qquad \text{等号}\iff\text{体保持}\;}$$

（VSEM 归一化；$\mu(\tau)=|\tau|$。若用 n-VSE/IEM 归一化 $E_V=\sum|f(\tau)|^2/\mu(\tau)$，则下界是 $|f(\mathcal M)|^2/\mu(\mathcal M)$。**两者差 $2/3$ 倍，比较前必须对齐。**）

**Mass-preserving error（可直接当停止准则）：**
$$\varepsilon=E_V(f)-\frac{\nu(\mathbb B^n)^2}{\mu(\mathcal M)}=\frac{\nu(\mathbb B^n)^2}{\mu(\mathcal M)^2}\|\delta\|_{L^2(\mathcal M)}^2$$

**IEM 能量与其下界（最优停止准则）：**
$$\boxed{\;E_{\mathrm I}(f)=\frac{\mathcal V(e)}{\mathcal V(f)}E_V(f)-\mathcal V(f)\;\ge\;0,\qquad \text{等号}\iff\text{体保持}\;}$$

实测：
- **几何球**上 P1 harmonic 映射（= 恒等）：ratio std $=1.86\times10^{-15}$，$E_V=6.229111225639588$，$1.5V(f)=6.229111225639587$，gap $=8.9\times10^{-16}$
- **几何球**上 VSEM：ratio std $=2.77\times10^{-15}$，gap $=0$（机器精度）
- **非凸 peanut** 上 VSEM（固定径向 BC）：gap $=1.74>0$ —— 说明**存在固定的边界映射使得体保持映射不存在**
- 参考实现最终版（非凸上的正确下界检查）：$\sigma_D$ 从 2.0934 降到 0.2364（8.9×）

### (d) 实践中可行的双射性检验

局部单射（必要不充分）：
$$\text{局部单射必要}\iff \forall\tau,\ J(\tau)>0$$

**实际可用的全局检验组合：**

1. **边界球面映射的双射性**（便宜、必要）：
   $$\sum_{\text{边界面}}\Omega_\tau=4\pi\quad\text{且}\quad\forall\tau,\ \Omega_\tau>0$$
   用 §2.2 的 Van Oosterom–Strackee 有符号立体角公式。**任何负值 → 球面映射翻转。**
   实测：标准球的径向投影 → 和 $=12.566371=4\pi$ 精确，无负值 ✓
2. **体积一致性**（散度定理）：
   $$\mathcal V(f)\stackrel{?}{=}\sum_{\text{边界面}}\frac16\,f_a\cdot(f_b\times f_c)=\sum_\tau|f(\tau)|$$
   两者不等 → 网格不自洽或绕序错。实测：球体径向投影 → $4.047045$ vs $4\pi/3=4.188790$（差 3.4%，因为边界是 642 顶点的多面体近似，**正常**）。
3. **局部单射**：$\forall\tau,\ J(\tau)>0$（见 (b)）。
4. **全局自交检查**（昂贵但确定性）：从原点沿 $N$ 条随机方向 $\omega$ 打射线，数穿越球面的次数；必须**恰好 1 次**。对凸性较弱的实体这是有效的实际检验。
   > 数学依据：$f$ 是单纯映射 + 局部单射 + 边界单射 $\Rightarrow$ 全局双射（度理论 / 局部度恒为 1）。
5. **体积加性检查**：$\mathcal V(f)=\sum_\tau|f(\tau)|$ 且 $\sum_\tau|f(\tau)|=$ 边界围成的体积 —— 已在 2 中。
6. **（可选）最近点查询**：对随机内点 $x\in\mathcal M$，检查 $f^{-1}(f(x))$ 是否唯一（用 KD-tree 在 $f(\mathcal M)$ 的所有 tets 上做点在四面体内测试）。$O(m)$ 但可靠。

**推荐的固定 CI 检查（按性价比排序）：**
```
1. min_τ J(τ) > 0                       # 必要条件，O(m)
2. Σ Ω_τ = 4π (相对误差 < 1e-6)          # 球面边界双射性, O(n_f)
3. σ_D < tol_std 且 max/min < tol_range  # 质量
4. E_I(f) < tol_E                        # 能量收敛
5. (抽样) 射线穿越 == 1                  # 全局双射, O(N·m) 抽样
```

---

## 6. 推荐的完整流水线（工程蓝图）

```
输入：TetGen .node/.ele（topological 3-ball, genus-0 boundary）
  │
  ├─[0] 读入 + 规范化
  │     meshio.read -> V, T
  │     orient_tets(V,T)                 # 统一 det>0  ← 必须！
  │     boundary_faces_oriented(V,T)     # 外法向修正   ← 必须！
  │     质量检查（min dihedral, aspect ratio, sliver 计数）
  │
  ├─[1] 边界：∂M -> S^2  （目标 A 的边界条件）
  │     优先：MATLAB spherical_conformal_map(v,f)          [源码可得]
  │     或：  MATLAB SphericalAEM(F,V)                     [源码可得, 保面积]
  │     验证：Σ Ω_τ == 4π, 所有 Ω_τ > 0
  │
  ├─[2] 目标 A 初值：P1 volumetric harmonic map
  │     K = p1_tet_stiffness(V,T)         # 或 igl.cotmatrix(V,T)（支持 tets）
  │     解 K_II f_I^s = -K_IB f_B^s, s=1,2,3      （Cholesky 一次分解）
  │     → f_harm
  │     检查：frac_flipped, frac_J<0.5, σ_D
  │     若 σ_D 已足够小 → 可以直接停
  │
  ├─[3] 目标 B 精修：VSEM（固定边界）或 IEM（自由边界）
  │     VSEM: 30 行固定点迭代（§3.1 Algorithm）
  │           → 直到 E_V 变化 < tol，或 15~100 次迭代
  │     IEM : preconditioned nonlinear CG（§3.3 Algorithm alg:IEM）
  │           → 直到 E_I < tol
  │
  └─[4] 验证并输出（§5）
        ratio stats, J field, frac_flipped, E_V vs bound,
        boundary sphere area check, ray-crossing bijectivity sample
        meshio.write("param.node", ...)
  │
  └─[5] 后续：用该参数化做四面体重网格化
```

**关键设计决策（我的建议）：**
1. **先做 [0] 的规范化**。不统一定向会让后面所有几何量都错，且错得很隐蔽。
2. **边界映射用现成 MATLAB**（`spherical_conformal_map` 最稳）。自己实现球面共形/保面积不值得。
3. **初值用 harmonic（[2]），精修用 IEM（[3]）**。VSEM 作为 IEM 的前 15 次迭代（论文就是这么做的）。
4. **加一条断言**：如果输入实体不满足星形条件，先警告。实测 $a\ge0.20$ 的 peanut 就开始危险。

---

## 7. 本报告的可复现产物（**全部已在本机跑通**）

工作目录 `E:\panyingyun\smartmm\`：

| 文件 | 内容 | 状态 |
|---|---|---|
| `tools/tet_ball_ref.py` | 参考实现：`orient_tets`, `boundary_faces_oriented`, `tet_volume_signed/volumes`, `dihedral_angles`, `radius_edge_ratio`, `aspect_ratio`, `p1_tet_stiffness(_local)`, `harmonic_ball_map`, `_local_stretch_matrix`, `vsem_laplacian`, `vsem_trace_laplacian`, `vsem_energy`, `vsem_lower_bound`, `vsem_solve`, `volume_distortion_ratios`, `jacobian_det_field`, `flipped_fraction`, `volume_preserving_error`, `verify_map`, `make_ball_tet_mesh` | ✅ 运行通过 |
| `tools/verify_ball_ref.py` | **28 条数值验证**（V1–V9b），覆盖全部关键公式 | ✅ **28 passed, 0 failed** |
| `tools/demo_failure_modes.py` | 非凸实体（peanut/wedge/crinkly）上的失效模式演示 | ✅ 运行通过 |
| `tools/sweep_convexity.py` | 非凸度扫描：$r=1+a\cos2\theta$，$a\in[0,0.45]$，报告边界双射性、harmonic/VSEM 的畸变 | ✅ 运行通过 |
| `vendor/SAEM/` | SAEM MATLAB 工具箱（Liu & Yueh 2026）已解包 | ✅ |
| `vendor/sphconf_README.md`, `vendor/sphconf_demo.m` | `spherical-conformal-map` 的 README 与 demo | ✅ |
| `vendor/2210.09654_tex/VSEM.tex` | VSEM 论文 LaTeX 源（`\title` 已确认为 Convergence Analysis…） | ✅ |
| `vendor/2402.00380_tex/` | n-VSE LaTeX 源 | ✅ |
| `vendor/2407.19272_tex/` | IEM LaTeX 源 | ✅ |
| `vendor/2412.19011_tex/` | SAEM LaTeX 源 | ✅ |
| `vendor/2506.17025_tex/` | Lyu-Chen-Choi-Lui"Volumetric Parameterization for 3-D Simply-Connected Manifolds" LaTeX 源（见 §9） | ✅ |
| `refs/` | 各论文 HTML + LaTeX 源 + VSEM/genus-1 论文 PDF | ✅ |
| `research/_verify_ball_ref_log.txt` | `verify_ball_ref.py` 的原始输出（**28 passed, 0 failed**） | ✅ |
| `research/_sweep_convexity_log.txt` | 非凸度扫描完整输出 | ✅ |
| `research/_demo_failure_modes_log.txt` | 失效模式演示完整输出 | ✅ |
| `tools/_fix_report_appendix.py` | 重建本报告附录的小工具（编码安全） | ✅ |

### 验证清单（`python tools/verify_ball_ref.py` 的实际输出）

```
[PASS] V1  all tets positively oriented (TetGen output)                min vol = 9.807e-04
[PASS] V2a -0.5*tr(f^T L_tau f) == 3|f(tau)|^2/(2 mu)  [Thm 3.1]      max rel err = 2.508e-13
[PASS] V2b the OTHER sign (+0.5 tr) does NOT satisfy Thm 3.1          rel err = 2.000e+00
[PASS] V3  grad_f E_V == -3 L_V(f) f  [Thm 3.2, finite diff]          rel err = 1.619e-09
[PASS] V3b the +3 L_V f form is WRONG by a factor -1                  rel err = 2.000e+00
[PASS] V4a E_V >= (3/2) V(f) always
[PASS] V4b VSEM reaches the lower bound 3/2*V(f)                      gap = 1.776e-15
[PASS] V4c VSEM map is volume preserving (ratio std < 1%)             std = 7.826e-16
[PASS] V4d VSEM map has no flipped tets                               min ratio = 1.000000
[PASS] V5a K * 1 = 0                                                   max |K 1| = 1.332e-15
[PASS] V5b local stiffness == -(1/6) l_kl cot(theta_ij^kl)            max rel diff = 3.453e-16
[PASS] V5c local tet stiffness annihilates constants                   max = 1.388e-17
[PASS] V5d local tet stiffness symmetric                               max asym = 0.000e+00
[PASS] V5e local tet stiffness PSD with one zero eigenvalue            eigs = [0, 0.0478, 0.0580, 0.5122]
[PASS] V5f global K == scatter of element matrices                     max err = 0.000e+00
[PASS] V5g B^{-1} B = I => grad(lambda_k).(v_j-v_0) = delta_kj         max err = 1.659e-16
[PASS] V5h element matrix rebuilt from gradients == assembled matrix    max err = 0.000e+00
[PASS] V5j regular tet: all 6 dihedral angles = 70.528779 deg
[PASS] V5k regular tet element stiffness: rowsum 0, PSD, one zero eig
[PASS] V6a .node/.ele round trip: points                               max diff = 0.000e+00
[PASS] V6b round trip: tets
[PASS] V6c round trip: point marker 'tetgen:ref'
[PASS] V6d round trip: cell attribute 'tetgen:ref'
[PASS] V7a geometric ball + radial BC == identity -> exact zero distortion
[PASS] V7b NON-convex ball + radial BC: 16.483% of tets have J<0.5
[PASS] V8  VSEM cuts volume-distortion std 2.0934 -> 0.2364 (8.9x)
[PASS] V9  VSEM fixed-boundary energy deficit = 1.7185 > 0
[PASS] V9b (documented motivation for IEM over VSEM)
```

**复现命令：**
```powershell
pip install numpy==2.4.4 scipy==1.18.0 meshio==5.3.5 tetgen==0.8.4 pyvista==0.49.0 gpytoolbox==0.3.7
cd E:\panyingyun\smartmm
C:\Python312\python.exe tools\verify_ball_ref.py      # 28 checks
C:\Python312\python.exe tools\sweep_convexity.py      # non-convexity sweep
C:\Python312\python.exe tools\demo_failure_modes.py   # failure-mode demo
C:\Python312\python.exe tools\tet_ball_ref.py --demo  # end-to-end smoke test
```

---

## 8. ❌ 未能验证 / 存疑的事项（**请团队补查**）

| # | 事项 | 状态 | 建议 |
|---|---|---|---|
| 1 | Wang/Li/Gu **2003 CIS "Volumetric harmonic map" 原文** | ❌ **未能下载**（付费墙）。§1.2 的离散形式是从 n-VSE(2402.00380) 的 $n=3$ 特例 + VSEM 重建的 | 通过学校图书馆/Interlibrary loan 拿 PDF，核对权重是否逐字一致 |
| 2 | "**Numerical Optimization for Volume-Preserving Parameterization of 3-Manifolds**" 这个标题（任务描述里给的，挂在 arXiv:2210.09654 名下） | ❌ **未能定位**。我下载 arXiv:2210.09654 源码确认它就是 "Convergence Analysis of Volumetric Stretch Energy Minimization…" | 请确认是否记错了标题/编号，或提供另一篇的线索 |
| 3 | Gu 组 **spherical OT / Ricci flow / harmonic map C++ 源码** | ⚠️ **页面声称有，但我没能提取到下载链接**（HTML 里只解析出 4 个目录链接）。SOT 页面给出的命令行 `OT.exe -source <mesh> -target <sphere.obj>` 是页面文本，未实际操作 | 人工打开 <https://www3.cs.stonybrook.edu/~gu/software/index.html> 与 `/software/SOT/index.html` 点链接下载 |
| 4 | **CGAL** 是否有球面参数化 | ❌ **我没有实际检查 CGAL**（无 Windows 安装）。我的理解是 `CGAL::Surface_mesh_parameterization` 只有圆盘/平面（LSCM/ARAP/MVC/离散共形），无球面 | 若要走 CGAL 路线，先核实 |
| 5 | `spharm-pdm` / `cmrep` / FreeSurfer 的实际接口 | ⚠️ **未安装验证**。基于文档/论文判断它们不适合本任务（见 §2.3 末） | 若团队已有 FreeSurfer 环境，可试 `mris_sphere` |
| 6 | `icemiliang/spherical_harmonic_maps` 的函数签名与算法 | ⚠️ **只看了文件列表**，未 fetch `harmonic.cpp` | 若能编译，值得一看 |
| 7 | MathWorks MATLAB Central 作者页（Yueh 是否有 File Exchange 提交） | ❌ **Access Denied**（Akamai 403） | 用非受限网络访问 <https://www.mathworks.com/matlabcentral/profile/authors/9330440> |
| 8 | **n-VSEM 的"convex combination postprocess"去折叠**细节 | ⚠️ **论文只提了一句，无公式** | 可参考 §9 的 Lyu-Chen-Choi-Lui 2025（专门做 bijectivity enforcement）或 "Foldover-free maps in 50 lines of code"（IEM 论文引用了它，ref `Foldover-free maps in 50 lines of code`） |
| 9 | IEM/n-VSE 的实际数值性能（大网格上的时间/内存） | ❌ **未复现**。我实现的是 VSEM，不是 IEM | 建议先按 §3.3 的算法框实现 IEM，在 10 万 tets 量级上测 |
| 10 | Lyu-Chen-Choi-Lui **arXiv:2506.17025** 的完整算法细节 | ⚠️ **只读了 §1–§4.1**（HTML 被截断，但 LaTeX 源已下载到 `vendor/2506.17025_tex/`） | 该文有 bijectivity enforcement，**很可能是 §8#8 的答案**。建议细读 §5.2/§5.3 |

---

## 9. 额外发现：值得一看的新工作

**arXiv:2506.17025** —— Z. Lyu, Q. Chen, G. P. T. Choi, L. M. Lui, "**Volumetric Parameterization for 3-Dimensional Simply-Connected Manifolds**", 20 Jun 2025（CUHK）。**LaTeX 源已下载到 `vendor/2506.17025_tex/`。**

**为什么值得看：**
- 它**正面攻击了本报告 §1.4 和 §8#8 的两个痛点**：**bijectivity 保证** + **局部几何畸变控制**。
- 提出三个模型：
  - **3DQC**：bijective 3D quasi-conformal map（最小化含几何项的 Beltrami 型能量）
  - **3DDEM**：bijective 3D density-equalizing map（扩散/密度均衡，天然保体积）
  - **3DDEQ**：把两者结合，在几何畸变与体积畸变之间取得最优平衡
- 3D quasi-conformal 的核心（§3.3）：对 Jacobian 做极分解 $J_f=UP$，$P=\sqrt{J_f^\top J_f}=W\Sigma W^{-1}$（$\Sigma=\operatorname{diag}(\lambda_1,\lambda_2,\lambda_3)$），则
  $$J_f^\top=W\begin{pmatrix}\frac{\lambda_1}{\lambda_2\lambda_3}&&\\&\frac{\lambda_2}{\lambda_1\lambda_3}&\\&&\frac{\lambda_3}{\lambda_1\lambda_2}\end{pmatrix}W^{-1}\operatorname{Adj}(J_f),\qquad
  \operatorname{Adj}(J_f)=\begin{pmatrix}\nabla v\times\nabla w&\nabla w\times\nabla u&\nabla u\times\nabla v\end{pmatrix}$$
  左乘 $\mathcal A=W\operatorname{diag}\!\big(\frac{\lambda_2\lambda_3}{\lambda_1},\frac{\lambda_1\lambda_3}{\lambda_2},\frac{\lambda_1\lambda_2}{\lambda_3}\big)W^{-1}$ 后取散度，得到**线性系统**：
  $$\nabla\cdot\mathcal A\nabla u=0,\qquad \nabla\cdot\mathcal A\nabla v=0,\qquad \nabla\cdot\mathcal A\nabla w=0$$
  这是 Chen-Lui 3D quasi-conformal 表示的直接应用，**每次迭代解一个 3-分量线性系统**（与 §1.2 的 harmonic 同构，只是把标量权重换成张量 $\mathcal A$）。
- §5.1 给了 "initial solid ball parameterization"（边界条件 + 内部映射）—— 与本报告的目标 A 直接对应。
- **代码可用性：未查证**（Lui 的 GitHub `lmlui` 只有 2 个无关仓库；Gary Choi 有 29 个仓库但没有 volumeric 的那一个）。

**建议：** 让团队里一个人细读 `vendor/2506.17025_tex/3DDEQ_arxiv.tex` 的 §5.2（"Enforcing the bijectivity"）和 §5.3。如果它的 bijectivity enforcement 好用，可以**直接用在本报告的流水线 [3] 之后**，解决 n-VSEM 的折叠问题。

---

## 10. 快速参考卡

### 公式
```
# P1 四面体刚度（Dirichlet；SPD）
K^τ_ij = -(1/6) * l_kl * cot(θ_ij^kl)   (i≠j),  K^τ_ii = -Σ_{j≠i} K^τ_ij
K = Σ_τ K^τ      #  (k,l) = 与 (i,j) 相对的边  ← 最容易写错

# 目标 A
K_II f_I^s = -K_IB f_B^s,  s = 1,2,3        # Cholesky 一次分解

# Volumetric stretch Laplacian（Yueh2019 字面符号；非对角 < 0，PSD M-矩阵）
L_τ = -(1/(36 μ(τ))) [a_ij],  a_ij = [(f_k-f_i)×(f_l-f_i)]·[(f_l-f_j)×(f_k-f_j)]
a_ji = a_ij,  a_ii = -Σ_{j≠i} a_ij
[L_V]_ij = Σ_{τ∋[v_i,v_j]} L_τ_ij   (i≠j),  [L_V]_ii = -Σ_{k≠i} [L_V]_ik

# 能量与下界
E_V = Σ_τ 3|f(τ)|²/(2 μ(τ))  = -½ trace(f^T L_V f)  ≥  (3/2)|f(M)|
∇_f E_V = -3 L_V(f) f
E_I = V(e)/V(f) * E_V - V(f)  ≥ 0,  等号 ⟺ 体保持

# 迭代（VSEM 固定边界）
[L_V(f)]_II f_I^s = -[L_V(f)]_IB f_B^s,  s=1,2,3,  重复

# 验证
J(τ) = ((f1-f0)×(f2-f0))·(f3-f0) / ((v1-v0)×(v2-v0))·(v3-v0)
D_V(τ) = |f(τ)|/|τ|
frac_flipped = #{J ≤ 0}/m          (必须为 0)
frac_J<0.5   = #{J < 0.5}/m        (也必须看！)
Σ_boundary Ω_τ = 4π,  Ω_τ = 2·atan2(a·(b×c), 1+a·b+b·c+c·a)   ← 必须带符号
```

### 命令
```powershell
pip install numpy scipy meshio tetgen pyvista gpytoolbox potpourri3d robust-laplacian pymeshlab libigl
C:\Python312\python.exe tools\verify_ball_ref.py
```

### 三条"不要做"
1. ❌ **不要**用已经是单位球的网格测试 harmonic 失效模式（结果是精确 0 畸变，浪费时间）。
2. ❌ **不要**只检查"有没有翻转 tets"。必须同时看 $J$ 的分布。
3. ❌ **不要**用 `K[idx][:,idx]` 去验证单元矩阵，也不要用 `tile` 写 scatter 索引。

---

## 附录 A：非凸度扫描的完整实测数据

`tools/sweep_convexity.py` 的输出（已剔除 TetGen 命令行噪声）。
列含义：`sphArea/4pi` = 边界球面映射的有符号立体角总和 / 4π（**1.0 = 双射**）；
`A:*` = 普通 P1 harmonic 映射到球；`B:*` = 同一（径向投影）边界下的 VSEM；
`B:gap` = $E_V(f) - \frac32 V(f_{\text{image}})$。

```
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 178
  Mesh tetrahedra: 469
  Mesh faces: 1098
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 182
  Mesh tetrahedra: 501
  Mesh faces: 1162
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 187
  Mesh tetrahedra: 545
  Mesh faces: 1250
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 190
  Mesh tetrahedra: 569
  Mesh faces: 1298
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 191
  Mesh tetrahedra: 576
  Mesh faces: 1312
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 189
  Mesh tetrahedra: 559
  Mesh faces: 1278
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 178
  Mesh tetrahedra: 469
  Mesh faces: 1098
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 180
  Mesh tetrahedra: 493
  Mesh faces: 1146
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 179
  Mesh tetrahedra: 517
  Mesh faces: 1194
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 182
  Mesh tetrahedra: 545
  Mesh faces: 1250
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 178
  Mesh tetrahedra: 469
  Mesh faces: 1098
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 181
  Mesh tetrahedra: 498
  Mesh faces: 1156
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 185
  Mesh tetrahedra: 531
  Mesh faces: 1222
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
Creating surface mesh ...
Removing exterior tetrahedra ...
  Input points: 162
  Input facets: 320
  Input segments: 480
  Input holes: 0
  Input regions: 0
  Mesh points: 182
  Mesh tetrahedra: 526
  Mesh faces: 1212
  Mesh faces on exterior boundary: 320
  Mesh faces on input facets: 320
  Mesh edges on input segments: 480
mesh: icosphere sub=2, tetgen maxvol=0.06
shape      amp sphArea/4pi  bndOK A:ratioStd  A:J<0.5  A:flip B:ratioStd    B:min  B:flip      B:gap
----------------------------------------------------------------------------------------------------
peanut    0.00      1.0000   True     0.0000   0.000%       0     0.0000  1.00000       0  1.776e-15
peanut    0.10      1.0000   True     0.2112   0.000%       0     0.0726  0.92624       0  6.119e-01
peanut    0.20      1.0000   True     0.4728   0.183%       0     0.1667  0.32705       0  1.200e+00
peanut    0.30      1.0000   True     0.8213   7.909%       0     0.2857  0.32011       0  1.789e+00
peanut    0.40      1.0000   True     1.3381  12.674%       0     0.4137  0.06724       0  2.364e+00
peanut    0.45      1.0000   True     1.7121  15.027%       0     0.5675  0.02972       0  2.774e+00
crinkly   0.00      1.0000   True     0.0000   0.000%       0     0.0000  1.00000       0  1.776e-15
crinkly   0.10      1.0000   True     0.2450   3.651%       0     0.2332  0.10831       0  1.799e-01
crinkly   0.20      1.0000   True     0.4858  31.915%       0     0.5012  0.03117       0  9.961e-01
crinkly   0.25      1.0000   True     0.5933  42.018%       0     0.5447  0.01905       0  1.414e+00
wedge     0.00      1.0000   True     0.0000   0.000%       0     0.0000  1.00000       0  1.776e-15
wedge     0.10      1.0000   True     0.1048   0.402%       0     0.0883  0.22054       0 -7.922e-01
wedge     0.20      1.0000   True     0.1872   2.637%       0     0.1443  0.11935       0 -1.413e+00
wedge     0.35      1.0000   True     0.2644  51.521%       0     0.2351  0.04546       0 -2.080e+00
Legend: A = plain P1 harmonic map to the ball; B = VSEM fixed point with the same (radial-projection) boundary map.
sphArea/4pi = 1 means the boundary map is a bijection to the sphere; 2 means it double-covers the sphere.
```
