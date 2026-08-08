# 问题 4 Spec：A/B 混合体系 P_conn≥90% 的最低成本组合

## Why

问题 4 在问题 3 基础上引入第二种导电介质 B（半径 200nm 的球，成本 0.05 元/μm³，远低于 A 的 1.05 元/μm³），要求寻找在 `P_conn ≥ 0.90`（95% Wilson 单侧置信下界判据）约束下使总成本 `C(N_A,N_B)` 最小的整数组合。这是把概率云理论从单一介质推广到双介质连接矩阵（λ_max 指标）的关键一步，也是 A 题优化环节的收官。

## What Changes

- **生成层**：新增 B 球生成 `generate_b_spheres(rng, n_b)`——球心均匀落在 `[-(5000-R_B), 5000-R_B]³`（球完全在盒内，不做切段/回绕，用户已确认）。
- **几何判据**：新增 A-B（点-线段距离 ≤ R_A+R_B+δ=231.8）、B-B（球心距 ≤ 2R_B+δ=401.8）、B-电极（球面到平面 ≤ δ）；A-A、A-电极沿用现有实现。暴露 `point_segment_distance`。
- **宽相位**：`aabb_candidates` 扩展支持线段+球混合输入（球 AABB=[c±R_B]），A-B/B-B 对同样走 AABB 候选+精确判距。
- **混合单次仿真**：`run_single_trial_mixed(rng, n_a, n_b, mode)`——A 圆柱片段（沿用 `generate_batch`）+ B 球（独立节点，不参与 PERIODIC 端点合并），统一建图、DSU、左右电极贯通判定；保留双边界模式。
- **概率云扩展**：新增 `estimate_q_ab(r)`（A 随机圆柱 + B 球球壳采样）、`q_bb(r)` 解析基准 `I(r≤401.8)`、`V_AB=4π∫r²q_AB dr`、`V_BB=4/3π(401.8)³`；A-B 用 capsule 近似 `V_AB^approx=πR_AB²L_A+4/3πR_AB³`（R_AB=231.8）做 sanity check。
- **连接矩阵与理论筛选**：`M=[[ρ_A V_AA, ρ_B V_AB],[ρ_A V_AB, ρ_B V_BB]]`，`λ_max(M)` 作为理论筛选指标（非精确判据），快速缩小二维搜索范围。
- **成本函数**：`C(N_A,N_B)=c_A·N_A·V_A(μm³)+c_B·N_B·V_B(μm³)`（1μm³=10⁹nm³）。
- **主脚本 `solve_problem4.py`** 三级流程：
  1. Level 1 理论筛选：λ_max 扫描 N_A×N_B，标定理论可行区；
  2. Level 2 MC 边界：固定 N_A，二分最小 N_B 使 95% 单侧置信下界 ≥ 0.90，得到可行边界 N_B_min(N_A)；
  3. Level 3 成本优化：沿边界比较成本取全局最小候选，对最优点及邻域高次数（4000）MC 复核确认。
- **可视化**：λ_max 热图、MC P_conn 热图、成本等高线 + P=0.9 可行边界（规划 §49 Figure 6/7/8）。
- **主流程运行**：实现并通过单元测试后，立即运行 Q4 主流程（用户已确认）。
- **文档**：`问题4_说明文档.md`、`版本说明.md`、根 `修改日志.md` 补录。

## Impact

- Affected specs：`implement-problem3-min-fill-fraction`（沿用二分+Wilson 判据+可复现种子机制，从一维 φ 推广到二维 (N_A,N_B)）。
- Affected code：
  - 新增：`src/cloud/mixed_cloud.py`（q_AB/q_BB/V_AB/V_BB）、`scripts/solve_problem4.py`、`tests/test_problem4.py`、`version_log/2.0.0/问题4_说明文档.md`
  - 修改：`src/generation/medium_generator.py`（球生成）、`src/geometry/cylinder_distance.py`（A-B/B-B/B-电极判据）、`src/geometry/segment_distance.py`（暴露点-线段距离）、`src/graph/spatial_index.py`（混合宽相位）、`src/simulation/single_trial.py`（混合单 trial）、`src/config.py`（B 相关阈值，仅补常量不改现有）、`version_log/2.0.0/版本说明.md`、根 `修改日志.md`

## ADDED Requirements

### Requirement: B 球生成
系统 SHALL 提供 `generate_b_spheres(rng, n_b)`，返回 (n_b,3) 球心数组；球心各分量在 `[-BOX_HALF+B_RADIUS, BOX_HALF-B_RADIUS]` 均匀采样（球完全在盒内）。

#### Scenario: 球完全在盒内
- **WHEN** 调用 `generate_b_spheres(rng, 100)`
- **THEN** 返回 100 个球心，全部满足 `|x|,|y|,|z| ≤ 4800`，且单元测试断言范围成立

### Requirement: 混合导通判据
系统 SHALL 提供 A-B、B-B、B-电极三类判据：
- A-B：轴线段到球心距离 ≤ R_A+R_B+δ（231.8）；
- B-B：球心距 ≤ 2R_B+δ（401.8）；
- B-电极：球面到 x=±5000 平面距离 ≤ δ，即 `|5000-|x_c|| ≤ R_B+δ`（201.8）。

#### Scenario: B-B 临界距离
- **WHEN** 两球心距 = 401.8（含），**THEN** 判定导通
- **WHEN** 两球心距 = 401.81，**THEN** 判定不导通

### Requirement: 混合宽相位与混合单次仿真
系统 SHALL 扩展候选对筛选以同时处理 A 圆柱线段与 B 球，并提供 `run_single_trial_mixed(rng, n_a, n_b, mode)`：B 球为独立节点（不参与 PERIODIC 端点合并），A 沿用现有切段/合并语义；返回与 Q2/Q3 一致的统计 dict（connected/node_count/edge_count/max_component_ratio/mean_degree 等）。双边界模式均支持。

#### Scenario: 纯 B 体系
- **WHEN** 调用 `run_single_trial_mixed(rng, 0, n_b, mode)`（n_b>0）
- **THEN** 正确判定纯 B 球网络的左右贯通，且不因缺少 A 而崩溃

### Requirement: 双介质概率云
系统 SHALL 提供 `estimate_q_ab(r, ...)`（A 圆柱随机姿态 + B 球球壳采样）与解析 `q_bb(r)=I(r≤401.8)`；计算 `V_AB=4π∫r²q_AB dr`、`V_BB=4/3π·401.8³`；`V_AB` 与 capsule 近似 `π·231.8²·5000+4/3π·231.8³` 相对偏差 ≤ 15%。

### Requirement: 连接矩阵与 λ_max 理论筛选
系统 SHALL 提供 `connection_matrix(n_a, n_b, v_aa, v_ab, v_bb)` 返回对称 2×2 矩阵（元素为 ρ·V，ρ=N/V₀），及 `lambda_max(m)`（np.linalg.eigvalsh）。λ_max 仅用于理论筛选，不作最终判据。

### Requirement: 成本函数
系统 SHALL 提供 `total_cost(n_a, n_b)`，按 `c_A·N_A·V_A(μm³)+c_B·N_B·V_B(μm³)` 计算，单位换算 1μm³=10⁹nm³；V_A=π·30²·5000≈1.4137e7 nm³，V_B=4/3π·200³≈3.3510e7 nm³。

#### Scenario: 成本单调性
- **WHEN** 固定 N_B，增大 N_A，**THEN** 成本严格增大（两者均为正成本）

### Requirement: 三级搜索主流程
系统 SHALL 在 `solve_problem4.py` 实现：
1. Level 1：λ_max 快速扫描 N_A×N_B（理论可行区，λ_max≥1 启发式）；
2. Level 2：对 N_A 网格点（默认步长 10，范围 0~360 可 CLI 调），二分最小 N_B（默认 trials=500，区间 [0, N_B_max]），判据为 `wilson_one_sided_lower(connected, trials) ≥ 0.90`，得到可行边界；
3. Level 3：沿边界比较成本，取全局最小候选；对候选及邻域（N_A±步长、N_B±固定窗）用 confirm_trials=4000 复核，确认最终 `(N_A*, N_B*)` 与最小成本。
N_B 上界默认 5000（纯 B 理论 k̄=1 ≈ 3681 球的余量）。

#### Scenario: 边界单调
- **WHEN** 某 (N_A, N_B) 满足判据，**THEN** 任何 N_B'≥N_B 同 N_A 也应满足（单调性校验，二分前提）

### Requirement: 可复现性
系统 SHALL 沿用 Q3 种子机制：主 rng（seed=42，CLI 可调）一次性派生全部评估点子种子池，固定评估顺序（Level 2 逐 N_A 点、逐二分迭代），串/并行结果一致；最终 JSON/CSV 记录 seed、config、V_AA/V_AB/V_BB、λ_max 扫描与边界、成本、确认记录。

### Requirement: 输出与图
系统 SHALL 输出 `results/problem4/problem4_result.csv` / `.json`（含理论筛选、可行边界、成本比较、确认记录），并生成三张图：λ_max 热图、MC P_conn 热图、成本等高线 + P=0.9 可行边界（中文标签、Microsoft YaHei 字体、k̄ 类组合字符用 mathtext）。

## MODIFIED Requirements

### Requirement: 单次仿真入口（扩展）
现有 `run_single_trial(rng, n_a, mode)` 语义不变（纯 A）；新增混合入口 `run_single_trial_mixed`。`n_a=0` 或 `n_b=0` 时退化为对方单介质体系，不得崩溃。

### Requirement: 判定标准（沿用）
问题 4 沿用 Q3 判定标准：**95% Wilson 单侧置信下界 `p_lower,95% ≥ 0.90`**，不允许仅以 `p̂ ≥ 0.90` 作为最终可行判据（规划 §56 原则 8）。

### Requirement: 边界模式（沿用）
保留 `PERIODIC_CONNECTED` / `WRAPPED_GEOMETRY_ONLY` 双模式（用户已确认），Q4 对两模式各出一份边界与最优解。
