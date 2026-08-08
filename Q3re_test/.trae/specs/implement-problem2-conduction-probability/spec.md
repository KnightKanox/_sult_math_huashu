# 问题 2：导通概率 Monte Carlo 仿真 Spec

## Why
Q1（版本 1.0.0）已确定性判定附件三组介质 A 是否导通。问题 2 要求：在随机填充介质 A 的前提下，给定体积分数 `φ∈{0.50%, 0.60%, 0.70%, 1.00%}`，求左右电极导通的概率 `P_conn`，并用概率云理论（`q_AA(r) → V_AA → k̄`）解释导通概率随填充率上升的机理。本次推进版本 **V2.0.0**。

## What Changes
- **新增生成层** `src/generation/`：随机位置 `U(Ω)`、球面均匀方向、5000nm 圆柱周期回绕切段，输出与附件一致的"在箱片段"格式
- **新增概率云层** `src/cloud/`：`q_AA(r)` Monte Carlo 估计、`V_AA = 4π∫r²q(r)dr`、平均连接度 `k̄ = ρ_A·V_AA`
- **新增仿真层** `src/simulation/`：单次微构体导通仿真 `single_trial`、置信区间 `confidence`（Wilson / Clopper-Pearson，95% 双侧 + 95% 单侧下界）
- **新增宽相位** `src/graph/spatial_index.py`：AABB 向量化候选对筛选（N=700 级下 O(N²) 判距不可行）
- **新增脚本** `scripts/build_aa_cloud.py`、`scripts/solve_problem2.py`
- **新增可视化** `src/visualization/probability_plot.py`（matplotlib，中文字体）
- **修改** `src/graph/connectivity.py`：PERIODIC 端点合并从 O(N²) 配对改为坐标哈希分组 O(N)，行为语义不变
- **新增输出** `results/problem2/`、`results/cloud/`、`results/figures/`
- `requirements.txt` 增加 `matplotlib`
- **版本推进**：新建 `version_log/2.0.0/`（复制 1.0.0 工程为基线，交付 Q2 代码、版本说明、问题2说明文档），更新工作区 `修改日志.md`

## Impact
- Affected specs: `implement-microstructure-conduction-simulation`（对应其 Task 3 概率云层、Task 4 Monte Carlo 层）
- Affected code: `src/graph/connectivity.py`（合并性能）、`src/config.py`（如需补充参数）、`scripts/`（新增两个脚本）
- 边界模式：**双模式对比**（`PERIODIC_CONNECTED` / `WRAPPED_GEOMETRY_ONLY`，用户已确认）
- Monte Carlo 次数：默认 **2000 trials/点**（用户已确认）
- 几何判据：延续 Q1 的 capsule 近似（轴距 ≤ 61.8nm；电极 ≤ 31.8nm）

## ADDED Requirements

### Requirement: 随机介质 A 生成器
系统 SHALL 提供固定 seed 的随机生成：中心 `c ~ U(Ω)`；方向球面均匀（`z=cosθ ~ U(-1,1)`、`φ ~ U(0,2π)`，禁止 `θ~U(0,π)`）；端点 `p=c-u·L_A/2`、`q=c+u·L_A/2`；按周期回绕规则把越界部分**解析切段**回绕为多段"在箱片段"，每片段 6 个坐标（与附件格式一致）。

#### Scenario: 切段正确性
- **WHEN** 生成一根沿 x 方向、跨过 `x=±5000` 的 5000nm 圆柱
- **THEN** 输出多段在箱片段，所有端点坐标均在 `[-5000,5000]³` 内；`PERIODIC_CONNECTED` 模式能将回绕重合的片段合并回同一节点

### Requirement: 概率云 q_AA(r)
系统 SHALL 用 Monte Carlo 估计 `q_AA(r)`：固定中心距 r（相对位置球面均匀、两圆柱方向独立球面均匀），采样 M 次判断两圆柱是否导通（轴距 ≤ 61.8nm），输出 `q̂`、样本数、95% CI。

#### Scenario: 远距离衰减
- **WHEN** `r > L_A + 61.8 ≈ 5061.8nm`
- **THEN** `q̂_AA(r) = 0`

### Requirement: 等效连接体积与平均连接度
系统 SHALL 计算 `V_AA = 4π∫₀^∞ r²q(r)dr`（`np.trapezoid`）与 `k̄ = ρ_A·V_AA`（`ρ_A = N_A/V₀`）。

### Requirement: 概率云框架解析校验（BB 基准）
系统 SHALL 用 B-B 解析解校验二体采样/积分框架：`q_BB(r) ≈ I(r ≤ 401.8)`，`V_BB ≈ 4/3π·401.8³`。

### Requirement: 单次 Monte Carlo 仿真
系统 SHALL 实现 `single_trial`：生成 N_A 根 A → 切段 → AABB 宽相位筛选候选对 → 精确判距建边 → DSU → 判定 L/R 导通；返回导通与否及统计量（节点数、边数、最大连通分量比例、平均度等）。

#### Scenario: 确定性基准
- **WHEN** `N_A = 0`
- **THEN** 一定不导通
- **WHEN** 构造一根横跨左右电极的圆柱
- **THEN** 一定导通

### Requirement: 置信区间
系统 SHALL 实现 Wilson 与 Clopper-Pearson 的 95% 双侧区间及 95% 单侧下界（单侧下界为 Q3 复用准备）。

### Requirement: solve_problem2 主脚本
系统 SHALL 提供 `solve_problem2.py`：`φ∈{0.50%,0.60%,0.70%,1.00%}` × 双边界模式 × trials（默认 2000，可配）→ 输出 CSV（`phi, N_A, boundary_mode, trials, connected_count, p_hat, ci_low, ci_high, rho, k̄`）与汇总 JSON（含 seed/config）。

#### Scenario: 可复现
- **WHEN** 相同 `--seed` 运行两次
- **THEN** 输出 CSV 完全一致

### Requirement: 可视化
系统 SHALL 绘制 `P_conn(φ)`（双模式两条曲线，含 CI 误差棒）与 `k̄(φ)`，中文字体正常显示。

### Requirement: 问题2说明文档与版本归档
系统 SHALL 输出 `问题2_说明文档.md`（建模假设、参数、结果表、与概率云理论的对照、结论），更新工作区 `修改日志.md`，新建 `version_log/2.0.0/版本说明.md`。

## MODIFIED Requirements

### Requirement: PERIODIC 合并性能
原 `_build_nodes` 的 PERIODIC 模式用双重循环配对"回绕重合端点"为 O(N²)；Q2 单次仿真有 707 根圆柱 × 多片段，每次 trial 均调用，O(N²) 不可行。改为按**回绕后端点坐标哈希分组**（O(N)），合并语义与结果不变。

## REMOVED Requirements
无
