# 问题 3：P_conn≥90% 最低 A 填充率（二分搜索）Spec

## Why

问题 2 已给出 4 个固定填充率下的导通概率，问题 3 是反向优化问题：求满足 `P_conn ≥ 0.90` 的最低介质 A 体积分数 φ_min。判定必须用 **95% 单侧置信下界 ≥ 0.90**（而非简单 `p̂ ≥ 0.90`），以保证统计保证。同时利用概率云理论的 `k̄ = ρ_A·V_AA` 给出理论预测（k̄=1 → φ≈0.5546%），与二分搜索结果对照，完成理论-仿真的定量互证。

## What Changes

- 新增 `scripts/solve_problem3.py`：二分搜索主脚本（本阶段**只生成代码，不执行**主流程）
  - 复用 Q2 基础设施：`run_single_trial`（仿真）、`wilson_one_sided_lower`（95% 单侧下界）、`n_a_from_phi`/`num_density`/`mean_degree`/`load_v_aa`（换算与理论量）
  - 二分搜索：区间 `[phi_low, phi_high]`（默认 0.20%–0.70%），中点跑 `trials` 次 MC，按 `p_lower,95% ≥ 0.90` 判定收缩，至 φ 绝对精度 `tol=0.0001`（0.01%），最大迭代保护 40 次
  - 最终确认：对收敛后 φ_min 及其上下邻点加大仿真量（`confirm_trials=4000`）复核
  - 双边界模式（默认 `--mode all`，延续 Q2 敏感性对比）
  - 理论对照：`φ_theory(k̄=1) = V_A / V_AA` 随结果一并输出
  - 输出：`results/problem3/problem3_result.csv`（φ_min、N_A、模式、p_hat、单侧下界、k̄、理论预测）+ `problem3_result.json`（seed/config/全部二分迭代记录）+ 可视化图（P_conn(φ) 曲线 + φ_min 竖线标注）
- 新增 `tests/test_problem3.py`：二分逻辑用 mock MC 验证收敛与判定边界（不跑真实仿真）
- 新增文档：`version_log/2.0.0/问题3_说明文档.md`（含运行命令与结果表占位，实验后填写）、`version_log/2.0.0/版本说明.md` 更新、工作区 `修改日志.md` 更新
- **不修改**任何既有源码模块（Q2 全部复用）

## Impact

- Affected specs: `implement-microstructure-conduction-simulation`（Q1–Q4 总体规划中的 Q3 部分）
- Affected code: `version_log/2.0.0/conductive_microstructure/`（新增 `scripts/solve_problem3.py`、`tests/test_problem3.py`、`results/problem3/`、文档）
- 复用不复制：从 `scripts/solve_problem2.py` **导入** `n_a_from_phi`、`load_v_aa`、`run_trials` 相关辅助；从 `src.simulation.confidence` 导入 `wilson_one_sided_lower`

## ADDED Requirements

### Requirement: 二分搜索主脚本

系统 SHALL 提供 `scripts/solve_problem3.py`，对双边界模式在 `[phi_low, phi_high]` 内二分搜索满足 `P_conn ≥ 0.90`（判定：95% Wilson 单侧下界 `p_lower,95% ≥ 0.90`）的最低 φ。

#### Scenario: 成功搜索
- **WHEN** 运行 `python scripts\solve_problem3.py`（默认参数）
- **THEN** 输出 φ_min（绝对精度 ≤0.0001）、对应 N_A、`p_hat`、`p_lower,95%`、`k̄`、理论预测 φ_theory；写 CSV/JSON/图；收敛后对 φ_min 邻点做 confirm_trials 复核

#### Scenario: 判定与收缩
- **WHEN** 中点 φ_mid 的 `p_lower,95% ≥ 0.90`
- **THEN** `phi_high = phi_mid`；否则 `phi_low = phi_mid`

#### Scenario: 理论对照
- **WHEN** 脚本加载 `results/cloud/cloud_AA.csv` 积分得 V_AA
- **THEN** 计算并输出 `φ_theory = V_A / V_AA`（k̄=1 的临界填充率理论值）用于与二分实测对比

#### Scenario: 本阶段不运行
- **WHEN** 实施本 spec
- **THEN** 只生成 `solve_problem3.py` 代码与测试/文档，**不执行**主流程 Monte Carlo（运行命令写入说明文档，留待用户批准后执行）

### Requirement: 单元测试（mock，不跑真实仿真）

系统 SHALL 提供 `tests/test_problem3.py`，用可注入的假 MC 函数验证二分收敛方向、`tol` 终止条件、`p_lower,95% ≥ 0.90` 判定边界及 φ↔N_A 换算。

#### Scenario: mock 收敛
- **WHEN** 假 MC 给定单调真值函数（如 p̂ = sigmoid(k̄)）
- **THEN** 二分在有限迭代内收敛到真值附近，且随 tol 减小结果单调收紧

#### Scenario: 换算正确
- **WHEN** 给定 φ 与已知 V_A、V_0
- **THEN** N_A = round(φ·V_0/V_A)，φ_theory = V_A/V_AA 与手算一致

## MODIFIED Requirements

无（不修改既有需求）。

## REMOVED Requirements

无。
