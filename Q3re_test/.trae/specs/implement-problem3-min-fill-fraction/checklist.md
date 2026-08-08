# Checklist — 问题 3：P_conn≥90% 最低 A 填充率

- [x] `scripts/solve_problem3.py` 已生成，`py_compile` 语法检查通过（主流程未运行）
- [x] 二分搜索核心正确：中点 MC → `p_lower,95% ≥ 0.90` 判定收缩，`tol=0.0001`（0.01% 精度）、最大迭代 40 保护
- [x] 判定采用 95% Wilson **单侧下界**（`p_lower,95% ≥ 0.90`），而非简单 `p̂ ≥ 0.90`
- [x] 最终确认：φ_min 上下邻点以 `confirm_trials=4000` 复核
- [x] 理论对照：由 `cloud_AA.csv` 积分 V_AA，输出 `φ_theory = V_A/V_AA`（k̄=1 → ≈0.5546%）与 `k̄(φ_min)`
- [x] 复用 Q2 基础设施（导入而非复制）：`run_single_trial`/`run_trials`、`wilson_one_sided_lower`、`n_a_from_phi`、`load_v_aa`
- [x] CLI 参数齐全：`--phi-low/--phi-high/--tol/--trials/--confirm-trials/--seed/--mode/--workers/--out-dir/--plot`
- [x] 输出 `problem3_result.csv` + `problem3_result.json`（含二分迭代记录）与可视化图（输出代码已实现；实际文件待运行主流程后生成）
- [x] `tests/test_problem3.py` 全部通过（mock 二分收敛、判定边界、换算、迭代保护；不跑真实仿真）
- [x] `问题3_说明文档.md`（含运行命令与结果表占位）、`版本说明.md`、`修改日志.md` 已更新
- [x] 主流程 Monte Carlo 未执行（符合本阶段"只生成脚本不运行"要求）
