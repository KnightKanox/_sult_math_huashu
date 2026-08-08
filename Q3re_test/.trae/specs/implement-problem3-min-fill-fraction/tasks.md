# Tasks — 问题 3：P_conn≥90% 最低 A 填充率（V2.0.0 追加）

> 主规范：`.trae/specs/implement-problem3-min-fill-fraction/spec.md`；理论依据：`MODELING_PLAN_概率云_微构体导电优化.md` §24–26。
> 实现原则：复用 Q2 基础设施（不复制代码）、95% 单侧下界判定（`p_lower,95% ≥ 0.90`）、二分至 0.01% 精度、双边界模式、固定 seed、每个函数前单行中文注释、**本阶段只生成脚本不运行主流程**、完成后更新 `修改日志.md`。

- [x] Task 1: `scripts/solve_problem3.py` 二分搜索主脚本
  - [x] SubTask 1.1: 从 `solve_problem2.py` 导入复用 `n_a_from_phi`、`load_v_aa`、`run_trial_batch`/`run_trials` 等；从 `src.simulation.confidence` 导入 `wilson_one_sided_lower`
  - [x] SubTask 1.2: 二分核心 `binary_search_min_phi(...)`：区间 [phi_low, phi_high]，中点 MC（默认 trials=1000）→ `p_lower,95% ≥ 0.90` 判定收缩，`tol=0.0001`、最大迭代 40；返回 φ_min 与全部迭代记录
  - [x] SubTask 1.3: 最终确认：对 φ_min 上下邻点用 `confirm_trials=4000` 复核单侧下界
  - [x] SubTask 1.4: 理论对照：由 cloud CSV 积分 V_AA，输出 `φ_theory = V_A/V_AA`（k̄=1 临界填充率）与 `k̄(φ_min)`
  - [x] SubTask 1.5: CLI 参数：`--phi-low 0.002 --phi-high 0.007 --tol 0.0001 --trials 1000 --confirm-trials 4000 --seed 42 --mode all --workers 8 --out-dir --plot/--no-plot`
  - [x] SubTask 1.6: 输出 `results/problem3/problem3_result.csv`（phi_min、N_A、boundary_mode、trials、p_hat、p_lower_95、kbar、phi_theory、phi_low_end、phi_high_end）+ `problem3_result.json`（seed/config/迭代记录/理论预测）+ 可视化图（P_conn(φ) 曲线 + φ_min 竖线，复用中文字体设置）
  - [x] SubTask 1.7: `python -m py_compile scripts\solve_problem3.py` 语法通过（不运行主流程）
- [x] Task 2: `tests/test_problem3.py` 单元测试（mock，不跑真实仿真）
  - [x] SubTask 2.1: 假 MC 单调真值函数下二分收敛（迭代次数 ≤ 上限、结果逼近真值）
  - [x] SubTask 2.2: `tol` 终止与最大迭代保护（不收敛时正确报错）
  - [x] SubTask 2.3: 判定边界：`p_lower,95% ≥ 0.90` 与 `p_hat` 换算正确（与 `wilson_one_sided_lower` 一致）
  - [x] SubTask 2.4: φ↔N_A 换算、`φ_theory = V_A/V_AA` 手算对照
  - [x] SubTask 2.5: 运行 `pytest tests\test_problem3.py` 全部通过（仅此测试，不跑主流程）
- [x] Task 3: 文档与版本收尾
  - [x] SubTask 3.1: 新建 `version_log/2.0.0/问题3_说明文档.md`：建模与方法（二分+单侧下界+理论对照）、CLI 运行命令、结果表占位（待运行后填写）
  - [x] SubTask 3.2: 更新 `version_log/2.0.0/版本说明.md`（追加 Q3 脚本就绪、实验待运行状态）
  - [x] SubTask 3.3: 更新工作区 `修改日志.md`（新增文件/修改记录）

# Task Dependencies
- [Task 2] 依赖 [Task 1]（导入脚本函数）
- [Task 3] 依赖 [Task 1]（脚本接口与参数已在代码中定型）

# 可并行任务
- 无（顺序执行）
