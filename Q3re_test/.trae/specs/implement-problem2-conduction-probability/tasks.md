# Tasks — 问题 2：导通概率 Monte Carlo 仿真（V2.0.0）

> 主规范：`.trae/specs/implement-problem2-conduction-probability/spec.md`；理论依据：`MODELING_PLAN_概率云_微构体导电优化.md` §22–23、§41–43、§50。
> 实现原则：正确优先、固定 seed、球面均匀方向、双边界模式、默认 2000 trials、95% CI、每个函数前单行中文注释、完成后更新 `修改日志.md`。

- [x] Task 1: V2.0.0 环境与基线
  - [x] SubTask 1.1: 新建 `version_log/2.0.0/`，复制 `version_log/1.0.0/conductive_microstructure/` 为基线
  - [x] SubTask 1.2: `requirements.txt` 增加 `matplotlib`，安装依赖并确认可导入
- [x] Task 2: 随机生成层（`src/generation/`）
  - [x] SubTask 2.1: `random_orientation.py`：球面均匀方向采样（z=cosθ~U(-1,1)）
  - [x] SubTask 2.2: `random_position.py`：盒内均匀位置 U(Ω)
  - [x] SubTask 2.3: `medium_generator.py`：圆柱生成 + 解析求越界切分时刻，把 5000nm 圆柱切成"在箱片段"（每片段 6 坐标，与附件格式一致）
  - [x] SubTask 2.4: `tests/test_generation.py`：方向矩检验（<u_x²>≈1/3 等）、切段端点均在盒内、跨边界圆柱切段数正确
- [x] Task 3: 概率云层（`src/cloud/`）
  - [x] SubTask 3.1: `aa_cloud.py`：q_AA(r) 估计（固定中心距，方向随机，M 次采样）
  - [x] SubTask 3.2: `effective_volume.py`：V_AA=4π∫r²q(r)dr（np.trapezoid）、k̄=ρ_A·V_AA
  - [x] SubTask 3.3: `scripts/build_aa_cloud.py`：r∈[0,5100]、默认 step 25、samples_per_r 默认 2000 → `results/cloud/cloud_AA.csv`（r_nm,q_hat,success_count,sample_count,ci_low,ci_high）
  - [x] SubTask 3.4: `tests/test_cloud.py`：q(0)=1、远距（>5061.8）=0、BB 解析基准（q_BB≈I(r≤401.8)、V_BB≈4/3π·401.8³）
- [x] Task 4: 仿真层（`src/simulation/` + 性能）
  - [x] SubTask 4.1: `src/graph/spatial_index.py`：AABB 向量化宽相位，返回候选片段对
  - [x] SubTask 4.2: 修改 `src/graph/connectivity.py`：PERIODIC 端点合并改为哈希分组 O(N)，并用 Q1 三组结果回归（结果不变）
  - [x] SubTask 4.3: `single_trial.py`：生成→切段→宽相位→精确判距→DSU→L/R 判定，返回导通与统计量
  - [x] SubTask 4.4: `confidence.py`：Wilson/Clopper-Pearson 95% 双侧 + 95% 单侧下界
  - [x] SubTask 4.5: `tests/test_simulation.py`：N_A=0 不导通、横跨圆柱导通、CI sanity、种子可复现
- [x] Task 5: solve_problem2 主脚本与可视化
  - [x] SubTask 5.1: `scripts/solve_problem2.py`：φ∈{0.50%,0.60%,0.70%,1.00%} × 双模式 × trials=2000 → `results/problem2/problem2_result.csv` + 汇总 JSON（seed/config）
  - [x] SubTask 5.2: `src/visualization/probability_plot.py`：P_conn(φ)（双模式+CI 误差棒）、k̄(φ)，中文字体
  - [x] SubTask 5.3: 运行全流程，生成 `results/problem2/`、`results/cloud/`、`results/figures/`
- [x] Task 6: 验证与结果整理
  - [x] SubTask 6.1: 收敛性实验：M∈{500,1000,2000} 对比 P̂ 与 CI 宽度（记录到说明文档）
  - [x] SubTask 6.2: 不同 seed（如 1–3）波动验证，确认可复现
  - [x] SubTask 6.3: 撰写 `version_log/2.0.0/问题2_说明文档.md`（建模假设、参数、结果表、概率云对照、结论）与 `version_log/2.0.0/版本说明.md`，更新工作区 `修改日志.md`

# Task Dependencies
- [Task 2] 依赖 [Task 1]
- [Task 3] 依赖 [Task 2]（方向生成）
- [Task 4] 依赖 [Task 2]
- [Task 5] 依赖 [Task 3][Task 4]
- [Task 6] 依赖 [Task 5]

# 可并行任务
- [Task 3] 与 [Task 4] 在 [Task 2] 完成后可并行开发
