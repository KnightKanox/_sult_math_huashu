# Tasks

> 主规范：`MODELING_PLAN_概率云_微构体导电优化.md`（§37 代码架构、§39–45 分阶段任务、§55 执行顺序）。
> 实现原则：正确优先、固定 seed、球面均匀方向、核心函数单元测试、CI 置信区间、Q3/Q4 用 95% 单侧置信下界、两种边界模式。每个函数前单行注释（中文），修改代码后同步更新工作区 `修改日志.md`。

- [ ] Task 1: 项目骨架与基础几何层
  - [ ] SubTask 1.1: 创建 `conductive_microstructure/` 目录骨架、`requirements.txt`、`src/config.py`（SimulationConfig、BoundaryMode 枚举、几何常数）
  - [ ] SubTask 1.2: 数据结构 `CylinderA`、`SphereB`（`src/geometry/primitives.py`）
  - [ ] SubTask 1.3: 随机位置生成（U(Ω)）、球面均匀方向采样（拒绝 θ~U(0,π) 方案，用 z=cosθ~U(-1,1)）（`src/generation/`）
  - [ ] SubTask 1.4: 线段-线段最短距离（`segment_distance.py`）
  - [ ] SubTask 1.5: A-A（capsule 近似，预留精确接口）、A-B、B-B、介质-电极距离判据（`cylinder_distance.py`、`sphere_distance.py`、`electrode_distance.py`）
  - [ ] SubTask 1.6: 周期边界回绕 + 两种边界模式（`periodic_boundary.py`）
  - [ ] SubTask 1.7: 单元测试（`tests/test_orientation.py`、`test_segment_distance.py`、`test_bb_distance.py`、`test_boundary.py`）并跑通 pytest
- [ ] Task 2: 图论层与问题 1（确定性连通判定）
  - [ ] SubTask 2.1: DSU（`src/graph/dsu.py`）与连通判定（`connectivity.py`）
  - [ ] SubTask 2.2: 从 `附件.xlsx` 三个分表读取端点数据（`组1`=13 根、`组2`=50 根、`组3`=536 根）
  - [ ] SubTask 2.3: `scripts/solve_problem1.py`：建图（虚拟电极 L/R）→ find(L)==find(R)，输出导通情况、总介质数、总边数、电极直连数、一条贯通路径
  - [ ] SubTask 2.4: 单元测试（`tests/test_dsu.py`）与 Q1 结果验证
- [ ] Task 3: 概率云层（理论模块）
  - [ ] SubTask 3.1: q_AA(r)、q_AB(r)、q_BB(r) 数值计算（`src/cloud/`：aa_cloud.py、ab_cloud.py、bb_cloud.py）
  - [ ] SubTask 3.2: 等效连接体积 V_eff=4π∫r²q(r)dr（np.trapz）与平均连接度 k̄（`effective_volume.py`）
  - [ ] SubTask 3.3: `scripts/build_aa_cloud.py`（参数 r_min/r_max/r_step/samples_per_r/seed，输出 results/cloud/cloud_AA.csv）
  - [ ] SubTask 3.4: BB 解析基准校验测试（q_BB≈I(r≤401.8)、V_BB≈4/3π·401.8³）（`tests/test_cloud_bb.py`）
- [ ] Task 4: Monte Carlo 仿真层与问题 2
  - [ ] SubTask 4.1: 单次仿真 `single_trial.py`（生成介质→边界处理→broad+narrow 判距→DSU→L/R 连通）
  - [ ] SubTask 4.2: 置信区间 Wilson/Clopper-Pearson（`confidence.py`），含 95% 单侧下界
  - [ ] SubTask 4.3: `scripts/solve_problem2.py`：φ∈{0.50%,0.60%,0.70%,1.00%} → CSV（phi,N_A,p_hat,ci_low,ci_high,rho,k̄）
  - [ ] SubTask 4.4: 绘制 P_conn(φ)、k̄(φ) 图（`visualization/probability_plot.py`）
  - [ ] SubTask 4.5: 单调性测试（`tests/test_monotonicity.py`）与不同 seed 波动验证（可复现性）
- [ ] Task 5: 问题 3（90% 最低 A 填充率）
  - [ ] SubTask 5.1: 单调性二分搜索（`src/optimization/problem3_search.py`），判据 95% 单侧下界 ≥ 0.90，精度 0.01%
  - [ ] SubTask 5.2: `scripts/solve_problem3.py`：输出最低 φ_A、N_A、p_hat、p_lower、trials，最终候选高次数复验
- [ ] Task 6: 问题 4（A/B 混合最低成本）
  - [ ] SubTask 6.1: 双类型连接矩阵 M 与 λmax 理论筛选（`src/cloud/effective_volume.py` 扩展 + `src/optimization/cost.py`）
  - [ ] SubTask 6.2: 固定 N_A 二分最小可行 N_B（Monte Carlo + 单侧下界），得 90% 可行边界（`problem4_search.py`）
  - [ ] SubTask 6.3: 成本比较取全局最低，最优点邻域高次数复验（`scripts/solve_problem4.py`）
  - [ ] SubTask 6.4: λmax 热图、P_conn 热图、成本等高线+可行边界图（`visualization/phase_plot.py`）
- [ ] Task 7: 性能优化（在正确 reference 基础上）
  - [ ] SubTask 7.1: 空间哈希/cell list broad phase（`src/graph/spatial_index.py`），候选对才做精确几何判断
  - [ ] SubTask 7.2: Monte Carlo 并行（concurrent.futures/multiprocessing），保持 seed 可复现
- [ ] Task 8: 全流程验证与结果整理
  - [ ] SubTask 8.1: 跑通 Q1–Q4 全部脚本，输出 results/tables 与 figures
  - [ ] SubTask 8.2: 实验日志 run_xxx.json（seed/config/boundary_mode/trials/CI 方法）
  - [ ] SubTask 8.3: 更新工作区 `修改日志.md` 汇总本次全部改动

# Task Dependencies
- [Task 2] depends on [Task 1]
- [Task 3] depends on [Task 1]
- [Task 4] depends on [Task 3]（k̄ 使用 V_AA）
- [Task 5] depends on [Task 4]
- [Task 6] depends on [Task 3][Task 4][Task 5]
- [Task 7] depends on [Task 4]（正确性先行）
- [Task 8] depends on [Task 2–7]

# 可并行任务
- [Task 2] 与 [Task 3] 在 [Task 1] 完成后可并行；
- [Task 7.2] 的并行化与 [Task 5]/[Task 6] 串行执行后再做（避免以正确性为代价）。
