# Tasks — 问题 4：A/B 混合体系最低成本组合

> 依赖主线：几何/生成扩展（T2）→ 混合单 trial（T3）；概率云扩展（T1）→ 连接矩阵（T4）；T3+T4 → 主脚本（T5）→ 测试（T6）→ 运行（T7）→ 文档（T8）。
> T1、T2 无相互依赖，可并行。

- [x] T1: 概率云层扩展（双介质）
  - 新增 `src/cloud/mixed_cloud.py`：
    - `estimate_q_ab(r, rng, samples_per_r)`：A 圆柱随机姿态 + B 球球壳采样，判据点-线段距离 ≤ R_A+R_B+δ；返回 q_AB(r)
    - `q_bb(r)`：解析 `I(r ≤ 2R_B+δ)`
    - `estimate_v_ab(r_max, dr, rng, samples_per_r)`：`4π∫r²q_AB dr`（np.trapezoid）
    - `v_bb_analytic()`：`4/3π(2R_B+δ)³`
    - `v_ab_capsule_approx()`：`πR_AB²L_A+4/3πR_AB³`，R_AB=R_A+R_B+δ
  - 验证：q_bb 与解析基准一致；V_BB 数值=解析；V_AB 与 capsule 近似相对偏差 ≤ 15%；可复现（同 seed 同结果）

- [x] T2: 生成层与几何判据扩展（可与 T1 并行）
  - `src/generation/medium_generator.py` 新增 `generate_b_spheres(rng, n_b)`：球心在 `[-4800,4800]³` 均匀（球完全在盒内）
  - `src/geometry/segment_distance.py` 暴露 `point_segment_distance(p, a, b)`（复用内部 `_point_segment_dist`）
  - `src/geometry/cylinder_distance.py` 新增：
    - `cylinder_sphere_connected(p, q, sphere_center, delta)`：点-线段距 ≤ R_A+R_B+δ
    - `spheres_connected(c1, c2, delta)`：中心距 ≤ 2R_B+δ
    - `sphere_electrode_connected(c, plane_x, delta)`：`|plane_x-|x_c|| ≤ R_B+δ`
  - `src/config.py` 补常量（不改现有）：`AB_SEG_THRESHOLD=231.8`、`BB_CENTER_THRESHOLD=401.8`、`BE_AXIS_THRESHOLD=201.8`、`B_SPHERE_BOX_HALF=4800.0`、成本常量 `COST_A_PER_UM3=1.05`、`COST_B_PER_UM3=0.05`、`NM3_PER_UM3=1e9`

- [x] T3: 混合宽相位与混合单次仿真
  - `src/graph/spatial_index.py` 新增 `aabb_candidates_mixed(seg_boxes, sphere_boxes, threshold)`（或等价统一入口）：A-A 对沿用现有 `aabb_candidates`；A-B/B-B 对按 AABB 扩边筛候选
  - `src/simulation/single_trial.py` 新增 `run_single_trial_mixed(rng, n_a, n_b, mode)`：
    - A：`generate_batch`（沿用切段/PERIODIC 合并语义）；B：`generate_b_spheres`（每球独立节点）
    - 电极：A-电极沿用；B-电极 `sphere_electrode_connected`
    - 连边：A-A（现有）、A-B（点-线段距）、B-B（中心距），均先 AABB 筛后精确判距
    - 返回与 Q2/Q3 一致的统计 dict；`n_a=0` 或 `n_b=0` 时退化为单介质不崩溃
  - 单元测试：纯 B 体系贯通、B-B 临界 401.8/401.81、A-B 临界 231.8/231.81、混合链式贯通、可复现

- [x] T4: 连接矩阵与 λ_max 理论筛选
  - `src/cloud/mixed_cloud.py`（或 `src/cloud/effective_volume.py` 扩展）新增：
    - `connection_matrix(n_a, n_b, v_aa, v_ab, v_bb)`：M=[[ρ_A V_AA, ρ_B V_AB],[ρ_A V_AB, ρ_B V_BB]]
    - `lambda_max(m)`：`np.linalg.eigvalsh(m).max()`
    - `total_cost(n_a, n_b)`：`c_A·N_A·V_A(μm³)+c_B·N_B·V_B(μm³)`（1μm³=10⁹nm³）
  - 单元测试：矩阵对称性、λ_max≥0、成本单调性、单位换算（V_A/V_B 的 nm³→μm³）

- [x] T5: 主脚本 `scripts/solve_problem4.py` 三级流程
  - CLI：`--n-a-max 360 --n-a-step 10 --n-b-max 5000 --trials 500 --confirm-trials 4000 --seed 42 --mode all --workers 8 --out-dir --plot/--no-plot --theory-only`
  - Level 1 理论筛选：λ_max 扫描 N_A×N_B 网格（0~N_A_max × 0~N_B_max），保存热图数据与理论可行区
  - Level 2 MC 边界：对每个 N_A 网格点二分最小 N_B（判据 `wilson_one_sided_lower ≥ 0.90`），记录各步 p̂/p_lower/动作；边界单调性校验
  - Level 3 成本优化：沿边界计算 `total_cost`，取全局最小候选；对候选及邻域（N_A±step、N_B±50）confirm_trials 复核；输出最终 (N_A*, N_B*)、成本、确认记录
  - 可复现：主 rng 一次性派生评估点子种子池，固定评估顺序（逐模式、逐 N_A、逐二分迭代）；`--workers` 并行且结果与串行一致
  - 输出：`results/problem4/problem4_result.csv` / `.json`
  - 图（`--plot`）：λ_max 热图（Figure 6）、MC P_conn 热图（Figure 7）、成本等高线 + P=0.9 可行边界（Figure 8）；中文标签 Microsoft YaHei，组合字符用 mathtext

- [x] T6: 单元测试 `tests/test_problem4.py`（不跑真实重 MC，mock/小样本）
  - 球生成范围、B-B/A-B/B-电极判据临界、混合单 trial（纯 B、纯 A、混合链式、可复现）、q_bb/V_BB 解析、V_AB 与 capsule 近似、连接矩阵/λ_max、total_cost 单位换算与单调性、Wilson 判定沿用
  - 全量 `python -m pytest tests -q` 通过（既有 42 + 新增）

- [x] T7: 运行 Q4 主流程（用户已确认"实现后立即运行"）
  - 先 `python scripts/solve_problem4.py --theory-only`（快）确认理论筛选与代码链路
  - 再运行完整三级流程（双模式、并行 workers），记录运行时长与输出
  - 解析结果：两模式的最优 (N_A*, N_B*) 与最低成本、λ_max 理论区 vs MC 边界对照

- [x] T8: 文档补录
  - 新建 `version_log/2.0.0/问题4_说明文档.md`（建模、成本函数、连接矩阵、搜索方法、参数表、结果、理论-仿真对照）
  - 更新 `version_log/2.0.0/版本说明.md`（Q4 状态与结果）
  - 根 `修改日志.md` 新增条目

# Task Dependencies
- T1 与 T2 互不依赖（可并行）
- T3 依赖 T2；T4 依赖 T1（数据准备）与 T2（成本常量）
- T5 依赖 T3、T4
- T6 依赖 T1-T4
- T7 依赖 T5、T6
- T8 依赖 T7（结果）与 T5（脚本）
