# 微构体导电介质仿真优化 Spec（概率云—随机几何—连续渗流—Monte Carlo）

## Why
2026 第七届"华数杯"A 题要求：分析边长为 10000nm 立方体微构体中随机填充导电介质（介质 A 直圆柱、介质 B 球体）的导通问题，回答四问——(1) 给定三组介质 A 坐标是否导通；(2) 给定体积分数下的导通概率；(3) P_conn≥90% 时 A 的最低体积分数；(4) P_conn≥90% 时 A/B 混合的最低成本组合。

现有工作区仅有题目 PDF、论文模板、附件 Excel 与建模方案文档 `MODELING_PLAN_概率云_微构体导电优化.md`，尚无任何可运行代码。本 Spec 依据该建模方案落地一个**正确优先、可复现**的 Python 仿真项目。

## What Changes
- 新建完整 Python 项目（`conductive_microstructure/`），按方案 §37 目录结构组织：`src/`（geometry、generation、graph、cloud、simulation、optimization、visualization）、`scripts/`、`tests/`、`results/`。
- 实现三维计算几何核心：球面均匀方向采样、线段-线段距离、A-A/A-B/B-B/介质-电极导通判据、周期边界截断规则（两种边界模式可切换）。
- 实现概率云理论模块：q_AA(r)/q_AB(r)/q_BB(r) 数值计算、等效连接体积 V_eff 积分、平均连接度 k̄、双类型连接矩阵与 λmax 指标。
- 实现完整 Monte Carlo 仿真器：单次导通判定（DSU）+ 导通概率点估计与 95% 置信区间（Wilson/Clopper-Pearson）。
- 实现问题 1–4 求解脚本：确定性图连通判定；指定体积分数概率计算；单调性二分搜索求 90% 最低填充率；A/B 两级搜索（理论筛选 + 边界搜索）求最低成本。
- 单元测试覆盖几何正确性、BB 解析基准（401.8）、概率云校验、单调性、可复现性（固定 seed）。
- 生成论文所需图表（概率云、4πr²q(r)、P-φ、k̄-φ、λmax 热图、P_conn 热图、成本等高线+可行边界）。
- 所有实验记录 seed/config/边界模式/试验次数等，输出 `results/logs/run_xxx.json` 保证可复现。

## Impact
- 受影响目录：`d:\000AAAitaem\math\26huashu\A_try\`（新增 `conductive_microstructure/` 项目目录，不删除、不回滚任何现有文件）。
- 依赖：Python 3（numpy、scipy、pandas、openpyxl、matplotlib、pytest），写入 `requirements.txt`。
- 不改动题目原始数据 `附件.xlsx`。

## ADDED Requirements

### Requirement: 三维几何与导通判据
系统 SHALL 实现介质 A（有限圆柱）与介质 B（球体）的三维表示、随机生成（位置 U(Ω)，方向球面均匀采样）、以及以下导通判据：
- A-A：两有限圆柱真实最短距离 ≤ 1.8nm（默认 capsule/spherocylinder 近似，预留精确有限圆柱距离修正接口）；
- A-B：球心到圆柱最短距离 ≤ R_A+R_B+δ = 231.8nm；
- B-B：球心距 ≤ 2R_B+δ = 401.8nm；
- 介质-电极：介质到平面 x=-5000（L）/ x=5000（R）最短距离 ≤ 1.8nm。

#### Scenario: 单元测试
- **WHEN** 运行 `pytest tests/` 且构造两端点/平行/分离等人工样例
- **THEN** 所有几何判据结果与解析期望一致；B-B 边界 401.8/401.81 判断正确

### Requirement: 周期边界截断（两种模式）
系统 SHALL 按题目"边界截断规则"实现周期回绕 `x' = ((x+5000) mod 10000) - 5000`（y/z 同理），并支持两种边界语义：
- `PERIODIC_CONNECTED`（默认）：越界回绕后仍视为同一电学连续导体（考虑周期镜像连接）；
- `WRAPPED_GEOMETRY_ONLY`：仅几何回绕，回绕部分不额外建立连接。
两种模式通过 `SimulationConfig.boundary_mode` 切换，用于敏感性分析。

#### Scenario: 越界圆柱
- **WHEN** 一根介质 A 两端为 (3500,y1,z1) 与 (6000,y2,z2)
- **THEN** 越界部分回绕至 x=-5000~-4000，且整体仍为一个节点（PERIODIC_CONNECTED）

### Requirement: 概率云与等效连接体积
系统 SHALL 数值计算径向连接概率云 q_AA(r)、q_AB(r)、q_BB(r)（对每个固定 r 做 M≈1e4~1e5 次随机试验），并积分得到等效连接体积 V_eff = 4π∫r²q(r)dr，进而计算平均连接度 k̄ = ρ·V_eff。

#### Scenario: BB 解析基准校验
- **WHEN** 对 BB 概率云做 Monte Carlo 计算
- **THEN** q_BB(r) 逼近阶跃函数 I(r≤401.8)，且 V_BB ≈ 4/3π·401.8³（相对误差在设定容差内）

### Requirement: Monte Carlo 导通概率与置信区间
系统 SHALL 对给定 (N_A,N_B) 重复 M 次独立试验，输出导通概率点估计 p̂ = K/M 与 95% 置信区间；Q3/Q4 的可行判据 SHALL 采用 95% **单侧置信下界** p_lower ≥ 0.90（不允许仅凭 p̂ ≥ 0.90 判定），置信方法支持 Wilson 与 Clopper-Pearson。

#### Scenario: 概率估计
- **WHEN** 运行 `solve_problem2.py --phi 0.005 --trials 1000 --seed 42`
- **THEN** 输出 phi、N_A、p_hat、ci_low、ci_high、rho、k̄ 至 CSV

### Requirement: 问题 1 确定性连通判定
系统 SHALL 从 `附件.xlsx` 三个分表（组1=13根、组2=50根、组3=536根）读取两端点坐标，建图（节点=介质+虚拟电极 L/R），用 DSU 判断 find(L)==find(R)，输出每组导通情况，并输出总介质数、总边数、左右电极直接连接数及一条实际贯通路径（若有）。

### Requirement: 问题 2 指定体积分数导通概率
系统 SHALL 对 φ_A ∈ {0.50%, 0.60%, 0.70%, 1.00%}，按 N_A = round(φ_A·V_0/V_A) 生成介质并计算导通概率与 95% CI，同时给出 k̄，并绘制 P_conn(φ) 与 k̄(φ) 图。

### Requirement: 问题 3 最低 90% 填充率
系统 SHALL 利用单调性 P(N_A+1) ≥ P(N_A) 对 φ_A 做二分搜索，以 95% 单侧置信下界 ≥ 0.90 为可行判据，求解最低体积分数（精确到百分号下小数点后两位，即 0.01%），并对最终候选点提高试验次数复验。

### Requirement: 问题 4 A/B 混合最低成本
系统 SHALL 求解 min C(N_A,N_B) = 1.05·N_A·V_A + 0.05·N_B·V_B（成本单位为元，体积换算 1μm³=10⁹nm³），约束 p_lower,95% ≥ 0.90，N_A,N_B ∈ Z≥0。采用三级搜索：λmax(连接矩阵) 理论筛选 → 固定 N_A 二分最小可行 N_B（Monte Carlo）→ 边界成本比较取全局最低，并对最优点邻域高次数复验。

#### Scenario: 成本最优
- **WHEN** 输出问题 4 结果
- **THEN** 给出最优 (N_A, N_B)、总成本、p_hat、p_lower、试验次数，并绘制成本等高线 + P=0.9 可行边界图

### Requirement: 可复现性
所有随机过程 SHALL 支持固定 seed；每次实验 SHALL 记录 random_seed、git_commit、simulation_config、boundary_mode、number_of_trials、confidence_method、geometry_method 到 `results/logs/run_xxx.json`。

## MODIFIED Requirements
（无——本项目从零开始，无既有需求被修改。）

## REMOVED Requirements
（无。）
