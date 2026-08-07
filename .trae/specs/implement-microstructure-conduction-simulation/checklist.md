# Checklist

- [ ] 项目骨架按方案 §37 目录结构创建，`requirements.txt` 与 `SimulationConfig` 定义齐全
- [ ] 随机方向采用球面均匀采样（z=cosθ~U(-1,1)）并通过均匀性测试
- [ ] 线段-线段距离、A-A/A-B/B-B/介质-电极距离判据实现正确并通过人工样例单元测试
- [ ] 周期边界截断规则实现正确，`PERIODIC_CONNECTED` / `WRAPPED_GEOMETRY_ONLY` 两模式可切换
- [ ] DSU 与连通判定正确（find(L)==find(R)）
- [ ] 问题 1：三个微构体导通判定输出正确，含总介质数/总边数/电极直连数/贯通路径
- [ ] 概率云：q_AA/q_AB/q_BB 数值计算实现，BB 概率云逼近阶跃函数 I(r≤401.8)
- [ ] 等效连接体积：V_BB 数值积分 ≈ 4/3π·401.8³（解析基准校验通过）
- [ ] Monte Carlo：输出点估计与 95% CI（Wilson/Clopper-Pearson），Q3/Q4 用 95% 单侧置信下界判定
- [ ] 问题 2：φ∈{0.50%,0.60%,0.70%,1.00%} 的导通概率表与 P_conn(φ)、k̄(φ) 图生成
- [ ] 问题 3：二分搜索得到 90% 最低 A 体积分数（精确到 0.01%），最终候选已复验
- [ ] 问题 4：λmax 理论筛选 + 边界搜索得到最低成本 (N_A, N_B)，最优点已复验
- [ ] 性能优化：空间索引 broad phase 与并行 MC 生效，未牺牲几何正确性
- [ ] 可复现性：固定 seed，实验结果记录 run_xxx.json（seed/config/boundary_mode/trials/CI 方法）
- [ ] 结果整理：results/tables、results/figures、results/logs 齐全，`修改日志.md` 已更新
