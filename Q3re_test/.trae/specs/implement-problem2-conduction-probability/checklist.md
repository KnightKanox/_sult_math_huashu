# Checklist — 问题 2：导通概率 Monte Carlo 仿真（V2.0.0）

- [x] `version_log/2.0.0/` 基线复制完成，matplotlib 可用
- [x] 随机方向球面均匀（矩检验通过），生成器切段后所有端点均在 `[-5000,5000]³` 内
- [x] 概率云框架通过 BB 解析基准校验：`q_BB(r)≈I(r≤401.8)`、`V_BB≈4/3π·401.8³`
- [x] `q_AA(0)=1`、`q_AA(r>5061.8)=0`，`cloud_AA.csv` 生成
- [x] `single_trial` 确定性基准通过：`N_A=0` 不导通；横跨左右电极的圆柱导通
- [x] `connectivity.py` PERIODIC 合并改为哈希后，Q1 三组结果回归一致
- [x] `confidence.py` Wilson/Clopper-Pearson 双侧 95% 与单侧下界正确
- [x] `solve_problem2.py` 输出 `4φ × 2 模式` 结果表，相同 seed 完全可复现
- [x] `P_conn(φ)`、`k̄(φ)` 图生成且中文正常
- [x] 收敛性验证：M 增大 CI 收窄、P̂ 稳定；不同 seed 波动合理
- [x] `results/problem2/`、`results/cloud/`、`results/figures/` 齐全
- [x] `问题2_说明文档.md`、`version_log/2.0.0/版本说明.md` 完成，`修改日志.md` 已更新
