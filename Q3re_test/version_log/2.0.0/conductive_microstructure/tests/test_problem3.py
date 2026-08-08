# 问题 3 单元测试：二分搜索收敛/容差/迭代上限、Wilson 单侧下界判定边界、φ↔N_A 与理论换算
# 全部为 mock 快速测试：注入假评估回调（p_hat=p_lower 单调递增、在 phi0 处恰跨 0.90），
# 不调用真实 Monte Carlo。判定标准即问题 3 的 p_lower(95%) >= 0.90。
import os
import sys

import numpy as np
import pytest

# 保证从项目根以 python -m pytest 运行时能 import scripts/ 下的模块：pytest 的 -m 运行
# 模式已把项目根加入 sys.path（与 tests/test_simulation.py 的 from src... 依赖一致），
# 这里按 solve_problem2.py 自身的 sys.path 处理方式补充 scripts 目录。
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))

from solve_problem2 import CYLINDER_VOLUME, n_a_from_phi  # noqa: E402
from solve_problem3 import (  # noqa: E402
    binary_search_min_phi,
    wilson_one_sided_lower,
)


# 构造单调递增的假评估回调：p_hat=p_lower，恰在 phi0 处等于 0.90（(φ-φ0)/0.002 斜率）
def make_fake_evaluator(phi0, trials=1000):
    """返回假 evaluate(phi, mode)：p_lower 随 φ 单调递增且在 phi0 处恰为 0.90。"""
    def evaluate(phi, mode):
        p = float(np.clip(0.90 + (phi - phi0) / 0.002, 0.0, 1.0))
        k = int(round(p * trials))
        return {"connected": k, "trials": trials, "p_hat": p, "p_lower": p}
    return evaluate


# 二分收敛：单调真值回调下 phi_min 逼近 phi0，迭代次数不超上限且末次区间宽度 <= tol
def test_bisection_converges_to_threshold():
    phi0 = 0.005
    evaluate = make_fake_evaluator(phi0)
    res = binary_search_min_phi(0.002, 0.007, evaluate, "periodic_connected", tol=1e-4)
    assert res["iterations"] <= 40
    assert res["phi_high_end"] - res["phi_low_end"] <= 1e-4
    assert abs(res["phi_min"] - phi0) <= 1e-4
    assert res["phi_min"] == res["phi_high_end"]


# tol 更小时结果更紧：两次不同 tol 运行，区间宽度更小且 phi_min 更贴近 phi0
def test_smaller_tol_tighter_interval():
    phi0 = 0.005
    r_loose = binary_search_min_phi(0.002, 0.007, make_fake_evaluator(phi0),
                                    "wrapped_geometry_only", tol=1e-3)
    r_tight = binary_search_min_phi(0.002, 0.007, make_fake_evaluator(phi0),
                                    "wrapped_geometry_only", tol=1e-6)
    w_loose = r_loose["phi_high_end"] - r_loose["phi_low_end"]
    w_tight = r_tight["phi_high_end"] - r_tight["phi_low_end"]
    assert w_tight < w_loose
    assert w_tight <= 1e-6
    assert abs(r_tight["phi_min"] - phi0) <= 1e-6


# 最大迭代保护：max_iter=1 且 tol 极小而区间非空时抛 RuntimeError
def test_max_iter_protection_raises():
    evaluate = make_fake_evaluator(0.005)
    with pytest.raises(RuntimeError):
        binary_search_min_phi(0.002, 0.007, evaluate, "periodic_connected",
                              tol=1e-12, max_iter=1)


# Wilson 单侧 95% 下界判定边界换算：950/1000 满足 >=0.90，850/1000 不满足，且下界 <= p_hat
def test_wilson_decision_boundary():
    assert wilson_one_sided_lower(950, 1000) >= 0.90
    assert wilson_one_sided_lower(850, 1000) < 0.90
    for k in (850, 950):
        assert wilson_one_sided_lower(k, 1000) <= k / 1000.0


# φ↔N_A 与理论换算：n_a_from_phi(0.005)==354，CYLINDER_VOLUME/2.5493e9≈0.005546（k̄=1 临界填充率）
def test_phi_to_n_a_and_theory():
    assert n_a_from_phi(0.005) == 354
    assert CYLINDER_VOLUME / 2.5493e9 == pytest.approx(0.005546, rel=1e-3)


# binary_search 的 history 记录完整：每步含 phi/p_hat/p_lower/action，区间单调收窄且末步收敛
def test_history_complete_and_final_interval():
    res = binary_search_min_phi(0.002, 0.007, make_fake_evaluator(0.005),
                                "periodic_connected", tol=1e-4)
    hist = res["history"]
    assert len(hist) == res["iterations"]
    for h in hist:
        assert {"phi", "n_a", "connected", "trials", "p_hat", "p_lower",
                "action", "phi_low", "phi_high"}.issubset(h.keys())
        assert h["action"] in ("high", "low")
        assert h["phi_low"] <= h["phi"] <= h["phi_high"]
    los = [h["phi_low"] for h in hist]
    his = [h["phi_high"] for h in hist]
    assert all(b >= a for a, b in zip(los, los[1:]))
    assert all(b <= a for a, b in zip(his, his[1:]))
    assert hist[-1]["phi_high"] - hist[-1]["phi_low"] <= 1e-4
    assert hist[-1]["phi_low"] == res["phi_low_end"]
    assert hist[-1]["phi_high"] == res["phi_high_end"]


# 防回归：真实 MC evaluator 只消耗 trials 个子种子（池行宽为 max(trials, confirm_trials)），
# 且 p_hat = connected/trials 正确（此前的 bug：整行子种子全部传入导致 connected > trials、sqrt 域错误）
def test_mc_evaluator_uses_only_trials_seeds(monkeypatch):
    import solve_problem3 as sp3
    calls = {}

    # 记录 run_trials 收到的种子长度，并返回全导通（connected == trials）
    def fake_run_trials(seeds, n_a, mode, workers):
        calls["n_seeds"] = len(seeds)
        return len(seeds)

    monkeypatch.setattr(sp3, "run_trials", fake_run_trials)
    pool = np.full((3, 20), 7)  # 池行宽 20 > trials=10
    ev = sp3.make_mc_evaluator(pool, trials=10, workers=1)
    rec = ev(0.005, sp3.BoundaryMode.PERIODIC_CONNECTED)
    assert calls["n_seeds"] == 10
    assert rec["connected"] == 10
    assert rec["trials"] == 10
    assert rec["p_hat"] == 1.0
    # p̂=1 时 Wilson 单侧下界 = 1/(1+z²/n) < 1（与函数本身一致）
    assert rec["p_lower"] == sp3.wilson_one_sided_lower(10, 10)
    assert rec["p_lower"] < 1.0
    # 每行子种子消耗一次，下一次调用使用下一行
    ev(0.006, sp3.BoundaryMode.PERIODIC_CONNECTED)
    assert calls["n_seeds"] == 10
