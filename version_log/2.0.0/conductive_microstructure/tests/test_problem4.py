# 问题 4 单元测试：B 球生成、A-B/B-B/B-电极判据临界、混合单 trial（run_single_trial_mixed）、
# 概率云（q_BB/V_BB/V_AB 与 capsule 近似）、连接矩阵 λ_max、total_cost 与 Wilson 判定
# 全部为快速测试：混合单 trial 用少量样本（n_a<=30、n_b<=50）或 monkeypatch 注入固定几何，
# V_AB 用小 dr=25、samples_per_r=800（约 1.5s），不运行主流程与重型 MC。
import numpy as np
import pytest

import src.simulation.single_trial as st
from src.cloud.mixed_cloud import (connection_matrix, estimate_v_ab,
                                   lambda_max, q_bb, total_cost,
                                   v_ab_capsule_approx, v_bb_analytic)
from src.config import (AB_SEG_THRESHOLD, A_LENGTH, A_RADIUS,
                        B_SPHERE_BOX_HALF, B_RADIUS, BB_CENTER_THRESHOLD,
                        BE_AXIS_THRESHOLD, BoundaryMode, COST_A_PER_UM3,
                        COST_B_PER_UM3, LEFT_PLANE_X, RIGHT_PLANE_X)
from src.generation.medium_generator import generate_b_spheres, make_cylinder_segments
from src.geometry.cylinder_distance import (cylinder_sphere_connected,
                                            sphere_electrode_connected,
                                            spheres_connected)
from src.simulation.confidence import wilson_one_sided_lower
from src.simulation.single_trial import run_single_trial, run_single_trial_mixed


# B 球生成：n_b 个球心全部落在 [-4800,4800]^3 内（球完全在盒内），n_b=0 返回 (0,3) 空数组
def test_generate_b_spheres_range_and_count():
    spheres = generate_b_spheres(np.random.default_rng(1), 100)
    assert spheres.shape == (100, 3)
    assert np.all(spheres >= -B_SPHERE_BOX_HALF)
    assert np.all(spheres <= B_SPHERE_BOX_HALF)
    empty = generate_b_spheres(np.random.default_rng(1), 0)
    assert empty.shape == (0, 3)


# B-B 判据临界：中心距恰为 401.8（2R_B+δ）导通，401.81 不导通
def test_bb_critical_threshold():
    assert spheres_connected((0.0, 0.0, 0.0), (BB_CENTER_THRESHOLD, 0.0, 0.0)) is True
    assert spheres_connected((0.0, 0.0, 0.0), (BB_CENTER_THRESHOLD + 0.01, 0.0, 0.0)) is False


# A-B 判据临界：球心到 A 轴线段距离恰为 231.8（R_A+R_B+δ）导通，231.81 不导通
def test_ab_critical_threshold():
    seg_p = np.array([-2500.0, 0.0, 0.0])
    seg_q = np.array([2500.0, 0.0, 0.0])
    assert cylinder_sphere_connected(seg_p, seg_q, (0.0, AB_SEG_THRESHOLD, 0.0)) is True
    assert cylinder_sphere_connected(seg_p, seg_q, (0.0, AB_SEG_THRESHOLD + 0.01, 0.0)) is False


# B-电极判据临界：球心距电极平面 ≤ 201.8（R_B+δ）导通、> 201.8 不导通，且左右平面按 x 分量方向区分
# 注：因 double 减法舍入（-5000-(-4798.2) 恰好 201.80000000000018>201.8），用 ±0.01 裕度逼近临界
def test_be_critical_threshold():
    c_left = (LEFT_PLANE_X + BE_AXIS_THRESHOLD - 0.01, 0.0, 0.0)
    c_right = (RIGHT_PLANE_X - BE_AXIS_THRESHOLD + 0.01, 0.0, 0.0)
    # 注：sphere_electrode_connected 返回 np.bool_，故用 == True/False（等价断言）而非 is
    assert sphere_electrode_connected(c_left, LEFT_PLANE_X) == True
    assert sphere_electrode_connected(
        (LEFT_PLANE_X + BE_AXIS_THRESHOLD + 0.01, 0.0, 0.0), LEFT_PLANE_X) == False
    assert sphere_electrode_connected(c_right, RIGHT_PLANE_X) == True
    assert sphere_electrode_connected(
        (RIGHT_PLANE_X - BE_AXIS_THRESHOLD - 0.01, 0.0, 0.0), RIGHT_PLANE_X) == False
    # x 分量方向性：靠近左电极的球不导通右电极，反之亦然
    assert sphere_electrode_connected(c_left, RIGHT_PLANE_X) == False
    assert sphere_electrode_connected(c_right, LEFT_PLANE_X) == False


# 混合单 trial：纯 B（n_a=0, n_b=50）不崩溃，返回字段齐全且 B 球每球一个节点
def test_mixed_pure_b_no_crash():
    res = run_single_trial_mixed(np.random.default_rng(4), 0, 50,
                                 BoundaryMode.WRAPPED_GEOMETRY_ONLY)
    expected = {"connected", "node_count", "segment_count", "edge_count",
                "left_node_count", "right_node_count", "max_component_ratio", "mean_degree"}
    assert set(res) == expected
    assert res["node_count"] == 50
    assert res["segment_count"] == 0
    assert 0.0 <= res["max_component_ratio"] <= 1.0


# 一致性回归：n_b=0 时 run_single_trial_mixed 与 run_single_trial 同种子逐字段一致
def test_mixed_pure_a_consistent_with_run_single_trial():
    mode = BoundaryMode.PERIODIC_CONNECTED
    a = run_single_trial(np.random.default_rng(3), 20, mode)
    b = run_single_trial_mixed(np.random.default_rng(3), 20, 0, mode)
    assert set(a) == set(b)
    for key in a:
        assert a[key] == b[key], f"纯A一致性失败: {key} 不一致"


# 构造混合链式贯通：A 圆柱触左电极 + B 球链 + A 圆柱触右电极（monkeypatch 注入固定几何），两种模式均导通
def test_mixed_chain_connected(monkeypatch):
    segs = np.concatenate([
        make_cylinder_segments(np.array([-2500.0, 0.0, 0.0]),
                               np.array([1.0, 0.0, 0.0])).reshape(1, 6),
        make_cylinder_segments(np.array([2499.5, 0.0, 0.0]),
                               np.array([1.0, 0.0, 0.0])).reshape(1, 6),
    ], axis=0)
    # 球1(150,150,0)距圆柱1 右端点 (0,0,0) 为 212.1<=231.8；球2(250,0,0) 距球1 为 180.3<=401.8、
    # 且落在圆柱2 轴线上（距离 0<=231.8），构成 电极L-A1-球1-球2-A2-电极R 链
    spheres = np.array([[150.0, 150.0, 0.0], [250.0, 0.0, 0.0]])
    monkeypatch.setattr(st, "generate_batch", lambda rng, n: segs)
    monkeypatch.setattr(st, "generate_b_spheres", lambda rng, n: spheres)
    for mode in BoundaryMode:
        res = st.run_single_trial_mixed(np.random.default_rng(0), 2, 2, mode)
        assert res["connected"] is True
        assert res["left_node_count"] >= 1
        assert res["right_node_count"] >= 1
        assert res["node_count"] == 4  # 2 个 A 节点 + 2 个 B 球节点


# 可复现性：相同 seed 两次 run_single_trial_mixed 所有返回字段完全一致
def test_mixed_reproducible():
    a = run_single_trial_mixed(np.random.default_rng(5), 20, 20,
                               BoundaryMode.PERIODIC_CONNECTED)
    b = run_single_trial_mixed(np.random.default_rng(5), 20, 20,
                               BoundaryMode.PERIODIC_CONNECTED)
    for key in a:
        assert a[key] == b[key], f"可复现失败: {key} 不一致"


# 概率云：q_BB 为阶跃函数，中心距 401.8 处 q=1、401.81 处 q=0（与方向无关）
def test_q_bb_step_function():
    assert q_bb(BB_CENTER_THRESHOLD) == 1.0
    assert q_bb(BB_CENTER_THRESHOLD + 0.01) == 0.0


# 概率云：V_BB 解析值 = 4/3·π·401.8³，与手算解析值相对误差 < 1e-6
def test_v_bb_analytic_value():
    v_true = 4.0 / 3.0 * np.pi * BB_CENTER_THRESHOLD ** 3
    assert abs(v_bb_analytic() - v_true) / v_true < 1e-6


# 概率云：V_AB 数值积分（dr=25、M=800 约 1.5s）与 capsule 近似理论值相对偏差 ≤ 15%（实测约 1%）
def test_v_ab_vs_capsule_approx():
    v_ab = estimate_v_ab(r_max=2731.8, dr=25.0,
                         rng=np.random.default_rng(2024), samples_per_r=800)
    cap = v_ab_capsule_approx()
    assert cap > 0.0
    assert abs(v_ab - cap) / cap <= 0.15


# 连接矩阵/λ_max：矩阵对称、元素非负、λ_max≥0，且 λ_max 随 n_a 增大不降（理论单调 sanity）
def test_connection_matrix_and_lambda_max():
    m = connection_matrix(100, 200, 1e10, 2e9, 3e9)
    assert m.shape == (2, 2)
    assert np.allclose(m, m.T)
    assert np.all(m >= 0.0)
    assert lambda_max(m) >= 0.0
    lam_lo = lambda_max(connection_matrix(50, 200, 1e10, 2e9, 3e9))
    lam_hi = lambda_max(connection_matrix(100, 200, 1e10, 2e9, 3e9))
    assert lam_hi >= lam_lo


# total_cost：nm³→μm³ 换算正确（1.05·V_A/1e9 与 0.05·V_B/1e9 与手算一致），且随 N_A/N_B 单调不减
def test_total_cost_units_and_monotone():
    v_a = np.pi * A_RADIUS ** 2 * A_LENGTH
    v_b = 4.0 / 3.0 * np.pi * B_RADIUS ** 3
    assert total_cost(1, 0) == pytest.approx(COST_A_PER_UM3 * v_a / 1e9)
    assert total_cost(0, 1) == pytest.approx(COST_B_PER_UM3 * v_b / 1e9)
    assert total_cost(3, 5) == pytest.approx(3.0 * total_cost(1, 0) + 5.0 * total_cost(0, 1))
    assert total_cost(10, 0) <= total_cost(11, 0)
    assert total_cost(0, 10) <= total_cost(0, 11)


# Wilson 判定沿用 Q3 判据：950/1000 下界 ≥ 0.90，850/1000 下界 < 0.90
def test_wilson_decision_boundary():
    assert wilson_one_sided_lower(950, 1000) >= 0.90
    assert wilson_one_sided_lower(850, 1000) < 0.90
