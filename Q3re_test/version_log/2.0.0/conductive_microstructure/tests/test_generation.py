# 随机生成层（问题 2）单元测试：方向/位置采样、跨边界切段、周期合并复现
import numpy as np

from src.config import BoundaryMode
from src.graph.connectivity import analyze_group
from src.generation.random_orientation import random_unit_vector
from src.generation.random_position import random_position_in_box
from src.generation.medium_generator import generate_batch, make_cylinder_segments


# 方向采样：单位模长 + 二阶矩均匀（<u_x²>=<u_y²>=<u_z²>≈1/3）+ 一阶/混合矩≈0
def test_orientation_norm_and_moments():
    rng = np.random.default_rng(42)
    us = np.array([random_unit_vector(rng) for _ in range(20000)])
    # 模长为 1
    norms = np.linalg.norm(us, axis=1)
    assert np.max(np.abs(norms - 1.0)) < 1e-9
    # 二阶矩各向同性
    m2 = np.mean(us ** 2, axis=0)
    assert np.all(np.abs(m2 - 1.0 / 3.0) < 0.02)
    # 一阶矩与交叉矩接近 0
    assert np.all(np.abs(np.mean(us, axis=0)) < 0.02)
    assert abs(np.mean(us[:, 0] * us[:, 1])) < 0.02


# 位置采样：全部落在 [-5000,5000]^3 内，且各维均值接近盒中心 0
def test_position_in_box():
    rng = np.random.default_rng(42)
    pts = np.array([random_position_in_box(rng) for _ in range(2000)])
    assert np.all(pts >= -5000.0)
    assert np.all(pts <= 5000.0)
    assert np.all(np.abs(np.mean(pts, axis=0)) < 200.0)


# 不越界圆柱：端点在 ±2500，应返回 1 段
def test_no_cross_returns_one_segment():
    segs = make_cylinder_segments(np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]))
    assert segs.shape == (1, 2, 3)


# 沿 x 轴跨边界圆柱（q 端 x=5500 越界 500nm）：应切成 2 段且全部端点回绕在盒内
def test_cross_x_returns_two_segments():
    segs = make_cylinder_segments(np.array([3000.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]))
    assert len(segs) == 2
    # 两段所有端点 x 坐标都在盒内
    xs = segs[..., 0].ravel()
    assert np.all(xs >= -5000.0 - 1e-9)
    assert np.all(xs <= 5000.0 + 1e-9)
    # 穿越点未回绕坐标为 x=5000，回绕后为 x=-5000（二者周期等价），
    # 一段应含该边界重合点（|x-5000|<1e-6 或回绕等价式 |x+5000|<1e-6）
    assert np.any(np.abs(xs - 5000.0) < 1e-6) or np.any(np.abs(xs + 5000.0) < 1e-6)
    # 另一段应含 x=-4500 附近（q 端 5500 回绕后），说明回绕发生
    assert np.any(np.abs(xs + 4500.0) < 1e-6)


# 批量生成 50 根：所有端点坐标均在盒内（含 1e-6 浮点容差）
def test_batch_endpoints_in_box():
    batch = generate_batch(np.random.default_rng(7), 50)
    assert batch.ndim == 2 and batch.shape[1] == 6 and batch.shape[0] >= 50
    assert np.all(batch >= -5000.0 - 1e-6)
    assert np.all(batch <= 5000.0 + 1e-6)


# 关键验证：同一根圆柱的跨边界片段经 PERIODIC_CONNECTED 合并回同一节点，node_count 应等于圆柱数
def test_periodic_merge_recovers_cylinder_count():
    batch = generate_batch(np.random.default_rng(7), 50)
    res = analyze_group(batch, BoundaryMode.PERIODIC_CONNECTED)
    assert res["node_count"] == 50


# 可复现性：相同种子两次生成结果完全一致
def test_reproducible_seed():
    a = generate_batch(np.random.default_rng(123), 20)
    b = generate_batch(np.random.default_rng(123), 20)
    assert np.array_equal(a, b)
