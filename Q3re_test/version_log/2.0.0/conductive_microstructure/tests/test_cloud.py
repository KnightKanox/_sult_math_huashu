# 概率云层（问题 3 Task 3）单元测试：AA 云采样、BB 框架校验、置信区间与等效体积
import numpy as np

from src.cloud.aa_cloud import estimate_q, estimate_q_aa
from src.cloud.effective_volume import effective_volume
from src.simulation.confidence import wilson_ci


# r=0：两圆柱中心重合、轴必交（轴距 0 ≤ 61.8，与方向无关），q_hat 应精确为 1.0
def test_estimate_q_aa_zero_distance():
    q_hat, success_count, sample_count, ci_low, ci_high = \
        estimate_q_aa(0.0, np.random.default_rng(0), samples_per_r=500)
    assert q_hat == 1.0
    assert success_count == 500


# r=5200 > 5061.8（同轴首尾相接极限）：任意方向两轴线段距离 ≥ 200 > 61.8，q_hat 应精确为 0.0
def test_estimate_q_aa_far_zero():
    q_hat, success_count, sample_count, ci_low, ci_high = \
        estimate_q_aa(5200.0, np.random.default_rng(0), samples_per_r=500)
    assert q_hat == 0.0
    assert success_count == 0


# BB 连接与方向无关：q_BB(r)=I(r≤401.8)，用通用 estimate_q 框架验证其正确性
def test_bb_step_function():
    rng = np.random.default_rng(0)
    for r, expected in ((400.0, 1.0), (410.0, 0.0)):
        # 用默认参数绑定外层 r，避免循环变量捕获问题
        def connected_func(r_hat, _r=r):
            return _r <= 401.8
        q_hat, success_count, sample_count, ci_low, ci_high = \
            estimate_q(r, rng, connected_func, 100)
        assert q_hat == expected


# 解析 BB 云（半径 401.8 的球）数值积分验证 effective_volume
# 注：r 步长 10nm 的梯形近似对阶跃 q 有约 2.4% 离散误差（r 截断于 400 而非 401.8），故容差取 0.03
def test_bb_volume_analytic():
    r_grid = np.arange(0, 501, 10)
    q_grid = (r_grid <= 401.8).astype(float)
    v_num = effective_volume(r_grid, q_grid)
    v_true = 4.0 / 3.0 * np.pi * 401.8 ** 3
    assert abs(v_num - v_true) / v_true < 0.03


# Wilson 区间合理性：wilson_ci(1000, 2000) 包含 0.5，区间关于 0.5 对称，半宽在 [0.018, 0.026]
def test_wilson_ci_sanity():
    lo, hi = wilson_ci(1000, 2000)
    assert lo < 0.5 < hi
    half_width = (hi - lo) / 2.0
    assert 0.018 <= half_width <= 0.026


# 可复现性：相同 seed 的两次估计结果完全一致
def test_reproducible_cloud():
    a = estimate_q_aa(500.0, np.random.default_rng(1), 200)
    b = estimate_q_aa(500.0, np.random.default_rng(1), 200)
    assert a == b
