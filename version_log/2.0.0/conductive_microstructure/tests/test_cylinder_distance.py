# 圆柱导通判据（胶囊近似）与电极判据的单元测试
import numpy as np

from src.geometry.cylinder_distance import cylinders_connected, electrode_connected
from src.config import AA_SEG_THRESHOLD, A_RADIUS, LEFT_PLANE_X, RIGHT_PLANE_X


# 两根平行圆柱：轴距 61.8（恰好阈值）应导通，61.81 应不导通
def test_aa_threshold_parallel():
    p1 = np.array([0.0, 0, 0]); q1 = np.array([100.0, 0, 0])
    p2 = np.array([0.0, 61.8, 0]); q2 = np.array([100.0, 61.8, 0])
    assert cylinders_connected(p1, q1, p2, q2, A_RADIUS, A_RADIUS)
    p2b = np.array([0.0, 61.81, 0]); q2b = np.array([100.0, 61.81, 0])
    assert not cylinders_connected(p1, q1, p2b, q2b, A_RADIUS, A_RADIUS)


# 交叉圆柱：轴线段相交，必然导通
def test_aa_crossing():
    p1 = np.array([-10.0, 0, 0]); q1 = np.array([10.0, 0, 0])
    p2 = np.array([0.0, -10, 0]); q2 = np.array([0.0, 10, 0])
    assert cylinders_connected(p1, q1, p2, q2, A_RADIUS, A_RADIUS)


# 左电极：端点贴 x=-5000 时导通；远离时不导通
def test_left_electrode():
    p = np.array([-5000.0, 0, 0]); q = np.array([-4900.0, 0, 0])
    assert electrode_connected(p, q, A_RADIUS, LEFT_PLANE_X)
    p2 = np.array([-4800.0, 0, 0]); q2 = np.array([-4700.0, 0, 0])
    assert not electrode_connected(p2, q2, A_RADIUS, LEFT_PLANE_X)


# 右电极：端点贴 x=5000 时导通
def test_right_electrode():
    p = np.array([4900.0, 0, 0]); q = np.array([5000.0, 0, 0])
    assert electrode_connected(p, q, A_RADIUS, RIGHT_PLANE_X)
