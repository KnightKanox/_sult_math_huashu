# 线段-线段最短距离算法的单元测试
import numpy as np
import pytest

from src.geometry.segment_distance import segment_segment_distance


# 两条平行线段：水平共面，间距应为 5.0
def test_parallel_segments_distance():
    d = segment_segment_distance(np.array([0.0, 0, 0]), np.array([10.0, 0, 0]),
                                 np.array([0.0, 5, 0]), np.array([10.0, 5, 0]))
    assert d == pytest.approx(5.0, abs=1e-9)


# 相交线段：距离应为 0
def test_intersecting_segments_distance():
    d = segment_segment_distance(np.array([-1.0, 0, 0]), np.array([1.0, 0, 0]),
                                 np.array([0.0, -1, 0]), np.array([0.0, 1, 0]))
    assert d == pytest.approx(0.0, abs=1e-9)


# 共线且首尾相接的线段：距离应为 0
def test_touching_endpoint_segments_distance():
    d = segment_segment_distance(np.array([0.0, 0, 0]), np.array([5.0, 0, 0]),
                                 np.array([5.0, 0, 0]), np.array([9.0, 0, 0]))
    assert d == pytest.approx(0.0, abs=1e-9)


# 异面（skew）线段：s1 沿 x 轴(0,0,0)→(1,0,0)，s2 在 x=0,y=1 平面内 z∈[1,2]
# 最近点为 (0,0,0) 与 (0,1,1)，距离 sqrt(2)
def test_skew_segments_distance():
    d = segment_segment_distance(np.array([0.0, 0, 0]), np.array([1.0, 0, 0]),
                                 np.array([0.0, 1, 1]), np.array([0.0, 1, 2]))
    assert d == pytest.approx(np.sqrt(2.0), abs=1e-9)


# 点(退化线段)到线段：距离为点到直线垂距
def test_point_to_segment_distance():
    d = segment_segment_distance(np.array([0.0, 0, 0]), np.array([0.0, 0, 0]),
                                 np.array([0.0, 3, 0]), np.array([3.0, 0, 0]))
    assert d == pytest.approx(3.0 / np.sqrt(2.0), abs=1e-9)


# 线段端点位于另一线段内部附近（垂足落在线段内）：距离应为 0
def test_endpoint_inside_projection():
    d = segment_segment_distance(np.array([2.0, 2, 0]), np.array([2.0, -2, 0]),
                                 np.array([0.0, 0, 0]), np.array([4.0, 0, 0]))
    assert d == pytest.approx(0.0, abs=1e-9)
