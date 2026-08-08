# 导通判据：介质-介质（胶囊近似）、介质-电极平面
import numpy as np

from .segment_distance import segment_segment_distance, point_segment_distance
from ..config import A_RADIUS, B_RADIUS, CONNECT_DELTA


# 判断两条有限圆柱（胶囊近似：轴线段向四周膨胀各自半径）是否导通
def cylinders_connected(p1, q1, p2, q2, r1, r2, delta=CONNECT_DELTA):
    """当两轴线段最短距离 <= r1+r2+delta 时，两圆柱视为电学导通。"""
    d = segment_segment_distance(np.asarray(p1, float), np.asarray(q1, float),
                                 np.asarray(p2, float), np.asarray(q2, float))
    return d <= r1 + r2 + delta


# 判断一根圆柱（胶囊近似）是否与垂直于X轴的电极平面导通
def electrode_connected(p, q, radius, plane_x, delta=CONNECT_DELTA):
    """圆柱(半径radius)到平面 x=plane_x 的距离 <= delta 即导通。

    轴到平面的距离 d_axis = max(0, x_min-plane_x, plane_x-x_max)，
    胶囊表面到平面距离 = max(0, d_axis - radius)。
    """
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    x_min = float(min(p[0], q[0]))
    x_max = float(max(p[0], q[0]))
    d_axis = max(0.0, x_min - plane_x, plane_x - x_max)
    d = max(0.0, d_axis - radius)
    return d <= delta


# 判断介质A圆柱（胶囊近似）与介质B球是否导通
def cylinder_sphere_connected(p, q, sphere_center, delta=CONNECT_DELTA):
    """A圆柱轴线段(p→q)到球心距离 <= A_RADIUS+B_RADIUS+delta 时视为导通。"""
    d = point_segment_distance(sphere_center, np.asarray(p, float), np.asarray(q, float))
    return d <= A_RADIUS + B_RADIUS + delta


# 判断两个介质B球是否导通
def spheres_connected(c1, c2, delta=CONNECT_DELTA):
    """两球心距离 <= 2*B_RADIUS+delta 时视为导通。"""
    d = float(np.linalg.norm(np.asarray(c1, float) - np.asarray(c2, float)))
    return d <= 2.0 * B_RADIUS + delta


# 判断介质B球是否与垂直于X轴的电极平面导通
def sphere_electrode_connected(c, plane_x, delta=CONNECT_DELTA):
    """球心 c 到平面 x=plane_x 的最近距离（max(0, |plane_x-c[0]|-B_RADIUS)）<= delta 即导通。"""
    c = np.asarray(c, float)
    d = max(0.0, abs(plane_x - c[0]) - B_RADIUS)
    return d <= delta
