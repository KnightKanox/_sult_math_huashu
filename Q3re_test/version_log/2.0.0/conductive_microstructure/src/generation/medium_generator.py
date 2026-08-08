# 介质A圆柱生成核心：随机方向/中心采样 + 周期回绕边界解析切段
# 微构体 Ω=[-5000,5000]^3；介质A为 5000nm 长、半径 30nm 的直圆柱。
# 圆柱随机放置可能越界，按周期回绕规则 x' = ((x+5000) mod 10000) - 5000 被"切开"为多段在箱片段，
# 每段端点均在盒内，行格式 (p_x,p_y,p_z,q_x,q_y,q_z) 与问题 1 附件一致。
import numpy as np

from .random_orientation import random_unit_vector
from .random_position import random_position_in_box
from ..config import BOX_HALF, A_LENGTH

# 介质A圆柱半长（nm）
A_HALF = A_LENGTH / 2.0


# 将坐标按周期边界回绕到盒内（对 numpy 数组按元素运算）
def _wrap(v):
    """等价于 connectivity._periodic_wrap：x' = ((x+5000) mod 10000) - 5000。"""
    return ((np.asarray(v) + BOX_HALF) % (2.0 * BOX_HALF)) - BOX_HALF


# 解析求解圆柱沿参数 t 运动时与六个盒面的全部穿越时刻
def _cross_times(c, u):
    """返回排序后的穿越时刻列表（t∈(-A_HALF, A_HALF)，端点处视为不穿越）。

    对每轴：若 |u_axis|>1e-12，解 t1=(BOX_HALF-c_axis)/u_axis、
    t2=(-BOX_HALF-c_axis)/u_axis；收集严格落在 (-A_HALF, A_HALF) 内、
    且与 A_HALF / -A_HALF 差 >= 1e-9（浮点边界容差）的值。
    """
    times = []
    for axis in range(3):
        ua = float(u[axis])
        if abs(ua) > 1e-12:
            for target in (BOX_HALF, -BOX_HALF):
                t = (target - float(c[axis])) / ua
                if -A_HALF < t < A_HALF and abs(t - A_HALF) >= 1e-9 and abs(t + A_HALF) >= 1e-9:
                    times.append(t)
    return sorted(times)


# 由显式中心与单位方向构造该圆柱的"在箱片段集"
def make_cylinder_segments(c, u):
    """给定中心 c 与单位方向 u，返回 (K,2,3) 片段数组。

    按穿越时刻把参数区间 [-A_HALF, A_HALF] 切成 K 段；每段端点先取
    c+t*u 再周期回绕进盒；若回绕后两点距离 <= 1e-9（退化零长）则丢弃。
    """
    c = np.asarray(c, dtype=float)
    u = np.asarray(u, dtype=float)
    ts = [-A_HALF] + _cross_times(c, u) + [A_HALF]
    segments = []
    for i in range(len(ts) - 1):
        a = _wrap(c + ts[i] * u)
        b = _wrap(c + ts[i + 1] * u)
        if np.linalg.norm(a - b) > 1e-9:
            segments.append([a, b])
    return np.array(segments)


# 随机生成一根介质A圆柱的"在箱片段集"
def generate_a_cylinder(rng):
    """随机方向（球面均匀）+ 随机中心（盒内均匀），返回 (K,2,3) 片段数组。

    不越界时 K=1；越界被边界"切开"时 K>=2，同一圆柱的相邻片段
    在回绕后共享重合端点（坐标差 < PERIODIC_MATCH_TOL=1e-6）。
    """
    u = random_unit_vector(rng)
    c = random_position_in_box(rng)
    return make_cylinder_segments(c, u)


# 批量生成 n 根圆柱，拼接为与问题 1 附件一致的行格式
def generate_batch(rng, n):
    """生成 n 根圆柱的全部在箱片段，返回 (M,6) 数组（每行 [p_x,p_y,p_z,q_x,q_y,q_z]）。

    每根圆柱按序调用 generate_a_cylinder（固定 RNG 消耗顺序，保证可复现）。
    """
    rows = []
    for _ in range(n):
        segs = generate_a_cylinder(rng)
        rows.append(segs.reshape(-1, 6))
    if not rows:
        return np.empty((0, 6))
    return np.concatenate(rows, axis=0)
