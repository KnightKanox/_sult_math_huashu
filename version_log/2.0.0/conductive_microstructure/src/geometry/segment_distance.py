# 计算三维空间中两条有限线段之间的最短距离（经典算法，含退化/平行稳健处理）
import numpy as np


# 计算点 p 到有限线段 [a,b] 的最短距离
def _point_segment_dist(p, a, b):
    """返回点到线段的最短欧氏距离。"""
    d = b - a
    ll = float(d @ d)
    if ll <= 1e-12:
        return float(np.linalg.norm(p - a))
    t = max(0.0, min(1.0, float((p - a) @ d) / ll))
    return float(np.linalg.norm(p - (a + t * d)))


# 公开接口：计算点 p 到有限线段 [a,b] 的最短距离
def point_segment_distance(p, a, b):
    """返回点 p 到线段 a→b 的最短欧氏距离。"""
    return _point_segment_dist(np.asarray(p, float), np.asarray(a, float), np.asarray(b, float))


def segment_segment_distance(p1, q1, p2, q2):
    """计算线段1(p1→q1)与线段2(p2→q2)的最短欧氏距离。

    算法：两线段参数化为 p1+s*d1 与 p2+t*d2，先解无约束最近参数(s,t)，
    再按边界裁剪；平行/退化情形枚举端点投影取最小，保证稳健。
    """
    p1 = np.asarray(p1, float)
    q1 = np.asarray(q1, float)
    p2 = np.asarray(p2, float)
    q2 = np.asarray(q2, float)
    d1 = q1 - p1
    d2 = q2 - p2
    r = p1 - p2
    a = float(d1 @ d1)
    e = float(d2 @ d2)
    # 两条线段均退化为点
    if a <= 1e-12 and e <= 1e-12:
        return float(np.linalg.norm(r))
    # 线段1退化为点
    if a <= 1e-12:
        return _point_segment_dist(p1, p2, q2)
    # 线段2退化为点
    if e <= 1e-12:
        return _point_segment_dist(p2, p1, q1)
    b = float(d1 @ d2)
    c = float(d1 @ r)
    f = float(d2 @ r)
    denom = a * e - b * b
    if denom > 1e-12:
        # 非平行：解二次方程并逐段裁剪
        s = (b * f - c * e) / denom
        s = max(0.0, min(1.0, s))
        t = (b * s + f) / e
        if t < 0.0:
            t = 0.0
            s = max(0.0, min(1.0, -c / a))
        elif t > 1.0:
            t = 1.0
            s = max(0.0, min(1.0, (b - c) / a))
        return float(np.linalg.norm((p1 + s * d1) - (p2 + t * d2)))
    # 平行：最小距离必然出现在端点投影组合，取四种组合最小值
    return min(
        _point_segment_dist(p2, p1, q1),
        _point_segment_dist(q2, p1, q1),
        _point_segment_dist(p1, p2, q2),
        _point_segment_dist(q1, p2, q2),
    )
