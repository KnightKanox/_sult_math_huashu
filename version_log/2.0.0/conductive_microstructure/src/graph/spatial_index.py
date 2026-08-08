# 宽相位筛选：用轴对齐包围盒（AABB）距离快速筛出可能导通的线段候选对
# 数学依据：线段间的真实距离 >= 其 AABB 间的最小距离，故 AABB 距离 <= 阈值
# 的候选集必然包含所有真实距离 <= 阈值的线段对（无漏判），再对候选做精确判距。
import numpy as np


# 宽相位筛选（A-A 快速版）：AABB 筛候选后，再用包含直线距离预筛，减少精确判距调用
def aabb_candidates_fast(segs, threshold):
    """输入 (N,2,3) 轴线段，返回 (K,2) int 数组：真实距离可能 <= threshold 的候选对 (i,j)，i<j。

    AABB 宽相位对长线段（跨越盒子大部分）会产生大量候选（投影重叠但实际距离很远，
    实测 360 根圆柱时 ~1.6 万对）。再叠加"包含直线距离"预筛：两线段距离 >= 其所在
    包含直线的距离，故直线距离 > threshold 的对必不可能导通，直接向量化剔除，
    把精确判距的调用量降到数百量级；平行/近平行时直线距离不适用，保守保留。
    """
    pairs = aabb_candidates(segs, threshold)
    if len(pairs) == 0:
        return pairs
    p = segs[pairs[:, 0], 0]
    q = segs[pairs[:, 0], 1]
    r = segs[pairs[:, 1], 0]
    s = segs[pairs[:, 1], 1]
    n = np.cross(q - p, s - r)  # 两包含直线方向叉积（法向）
    n2 = np.einsum("ij,ij->i", n, n)
    ok = n2 > 1e-12  # 非平行（|n|² 足够大）时才可用直线距离
    line_dist2 = np.zeros(len(pairs), dtype=float)
    if np.any(ok):
        w = p[ok] - r[ok]
        num = np.abs(np.einsum("ij,ij->i", w, n[ok]))
        line_dist2[ok] = num * num / n2[ok]
    keep = (~ok) | (line_dist2 <= threshold * threshold)
    return pairs[keep]


# 对全部轴线段计算 AABB 间最小距离，返回距离 <= threshold 的下三角候选对
def aabb_candidates(segs, threshold):
    """输入 (N,2,3) 轴线段数组，返回 (K,2) int 数组：AABB 距离 <= threshold 的候选对 (i,j)，i<j。

    逐轴构造 (N,N) 区间距离矩阵并累加平方，避免一次性构造 (N,N,3) 大数组；
    N≈2500 时峰值内存为单个 (N,N) 矩阵量级，可接受。
    """
    segs = np.asarray(segs, dtype=float)
    n = len(segs)
    if n == 0:
        return np.empty((0, 2), dtype=np.int64)
    # 每根线段的 AABB：各轴 min/max
    lo = segs.min(axis=1)  # (N,3)
    hi = segs.max(axis=1)  # (N,3)
    sq = np.zeros((n, n), dtype=float)
    for ax in range(3):
        # 两区间沿该轴的间隔距离：max(0, lo_i - hi_j, lo_j - hi_i)
        d = lo[:, ax][:, None] - hi[None, :, ax]
        np.maximum(d, lo[None, :, ax] - hi[:, ax][:, None], out=d)
        np.maximum(d, 0.0, out=d)
        sq += d * d
    ii, jj = np.nonzero(sq <= threshold * threshold)
    keep = ii < jj
    if not np.any(keep):
        return np.empty((0, 2), dtype=np.int64)
    return np.column_stack((ii[keep], jj[keep]))


# 宽相位筛选（A/B 混合）：线段 AABB 与球心点 AABB 间距 <= threshold 的候选对
def aabb_candidates_mixed(segs, sphere_centers, threshold):
    """输入 (M,2,3) 轴线段与 (N,3) 球心数组，返回 (K,2) int 数组：候选对 (i,j)，第 0 列为线段索引、第 1 列为球索引。

    线段 AABB 为 lo=segs.min(axis=1)、hi=segs.max(axis=1)，B 球退化为点（lo=hi=球心）；
    逐轴构造 (M,N) 间隔距离矩阵并累加平方，避免一次性构造 (M,N,3) 大数组。
    点-线段真实距离 >= 点-AABB 距离，故所有真实距离 <= threshold 的对必被筛出（无漏判）。
    """
    segs = np.asarray(segs, dtype=float)
    sphere_centers = np.asarray(sphere_centers, dtype=float)
    m = len(segs)
    n = len(sphere_centers)
    if m == 0 or n == 0:
        return np.empty((0, 2), dtype=np.int64)
    lo = segs.min(axis=1)  # (M,3)
    hi = segs.max(axis=1)  # (M,3)
    sq = np.zeros((m, n), dtype=float)
    for ax in range(3):
        # 线段区间 [lo_i, hi_i] 与球心点 c_j 的间隔距离：max(0, lo_i-c_j, c_j-hi_i)
        d = lo[:, ax][:, None] - sphere_centers[None, :, ax]
        np.maximum(d, sphere_centers[None, :, ax] - hi[:, ax][:, None], out=d)
        np.maximum(d, 0.0, out=d)
        sq += d * d
    ii, jj = np.nonzero(sq <= threshold * threshold)
    if len(ii) == 0:
        return np.empty((0, 2), dtype=np.int64)
    return np.column_stack((ii, jj))


# 宽相位筛选（B 球）：球心点对距离 <= threshold 的下三角候选对（cKDTree 实现，O(N log N)）
def ball_candidates(centers, threshold):
    """输入 (N,3) 球心数组，返回 (K,2) int 数组：中心距 <= threshold 的候选对 (i,j)，i<j。

    用 scipy.spatial.cKDTree.query_pairs 精确枚举距离 <= threshold 的球心对，
    替代原 O(N²) 距离矩阵构造（N=5000 时 25M 元素 ~200MB 内存、~1.3s），
    使 N_B 高达数千的混合仿真单次耗时从秒级降到毫秒级；无候选时返回 (0,2) 空数组。
    """
    centers = np.asarray(centers, dtype=float)
    n = len(centers)
    if n < 2:
        return np.empty((0, 2), dtype=np.int64)
    from scipy.spatial import cKDTree
    pairs = cKDTree(centers).query_pairs(float(threshold))
    if not pairs:
        return np.empty((0, 2), dtype=np.int64)
    return np.asarray(sorted(pairs), dtype=np.int64)
