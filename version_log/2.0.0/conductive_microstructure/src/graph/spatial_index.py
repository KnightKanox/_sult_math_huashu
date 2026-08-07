# 宽相位筛选：用轴对齐包围盒（AABB）距离快速筛出可能导通的线段候选对
# 数学依据：线段间的真实距离 >= 其 AABB 间的最小距离，故 AABB 距离 <= 阈值
# 的候选集必然包含所有真实距离 <= 阈值的线段对（无漏判），再对候选做精确判距。
import numpy as np


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
