# 置信区间基础：二项比例的 Wilson 分数区间、Clopper-Pearson 精确区间与 Wilson 单侧下界
# 供概率云（q_AA 采样比例）与蒙特卡洛判定共用；k 为成功数，n 为总试验数
import math

from scipy.stats import beta, norm


# 双侧 Wilson 分数区间（默认 95%，z=1.96）
def wilson_ci(k, n, z=1.96):
    """Wilson 分数区间（双侧），返回 (low, high)。

    p̂=0 或 p̂=1 时公式自然给出 0 / 1 边界，再做 clip 保证落在 [0,1]。
    """
    if n <= 0:
        raise ValueError("n 必须为正整数（总试验数）")
    p = float(k) / float(n)
    denom = 1.0 + z * z / float(n)
    center = p + z * z / (2.0 * float(n))
    half = z * math.sqrt(p * (1.0 - p) / float(n)
                         + z * z / (4.0 * float(n) * float(n)))
    low = max(0.0, (center - half) / denom)
    high = min(1.0, (center + half) / denom)
    return (low, high)


# 双侧 Clopper-Pearson 精确区间（默认 95%）
def clopper_pearson_ci(k, n, alpha=0.05):
    """Clopper-Pearson 精确区间（双侧），返回 (low, high)。

    用 Beta 分布分位：low=B(α/2; k, n-k+1)、high=B(1-α/2; k+1, n-k)；
    k=0 时 low=0，k=n 时 high=1。
    """
    if n <= 0:
        raise ValueError("n 必须为正整数（总试验数）")
    if k == 0:
        low = 0.0
    else:
        low = float(beta.ppf(alpha / 2.0, k, n - k + 1))
    if k == n:
        high = 1.0
    else:
        high = float(beta.ppf(1.0 - alpha / 2.0, k + 1, n - k))
    return (low, high)


# Wilson 单侧置信下界（用于问题 3 判定 P≥0.90 一类假设）
def wilson_one_sided_lower(k, n, confidence=0.95):
    """Wilson 单侧置信下界，返回下界值。

    z 取标准正态的 confidence 分位（如 0.95 → 1.645），公式与 Wilson 双侧
    下界相同但 z 为单侧分位。
    """
    if n <= 0:
        raise ValueError("n 必须为正整数（总试验数）")
    z = float(norm.ppf(confidence))
    p = float(k) / float(n)
    denom = 1.0 + z * z / float(n)
    center = p + z * z / (2.0 * float(n))
    half = z * math.sqrt(p * (1.0 - p) / float(n)
                         + z * z / (4.0 * float(n) * float(n)))
    return max(0.0, (center - half) / denom)
