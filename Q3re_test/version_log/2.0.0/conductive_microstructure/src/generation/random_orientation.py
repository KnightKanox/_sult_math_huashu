# 随机方向采样：球面均匀单位向量（z 均匀 + 方位角均匀，避免极区聚集）
import numpy as np


# 采样一个球面均匀分布的 (3,) 单位方向向量
def random_unit_vector(rng):
    """用 z=U(-1,1)、phi=U(0,2π) 采样球面均匀单位向量（形状 (3,)）。

    此方法保证立体角均匀；禁止 θ~U(0,π) 方式（会在两极聚集）。
    rng 必须为 numpy.random.Generator（调用方传入）。
    """
    z = rng.uniform(-1.0, 1.0)
    phi = rng.uniform(0.0, 2.0 * np.pi)
    s = np.sqrt(1.0 - z * z)
    return np.array([s * np.cos(phi), s * np.sin(phi), z])
