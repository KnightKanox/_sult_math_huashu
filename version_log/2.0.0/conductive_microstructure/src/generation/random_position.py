# 随机位置采样：在立方盒 [-box_half, box_half]^3 内均匀取点
import numpy as np


# 在立方盒内采样一个均匀分布的三维位置（形状 (3,)）
def random_position_in_box(rng, box_half=5000.0):
    """返回盒内均匀随机位置，每维独立服从 [-box_half, box_half] 均匀分布。

    默认 box_half 取全局 BOX_HALF=5000.0（纳米）。
    rng 必须为 numpy.random.Generator（调用方传入）。
    """
    return rng.uniform(-box_half, box_half, size=3)
