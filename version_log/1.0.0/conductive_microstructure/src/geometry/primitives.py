# 三维几何基元：有限线段（介质A轴）、节点表示
from dataclasses import dataclass
from typing import List

import numpy as np


# 节点：由一行或多行（同一导体的跨边界片段）合并而成，保存其全部轴线段
@dataclass
class Node:
    # 节点编号
    index: int
    # 该节点包含的原始行索引列表
    row_indices: List[int]
    # 轴线段端点数组，形状 (K, 2, 3)，[段k] = [p_k, q_k]
    segments: np.ndarray

    # 返回该节点全部端点中 x 的最小值
    def min_x(self) -> float:
        return float(self.segments[:, :, 0].min())

    # 返回该节点全部端点中 x 的最大值
    def max_x(self) -> float:
        return float(self.segments[:, :, 0].max())
