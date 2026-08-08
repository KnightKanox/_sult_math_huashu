# 连通性分析（含两种边界模式与跨边界合并）的单元测试
import numpy as np

from src.config import BoundaryMode
from src.graph.connectivity import analyze_group


# 三根介质A构成左右贯通的链：应判定导通（两种模式均可，贯通路径非空）
def test_chain_connected():
    # 链1：左电极(x=-5000) -> (0,0,0) ；链2：接通两段；链3：(0,0,0)->右电极(x=5000)
    e = np.array([
        [-5000.0, 0, 0, -1000.0, 0, 0],
        [-1000.0, 0, 0, 1000.0, 0, 0],
        [1000.0, 0, 0, 5000.0, 0, 0],
    ])
    for mode in BoundaryMode:
        res = analyze_group(e, mode)
        assert res["connected"] is True
        assert len(res["path_node_ids"]) >= 1
    # 行独立模式下路径应经过三个节点
    wo = analyze_group(e, BoundaryMode.WRAPPED_GEOMETRY_ONLY)
    assert len(wo["path_node_ids"]) == 3


# 左右两段中间断开（且无跨边界重合端点）：两种模式均应判定不导通
def test_chain_broken():
    e = np.array([
        [-5000.0, 0, 0, -3000.0, 0, 0],
        [3500.0, 0, 100, 5000.0, 0, 100],
    ])
    for mode in BoundaryMode:
        res = analyze_group(e, mode)
        assert res["connected"] is False


# 跨边界片段：同一圆柱在 x=-5000 与 x=5000 处有重合端点，
# PERIODIC_CONNECTED 应合并为一节点并导通，WRAPPED_GEOMETRY_ONLY 应不合并
def test_periodic_merge_diff():
    e = np.array([
        [-5000.0, 10.0, 20.0, -2000.0, 15.0, 25.0],
        [2000.0, 15.0, 25.0, 5000.0, 10.0, 20.0],
    ])
    pc = analyze_group(e, BoundaryMode.PERIODIC_CONNECTED)
    wo = analyze_group(e, BoundaryMode.WRAPPED_GEOMETRY_ONLY)
    assert pc["node_count"] == 1
    assert pc["connected"] is True
    assert wo["node_count"] == 2
    assert wo["connected"] is False
