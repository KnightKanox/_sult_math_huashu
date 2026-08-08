# 单次仿真（问题 2）：随机生成的圆柱片段 -> 建图 -> 左右电极贯通判定
# 导体连接使用 AABB 宽相位筛候选对 + 精确判距，避免 O(M^2) 全对判距，
# 支撑 N_A≈700（M 可达 1500-2500）下数千次 Monte Carlo 仿真的性能需求。
import numpy as np

from ..config import (AA_SEG_THRESHOLD, A_RADIUS, BoundaryMode,
                      LEFT_PLANE_X, RIGHT_PLANE_X)
from ..generation.medium_generator import generate_batch
from ..geometry.cylinder_distance import electrode_connected
from ..geometry.segment_distance import segment_segment_distance
from ..graph.connectivity import _build_nodes
from ..graph.dsu import DSU
from ..graph.spatial_index import aabb_candidates


# 由"在箱片段"行数组完成一次完整建图与左右电极贯通判定
def single_trial_from_segments(segments, mode):
    """输入 (M,6) 片段数组（行=[p,q]），返回连通性统计 dict。

    步骤：_build_nodes 建节点（PERIODIC 模式自动合并跨边界片段）-> 节点级 DSU ->
    电极导通收集 left/right 节点集 -> AABB 宽相位筛候选对并精确判距连边 ->
    判定左右电极是否处于同一连通分量；统计 node/segment/edge 数与最大分量占比。
    """
    segments = np.asarray(segments, dtype=float)
    nodes, _merged_pairs = _build_nodes(segments, mode)
    m = len(nodes)
    node_dsu = DSU(m)
    # 电极连接：节点内任一轴线段与左/右电极平面导通即登记
    left_nodes, right_nodes = set(), set()
    for node in nodes:
        for seg in node.segments:
            if electrode_connected(seg[0], seg[1], radius=A_RADIUS, plane_x=LEFT_PLANE_X):
                left_nodes.add(node.index)
            if electrode_connected(seg[0], seg[1], radius=A_RADIUS, plane_x=RIGHT_PLANE_X):
                right_nodes.add(node.index)
    # 导体连接：拼接全部轴线段 (N,2,3) 并记录 seg -> node 映射
    all_segs_list, seg_to_node = [], []
    for node in nodes:
        for k in range(len(node.segments)):
            all_segs_list.append(node.segments[k])
            seg_to_node.append(node.index)
    all_segs = np.array(all_segs_list) if all_segs_list else np.empty((0, 2, 3))
    seg_to_node = np.array(seg_to_node, dtype=np.int64)
    edge_count = 0
    for i, j in aabb_candidates(all_segs, AA_SEG_THRESHOLD):
        # 跳过同一节点内的片段对，仅统计跨节点成功连边
        if seg_to_node[i] != seg_to_node[j]:
            d = segment_segment_distance(all_segs[i, 0], all_segs[i, 1],
                                         all_segs[j, 0], all_segs[j, 1])
            if d <= AA_SEG_THRESHOLD:
                ni, nj = seg_to_node[i], seg_to_node[j]
                if node_dsu.find(ni) != node_dsu.find(nj):
                    node_dsu.union(ni, nj)
                    edge_count += 1
    # 左右电极贯通判定 + 最大连通分量占比
    connected = any(node_dsu.find(ln) == node_dsu.find(rn)
                    for ln in left_nodes for rn in right_nodes)
    comp_sizes = {}
    for i in range(m):
        r = node_dsu.find(i)
        comp_sizes[r] = comp_sizes.get(r, 0) + 1
    max_comp = max(comp_sizes.values()) if comp_sizes else 0
    return {
        "connected": bool(connected),
        "node_count": m,
        "segment_count": len(segments),
        "edge_count": edge_count,
        "left_node_count": len(left_nodes),
        "right_node_count": len(right_nodes),
        "max_component_ratio": float(max_comp / m) if m > 0 else 0.0,
        "mean_degree": float(2.0 * edge_count / m) if m > 0 else 0.0,
    }


# 一次完整 Monte Carlo 单次仿真：随机生成 N_A 根圆柱后做贯通判定
def run_single_trial(rng, n_a, mode):
    """随机生成 n_a 根介质A圆柱并判定左右电极是否贯通，返回单次仿真统计 dict。

    generate_batch 输出 (M,6)；n_a=0 时返回空数组，直接给出全 0 统计。
    """
    segments = generate_batch(rng, n_a)
    if len(segments) == 0:
        return {
            "connected": False,
            "node_count": 0,
            "segment_count": 0,
            "edge_count": 0,
            "left_node_count": 0,
            "right_node_count": 0,
            "max_component_ratio": 0.0,
            "mean_degree": 0.0,
        }
    return single_trial_from_segments(segments, mode)
