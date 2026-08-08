# 单次仿真（问题 2）：随机生成的圆柱片段 -> 建图 -> 左右电极贯通判定
# 导体连接使用 AABB 宽相位筛候选对 + 精确判距，避免 O(M^2) 全对判距，
# 支撑 N_A≈700（M 可达 1500-2500）下数千次 Monte Carlo 仿真的性能需求。
import numpy as np

from ..config import (AA_SEG_THRESHOLD, AB_SEG_THRESHOLD, A_RADIUS,
                      BB_CENTER_THRESHOLD, BE_AXIS_THRESHOLD, BoundaryMode,
                      LEFT_PLANE_X, RIGHT_PLANE_X)
from ..generation.medium_generator import generate_b_spheres, generate_batch
from ..geometry.cylinder_distance import electrode_connected
from ..geometry.segment_distance import point_segment_distance, segment_segment_distance
from ..graph.connectivity import _build_nodes
from ..graph.dsu import DSU
from ..graph.spatial_index import (aabb_candidates, aabb_candidates_fast,
                                   aabb_candidates_mixed, ball_candidates)


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


# 问题4：一次完整 Monte Carlo 单次仿真：随机生成 n_a 根 A 圆柱 + n_b 个 B 球后做贯通判定
def run_single_trial_mixed(rng, n_a, n_b, mode):
    """随机生成 n_a 根介质A圆柱与 n_b 个介质B球并判定左右电极是否贯通，返回单次仿真统计 dict。

    A 节点由 _build_nodes 按边界模式构建（PERIODIC 模式自动合并跨边界片段），
    B 球每球一个节点（索引偏移 m_a）；三种导体连边（A-A/A-B/B-B）均先经
    AABB 宽相位筛候选、再精确判距；n_a=0 或 n_b=0 均正常返回。
    返回字段与 single_trial_from_segments 完全一致。
    """
    segments = generate_batch(rng, n_a)  # (M,6)
    spheres = generate_b_spheres(rng, n_b)  # (N,3)
    # 1) A 部分：构建 A 节点并拼接全部轴线段与 seg->node 映射
    if len(segments) > 0:
        nodes_a, _merged_pairs = _build_nodes(segments, mode)
        m_a = len(nodes_a)
        all_segs_list, seg_to_node_list = [], []
        for node in nodes_a:
            for k in range(len(node.segments)):
                all_segs_list.append(node.segments[k])
                seg_to_node_list.append(node.index)
        all_segs = np.array(all_segs_list, dtype=float)
        seg_to_node = np.array(seg_to_node_list, dtype=np.int64)
    else:
        nodes_a, m_a = [], 0
        all_segs = np.empty((0, 2, 3))
        seg_to_node = np.empty((0,), dtype=np.int64)
    # 2) 总节点数：A 节点在前、B 球节点在后，建立 DSU
    m = m_a + len(spheres)
    node_dsu = DSU(m)
    # 3) 电极连接：A 圆柱节点按轴线段、B 球按球心分别登记左右电极节点
    left_nodes, right_nodes = set(), set()
    for node in nodes_a:
        for seg in node.segments:
            if electrode_connected(seg[0], seg[1], radius=A_RADIUS, plane_x=LEFT_PLANE_X):
                left_nodes.add(node.index)
            if electrode_connected(seg[0], seg[1], radius=A_RADIUS, plane_x=RIGHT_PLANE_X):
                right_nodes.add(node.index)
    # B 球电极判定向量化：|plane_x - c_x| <= R_B+δ（BE_AXIS_THRESHOLD）即导通
    if len(spheres):
        left_mask = np.abs(spheres[:, 0] - LEFT_PLANE_X) <= BE_AXIS_THRESHOLD
        right_mask = np.abs(spheres[:, 0] - RIGHT_PLANE_X) <= BE_AXIS_THRESHOLD
        left_nodes.update((m_a + np.nonzero(left_mask)[0]).tolist())
        right_nodes.update((m_a + np.nonzero(right_mask)[0]).tolist())
    # 4) 导体连接：AABB 宽相位筛候选 + 精确判距，统计跨节点成功连边数
    edge_count = 0
    # A-A：轴线段距离 <= AA_SEG_THRESHOLD，跳过同一节点内的片段对（快速宽相位：AABB+包含直线预筛）
    for i, j in aabb_candidates_fast(all_segs, AA_SEG_THRESHOLD):
        if seg_to_node[i] != seg_to_node[j]:
            d = segment_segment_distance(all_segs[i, 0], all_segs[i, 1],
                                         all_segs[j, 0], all_segs[j, 1])
            if d <= AA_SEG_THRESHOLD:
                ni, nj = seg_to_node[i], seg_to_node[j]
                if node_dsu.find(ni) != node_dsu.find(nj):
                    node_dsu.union(ni, nj)
                    edge_count += 1
    # A-B：球心到轴线段距离 <= AB_SEG_THRESHOLD
    for i, j in aabb_candidates_mixed(all_segs, spheres, AB_SEG_THRESHOLD):
        d = point_segment_distance(spheres[j], all_segs[i, 0], all_segs[i, 1])
        if d <= AB_SEG_THRESHOLD:
            ni, nj = seg_to_node[i], m_a + j
            if node_dsu.find(ni) != node_dsu.find(nj):
                node_dsu.union(ni, nj)
                edge_count += 1
    # B-B：球心距 <= BB_CENTER_THRESHOLD（ball_candidates 用 cKDTree 精确枚举，无需再复核）
    for i, j in ball_candidates(spheres, BB_CENTER_THRESHOLD):
        ni, nj = m_a + i, m_a + j
        if node_dsu.find(ni) != node_dsu.find(nj):
            node_dsu.union(ni, nj)
            edge_count += 1
    # 5) 左右电极贯通判定 + 最大连通分量占比
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
