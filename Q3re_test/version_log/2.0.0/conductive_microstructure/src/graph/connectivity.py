# 图论：建图（节点=导体、虚拟电极L/R）、连通判定、贯通路径提取
from collections import deque

import numpy as np

from .dsu import DSU
from ..config import (AA_SEG_THRESHOLD, A_RADIUS, BoundaryMode,
                      LEFT_PLANE_X, RIGHT_PLANE_X, PERIODIC_MATCH_TOL)
from ..geometry.cylinder_distance import cylinders_connected, electrode_connected
from ..geometry.primitives import Node
from ..geometry.segment_distance import segment_segment_distance


# 将一段端点坐标按周期边界回绕到盒内（x' = ((x+5000) mod 10000) - 5000）
def _periodic_wrap(v):
    return ((v + 5000.0) % 10000.0) - 5000.0


# 在 PERIODIC_CONNECTED 模式下，将跨边界重合片段的行合并为同一导体节点
def _build_nodes(endpoints, mode):
    """按边界模式把原始行组织为导体节点列表。

    PERIODIC_CONNECTED：某行端点与另一行端点经周期回绕后重合（且方向一致），
    则视为同一根圆柱的跨边界片段，合并为一个节点；
    WRAPPED_GEOMETRY_ONLY：每一行独立为一个节点。
    返回 (nodes, merged_pairs)；merged_pairs 为被合并的行号对列表。
    """
    n = len(endpoints)
    p = endpoints[:, :3]
    q = endpoints[:, 3:]
    merged_pairs = []
    if mode == BoundaryMode.WRAPPED_GEOMETRY_ONLY:
        return [Node(i, [i], np.array([[p[i], q[i]]])) for i in range(n)], merged_pairs
    # PERIODIC_CONNECTED：先做合并（哈希分组 O(M)，替代原 O(M^2) 四重循环配对）
    # 对每行两端点回绕坐标 round 到 6 位小数作为键，同一键下的行视为在该端点
    # 重合的跨边界片段，全部 DSU union；round 网格 0.5μm 远小于实际重合误差（~1e-12）
    dsu = DSU(n)
    key_rows = {}
    for i in range(n):
        for v in (p[i], q[i]):
            key = tuple(float(x) for x in np.round(_periodic_wrap(v), 6))
            key_rows.setdefault(key, []).append(i)
    for rows in key_rows.values():
        if len(rows) > 1:
            root = rows[0]
            for r in rows[1:]:
                dsu.union(root, r)
    groups = {}
    for i in range(n):
        groups.setdefault(dsu.find(i), []).append(i)
    nodes = []
    for idx, row_ids in enumerate(groups.values()):
        if len(row_ids) > 1:
            merged_pairs.append(row_ids)
        segs = np.array([np.array([p[i], q[i]]) for i in row_ids])
        nodes.append(Node(idx, row_ids, segs))
    return nodes, merged_pairs


# 分析一个微构体组：返回连通性统计与一条贯通路径（若有）
def analyze_group(endpoints, mode):
    """给定一组介质A端点数组（N行×6列），在指定边界模式下判定左右电极是否导通。

    返回 dict：connected、node_count、segment_count、edge_count、
    left_node_ids、right_node_ids、path_node_ids、path_row_ids、merged_pairs。
    """
    endpoints = np.asarray(endpoints, dtype=float)
    nodes, merged_pairs = _build_nodes(endpoints, mode)
    m = len(nodes)
    # 1) 电极连接：节点内任一轴线段与左/右电极平面导通
    left_nodes, right_nodes = set(), set()
    for node in nodes:
        for seg in node.segments:
            if electrode_connected(seg[0], seg[1], radius=A_RADIUS, plane_x=LEFT_PLANE_X):
                left_nodes.add(node.index)
            if electrode_connected(seg[0], seg[1], radius=A_RADIUS, plane_x=RIGHT_PLANE_X):
                right_nodes.add(node.index)
    # 2) 导体-导体连接：跨节点任一片段对，轴线段距离 <= 2*R_A+δ
    adjacency = {i: set() for i in range(m)}
    edge_count = 0
    for a in range(m):
        for b in range(a + 1, m):
            best = np.inf
            for sa in nodes[a].segments:
                for sb in nodes[b].segments:
                    d = segment_segment_distance(sa[0], sa[1], sb[0], sb[1])
                    if d < best:
                        best = d
            if best <= AA_SEG_THRESHOLD:
                adjacency[a].add(b)
                adjacency[b].add(a)
                edge_count += 1
    # 3) DSU 判断左右导通 + BFS 提取一条贯通路径
    dsu = DSU(m)
    for a in range(m):
        for b in adjacency[a]:
            dsu.union(a, b)
    connected = any(dsu.find(ln) == dsu.find(rn) for ln in left_nodes for rn in right_nodes)
    path_node_ids = _find_path(left_nodes, right_nodes, adjacency) if connected else []
    path_row_ids = []
    for ni in path_node_ids:
        path_row_ids.append(nodes[ni].row_indices)
    return {
        "connected": connected,
        "node_count": m,
        "segment_count": len(endpoints),
        "edge_count": edge_count,
        "left_node_ids": sorted(left_nodes),
        "right_node_ids": sorted(right_nodes),
        "path_node_ids": path_node_ids,
        "path_row_ids": path_row_ids,
        "merged_pairs": merged_pairs,
    }


# 用 BFS 从左电极可达节点出发，找到一条到达右电极节点的节点序列（多源搜索）
def _find_path(left_nodes, right_nodes, adjacency):
    """多源 BFS：所有左电极直连节点同时入队，最先到达右节点时回溯重建路径。"""
    prev = {ln: None for ln in left_nodes}
    queue = deque(sorted(prev))
    target = None
    while queue:
        u = queue.popleft()
        if u in right_nodes:
            target = u
            break
        for v in sorted(adjacency[u]):
            if v not in prev:
                prev[v] = u
                queue.append(v)
    if target is None:
        return []
    path = []
    cur = target
    while cur is not None:
        path.append(cur)
        cur = prev[cur]
    return path[::-1]
