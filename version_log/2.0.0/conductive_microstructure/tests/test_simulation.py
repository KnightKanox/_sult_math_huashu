# 仿真层（问题 2）单元测试：单次仿真判定、宽相位候选不漏判、周期合并哈希化回归
import numpy as np

from src.config import BoundaryMode
from src.geometry.segment_distance import segment_segment_distance
from src.generation.medium_generator import make_cylinder_segments
from src.graph.connectivity import analyze_group
from src.graph.spatial_index import aabb_candidates
from src.io.attachment_reader import read_attachment
from src.simulation.single_trial import run_single_trial, single_trial_from_segments


# 空输入：无圆柱时不导通，且各统计量为 0
def test_zero_cylinders_not_connected():
    res = single_trial_from_segments(np.empty((0, 6)), BoundaryMode.WRAPPED_GEOMETRY_ONLY)
    assert res["connected"] is False
    assert res["node_count"] == 0
    assert res["segment_count"] == 0
    assert res["edge_count"] == 0
    assert res["max_component_ratio"] == 0.0
    assert res["mean_degree"] == 0.0


# 两根圆柱首尾相接构成左右贯通链：两种边界模式均应导通
def test_two_cylinders_chain_connected():
    # 圆柱1：中心(-2500,0,0)、方向(1,0,0) -> 盒内片段 [(-5000,0,0),(0,0,0)]，左端点接通左电极
    s1 = make_cylinder_segments((-2500.0, 0.0, 0.0), (1.0, 0.0, 0.0)).reshape(1, 6)
    # 圆柱2：中心(2499.5,0,0)、方向(1,0,0) -> 盒内片段 [(-0.5,0,0),(4999.5,0,0)]，
    # 右端点距右电极面 0.5nm（<=31.8 导通），左端点与圆柱1 端点 (0,0,0) 轴距 0.5nm（<=61.8 导通）。
    # 注：中心不能取 2500，否则右端点恰为 x=5000 会被周期回绕成 x=-5000（生成器端点值域 [-5000,5000)）
    s2 = make_cylinder_segments((2499.5, 0.0, 0.0), (1.0, 0.0, 0.0)).reshape(1, 6)
    segs = np.concatenate([s1, s2], axis=0)
    assert segs.shape == (2, 6)
    for mode in BoundaryMode:
        res = single_trial_from_segments(segs, mode)
        assert res["connected"] is True
    # WRAPPED 模式两行各为一个节点，一条链边，最大分量占比 1.0
    wo = single_trial_from_segments(segs, BoundaryMode.WRAPPED_GEOMETRY_ONLY)
    assert wo["node_count"] == 2
    assert wo["edge_count"] == 1
    assert wo["max_component_ratio"] == 1.0


# 关键回归：_build_nodes 哈希化（O(M)）后，附件三组数据的 node_count/连通结果必须与
# 原 O(M^2) 版本完全一致（组1=9、组2=39、组3=357，且全部贯通）
def test_periodic_hash_regression():
    data = read_attachment(r"d:\000AAAitaem\math\26huashu\A_try\附件.xlsx")
    expected = {"组1": (9, True), "组2": (39, True), "组3": (357, True)}
    for name, (exp_nodes, exp_conn) in expected.items():
        res = analyze_group(data[name], BoundaryMode.PERIODIC_CONNECTED)
        assert res["node_count"] == exp_nodes, f"{name}: node_count 回归失败"
        assert res["connected"] is exp_conn, f"{name}: connected 回归失败"


# 可复现性：相同种子两次单次仿真，连通性/节点数/连边数完全一致
def test_run_single_trial_reproducible():
    a = run_single_trial(np.random.default_rng(5), 20, BoundaryMode.PERIODIC_CONNECTED)
    b = run_single_trial(np.random.default_rng(5), 20, BoundaryMode.PERIODIC_CONNECTED)
    for key in ("connected", "node_count", "edge_count"):
        assert a[key] == b[key], f"{key}: 两次仿真结果不一致"


# 宽相位筛不漏判：候选对必须包含全部"真实距离 <= 阈值"的线段对（以 O(N^2) 为基准）
def test_aabb_candidates_no_miss():
    rng = np.random.default_rng(11)
    n = 50
    centers = rng.uniform(-4500.0, 4500.0, (n, 3))
    dirs = rng.normal(size=(n, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    half = 100.0
    segs = np.stack([centers - half * dirs, centers + half * dirs], axis=1)  # (50,2,3)
    # 附加一条与第 0 条轴线平行、法向距离 30nm（< 61.8）的线段，保证存在真实连接
    u0 = dirs[0]
    aux = np.array([0.0, 0.0, 1.0]) if abs(u0[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    off = aux - (aux @ u0) * u0
    off /= np.linalg.norm(off)
    c_extra = centers[0] + 30.0 * off
    extra = np.stack([c_extra - half * u0, c_extra + half * u0], axis=0)  # (2,3)
    segs = np.concatenate([segs, extra[None, :, :]], axis=0)  # (51,2,3)
    threshold = 61.8
    cand = aabb_candidates(segs, threshold)
    cand_set = set(map(tuple, cand.tolist()))
    # ground truth：O(N^2) 精确判距
    truth = set()
    ns = len(segs)
    for i in range(ns):
        for j in range(i + 1, ns):
            d = segment_segment_distance(segs[i, 0], segs[i, 1], segs[j, 0], segs[j, 1])
            if d <= threshold:
                truth.add((i, j))
    assert len(truth) >= 1, "测试构造应至少有一条真实连接"
    for pair in truth:
        assert pair in cand_set, f"漏判候选对 {pair}（真实距离 <= {threshold}）"


# 空输入：不报错，返回空 (0,2) 数组
def test_aabb_candidates_empty_input():
    out = aabb_candidates(np.empty((0, 2, 3)), 61.8)
    assert out.shape == (0, 2)
    assert np.issubdtype(out.dtype, np.integer)
