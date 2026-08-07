# 并查集单元测试
from src.graph.dsu import DSU


# 基础合并与查找
def test_basic_union_find():
    d = DSU(5)
    d.union(0, 1)
    d.union(1, 2)
    assert d.find(0) == d.find(2)
    assert d.find(0) != d.find(3)


# 传递闭包：0-1-2-3-4 全连通
def test_transitive_connectivity():
    d = DSU(5)
    for a, b in [(0, 1), (1, 2), (2, 3), (3, 4)]:
        d.union(a, b)
    assert all(d.find(0) == d.find(i) for i in range(5))


# 独立集合不受影响
def test_independent_sets():
    d = DSU(6)
    d.union(0, 1)
    d.union(3, 4)
    assert d.find(0) == d.find(1)
    assert d.find(3) == d.find(4)
    assert d.find(0) != d.find(3)
    assert d.find(2) != d.find(3)
    assert d.find(5) == 5
