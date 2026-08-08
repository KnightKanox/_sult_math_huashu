# 并查集（DSU）：用于节点合并与连通性判断
class DSU:
    # 初始化 n 个独立集合
    def __init__(self, n):
        self.parent = list(range(n))

    # 查找元素 x 的根，带路径压缩
    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    # 合并元素 a、b 所在集合
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra
