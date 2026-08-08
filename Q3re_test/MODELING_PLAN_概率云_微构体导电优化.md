# A题：微构体中填充导电介质的仿真优化
## 概率云—随机几何—连续渗流—Monte Carlo 建模与代码实现方案

> 用途：将当前建模讨论整理为可直接交给 Codex 继续实现的项目文档。  
> 建议文件名：`MODELING_PLAN.md`

---

# 0. 项目目标

本题研究一个边长为 \(10000\,\mathrm{nm}\) 的立方体微构体中的导电介质随机填充问题。

微构体内部默认充满绝缘溶剂，导电介质由以下两类构成：

- **介质 A**：直圆柱体
  - 长度：\(L_A=5000\,\mathrm{nm}\)
  - 半径：\(R_A=30\,\mathrm{nm}\)
- **介质 B**：球体
  - 半径：\(R_B=200\,\mathrm{nm}\)

选取垂直于 \(X\) 轴的左右两个面作为带电表面：

\[
x=-5000,\qquad x=5000.
\]

当两个导电介质之间，或导电介质与带电表面之间的**最短距离**

\[
d\le \delta,\qquad \delta=1.8\,\mathrm{nm}
\]

时，认为二者导通。

如果左右两个带电表面之间存在至少一条由导电介质组成的完整连接路径，则认为整个微构体导通。

---

# 1. 整体建模思想

本项目不采用“只做 Monte Carlo 暴力仿真”的方式，而采用：

\[
\boxed{
\text{几何建模}
\rightarrow
\text{局部连接概率云}
\rightarrow
\text{随机几何图}
\rightarrow
\text{连续渗流}
\rightarrow
\text{Monte Carlo 验证}
\rightarrow
\text{临界填充率 / 成本优化}
}
\]

核心思想是：

1. 用**计算几何**定义任意两介质是否导通；
2. 用**连接概率云**描述一个导体周围出现“可与其连接的另一个导体”的空间概率分布；
3. 将局部概率云压缩成**等效连接体积**；
4. 用随机几何图和连续渗流理论解释从局部连接到宏观贯通的转变；
5. 用完整三维 Monte Carlo 仿真得到有限尺寸微构体中的实际导通概率；
6. 最终解决：
   - 问题 1：给定结构是否导通；
   - 问题 2：给定填充率时的导通概率；
   - 问题 3：\(P_{\mathrm{conn}}\ge 90\%\) 时 A 的最低填充率；
   - 问题 4：\(P_{\mathrm{conn}}\ge 90\%\) 时 A/B 的最低成本组合。

---

# 2. 模型层次

建议把整个模型划分为五层。

## 2.1 几何层

研究：

\[
d(S_i,S_j)
\]

其中 \(S_i,S_j\) 表示两个导体在三维空间中占据的几何集合。

目标：

\[
d(S_i,S_j)\le 1.8
\]

则建立导通关系。

---

## 2.2 局部概率层

建立：

\[
q(\mathbf r,\mathbf u)
\]

描述固定一个导体后，另一随机姿态导体位于相对位置 \(\mathbf r\) 时与其连接的概率。

进一步对方向进行各向同性平均，得到：

\[
q(r).
\]

该函数称为：

> **径向连接概率云（Radial Connectivity Probability Cloud）**

---

## 2.3 网络层

每个导电介质视作图节点。

若两介质满足导通条件，则连接一条边。

再增加左右带电表面两个虚拟节点：

\[
L,\quad R.
\]

得到随机几何图：

\[
G=(V,E).
\]

---

## 2.4 宏观渗流层

判断：

\[
L\leftrightarrow R.
\]

若左右表面属于同一连通分量，则微构体导通。

随着填充率增加，网络会从稀疏不连通状态逐渐跨越渗流临界区域，最终形成贯通簇。

---

## 2.5 优化层

利用：

\[
P_{\mathrm{conn}}(N+1)\ge P_{\mathrm{conn}}(N)
\]

的单调性，进行：

- 二分搜索；
- 置信区间判断；
- A/B 二维可行边界搜索；
- 成本最小化。

---

# 3. 坐标系与随机变量

微构体定义为：

\[
\Omega=[-5000,5000]^3.
\]

总体积：

\[
V_0=(10000)^3=10^{12}\,\mathrm{nm}^3.
\]

---

## 3.1 介质 A

一根介质 A 表示为：

\[
A_i=(\mathbf c_i,\mathbf u_i)
\]

其中：

\[
\mathbf c_i=(x_i,y_i,z_i)
\]

为中心坐标，

\[
\mathbf u_i=(u_x,u_y,u_z),\qquad \|\mathbf u_i\|=1
\]

为方向单位向量。

轴线参数方程：

\[
\mathbf l_i(t)=\mathbf c_i+t\mathbf u_i,
\qquad
-\frac{L_A}{2}\le t\le \frac{L_A}{2}.
\]

---

## 3.2 介质 B

球体 B 表示为：

\[
B_j=(\mathbf c_j,R_B).
\]

---

## 3.3 随机位置

默认建模假设：

\[
\mathbf c_i\sim U(\Omega).
\]

---

## 3.4 随机方向

必须保证三维球面方向真正均匀。

推荐：

\[
\phi\sim U(0,2\pi),
\]

\[
z=\cos\theta\sim U(-1,1),
\]

\[
\theta=\arccos z.
\]

则：

\[
\mathbf u=
\left(
\sqrt{1-z^2}\cos\phi,
\sqrt{1-z^2}\sin\phi,
z
\right).
\]

禁止直接使用：

\[
\theta\sim U(0,\pi)
\]

因为这样不会得到球面均匀方向。

---

# 4. 边界处理

题目给出“边界截断规则”：

若导体任一部分越出微构体边界，则将越界部分沿对应方向反向平移一个微构体边长，使其从对侧重新进入。

这在几何上接近周期边界：

\[
x'=
((x+5000)\bmod 10000)-5000.
\]

\(y,z\) 同理。

---

## 4.1 当前必须保留的建模风险

由于 \(X\) 方向同时是左右带电面的方向，需特别明确：

> 一根介质 A 如果跨越 \(x=5000\)，其越界部分被映射到 \(x=-5000\) 后，映射后的两部分是否仍然视为同一个电学连续导体？

这将显著影响结果。

代码中建议设计两种边界模式：

```text
BoundaryMode.PERIODIC_CONNECTED
BoundaryMode.WRAPPED_GEOMETRY_ONLY
```

方便后续敏感性分析。

---

# 5. 几何导通判据

统一阈值：

\[
\delta=1.8\,\mathrm{nm}.
\]

---

# 6. A-A 导通

## 6.1 推荐正式模型

计算两个有限圆柱集合的真实最短距离：

\[
d_{AA}
=
\min_{\mathbf x\in A_i,\mathbf y\in A_j}
\|\mathbf x-\mathbf y\|.
\]

若：

\[
d_{AA}\le \delta
\]

则：

\[
A_i\leftrightarrow A_j.
\]

---

## 6.2 可用于高速筛选的近似模型

先计算两个圆柱轴线段之间的最短距离：

\[
d_{\mathrm{seg}}.
\]

由于两个圆柱半径均为 \(R_A\)，近似：

\[
d_{AA}
\approx
\max(0,d_{\mathrm{seg}}-2R_A).
\]

因此：

\[
d_{AA}\le1.8
\]

等价近似为：

\[
d_{\mathrm{seg}}
\le
2R_A+1.8
=
61.8\,\mathrm{nm}.
\]

注意：

该方法把圆柱近似成了 capsule / spherocylinder，在端面附近存在误差。

建议采用：

> 宽相位筛选 + 精确窄相位判断。

即：

```text
broad phase:
    轴线段距离 > 某安全阈值 -> 一定不连接

narrow phase:
    对候选对使用更精确的有限圆柱距离计算
```

---

# 7. A-B 导通

球心：

\[
\mathbf c_B.
\]

先求球心到有限圆柱 A 的最短距离：

\[
d(\mathbf c_B,A).
\]

若：

\[
d(\mathbf c_B,A)
\le
R_B+\delta
\]

则 A-B 导通。

也可先求球心到轴线段距离 \(d_{\mathrm{seg}}\)，再做：

\[
d_{\mathrm{seg}}
\le
R_A+R_B+\delta.
\]

其中：

\[
R_A+R_B+\delta
=
231.8\,\mathrm{nm}.
\]

---

# 8. B-B 导通

两球球心：

\[
\mathbf c_i,\mathbf c_j.
\]

若：

\[
\|\mathbf c_i-\mathbf c_j\|
\le
2R_B+\delta
\]

则导通。

即：

\[
\boxed{
\|\mathbf c_i-\mathbf c_j\|
\le401.8\,\mathrm{nm}
}
\]

---

# 9. 介质与带电表面的导通

左右电极：

\[
\Pi_L:x=-5000,
\qquad
\Pi_R:x=5000.
\]

对任意介质 \(S_i\)：

\[
d(S_i,\Pi_L)\le1.8
\Rightarrow
L\leftrightarrow i.
\]

\[
d(S_i,\Pi_R)\le1.8
\Rightarrow
R\leftrightarrow i.
\]

---

# 10. 问题 1：确定性图连通模型

问题 1 给出三组介质 A 的端点数据。

每根 A 的输入为：

\[
(\mathbf p_i,\mathbf q_i).
\]

则：

\[
\mathbf c_i=\frac{\mathbf p_i+\mathbf q_i}{2}
\]

\[
\mathbf u_i=
\frac{\mathbf q_i-\mathbf p_i}
{\|\mathbf q_i-\mathbf p_i\|}.
\]

---

## 10.1 建图

节点：

\[
V=
\{A_1,\dots,A_N,L,R\}.
\]

边：

\[
(A_i,A_j)\in E
\iff
d(A_i,A_j)\le1.8.
\]

以及介质—电极边。

---

## 10.2 判断

采用：

- BFS；
- DFS；
- Union-Find / DSU。

推荐 DSU，因为 Q2-Q4 中大量重复使用。

判断：

```text
find(L) == find(R)
```

即可。

---

# 11. 概率云模型：核心理论

## 11.1 为什么不定义“物质概率云”

原始想法：

\[
P(\mathbf x\text{ 被 A 占据})
\]

只描述空间占据概率。

但题目真正关心的是：

\[
d(A_i,A_j)\le1.8.
\]

因此需要定义“连接概率”，而不是单纯“物质存在概率”。

---

# 12. 姿态相关连接概率场

固定参考介质：

\[
A_0=(\mathbf 0,\mathbf u_0).
\]

另一个 A 的中心位于：

\[
\mathbf r.
\]

其方向：

\[
\mathbf u_1\sim U(S^2).
\]

定义：

\[
\boxed{
q(\mathbf r\mid\mathbf u_0)
=
P_{\mathbf u_1}
\left(
d(A_0,A_1)\le1.8
\right)
}
\]

含义：

> 如果另一根 A 的中心位于相对位置 \(\mathbf r\)，它随机取向时，与中心 A 导通的概率。

由于 A 是细长圆柱：

\[
q(\mathbf r\mid\mathbf u_0)
\]

一般是各向异性的。

---

# 13. 径向连接概率云

由于参考 A 自身方向也随机，可以进一步平均：

\[
\boxed{
q(r)
=
P
\left(
A_0\leftrightarrow A_1
\mid
\|\mathbf c_1-\mathbf c_0\|=r
\right)
}
\]

其中：

\[
r=\|\mathbf r\|.
\]

这就是最终用于理论分析的：

> **径向连接概率云**

---

# 14. 数值求概率云

不强行追求 \(q(r)\) 的完整闭式解析解。

对每个固定 \(r\)：

1. 随机生成球面方向：
   \[
   \hat{\mathbf r}
   \]
2. 令：
   \[
   \mathbf c_1=r\hat{\mathbf r}
   \]
3. 随机生成：
   \[
   \mathbf u_0,\mathbf u_1
   \]
4. 判断：
   \[
   d(A_0,A_1)\le1.8
   \]
5. 进行 \(M\) 次试验。

定义：

\[
X_k=
\begin{cases}
1,&\text{连接}\\
0,&\text{不连接}
\end{cases}
\]

则：

\[
\boxed{
\hat q(r)=\frac1M\sum_{k=1}^{M}X_k
}
\]

---

# 15. 建议的概率云实验

例如：

```text
r = 0, 25, 50, ..., r_max
```

每个距离：

```text
M_cloud = 1e4 ~ 1e5
```

得到数据：

```text
r_nm, q_hat, ci_low, ci_high
```

保存为：

```text
results/cloud_AA.csv
```

然后绘制：

```text
q_AA(r)
```

---

# 16. 等效连接体积

球壳：

\[
[r,r+dr]
\]

体积：

\[
dV=4\pi r^2dr.
\]

若 A 的空间数密度：

\[
\rho_A=\frac{N_A}{V_0},
\]

则球壳中平均包含：

\[
\rho_A4\pi r^2dr
\]

根 A。

其中比例：

\[
q(r)
\]

能够与中心 A 连接。

因此：

\[
dN_{\mathrm{conn}}
=
4\pi\rho_A r^2q(r)dr.
\]

积分：

\[
\bar k
=
4\pi\rho_A
\int_0^\infty
r^2q(r)\,dr.
\]

定义：

\[
\boxed{
V_{\mathrm{eff}}
=
4\pi
\int_0^\infty r^2q(r)\,dr
}
\]

称为：

> **等效连接体积（Effective Connectivity Volume）**

于是：

\[
\boxed{
\bar k
=
\rho_A V_{\mathrm{eff}}
}
\]

---

# 17. \(\bar k\) 的意义

\[
\bar k
\]

可以解释为：

> 随机空间中，一根介质 A 平均能够直接连接的其他介质数量。

填充率增加时：

\[
\phi\uparrow
\Rightarrow
N_A\uparrow
\Rightarrow
\rho_A\uparrow
\Rightarrow
\bar k\uparrow.
\]

因此：

\[
\bar k
\]

是连接局部概率云和宏观渗流的重要桥梁。

---

# 18. 从局部连接到宏观渗流

不能直接使用：

\[
\bar k>1
\]

作为“左右导通”的充分必要条件。

原因：

- 有有限尺寸效应；
- 图中的边存在空间相关性；
- 左右表面贯通不同于“存在巨型连通分量”；
- 电极连接也具有边界效应。

因此概率云理论负责：

1. 解释导通概率随填充率上升的机理；
2. 预测临界区间；
3. 为后续 Monte Carlo 搜索缩小范围。

最终实际概率仍由完整微构体仿真得到。

---

# 19. 完整 Monte Carlo 仿真器

对于给定：

\[
(N_A,N_B)
\]

进行一次仿真：

```text
1. 生成 N_A 根随机 A
2. 生成 N_B 个随机 B
3. 应用边界处理
4. 计算所有潜在连接关系
5. 构造图或 DSU
6. 连接左右电极虚拟节点
7. 判断 L 与 R 是否同属一个连通分量
8. 返回 True / False
```

---

# 20. Monte Carlo 导通概率

重复：

\[
M
\]

次。

定义：

\[
X_m=
\begin{cases}
1,&\text{第 }m\text{ 次导通}\\
0,&\text{否则}
\end{cases}
\]

则：

\[
\boxed{
\hat P_{\mathrm{conn}}
=
\frac1M
\sum_{m=1}^{M}X_m
}
\]

---

# 21. 置信区间

因为：

\[
X_m\sim Bernoulli(p)
\]

所以：

\[
K=\sum X_m\sim Binomial(M,p).
\]

除点估计：

\[
\hat p=\frac KM
\]

外，还必须输出：

```text
confidence_interval_low
confidence_interval_high
```

第三、四问建议主要使用：

> 95% 单侧置信下界。

若要求：

\[
P_{\mathrm{conn}}\ge0.9
\]

则推荐判定：

\[
\boxed{
p_{\mathrm{lower},95\%}\ge0.90
}
\]

而不仅仅：

\[
\hat p\ge0.90.
\]

---

# 22. 问题 2：给定体积分数计算导通概率

一个 A 的体积：

\[
V_A
=
\pi R_A^2L_A.
\]

即：

\[
V_A
=
\pi\cdot30^2\cdot5000.
\]

给定体积分数：

\[
\phi_A
\]

对应数量：

\[
\boxed{
N_A
=
\operatorname{round}
\left(
\frac{\phi_A V_0}{V_A}
\right)
}
\]

需要计算：

\[
0.50\%,\quad0.60\%,\quad0.70\%,\quad1.00\%.
\]

---

# 23. 问题 2 输出

建议输出表：

| \(\phi_A\) | \(N_A\) | \(\hat P\) | 95% CI | \(\bar k\) |
|---:|---:|---:|---:|---:|
| 0.50% | | | | |
| 0.60% | | | | |
| 0.70% | | | | |
| 1.00% | | | | |

同时绘制：

### 图 A

\[
P_{\mathrm{conn}}(\phi)
\]

### 图 B

\[
\bar k(\phi)
\]

如果：

\[
P_{\mathrm{conn}}
\]

发生明显跃迁的位置与：

\[
\bar k
\]

进入临界区间的位置接近，则说明概率云理论能有效解释宏观导通趋势。

---

# 24. 问题 3：90% 导通最低 A 填充率

优化问题：

\[
\boxed{
\min_{\phi_A}\phi_A
}
\]

满足：

\[
P_{\mathrm{conn}}(\phi_A)\ge0.90.
\]

---

# 25. 单调性

增加新的导体不会删除原来的导通边。

因此理论上：

\[
P_{\mathrm{conn}}(N_A+1)
\ge
P_{\mathrm{conn}}(N_A).
\]

所以可以采用二分搜索，而不是全区间暴力扫描。

---

# 26. 问题 3 求解策略

```text
1. 使用概率云 / Q2 初步定位临界区间
2. 设置 phi_low, phi_high
3. 取 phi_mid
4. Monte Carlo
5. 计算 95% 单侧下界
6. 若 lower_bound >= 0.90:
       high = mid
   else:
       low = mid
7. 搜索至满足 0.01% 精度要求
8. 对最终候选点增加仿真次数进行确认
```

---

# 27. 问题 4：A+B 混合体系

建立三种局部连接概率：

\[
q_{AA}(r),
\qquad
q_{AB}(r),
\qquad
q_{BB}(r).
\]

---

# 28. 三类等效连接体积

定义：

\[
V_{AA}
=
4\pi\int_0^\infty r^2q_{AA}(r)\,dr
\]

\[
V_{AB}
=
4\pi\int_0^\infty r^2q_{AB}(r)\,dr
\]

\[
V_{BB}
=
4\pi\int_0^\infty r^2q_{BB}(r)\,dr.
\]

其中 B-B 可作为解析基准。

---

# 29. B-B 理论基准

B-B 导通条件：

\[
d_{\mathrm{center}}
\le
2R_B+\delta
=
401.8.
\]

因此理想无限空间中：

\[
q_{BB}(r)
=
\begin{cases}
1,&r\le401.8\\
0,&r>401.8
\end{cases}
\]

从而：

\[
\boxed{
V_{BB}
=
\frac43\pi(401.8)^3
}
\]

这可以用于校验概率云数值积分代码。

---

# 30. A-B 近似理论基准

粗略近似连接半径：

\[
R_{AB}=R_A+R_B+\delta
=
231.8.
\]

对于长圆柱，可以近似把连接区域看成一个 capsule：

\[
V_{AB}^{approx}
=
\pi R_{AB}^2L_A
+
\frac43\pi R_{AB}^3.
\]

该式可作为 Monte Carlo 得到 \(V_{AB}\) 的理论 sanity check。

---

# 31. 双类型连接矩阵

定义数密度：

\[
\rho_A=\frac{N_A}{V_0},
\qquad
\rho_B=\frac{N_B}{V_0}.
\]

建立：

\[
\boxed{
M=
\begin{pmatrix}
\rho_AV_{AA} & \rho_BV_{AB}\\
\rho_AV_{AB} & \rho_BV_{BB}
\end{pmatrix}
}
\]

矩阵元素表示不同介质类型之间的平均局部连接能力。

---

# 32. 最大特征值指标

计算：

\[
\lambda_{\max}(M).
\]

将其作为：

> **混合体系宏观渗流能力指标**

用途：

- 快速比较不同 \(N_A,N_B\) 组合；
- 预测可行区边界；
- 缩小 Monte Carlo 搜索范围。

注意：

\[
\lambda_{\max}
\]

不是最终 90% 导通概率的精确判据。

它是理论筛选指标。

---

# 33. 问题 4 成本函数

A 成本：

\[
c_A=1.05\;\mathrm{元}/\mu m^3.
\]

B 成本：

\[
c_B=0.05\;\mathrm{元}/\mu m^3.
\]

注意单位换算：

\[
1\,\mu m^3=10^9\,nm^3.
\]

于是：

\[
V_A^{(\mu m^3)}
=
\frac{V_A^{(nm^3)}}{10^9}
\]

\[
V_B^{(\mu m^3)}
=
\frac{\frac43\pi R_B^3}{10^9}.
\]

总成本：

\[
\boxed{
C(N_A,N_B)
=
c_AN_AV_A
+
c_BN_BV_B
}
\]

---

# 34. 问题 4 优化模型

\[
\boxed{
\begin{aligned}
\min_{N_A,N_B}\quad&
C(N_A,N_B)
\\
\mathrm{s.t.}\quad&
P_{\mathrm{conn}}(N_A,N_B)\ge0.90
\\
&
N_A,N_B\in\mathbb Z_{\ge0}.
\end{aligned}
}
\]

这是一个：

> **机会约束整数优化问题**

---

# 35. 问题 4 搜索方法

不能全二维暴力扫描。

利用：

\[
P(N_A+1,N_B)\ge P(N_A,N_B)
\]

和：

\[
P(N_A,N_B+1)\ge P(N_A,N_B).
\]

因此可行域具有单调性。

建议：

```text
1. 用连接矩阵 lambda_max 快速扫描理论可行区域
2. 对若干 N_A：
       二分寻找最小 N_B
       使 Monte Carlo 的 95% 单侧置信下界 >= 0.90
3. 得到一条“90% 可行边界”
4. 只比较边界上的成本
5. 找到全局最低成本候选
6. 对最优点及其邻域进行高次数 Monte Carlo 验证
```

---

# 36. 理论与仿真的职责分工

| 模块 | 作用 |
|---|---|
| 计算几何 | 给出真实局部导通判据 |
| 概率云 | 描述随机局部连接规律 |
| 等效连接体积 | 将空间概率场压缩为单个连接尺度 |
| 平均连接度 | 描述节点平均局部连接能力 |
| 连续渗流 | 解释从局部连接到宏观贯通的相变 |
| 随机几何图 | 表示单次微构体内部的连接网络 |
| Monte Carlo | 计算真实有限结构的导通概率 |
| 置信区间 | 判断概率估计是否足够可信 |
| 二分搜索 | 求临界填充率 |
| 多类型连接矩阵 | 描述 A/B 混合体系 |
| 整数优化 | 求最小成本 |

---

# 37. 推荐代码架构

建议项目：

```text
conductive_microstructure/
│
├── README.md
├── MODELING_PLAN.md
├── requirements.txt
│
├── data/
│   ├── raw/
│   │   └── problem1.xlsx
│   └── processed/
│
├── src/
│   ├── __init__.py
│   │
│   ├── config.py
│   │
│   ├── geometry/
│   │   ├── __init__.py
│   │   ├── primitives.py
│   │   ├── segment_distance.py
│   │   ├── cylinder_distance.py
│   │   ├── sphere_distance.py
│   │   ├── electrode_distance.py
│   │   └── periodic_boundary.py
│   │
│   ├── generation/
│   │   ├── __init__.py
│   │   ├── random_position.py
│   │   ├── random_orientation.py
│   │   └── medium_generator.py
│   │
│   ├── graph/
│   │   ├── __init__.py
│   │   ├── dsu.py
│   │   ├── connectivity.py
│   │   └── spatial_index.py
│   │
│   ├── cloud/
│   │   ├── __init__.py
│   │   ├── aa_cloud.py
│   │   ├── ab_cloud.py
│   │   ├── bb_cloud.py
│   │   └── effective_volume.py
│   │
│   ├── simulation/
│   │   ├── __init__.py
│   │   ├── single_trial.py
│   │   ├── monte_carlo.py
│   │   └── confidence.py
│   │
│   ├── optimization/
│   │   ├── __init__.py
│   │   ├── problem3_search.py
│   │   ├── problem4_search.py
│   │   └── cost.py
│   │
│   └── visualization/
│       ├── cloud_plot.py
│       ├── probability_plot.py
│       └── phase_plot.py
│
├── scripts/
│   ├── solve_problem1.py
│   ├── build_aa_cloud.py
│   ├── solve_problem2.py
│   ├── solve_problem3.py
│   ├── solve_problem4.py
│   └── benchmark.py
│
├── tests/
│   ├── test_orientation.py
│   ├── test_segment_distance.py
│   ├── test_bb_distance.py
│   ├── test_boundary.py
│   ├── test_dsu.py
│   ├── test_cloud_bb.py
│   └── test_monotonicity.py
│
└── results/
    ├── tables/
    ├── cloud/
    ├── figures/
    └── logs/
```

---

# 38. 推荐核心数据结构

## 38.1 CylinderA

```python
@dataclass
class CylinderA:
    center: np.ndarray
    direction: np.ndarray
    length: float = 5000.0
    radius: float = 30.0
```

---

## 38.2 SphereB

```python
@dataclass
class SphereB:
    center: np.ndarray
    radius: float = 200.0
```

---

## 38.3 SimulationConfig

```python
@dataclass
class SimulationConfig:
    box_size: float = 10000.0
    connect_threshold: float = 1.8

    a_length: float = 5000.0
    a_radius: float = 30.0
    b_radius: float = 200.0

    boundary_mode: str = "periodic_connected"

    seed: int = 42
```

---

# 39. Codex 第一阶段任务：先不要优化性能

目标：

> 建立一个**正确但可以较慢**的 reference implementation。

优先完成：

```text
[ ] 随机方向生成器
[ ] A / B 数据结构
[ ] 线段-线段距离
[ ] A-A 近似距离
[ ] A-B 距离
[ ] B-B 距离
[ ] 导体-左右平面距离
[ ] 周期边界处理
[ ] DSU
[ ] 单次微构体导通判断
[ ] Monte Carlo 外层循环
```

这一阶段优先保证几何逻辑正确。

---

# 40. Codex 第二阶段任务：问题 1

实现：

```text
scripts/solve_problem1.py
```

需求：

```text
输入：
    附件三个 sheet

输出：
    组1：导通 / 不导通
    组2：导通 / 不导通
    组3：导通 / 不导通

附加输出：
    - 总介质数量
    - 总连接边数
    - 左电极直接连接节点数
    - 右电极直接连接节点数
    - 一条实际贯通路径（若存在）
```

“实际贯通路径”用于论文可视化和人工验证。

---

# 41. Codex 第三阶段任务：概率云

首先完成：

```text
build_aa_cloud.py
```

参数：

```text
r_min
r_max
r_step
samples_per_r
seed
```

输出：

```text
r_nm
q_hat
success_count
sample_count
ci_low
ci_high
```

再实现：

\[
V_{AA}
=
4\pi\int r^2q(r)\,dr
\]

数值积分推荐：

```python
np.trapz(...)
```

或更高精度积分。

---

# 42. 概率云必须做的验证

## Test 1：远距离

当：

\[
r
\]

大于两根 A 最大可能接触的距离后：

\[
q_{AA}(r)=0.
\]

---

## Test 2：B-B

Monte Carlo 的：

\[
q_{BB}(r)
\]

必须接近阶跃函数：

\[
q_{BB}(r)
=
I(r\le401.8).
\]

---

## Test 3：积分

数值得到：

\[
V_{BB}^{MC}
\]

应接近解析：

\[
\frac43\pi(401.8)^3.
\]

---

# 43. Codex 第四阶段任务：问题 2

实现：

```text
scripts/solve_problem2.py
```

输入：

```text
volume_fraction = [0.005, 0.006, 0.007, 0.010]
trials
seed
```

输出 CSV：

```text
phi
N_A
connected
trials
p_hat
ci_low
ci_high
rho
mean_degree_theory
```

---

# 44. Codex 第五阶段任务：问题 3

实现：

```text
scripts/solve_problem3.py
```

算法：

```text
输入：
    lower_phi
    upper_phi
    target_probability = 0.90
    accuracy_phi = 0.0001   # 对应百分数小数点后两位时再确认单位
    confidence = 0.95

执行：
    二分搜索
    每个点运行 Monte Carlo
    使用单侧置信下界判定

输出：
    最低填充率候选
    对应 N_A
    p_hat
    lower_bound
    trials
```

最后候选需要单独高精度复验。

---

# 45. Codex 第六阶段任务：问题 4

实现：

```text
scripts/solve_problem4.py
```

两级搜索：

## Level 1：理论筛选

计算：

\[
\lambda_{\max}(M)
\]

快速扫描：

```text
N_A × N_B
```

区域。

---

## Level 2：Monte Carlo 精确搜索

固定 \(N_A\)：

```text
binary search minimum N_B
```

满足：

\[
p_{lower}\ge0.9.
\]

保存每个 \(N_A\) 对应的最低可行 \(N_B\)。

---

## Level 3：成本比较

计算：

\[
C(N_A,N_B)
\]

找：

\[
\arg\min C.
\]

---

# 46. 性能优化阶段

当 reference implementation 正确后再优化。

---

## 46.1 不允许直接 \(O(N^2)\) 暴力到最后

对于 \(N\sim500-1000\)，单次 \(O(N^2)\) 尚可测试，但大量 Monte Carlo 会很慢。

后续需要：

- uniform grid；
- spatial hashing；
- cell list；
- KD-tree（仅作候选筛选）；
- bounding sphere / AABB broad phase。

---

## 46.2 推荐空间哈希

将盒子划分为 cell。

每根导体记录其 bounding box 覆盖的 cell。

只对存在公共或邻近 cell 的导体对做精确几何判断。

---

## 46.3 并行

Monte Carlo 各 trial 独立：

```text
trial_1
trial_2
...
trial_M
```

天然适合：

- multiprocessing；
- joblib；
- concurrent.futures；
- numba 并行；
- 后期 C++。

---

# 47. 可复现实验原则

所有实验必须记录：

```text
random_seed
git_commit
simulation_config
boundary_mode
number_of_trials
confidence_method
geometry_method
```

建议每次实验自动生成：

```text
results/logs/run_xxx.json
```

---

# 48. 推荐统计量

除导通概率外，建议同时保存：

```text
largest_component_size
number_of_components
mean_degree
max_degree
left_electrode_degree
right_electrode_degree
shortest_path_length_if_connected
```

用于解释渗流相变。

---

# 49. 推荐论文图

至少生成以下图：

## Figure 1
介质 A 的局部连接概率云：

\[
q_{AA}(r).
\]

## Figure 2
概率云积分：

\[
4\pi r^2q(r).
\]

该曲线的面积即：

\[
V_{\mathrm{eff}}.
\]

## Figure 3

\[
P_{\mathrm{conn}}-\phi.
\]

## Figure 4

\[
\bar k-\phi.
\]

## Figure 5
最大连通分量比例：

\[
S_{\max}/N-\phi.
\]

## Figure 6
A/B 组合的：

\[
\lambda_{\max}
\]

热图。

## Figure 7
A/B 组合的 Monte Carlo：

\[
P_{\mathrm{conn}}
\]

热图。

## Figure 8
成本等高线 + \(P=0.9\) 可行边界。

---

# 50. 模型验证思路

## 50.1 几何验证

构造人工样例：

```text
两个明显分离的 A
两个明显相交的 A
端面对端
平行侧面
A 与左墙接触
B-B 正好距离 401.8
B-B 距离 401.81
```

---

## 50.2 统计验证

不同 seed：

```text
seed = 1, 2, 3, ..., 10
```

比较导通概率波动。

---

## 50.3 收敛性验证

绘制：

\[
\hat P(M)
\]

随：

\[
M
\]

变化。

例如：

```text
M = 100
500
1000
2000
5000
10000
```

---

## 50.4 理论—仿真验证

比较：

\[
\bar k(\phi)
\]

与：

\[
P_{\mathrm{conn}}(\phi).
\]

理论要求：

> \(\bar k\) 增大并进入临界区时，宏观导通概率应同步出现明显跃迁趋势。

不要求二者存在简单一一对应闭式关系。

---

# 51. 模型假设必须明确写入论文

推荐假设：

1. 所有介质位置独立均匀分布；
2. 介质 A 姿态服从球面各向同性分布；
3. 介质之间允许重叠和贯穿；
4. 不考虑重力、运动和动态重排；
5. 最短距离不超过 \(1.8\,\mathrm{nm}\) 即视为理想导通；
6. 同一介质内部视为整体等势导体；
7. 除题目指定边界处理外，不引入额外表面效应；
8. Monte Carlo 试验之间相互独立。

其中第 1、2 条属于为“导通概率”建立明确概率测度所作的建模假设。

---

# 52. 不应在论文中做出的过强结论

避免写：

> \(\bar k>1\) 时微构体一定导通。

避免写：

> 概率云可以直接解析得到左右贯通概率。

避免写：

> \(\lambda_{\max}>1\) 等价于 90% 导通。

正确表述应是：

> 概率云与连接矩阵用于描述局部连接能力并预测渗流趋势；有限尺寸微构体的实际左右贯通概率最终由 Monte Carlo 仿真确定。

---

# 53. 最重要的理论逻辑

最终整篇论文要反复保持这条链：

\[
\boxed{
\text{随机位置 + 随机姿态}
}
\]

\[
\Downarrow
\]

\[
\boxed{
q(\mathbf r,\mathbf u)
}
\]

姿态相关连接概率场

\[
\Downarrow
\]

\[
\boxed{
q(r)
}
\]

径向连接概率云

\[
\Downarrow
\]

\[
\boxed{
V_{\mathrm{eff}}
=
4\pi\int r^2q(r)dr
}
\]

等效连接体积

\[
\Downarrow
\]

\[
\boxed{
\bar k
=
\rho V_{\mathrm{eff}}
}
\]

平均连接能力

\[
\Downarrow
\]

\[
\boxed{
\text{随机几何图 / 连续渗流}
}
\]

\[
\Downarrow
\]

\[
\boxed{
P_{\mathrm{conn}}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{Monte Carlo 验证}
}
\]

\[
\Downarrow
\]

\[
\boxed{
\text{临界填充率 / 最低成本优化}
}
\]

---

# 54. 四问最终映射

| 问题 | 数学模型 | 数值方法 |
|---|---|---|
| Q1 | 计算几何 + 图论 | DSU/BFS |
| Q2 | 概率云 + 连续渗流 | Monte Carlo |
| Q3 | 单调概率约束 | 二分 + Monte Carlo + CI |
| Q4 | 多类型概率云 + 连接矩阵 + 机会约束优化 | 边界搜索 + Monte Carlo |

---

# 55. Codex 推荐执行顺序

不要让 Codex 一次性写完整项目。

推荐逐阶段执行。

---

## Phase 1：基础几何

要求 Codex：

```text
严格按照 MODELING_PLAN.md，实现 src/geometry 与 tests/ 中的基础几何函数。
暂时不要实现 Monte Carlo 与优化。
优先保证单元测试与几何正确性。
```

---

## Phase 2：Q1

```text
实现问题1的数据读取、建图、DSU连通判断。
输出三组结果和一条贯通路径。
```

---

## Phase 3：概率云

```text
实现 AA / AB / BB 二体连接概率云仿真。
必须先通过 BB 解析解校验。
```

---

## Phase 4：Q2

```text
实现完整微构体 Monte Carlo。
先对题目四个指定体积分数计算概率。
输出置信区间。
```

---

## Phase 5：Q3

```text
利用单调性和二分搜索寻找 90% 导通最低 A 填充率。
使用单侧置信下界，而非仅使用点估计。
```

---

## Phase 6：Q4

```text
实现双介质概率云、连接矩阵、理论筛选以及最低成本边界搜索。
```

---

## Phase 7：性能优化

```text
在已有正确 reference implementation 基础上加入空间索引、
并行 Monte Carlo 和向量化。
不得以牺牲几何正确性为代价。
```

---

# 56. 可以直接给 Codex 的总 Prompt

```text
你现在要基于仓库中的 MODELING_PLAN.md 实现“微构体导电介质仿真优化”项目。

原则：

1. MODELING_PLAN.md 是当前模型设计的唯一主规范。
2. 不要一次性实现所有功能，按照 Phase 1 -> Phase 7 分阶段完成。
3. 当前优先保证计算几何正确性，不要提前进行过度性能优化。
4. 所有随机过程必须支持固定 seed。
5. 三维方向必须采用球面均匀采样，不允许 theta~Uniform(0,pi)。
6. 所有核心函数必须有单元测试。
7. Monte Carlo 必须同时返回点估计与置信区间。
8. 问题3、问题4中，不允许仅以 p_hat >= 0.9 作为最终可行判据；
   应实现 95% 单侧置信下界。
9. 保留两种边界解释模式，便于后续敏感性分析。
10. 每完成一个 Phase，先运行测试并汇报：
    - 新增文件
    - 核心接口
    - 测试结果
    - 尚未解决的问题
    再继续下一阶段。

现在只执行 Phase 1。
```

---

# 57. 当前尚未解决、后续必须重点讨论的问题

## 高优先级

### 1. 边界映射后的同一导体是否保持跨边界电学连续？

这是当前最重要的题意解释问题。

### 2. A-A 精确有限圆柱距离

需要决定：

- 采用 capsule 近似；
- 还是实现有限圆柱真实几何距离。

推荐：

> 先用 capsule reference model 跑通全流程，再实现精确距离作为修正并比较两者结果。

### 3. 连接概率云是否需要考虑有限盒子边界？

建议分两种：

- 理论云：无限均匀空间；
- 实际 Monte Carlo：真实有限盒子 + 题目边界规则。

这样理论结构最清晰。

### 4. 90% 判据的置信区间方法

推荐后续比较：

- Wilson；
- Clopper-Pearson；
- Beta posterior credible interval（仅作为补充）。

最终优先采用频率学派单侧区间。

---

# 58. 当前项目核心结论

本项目的理论基础不是“模仿电子云”本身。

“电子概率云”仅作为直观类比。

正式数学表述应为：

> 基于随机几何建立条件连接概率场，通过姿态平均得到径向连接概率云；进一步由连接概率云积分得到等效连接体积，以平均局部连接度刻画随机导体网络的连接能力，并结合连续渗流理论解释宏观贯通行为。对于有限尺寸、周期边界和多体相关性等难以解析处理的因素，采用 Monte Carlo 三维仿真计算实际导通概率，并进一步完成临界填充率与最低成本优化。

---

# 59. 一句话项目主线

\[
\boxed{
\text{局部几何连接}
\rightarrow
\text{概率云}
\rightarrow
\text{等效连接体积}
\rightarrow
\text{渗流}
\rightarrow
\text{Monte Carlo}
\rightarrow
\text{优化}
}
\]

该主线应作为后续代码、实验、论文和答辩统一的逻辑框架。
