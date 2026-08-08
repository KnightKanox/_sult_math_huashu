# 概率云层（混合介质）：A 圆柱与 B 球二体连接概率的蒙特卡洛采样估计与等效连接体积
# A 为直圆柱（长 A_LENGTH、半径 A_RADIUS），B 为球（半径 B_RADIUS）；
# A-B 判据：中心轴线段到球心距离 ≤ R_AB = R_A+R_B+δ；B-B 判据：球心距 ≤ 2R_B+δ；
# V_AB 由 4π∫r²·q_AB(r)dr 数值积分，另有 V_BB 解析值与 V_AB 胶囊近似理论值
import numpy as np

from src.config import (A_LENGTH, A_RADIUS, B_RADIUS, CONNECT_DELTA,
                        COST_A_PER_UM3, COST_B_PER_UM3)
from src.geometry.segment_distance import point_segment_distance
from src.generation.random_orientation import random_unit_vector

# 模块级几何/成本常量：A 轴半长、A-B 与 B-B 连接距离阈值、盒体积与单介质体积（nm³）
A_HALF = A_LENGTH / 2.0
R_AB = A_RADIUS + B_RADIUS + CONNECT_DELTA
BB_THRESHOLD = 2.0 * B_RADIUS + CONNECT_DELTA
BOX_VOLUME = 1e12
CYLINDER_VOLUME = np.pi * A_RADIUS ** 2 * A_LENGTH
BALL_VOLUME = 4.0 / 3.0 * np.pi * B_RADIUS ** 3


# 返回 A-B 可发生连接的最大中心距：A 轴半长与 R_AB 之和（r 超过该值 q_AB=0）
def ab_max_contact_center_distance():
    """返回 A 圆柱与 B 球可连接的最大中心距（= A_LENGTH/2 + R_AB，nm）。

    A 轴线段过原点且长度 A_LENGTH，线段上任意点到球心距离 ≥ r - A_HALF；
    当 r > A_HALF + R_AB 时该距离恒 > R_AB，故 q_AB=0。
    """
    return A_HALF + R_AB


# 估计 A-B 在相对距离 r 处的导通概率：球面均匀采样 A 轴方向与 B 球心方向并统计判据命中比例
def estimate_q_ab(r, rng, samples_per_r=2000):
    """估计 A 圆柱与 B 球在中心距 r（nm）下的连接概率 q_AB(r)。

    A 圆柱以原点为中心、方向 u 球面均匀，轴线段为 [-A_HALF*u, A_HALF*u]；
    B 球心位于 r*r_hat（r_hat 球面均匀）；判据为
    point_segment_distance(r*r_hat, -A_HALF*u, A_HALF*u) <= R_AB。
    对 samples_per_r 个样本计数导通比例，返回 q_hat（float，∈[0,1]）。
    """
    success_count = 0
    for _ in range(samples_per_r):
        u = random_unit_vector(rng)
        r_hat = random_unit_vector(rng)
        d = point_segment_distance(r * r_hat, -A_HALF * u, A_HALF * u)
        if d <= R_AB:
            success_count += 1
    return success_count / float(samples_per_r)


# 解析计算 B-B 导通概率：球心距不超过 2*B_RADIUS+δ 时恒为 1，否则为 0（与方向无关）
def q_bb(r):
    """返回 B-B 导通概率 q_BB(r)（解析）= 1.0 if r <= BB_THRESHOLD else 0.0。"""
    return 1.0 if r <= BB_THRESHOLD else 0.0


# 数值积分估计 A-B 等效连接体积：在 [0, r_max] 网格逐点采样 q_AB 后按 4π∫r²q dr 梯形积分
def estimate_v_ab(r_max=None, dr=25.0, rng=None, samples_per_r=2000):
    """估计 A-B 等效连接体积 V_AB = 4π∫₀^r_max r²·q_AB(r) dr（nm³）。

    r 从 0 到 r_max 以步长 dr 取网格（与 build_aa_cloud 相同约定），
    每点调用 estimate_q_ab 采样 samples_per_r 次估计 q_AB(r)，最后 np.trapezoid 积分。
    r_max 默认取 ab_max_contact_center_distance()；rng 为 None 时新建默认 Generator。
    """
    if r_max is None:
        r_max = ab_max_contact_center_distance()
    if rng is None:
        rng = np.random.default_rng()
    r_grid = np.arange(0.0, r_max + dr * 0.5, dr)
    q_grid = np.array([estimate_q_ab(float(r), rng, samples_per_r) for r in r_grid])
    return 4.0 * np.pi * np.trapezoid(r_grid ** 2 * q_grid, r_grid)


# 解析计算 B-B 等效连接体积：半径 BB_THRESHOLD 的球体体积 4/3·π·BB_THRESHOLD³
def v_bb_analytic():
    """返回 B-B 等效连接体积 V_BB = 4/3·π·(2*B_RADIUS+CONNECT_DELTA)³（nm³）。"""
    return 4.0 / 3.0 * np.pi * BB_THRESHOLD ** 3


# 解析计算 A-B 等效连接体积的胶囊近似：半径 R_AB 的胶囊（圆柱段 + 两端半球）
def v_ab_capsule_approx():
    """返回 A-B capsule 近似理论体积 = π·R_AB²·A_LENGTH + 4/3·π·R_AB³（nm³）。"""
    return np.pi * R_AB ** 2 * A_LENGTH + 4.0 / 3.0 * np.pi * R_AB ** 3


# 构造双介质的对称连接矩阵：对角元为各自 ρ·V 项，交叉元为 ρ_B·V_AB（供求最大特征值）
def connection_matrix(n_a, n_b, v_aa, v_ab, v_bb):
    """返回对称 2×2 连接矩阵 M（numpy 数组）。

    M[0,0]=ρ_A·V_AA、M[0,1]=M[1,0]=ρ_B·V_AB、M[1,1]=ρ_B·V_BB，
    其中 ρ_A=N_A/V₀、ρ_B=N_B/V₀，V₀=BOX_VOLUME=1e12 nm³。
    """
    rho_a = n_a / BOX_VOLUME
    rho_b = n_b / BOX_VOLUME
    return np.array([[rho_a * v_aa, rho_b * v_ab],
                     [rho_b * v_ab, rho_b * v_bb]])


# 返回对称矩阵的最大特征值（升序特征值数组末位，用于平均连接度/渗流判定）
def lambda_max(m):
    """返回对称矩阵 m 的最大特征值（np.linalg.eigvalsh(m).max()）。"""
    return float(np.linalg.eigvalsh(m).max())


# 计算双介质成本函数：各自体积 × 数量 × 单价（nm³→μm³ 除以 1e9），单位元
def total_cost(n_a, n_b):
    """返回总成本 C（元）= COST_A_PER_UM3·N_A·V_A/1e9 + COST_B_PER_UM3·N_B·V_B/1e9。

    V_A=π·A_RADIUS²·A_LENGTH（单根圆柱体积）、V_B=4/3·π·B_RADIUS³（单球体积），
    1e9 为 nm³→μm³ 换算系数。
    """
    return (COST_A_PER_UM3 * n_a * CYLINDER_VOLUME / 1e9
            + COST_B_PER_UM3 * n_b * BALL_VOLUME / 1e9)
