# 概率云层（问题 3 Task 3）：二体连接概率的蒙特卡洛采样估计与等效连接体积
# 固定两介质中心距 r，方向独立球面均匀采样，统计"连接"比例 q_AA(r)，
# 再用 4π∫r²q(r)dr 得到等效连接体积 V_AA（平均连接度 k̄ = ρ_A·V_AA）
from src.config import A_LENGTH, AA_SEG_THRESHOLD
from src.geometry.segment_distance import segment_segment_distance
from src.generation.random_orientation import random_unit_vector
from src.simulation.confidence import wilson_ci


# 两圆柱最大可能连接的中心距：同轴首尾相接（L_A + δ），r 超过该值 q=0
def aa_max_contact_center_distance():
    """返回介质 A 两圆柱可发生连接的最大中心距（= A_LENGTH + AA_SEG_THRESHOLD，nm）。"""
    return A_LENGTH + AA_SEG_THRESHOLD


# 通用二体连接概率估计器：对 samples_per_r 个球面均匀方向采样并统计连接比例
def estimate_q(r, rng, connected_func, samples_per_r):
    """估计给定中心距 r 下的二体连接概率。

    connected_func(r_hat) 输入相对方向单位向量 r_hat（形状 (3,)），返回 bool
    （该次采样是否连接），第二根介质中心位于 r*r_hat，第一根中心在原点。
    返回 (q_hat, success_count, sample_count, ci_low, ci_high)，
    其中 CI 为 Wilson 双侧 95%（z=1.96）。
    """
    success_count = 0
    for _ in range(samples_per_r):
        r_hat = random_unit_vector(rng)
        if connected_func(r_hat):
            success_count += 1
    q_hat = success_count / float(samples_per_r)
    ci_low, ci_high = wilson_ci(success_count, samples_per_r, z=1.96)
    return (q_hat, success_count, samples_per_r, ci_low, ci_high)


# 制造介质 A 的连接判定闭包：每次采样生成两圆柱的独立球面均匀方向并判定轴距
def make_aa_connected_func(rng, r):
    """返回 connected_func(r_hat)，闭包捕获中心距 r 与随机源 rng。

    第一根圆柱以原点为中心、方向 u0；第二根以 r*r_hat 为中心、方向 u1；
    两轴线段长度均为 A_LENGTH，连接判据为轴线段最短距离 ≤ AA_SEG_THRESHOLD。
    """
    a_half = A_LENGTH / 2.0

    def connected_func(r_hat):
        u0 = random_unit_vector(rng)
        u1 = random_unit_vector(rng)
        p1 = -u0 * a_half
        q1 = u0 * a_half
        center2 = r * r_hat
        p2 = center2 - u1 * a_half
        q2 = center2 + u1 * a_half
        return segment_segment_distance(p1, q1, p2, q2) <= AA_SEG_THRESHOLD

    return connected_func


# 介质 A 连接概率 q_AA(r) 的包装接口：组合闭包制造与通用估计器
def estimate_q_aa(r, rng, samples_per_r=2000):
    """估计介质 A 两圆柱在中心距 r（nm）下的连接概率。

    返回与 estimate_q 相同结构：(q_hat, success_count, sample_count, ci_low, ci_high)。
    """
    connected_func = make_aa_connected_func(rng, r)
    return estimate_q(r, rng, connected_func, samples_per_r)
