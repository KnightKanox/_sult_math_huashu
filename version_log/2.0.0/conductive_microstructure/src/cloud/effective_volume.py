# 概率云层：等效连接体积与平均连接度（由 q(r) 曲线 4π∫r²q dr 积分得到）
import numpy as np


# 由连接概率曲线积分得到等效连接体积 V = 4π∫₀^∞ r²·q(r) dr
def effective_volume(r_grid, q_grid):
    """计算等效连接体积 V = 4π * trapezoid(r²·q, r)（nm³）。

    r_grid、q_grid 为同长一维 numpy 数组，r 单位 nm，q 为无量纲概率。
    """
    return 4.0 * np.pi * np.trapezoid(r_grid ** 2 * q_grid, r_grid)


# 平均连接度：体积数密度 × 等效连接体积
def mean_degree(rho, v_eff):
    """返回平均连接度 k̄ = rho * v_eff（rho 单位 nm⁻³，v_eff 单位 nm³）。"""
    return rho * v_eff


# 体积数密度：数量除以盒体积（nm⁻³）
def num_density(n_a, box_volume=1e12):
    """返回数量为 n_a 的介质数密度 ρ = n_a / box_volume（nm⁻³）。"""
    return n_a / box_volume
