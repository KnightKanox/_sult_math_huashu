# 构建介质 A 概率云：在 r 网格上采样 q_AA(r)，输出 CSV 并积分得到 V_AA（问题 3）
import argparse
import csv
import os
import sys

import numpy as np

# 允许从 scripts/ 直接运行（把项目根目录加入 sys.path）
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.cloud.aa_cloud import estimate_q_aa
from src.cloud.effective_volume import effective_volume, mean_degree, num_density


# 解析命令行参数：r 网格范围/步长、每点采样数、随机种子与输出路径
def parse_args():
    parser = argparse.ArgumentParser(description="构建介质 A 概率云（q_AA(r) 采样并积分 V_AA）")
    parser.add_argument("--r-min", type=float, default=0.0, help="r 网格起点（nm）")
    parser.add_argument("--r-max", type=float, default=5100.0, help="r 网格终点（nm）")
    parser.add_argument("--r-step", type=float, default=25.0, help="r 网格步长（nm）")
    parser.add_argument("--samples-per-r", type=int, default=2000, help="每个 r 点的方向采样数")
    parser.add_argument("--seed", type=int, default=42, help="随机种子")
    parser.add_argument("--n-a", type=int, default=None,
                        help="介质 A 数量 N_A（提供时额外打印平均连接度 ρ_A·V_AA）")
    parser.add_argument("--out", default=os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "results", "cloud", "cloud_AA.csv")),
        help="输出 CSV 路径")
    return parser.parse_args()


# 主流程：生成 r 网格，逐点采样 q_AA，写 CSV 并打印 V_AA（与可选的 ρ_A·V_AA）
def main():
    args = parse_args()
    r_grid = np.arange(args.r_min, args.r_max + args.r_step * 0.5, args.r_step)
    rng = np.random.default_rng(args.seed)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    rows = []
    for r in r_grid:
        q_hat, success_count, sample_count, ci_low, ci_high = \
            estimate_q_aa(float(r), rng, samples_per_r=args.samples_per_r)
        rows.append((float(r), q_hat, success_count, sample_count, ci_low, ci_high))

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["r_nm", "q_hat", "success_count", "sample_count",
                         "ci_low", "ci_high"])
        for row in rows:
            writer.writerow([
                f"{row[0]:.6f}", f"{row[1]:.10f}", f"{row[2]:d}", f"{row[3]:d}",
                f"{row[4]:.10f}", f"{row[5]:.10f}",
            ])

    r_arr = np.array([row[0] for row in rows], float)
    q_arr = np.array([row[1] for row in rows], float)
    v_aa = effective_volume(r_arr, q_arr)
    print(f"r 网格: {r_arr[0]:.1f} ~ {r_arr[-1]:.1f} nm, 共 {len(r_arr)} 点")
    print(f"V_AA = 4π∫r²·q_AA(r)dr ≈ {v_aa:.3e} nm³")
    if args.n_a is not None:
        rho_a = num_density(args.n_a)
        k_bar = mean_degree(rho_a, v_aa)
        print(f"N_A={args.n_a}: ρ_A={rho_a:.6e} nm⁻³, 平均连接度 k̄ = ρ_A·V_AA ≈ {k_bar:.4f}")
    print(f"结果已保存: {args.out}")


if __name__ == "__main__":
    main()
