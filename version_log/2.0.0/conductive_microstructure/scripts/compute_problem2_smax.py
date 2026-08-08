# 计算问题 2 中各体积分数下的最大连通分量占比 S_max/N，并输出 CSV 与图
import argparse
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np

# 允许从 scripts/ 直接运行（把项目根目录加入 sys.path）
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.config import A_LENGTH, A_RADIUS, BoundaryMode
from src.simulation.single_trial import run_single_trial
from src.visualization.probability_plot import plot_smax

# 盒体积 V0（nm^3）
BOX_VOLUME = 1e12
# 单根介质 A 圆柱体积
CYLINDER_VOLUME = np.pi * A_RADIUS ** 2 * A_LENGTH


# 解析命令行参数：体积分数、trial 数、随机种子、边界模式、并行 worker 数与输出路径
def parse_args():
    parser = argparse.ArgumentParser(
        description="问题2：计算最大连通分量占比 S_max/N 随体积分数变化")
    parser.add_argument("--phi", default="0.005,0.006,0.007,0.010",
                        help="逗号分隔的十进制体积分数（默认 0.005,0.006,0.007,0.010）")
    parser.add_argument("--trials", type=int, default=2000,
                        help="每组 Monte Carlo trial 数（默认 2000）")
    parser.add_argument("--seed", type=int, default=42, help="随机种子（默认 42）")
    parser.add_argument("--mode", default="all",
                        choices=["all", "periodic_connected", "wrapped_geometry_only"],
                        help="边界模式（默认 all=两种都跑）")
    parser.add_argument("--workers", type=int, default=8,
                        help="并行 worker 数（默认 8）")
    parser.add_argument("--out-csv", default=None,
                        help="输出 CSV（默认 <项目根>/results/problem2/problem2_smax.csv）")
    parser.add_argument("--out-dir", default=None,
                        help="图片输出目录（默认 <项目根>/results/figures）")
    return parser.parse_args()


# 由体积分数 φ 计算介质 A 数量 N_A = round(φ·V0/V_A)
def n_a_from_phi(phi):
    """返回 N_A = round(phi * V0 / V_A)（四舍五入取整）。"""
    return int(round(phi * BOX_VOLUME / CYLINDER_VOLUME))


# 对一批子种子运行单次仿真并汇总最大连通分量占比之和
def run_trial_batch_smax(seeds, n_a, mode):
    """对 seeds 中每个子种子各做一次 run_single_trial，返回 max_component_ratio 累计和。"""
    total = 0.0
    for s in seeds:
        total += run_single_trial(np.random.default_rng(int(s)), n_a, mode)["max_component_ratio"]
    return total


# 对给定 trial 种子序列与 (n_a, mode) 运行全部 trial，支持并行，返回平均最大连通分量占比
def run_trials_smax(trial_seeds, n_a, mode, workers):
    """串行或按批次并行运行全部 trial，返回平均 max_component_ratio。"""
    if workers <= 1 or len(trial_seeds) <= 1:
        return float(run_trial_batch_smax(trial_seeds, n_a, mode) / len(trial_seeds))
    batches = np.array_split(np.asarray(trial_seeds), workers)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        totals = list(ex.map(partial(run_trial_batch_smax, n_a=n_a, mode=mode), batches))
    return float(sum(totals) / len(trial_seeds))


# 写入问题 2 的 S_max/N 结果 CSV
def write_csv(rows, csv_path):
    """按列 phi,N_A,boundary_mode,trials,mean_max_component_ratio 写 CSV。"""
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["phi", "N_A", "boundary_mode", "trials", "mean_max_component_ratio"])
        for row in rows:
            writer.writerow([
                f"{row['phi']:.6f}",
                f"{row['N_A']:d}",
                row["boundary_mode"],
                f"{row['trials']:d}",
                f"{row['mean_max_component_ratio']:.6f}",
            ])


# 主流程：逐组计算平均最大连通分量占比，写 CSV，并绘制 S_max/N 图
def main():
    args = parse_args()
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    out_csv = (os.path.abspath(args.out_csv) if args.out_csv
               else os.path.join(project_root, "results", "problem2", "problem2_smax.csv"))
    out_dir = (os.path.abspath(args.out_dir) if args.out_dir
               else os.path.join(project_root, "results", "figures"))
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    phis = [float(x) for x in args.phi.split(",")]
    modes = ([BoundaryMode.PERIODIC_CONNECTED, BoundaryMode.WRAPPED_GEOMETRY_ONLY]
             if args.mode == "all" else [BoundaryMode(args.mode)])
    n_groups = len(phis) * len(modes)
    main_rng = np.random.default_rng(args.seed)
    all_trial_seeds = main_rng.integers(0, 2 ** 31, size=(n_groups, args.trials))

    rows = []
    total_t0 = time.perf_counter()
    group_idx = 0
    for phi in phis:
        n_a = n_a_from_phi(phi)
        for mode in modes:
            t0 = time.perf_counter()
            trial_seeds = all_trial_seeds[group_idx]
            group_idx += 1
            mean_ratio = run_trials_smax(trial_seeds, n_a, mode, args.workers)
            row = {
                "phi": phi,
                "N_A": n_a,
                "boundary_mode": mode.value,
                "trials": args.trials,
                "mean_max_component_ratio": mean_ratio,
            }
            rows.append(row)
            print(f"[{time.perf_counter() - t0:7.1f}s] phi={phi * 100:.2f}% N_A={n_a} "
                  f"mode={mode.value}: mean_max_component_ratio={mean_ratio:.4f}",
                  flush=True)

    write_csv(rows, out_csv)
    png = plot_smax(out_csv, out_dir)
    total_time = time.perf_counter() - total_t0
    print(f"\n总耗时: {total_time:.1f}s ({total_time / 60.0:.1f} min)")
    print(f"CSV 已保存: {out_csv}")
    print(f"图已保存: {png}")


if __name__ == "__main__":
    main()
