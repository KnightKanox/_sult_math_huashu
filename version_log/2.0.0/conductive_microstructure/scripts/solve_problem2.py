# 求解问题2：随机填充介质 A 的导通概率（Monte Carlo）
# 对 4 组体积分数 φ × 2 种边界模式，各做 trials 次独立随机仿真，统计左右电极
# 贯通概率 p_hat 与 95% Wilson 置信区间，并用概率云理论平均连接度 k̄ = ρ_A·V_AA
# 对照。每个 trial 使用独立子种子（由主 rng 一次性派生），保证串行/并行/重复
# 运行结果完全一致；--workers>1 时用 ProcessPoolExecutor 按 trial 批次并行。
import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np

# 允许从 scripts/ 直接运行（把项目根目录加入 sys.path）
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.config import A_LENGTH, A_RADIUS, BoundaryMode
from src.cloud.effective_volume import effective_volume, mean_degree, num_density
from src.simulation.confidence import wilson_ci
from src.simulation.single_trial import run_single_trial

# 盒体积 V₀（nm³）
BOX_VOLUME = 1e12
# 单根介质 A 圆柱体积 V_A = π·R²·L（nm³）
CYLINDER_VOLUME = np.pi * A_RADIUS ** 2 * A_LENGTH


# 解析命令行参数：体积分数、trial 数、随机种子、边界模式、并行 worker 数与输出目录
def parse_args():
    parser = argparse.ArgumentParser(
        description="问题2：随机填充介质A导通概率 Monte Carlo（4组体积分数×2种边界模式）")
    parser.add_argument("--phi", default="0.005,0.006,0.007,0.010",
                        help="逗号分隔的十进制体积分数（默认 0.005,0.006,0.007,0.010）")
    parser.add_argument("--trials", type=int, default=2000,
                        help="每组 Monte Carlo trial 数（默认 2000）")
    parser.add_argument("--seed", type=int, default=42, help="随机种子（默认 42）")
    parser.add_argument("--mode", default="all",
                        choices=["all", "periodic_connected", "wrapped_geometry_only"],
                        help="边界模式（默认 all=两种都跑）")
    parser.add_argument("--workers", type=int, default=1,
                        help="并行 worker 数（>1 时用 ProcessPoolExecutor 按 trial 批次并行，默认 1）")
    parser.add_argument("--out-dir", default=None,
                        help="输出目录（默认 <项目根>/results/problem2）")
    parser.add_argument("--plot", dest="plot", action="store_true",
                        help="运行结束后自动绘制可视化图")
    parser.add_argument("--no-plot", dest="plot", action="store_false", help="不绘制可视化图")
    parser.set_defaults(plot=True)
    return parser.parse_args()


# 由体积分数 φ 计算介质 A 数量 N_A = round(φ·V₀/V_A)
def n_a_from_phi(phi):
    """返回 N_A = round(phi * V0 / V_A)（四舍五入取整）。"""
    return int(round(phi * BOX_VOLUME / CYLINDER_VOLUME))


# 读取概率云 CSV 并积分得到等效连接体积 V_AA（nm³）
def load_v_aa(cloud_csv):
    """读取 cloud_AA.csv（列 r_nm,q_hat,...）并用 4π∫r²q dr 积分 V_AA；文件缺失返回 None。"""
    if not os.path.isfile(cloud_csv):
        return None
    r, q = [], []
    with open(cloud_csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            r.append(float(row["r_nm"]))
            q.append(float(row["q_hat"]))
    return effective_volume(np.array(r), np.array(q))


# 运行一批独立 trial（每个 trial 用独立子种子派生），返回贯通次数
def run_trial_batch(seeds, n_a, mode):
    """对 seeds 中每个子种子各做一次 run_single_trial，返回 connected 的累计次数。"""
    connected = 0
    for s in seeds:
        if run_single_trial(np.random.default_rng(int(s)), n_a, mode)["connected"]:
            connected += 1
    return connected


# 对给定 trial 种子序列与 (n_a, mode) 运行全部 trial，支持并行，返回贯通次数
def run_trials(trial_seeds, n_a, mode, workers):
    """串行或按批次并行运行全部 trial，返回贯通次数（结果与 workers 取值无关）。"""
    if workers <= 1 or len(trial_seeds) <= 1:
        return run_trial_batch(trial_seeds, n_a, mode)
    batches = np.array_split(np.asarray(trial_seeds), workers)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        counts = list(ex.map(partial(run_trial_batch, n_a=n_a, mode=mode), batches))
    return int(sum(counts))


# 写入问题 2 结果 CSV
def write_csv(rows, csv_path):
    """按列 phi,N_A,boundary_mode,trials,connected_count,p_hat,ci_low,ci_high,rho,mean_degree_theory 写 CSV。"""
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["phi", "N_A", "boundary_mode", "trials", "connected_count",
                         "p_hat", "ci_low", "ci_high", "rho", "mean_degree_theory"])
        for row in rows:
            writer.writerow([
                f"{row['phi']:.6f}", f"{row['N_A']:d}", row["boundary_mode"],
                f"{row['trials']:d}", f"{row['connected_count']:d}",
                f"{row['p_hat']:.6f}", f"{row['ci_low']:.6f}", f"{row['ci_high']:.6f}",
                f"{row['rho']:.6e}", f"{row['mean_degree_theory']:.6f}",
            ])


# 主流程：读取云数据、逐组运行 Monte Carlo、写 CSV/JSON、打印进度并可选绘图
def main():
    args = parse_args()
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    cloud_csv = os.path.join(project_root, "results", "cloud", "cloud_AA.csv")
    v_aa = load_v_aa(cloud_csv)
    if v_aa is None:
        print(f"缺少概率云数据: {cloud_csv}")
        print("请先运行: python scripts\\build_aa_cloud.py")
        sys.exit(1)
    out_dir = (os.path.abspath(args.out_dir) if args.out_dir
               else os.path.join(project_root, "results", "problem2"))
    os.makedirs(out_dir, exist_ok=True)

    phis = [float(x) for x in args.phi.split(",")]
    modes = ([BoundaryMode.PERIODIC_CONNECTED, BoundaryMode.WRAPPED_GEOMETRY_ONLY]
             if args.mode == "all" else [BoundaryMode(args.mode)])

    # 主 rng 一次性派生全部 (phi, mode) 组的 trial 子种子，保证任意顺序/并行均可复现
    n_groups = len(phis) * len(modes)
    main_rng = np.random.default_rng(args.seed)
    all_trial_seeds = main_rng.integers(0, 2 ** 31, size=(n_groups, args.trials))

    rows, groups_json = [], []
    total_t0 = time.perf_counter()
    group_idx = 0
    for phi in phis:
        n_a = n_a_from_phi(phi)
        rho = num_density(n_a, BOX_VOLUME)
        k_bar = mean_degree(rho, v_aa)
        for mode in modes:
            t0 = time.perf_counter()
            trial_seeds = all_trial_seeds[group_idx]
            group_idx += 1
            connected = run_trials(trial_seeds, n_a, mode, args.workers)
            p_hat = connected / args.trials
            ci_low, ci_high = wilson_ci(connected, args.trials)
            row = {
                "phi": phi, "N_A": n_a, "boundary_mode": mode.value,
                "trials": args.trials, "connected_count": connected,
                "p_hat": p_hat, "ci_low": ci_low, "ci_high": ci_high,
                "rho": rho, "mean_degree_theory": k_bar,
            }
            rows.append(row)
            groups_json.append(row)
            print(f"[{time.perf_counter() - t0:7.1f}s] phi={phi * 100:.2f}% N_A={n_a} "
                  f"mode={mode.value}: p_hat={p_hat:.4f} "
                  f"(95% CI {ci_low:.4f}~{ci_high:.4f}), mean_degree_theory={k_bar:.4f}", flush=True)

    csv_path = os.path.join(out_dir, "problem2_result.csv")
    write_csv(rows, csv_path)

    # 汇总统计：各模式下 p_hat 均值与 k̄ 均值
    summary = {}
    for mode in modes:
        vals = [r["p_hat"] for r in rows if r["boundary_mode"] == mode.value]
        summary["mean_p_hat_" + mode.value] = float(np.mean(vals)) if vals else None
    summary["mean_mean_degree_theory"] = float(np.mean([r["mean_degree_theory"] for r in rows]))
    result = {
        "seed": args.seed, "trials": args.trials, "workers": args.workers,
        "mode": args.mode, "V_AA": v_aa, "groups": groups_json, "summary": summary,
    }
    json_path = os.path.join(out_dir, "problem2_result.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    total_time = time.perf_counter() - total_t0
    print(f"\n总耗时: {total_time:.1f}s ({total_time / 60.0:.1f} min)")
    print(f"CSV 已保存: {csv_path}")
    print(f"JSON 已保存: {json_path}")

    if args.plot:
        try:
            from src.visualization.probability_plot import plot_kbar, plot_p_conn
            fig_dir = os.path.join(project_root, "results", "figures")
            png1 = plot_p_conn(csv_path, fig_dir)
            png2 = plot_kbar(csv_path, fig_dir)
            print(f"图已保存: {png1}\n{png2}")
        except Exception as exc:  # 绘图失败不影响结果文件
            print(f"绘图失败（结果文件已保存）: {exc}")


if __name__ == "__main__":
    main()
