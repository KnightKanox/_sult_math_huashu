# 求解问题3：二分搜索 P_conn≥90%（95% Wilson 单侧下界）的最低 A 填充率 φ_min
# 判定标准：evaluate(φ)["p_lower"] >= 0.90，其中 p_lower = wilson_one_sided_lower(connected, trials)，
# 即 95% 单侧置信下界 ≥ 0.90 才认为该填充率满足"导通概率 ≥ 90%"（理论对照 k̄=ρ·V_AA≈1，φ_theory≈0.005546）。
# 对每个边界模式在 [φ_low, φ_high] 上二分（含真实 Monte Carlo 评估），收敛后取 φ_high 端
# （保守上界）作为 φ_min，并用更多 trial 数对 φ_low_end/φ_high_end 两个端点做最终确认。
# 注意：本脚本主流程不在开发阶段运行（真实 MC 耗时较长），开发期只执行
#   python -m py_compile scripts\solve_problem3.py
#   python -m pytest tests\test_problem3.py -q
import argparse
import csv
import json
import os
import sys

import numpy as np

# 允许从 scripts/ 直接运行（把项目根目录加入 sys.path）
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.cloud.effective_volume import mean_degree, num_density
from src.config import BoundaryMode
from src.simulation.confidence import wilson_one_sided_lower
from solve_problem2 import (BOX_VOLUME, CYLINDER_VOLUME, load_v_aa,
                            n_a_from_phi, run_trials)

# 二分最大迭代次数（binary_search_min_phi 默认参数，不暴露为 CLI 选项）
MAX_ITER = 40
# 判定阈值：Wilson 单侧 95% 置信下界需达到的目标导通概率
TARGET_P = 0.90


# 解析命令行参数：二分区间与容差、MC trial 数、种子、边界模式、并行数与输出目录
def parse_args():
    parser = argparse.ArgumentParser(
        description="问题3：二分搜索 P_conn≥90%（95% Wilson 单侧下界）的最低 A 填充率")
    parser.add_argument("--phi-low", type=float, default=0.002,
                        help="二分下界 φ（十进制小数，默认 0.002=0.2%）")
    parser.add_argument("--phi-high", type=float, default=0.007,
                        help="二分上界 φ（十进制小数，默认 0.007=0.7%）")
    parser.add_argument("--tol", type=float, default=0.0001,
                        help="二分收敛容差（默认 0.0001=0.01% 绝对精度）")
    parser.add_argument("--trials", type=int, default=1000,
                        help="二分中点每次评估的 Monte Carlo trial 数（默认 1000）")
    parser.add_argument("--confirm-trials", type=int, default=4000,
                        help="最终确认点的 Monte Carlo trial 数（默认 4000）")
    parser.add_argument("--seed", type=int, default=42, help="随机种子（默认 42）")
    parser.add_argument("--mode", default="all",
                        choices=["all", "periodic_connected", "wrapped_geometry_only"],
                        help="边界模式（默认 all=两种都跑）")
    parser.add_argument("--workers", type=int, default=8,
                        help="并行 worker 数（>1 时按 trial 批次并行，默认 8）")
    parser.add_argument("--out-dir", default=None,
                        help="输出目录（默认 <项目根>/results/problem3）")
    parser.add_argument("--plot", dest="plot", action="store_true",
                        help="运行结束后绘制二分搜索过程图")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                        help="不绘制可视化图")
    parser.set_defaults(plot=True)
    return parser.parse_args()


# 二分搜索满足 P_conn≥90% 的最低填充率：evaluate 为可注入评估回调（返回 connected/trials/p_hat/p_lower）
def binary_search_min_phi(phi_low, phi_high, evaluate, mode, tol, max_iter=MAX_ITER):
    """在 [phi_low, phi_high] 上二分搜索最低 φ，使 evaluate(φ)["p_lower"] >= 0.90。

    p_lower>=0.90 时把上界收到中点（phi_high=phi_mid），否则抬高下界（phi_low=phi_mid）；
    当 phi_high-phi_low <= tol 时终止并令 phi_min = phi_high（保守上界）；
    若达到 max_iter 仍未达 tol 则抛 RuntimeError。
    返回 {phi_min, phi_low_end, phi_high_end, iterations, history}，
    history 每步含 iter/phi/n_a/connected/trials/p_hat/p_lower/action/phi_low/phi_high。
    """
    if not (0.0 < phi_low < phi_high):
        raise ValueError(f"需要 0 < phi_low < phi_high，实际 ({phi_low}, {phi_high})")
    if tol <= 0.0:
        raise ValueError("tol 必须为正数")
    if max_iter < 1:
        raise ValueError("max_iter 必须 >= 1")
    history = []
    lo, hi = float(phi_low), float(phi_high)
    for it in range(1, max_iter + 1):
        mid = 0.5 * (lo + hi)
        ev = evaluate(mid, mode)
        if ev["p_lower"] >= TARGET_P:
            hi, action = mid, "high"
        else:
            lo, action = mid, "low"
        history.append({
            "iter": it,
            "phi": mid,
            "n_a": n_a_from_phi(mid),
            "connected": ev["connected"],
            "trials": ev["trials"],
            "p_hat": ev["p_hat"],
            "p_lower": ev["p_lower"],
            "action": action,
            "phi_low": lo,
            "phi_high": hi,
        })
        if hi - lo <= tol:
            return {
                "phi_min": hi,
                "phi_low_end": lo,
                "phi_high_end": hi,
                "iterations": it,
                "history": history,
            }
    raise RuntimeError(
        f"二分搜索 {max_iter} 次迭代内未收敛到 tol={tol}，当前区间 "
        f"[{lo:.6f}, {hi:.6f}] 宽度 {hi - lo:.2e}（请扩大区间或增大 max_iter）")


# 构造真实 MC 评估回调：按调用顺序从预派生种子池取用一组子种子跑 run_trials 并算 Wilson 下界
def make_mc_evaluator(trial_seed_pool, trials, workers):
    """返回 evaluate(phi, mode) -> dict（connected/trials/p_hat/p_lower）。

    闭包内部维护行指针，每次调用消耗 seed_pool 的下一行子种子；子种子由主 rng 在
    main 中一次性派生，评估顺序固定为"逐模式、按二分迭代逐点推进"，因此 seed 相同
    时串行/并行/重复运行结果完全可复现。
    """
    state = {"idx": 0}

    def evaluate(phi, mode):
        # 每行子种子池宽度为 max(trials, confirm_trials)，这里只取前 trials 个用于二分评估
        seeds = trial_seed_pool[state["idx"]][:trials]
        state["idx"] += 1
        connected = run_trials(seeds, n_a_from_phi(phi), mode, workers)
        return {
            "connected": connected,
            "trials": trials,
            "p_hat": float(connected) / trials,
            "p_lower": wilson_one_sided_lower(connected, trials),
        }

    return evaluate


# 对单个填充率点用给定子种子与 trial 数做 MC 复核，返回确认记录 dict
def confirm_point(phi, seeds, mode, trials, workers):
    """对 phi 用 seeds 跑 trials 次 MC，返回 {phi, boundary_mode, trials, connected_count, p_hat, p_lower_95}。"""
    connected = run_trials(seeds, n_a_from_phi(phi), mode, workers)
    return {
        "phi": phi,
        "boundary_mode": mode.value,
        "trials": trials,
        "connected_count": connected,
        "p_hat": float(connected) / trials,
        "p_lower_95": wilson_one_sided_lower(connected, trials),
    }


# 绘制二分搜索过程图：φ(%) 为 x 轴，画 0.90 水平虚线、理论 φ 竖线与各模式迭代点/φ_min 竖线
def plot_bisection(result, phi_theory, out_dir):
    """按问题 3 结果 dict 画二分搜索图，保存 problem3_bisection.png，返回路径（失败抛异常）。"""
    import matplotlib
    matplotlib.use("Agg")  # 无显示环境使用 Agg 后端（须在导入 pyplot 前设置）
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    ax.axhline(0.90, color="gray", linestyle="--", linewidth=1.0, label="判定线 90%")
    x_theory = phi_theory * 100.0
    ax.axvline(x_theory, color="purple", linestyle="--", linewidth=1.0,
               label=f"理论 $\\bar{{k}}$=1（φ={x_theory:.3f}%）")
    colors = {"periodic_connected": "#1f77b4", "wrapped_geometry_only": "#d62728"}
    for mode_key, mdata in result["modes"].items():
        color = colors.get(mode_key, "#333333")
        xs = [h["phi"] * 100.0 for h in mdata["history"]]
        ax.plot(xs, [h["p_hat"] for h in mdata["history"]],
                marker="o", linestyle="-", color=color, label=f"{mode_key} p_hat")
        ax.plot(xs, [h["p_lower"] for h in mdata["history"]],
                marker="s", linestyle="--", color=color, alpha=0.7,
                label=f"{mode_key} p_lower(95%)")
        ax.axvline(mdata["phi_min"] * 100.0, color=color, linestyle="-", linewidth=1.2,
                   label=f"{mode_key} φ_min={mdata['phi_min'] * 100.0:.4f}%")
    ax.set_xlabel("体积分数 φ (%)")
    ax.set_ylabel("导通概率 P_conn")
    ax.set_title("问题 3 二分搜索：P_conn≥90% 的最低填充率 φ_min")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    png = os.path.join(out_dir, "problem3_bisection.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 写入问题 3 汇总 CSV：每 mode 一行（p_hat/p_lower_95 取 φ_min 即 phi_high_end 处的二分评估值）
def write_result_csv(rows, csv_path):
    """按列 phi_min,N_A,boundary_mode,trials,p_hat,p_lower_95,kbar,phi_theory,phi_low_end,phi_high_end 写 CSV。"""
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["phi_min", "N_A", "boundary_mode", "trials", "p_hat",
                         "p_lower_95", "kbar", "phi_theory", "phi_low_end", "phi_high_end"])
        for row in rows:
            writer.writerow([
                f"{row['phi_min']:.6f}", f"{row['N_A']:d}", row["boundary_mode"],
                f"{row['trials']:d}",
                "" if row["p_hat"] is None else f"{row['p_hat']:.6f}",
                "" if row["p_lower_95"] is None else f"{row['p_lower_95']:.6f}",
                f"{row['kbar']:.6f}", f"{row['phi_theory']:.6f}",
                f"{row['phi_low_end']:.6f}", f"{row['phi_high_end']:.6f}",
            ])


# 写入问题 3 确认点 CSV：φ_low_end/φ_high_end 各 mode 每行
def write_confirm_csv(confirm_rows, csv_path):
    """按列 phi,boundary_mode,trials,connected_count,p_hat,p_lower_95 写 CSV。"""
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["phi", "boundary_mode", "trials", "connected_count",
                         "p_hat", "p_lower_95"])
        for row in confirm_rows:
            writer.writerow([
                f"{row['phi']:.6f}", row["boundary_mode"], f"{row['trials']:d}",
                f"{row['connected_count']:d}", f"{row['p_hat']:.6f}",
                f"{row['p_lower_95']:.6f}",
            ])


# 主流程：解析参数 → 加载 V_AA → 各模式二分（真实 MC）→ 确认 → 写 CSV/JSON → 可选绘图 → 打印汇总
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
               else os.path.join(project_root, "results", "problem3"))
    os.makedirs(out_dir, exist_ok=True)

    # 双模式先 periodic 后 wrapped，与 solve_problem2 的枚举顺序一致
    modes = ([BoundaryMode.PERIODIC_CONNECTED, BoundaryMode.WRAPPED_GEOMETRY_ONLY]
             if args.mode == "all" else [BoundaryMode(args.mode)])
    phi_theory = float(CYLINDER_VOLUME / v_aa)

    # 可复现性方案：主 rng 一次性派生全部评估点（二分 + 确认）的子种子。
    # 评估顺序固定为"逐模式、按二分迭代逐点推进"（每模式一个独立种子池，双模式先 periodic
    # 后 wrapped）；每个评估点派生 max(trials, confirm_trials) 个子种子，二分用前 trials 个、
    # 确认用前 confirm_trials 个。池行数取 2*MAX_ITER + 2：二分最多 MAX_ITER 个点 + 确认 2 个点。
    pool_width = max(args.trials, args.confirm_trials)
    n_points = 2 * MAX_ITER + 2
    main_rng = np.random.default_rng(args.seed)
    seed_pools = main_rng.integers(0, 2 ** 31, size=(len(modes), n_points, pool_width))

    result_modes = {}
    csv_rows, confirm_rows = [], []
    for m_idx, mode in enumerate(modes):
        evaluator = make_mc_evaluator(seed_pools[m_idx], args.trials, args.workers)
        bs = binary_search_min_phi(args.phi_low, args.phi_high, evaluator, mode,
                                   args.tol, max_iter=MAX_ITER)
        # φ_min 即 phi_high_end，对应最后一次 action=="high" 的二分评估（p_hat/p_lower 取该点值）
        last_high = next((h for h in reversed(bs["history"]) if h["action"] == "high"), None)
        p_hat = last_high["p_hat"] if last_high else None
        p_lower = last_high["p_lower"] if last_high else None
        n_a = n_a_from_phi(bs["phi_min"])
        kbar = mean_degree(num_density(n_a, BOX_VOLUME), v_aa)

        # 最终确认：对 phi_low_end 与 phi_high_end 各用 confirm_trials 次复核
        confirms = []
        for j, phi in enumerate([bs["phi_low_end"], bs["phi_high_end"]]):
            rec = confirm_point(phi, seed_pools[m_idx, MAX_ITER + j, :args.confirm_trials],
                                mode, args.confirm_trials, args.workers)
            confirms.append(rec)
            confirm_rows.append(rec)

        csv_rows.append({
            "phi_min": bs["phi_min"], "N_A": n_a, "boundary_mode": mode.value,
            "trials": args.trials, "p_hat": p_hat, "p_lower_95": p_lower,
            "kbar": kbar, "phi_theory": phi_theory,
            "phi_low_end": bs["phi_low_end"], "phi_high_end": bs["phi_high_end"],
        })
        result_modes[mode.value] = {
            "phi_min": bs["phi_min"], "N_A": n_a,
            "p_hat": p_hat, "p_lower_95": p_lower,
            "kbar": kbar, "phi_low_end": bs["phi_low_end"],
            "phi_high_end": bs["phi_high_end"],
            "iterations": bs["iterations"], "history": bs["history"],
            "confirm": confirms,
        }
        print(f"mode={mode.value}: φ_min={bs['phi_min'] * 100:.4f}% N_A={n_a} "
              f"p_hat={p_hat:.4f} p_lower(95%)={p_lower:.4f} kbar={kbar:.4f} "
              f"φ_theory={phi_theory * 100:.4f}% iterations={bs['iterations']}", flush=True)

    csv_path = os.path.join(out_dir, "problem3_result.csv")
    confirm_csv_path = os.path.join(out_dir, "problem3_confirm.csv")
    json_path = os.path.join(out_dir, "problem3_result.json")
    write_result_csv(csv_rows, csv_path)
    write_confirm_csv(confirm_rows, confirm_csv_path)

    result = {
        "seed": args.seed,
        "config": vars(args),
        "V_AA": v_aa,
        "phi_theory": phi_theory,
        "modes": result_modes,
    }
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    if args.plot:
        try:
            png = plot_bisection(result, phi_theory, out_dir)
            print(f"图已保存: {png}")
        except Exception as exc:  # 绘图失败不影响结果文件
            print(f"绘图失败（结果文件已保存）: {exc}")

    print(f"\nCSV 已保存: {csv_path}")
    print(f"确认 CSV 已保存: {confirm_csv_path}")
    print(f"JSON 已保存: {json_path}")


if __name__ == "__main__":
    main()
