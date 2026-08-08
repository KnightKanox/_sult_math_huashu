# 求解问题4：A/B 混合体系 P_conn≥90%（95% Wilson 单侧置信下界）的最低成本组合 (N_A*, N_B*)
# 三级流程：Level 1 λ_max 理论筛选（N_A×N_B 网格，λ_max≥1 为理论可行区；--theory-only 时到此为止）→
# Level 2 对每个 N_A 网格点整数二分最小 N_B（判据 wilson_one_sided_lower ≥ 0.90）得到可行边界
# N_B_min(N_A) 并做单调性校验 → Level 2.5 MC 粗网格热图数据（供 Figure 7）→
# Level 3 沿边界比较 total_cost 取全局最小候选，对候选及其 8 邻域用 confirm_trials 复核得到最优解。
# 成本 C = 1.05·N_A·V_A/1e9 + 0.05·N_B·V_B/1e9（元）；可复现机制沿用 Q3：主 rng 一次性派生
# bisect/confirm/mc_grid 三个子种子池，固定评估顺序（逐模式、逐 N_A、逐二分迭代），串/并行结果一致。
# 注意：run_single_trial_mixed 由 T3 任务提供，本模块在其落地前仍可被 import（函数内延迟导入）。
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

from src.cloud.mixed_cloud import (ab_max_contact_center_distance,
                                   connection_matrix, estimate_q_ab,
                                   estimate_v_ab, lambda_max, total_cost,
                                   v_ab_capsule_approx, v_bb_analytic)
from src.config import BoundaryMode
from src.simulation.confidence import wilson_one_sided_lower
from solve_problem2 import load_v_aa

# 判定阈值：Wilson 单侧 95% 置信下界需达到的目标导通概率
TARGET_P = 0.90
# 整数二分子种子池行宽（20 次评估足够覆盖 [0,5000] 的二分与边界检查）
MAX_BISECT_ITER = 20
# Level 1 理论筛选网格步长（N_A 步长 10、N_B 步长 100，固定不暴露为 CLI）
L1_NA_STEP = 10
L1_NB_STEP = 100
# Level 3 确认邻域：N_B 固定窗半宽（N_A 方向取 ±N_A 步长）
CONFIRM_DB = 50
# Level 3 确认点最多数量（候选 + 8 邻域 = 9）
CONFIRM_MAX_POINTS = 9


# 解析命令行参数：搜索范围、MC trial 数、种子、边界模式、并行数与输出选项
def parse_args():
    parser = argparse.ArgumentParser(
        description="问题4：A/B混合体系 P_conn≥90%（95% Wilson 单侧下界）的最低成本组合")
    parser.add_argument("--n-a-max", type=int, default=360,
                        help="N_A 搜索上限（默认 360）")
    parser.add_argument("--n-a-step", type=int, default=10,
                        help="Level 2/3 的 N_A 网格步长（默认 10）")
    parser.add_argument("--n-b-max", type=int, default=5000,
                        help="N_B 搜索上限（默认 5000）")
    parser.add_argument("--trials", type=int, default=500,
                        help="二分评估每次 Monte Carlo trial 数（默认 500）")
    parser.add_argument("--confirm-trials", type=int, default=4000,
                        help="Level 3 确认点 Monte Carlo trial 数（默认 4000）")
    parser.add_argument("--seed", type=int, default=42, help="随机种子（默认 42）")
    parser.add_argument("--mode", default="all",
                        choices=["all", "periodic_connected", "wrapped_geometry_only"],
                        help="边界模式（默认 all=两种都跑）")
    parser.add_argument("--workers", type=int, default=8,
                        help="并行 worker 数（>1 时按 trial 批次并行，默认 8）")
    parser.add_argument("--out-dir", default=None,
                        help="输出目录（默认 <项目根>/results/problem4）")
    parser.add_argument("--plot", dest="plot", action="store_true",
                        help="运行结束后绘制可视化图")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                        help="不绘制可视化图")
    parser.set_defaults(plot=True)
    parser.add_argument("--theory-only", action="store_true",
                        help="只做 V_AB/V_BB 计算 + Level 1 λ_max 理论筛选 + 输出与图，不做任何 MC")
    parser.add_argument("--ab-samples", type=int, default=500,
                        help="q_AB(r) 每个距离点采样数（默认 500）")
    parser.add_argument("--r-dr", type=float, default=25.0,
                        help="V_AB 数值积分的 r 步长 nm（默认 25）")
    parser.add_argument("--mc-grid-na-step", type=int, default=60,
                        help="MC 粗网格 N_A 步长（默认 60）")
    parser.add_argument("--mc-grid-nb-step", type=int, default=1000,
                        help="MC 粗网格 N_B 步长（默认 1000）")
    parser.add_argument("--mc-grid-trials", type=int, default=300,
                        help="MC 粗网格每点 trial 数（默认 300）")
    return parser.parse_args()


# 对 N_A×N_B 网格逐点计算连接矩阵 λ_max，返回网格与 λ_max 矩阵（供 Level 1 理论筛选）
def lambda_max_grid(n_a_pts, n_b_pts, v_aa, v_ab, v_bb):
    """在 n_a_pts×n_b_pts 网格上计算 λ_max(connection_matrix(...))。

    返回 (n_a_grid, n_b_grid, lam_grid)，lam_grid 形状 (len(n_a_pts), len(n_b_pts))。
    """
    n_a_grid = np.asarray(n_a_pts)
    n_b_grid = np.asarray(n_b_pts)
    lam_grid = np.zeros((len(n_a_grid), len(n_b_grid)))
    for i, na in enumerate(n_a_grid):
        for j, nb in enumerate(n_b_grid):
            lam_grid[i, j] = lambda_max(
                connection_matrix(int(na), int(nb), v_aa, v_ab, v_bb))
    return n_a_grid, n_b_grid, lam_grid


# 对一批子种子各做一次混合单次仿真（run_single_trial_mixed），返回贯通次数
def run_trial_batch_mixed(seeds, n_a, n_b, mode):
    """对 seeds 中每个子种子各做一次 run_single_trial_mixed，返回 connected 的累计次数。

    run_single_trial_mixed 由 T3 任务在 src/simulation/single_trial.py 中新增，
    本函数内延迟导入以保证模块在 T3 落地前可被安全 import。
    """
    from src.simulation.single_trial import run_single_trial_mixed
    connected = 0
    for s in seeds:
        if run_single_trial_mixed(np.random.default_rng(int(s)), n_a, n_b, mode)["connected"]:
            connected += 1
    return connected


# 对给定 trial 种子序列与 (n_a, n_b, mode) 运行全部混合 trial，支持并行，返回贯通次数
def run_trials_mixed(trial_seeds, n_a, n_b, mode, workers, executor=None):
    """串行或按批次并行运行全部混合 trial，返回贯通次数（结果与 workers 取值无关）。

    executor 非 None 时复用外部传入的共享 ProcessPoolExecutor（由 main 创建一次，
    避免每次二分评估重建 8 个进程的启动开销）；否则每次调用自建并关闭进程池。
    """
    if workers <= 1 or len(trial_seeds) <= 1:
        return run_trial_batch_mixed(trial_seeds, n_a, n_b, mode)
    batches = np.array_split(np.asarray(trial_seeds), workers)
    task = partial(run_trial_batch_mixed, n_a=n_a, n_b=n_b, mode=mode)
    if executor is not None:
        return int(sum(executor.map(task, batches)))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        counts = list(ex.map(task, batches))
    return int(sum(counts))


# 构造固定 N_A 的二分评估回调：按调用次数依次取该 (mode, N_A) 点子种子池行跑混合 MC 并算 Wilson 下界
def make_nb_evaluator(seed_pool, n_a, mode, trials, workers, executor=None):
    """返回 evaluate(n_b) -> dict（connected/trials/p_hat/p_lower）。

    seed_pool 为 (max_iter, trials) 形状的该 (mode, N_A) 点子种子池；闭包内维护行指针，
    每次调用消耗下一行子种子，评估顺序固定为二分调用顺序，保证同 seed 串/并行结果一致；
    executor 非 None 时复用共享进程池（见 run_trials_mixed）。
    """
    state = {"it": 0}

    def evaluate(n_b):
        seeds = seed_pool[state["it"]][:trials]
        state["it"] += 1
        connected = run_trials_mixed(seeds, n_a, n_b, mode, workers, executor)
        return {
            "connected": connected,
            "trials": trials,
            "p_hat": float(connected) / trials,
            "p_lower": wilson_one_sided_lower(connected, trials),
        }

    return evaluate


# 整数二分搜索满足 evaluate(n_b)["p_lower"]>=0.90 的最小 N_B，返回 {n_b_min, iterations, history}
def binary_search_min_nb(n_b_low, n_b_high, evaluate, max_iter=MAX_BISECT_ITER):
    """在 [n_b_low, n_b_high] 上整数二分最小 N_B，使 evaluate(n_b)["p_lower"] >= 0.90。

    边界情形：evaluate(n_b_low) 已满足 → n_b_min=n_b_low；evaluate(n_b_high) 不满足 →
    该 N_A 不可行（n_b_min=None）；否则 while hi-lo>1 收窄（lo=不满足/未知、hi=满足），
    终止时 hi 即最小可行 N_B。history 每步含 iter/n_b/connected/trials/p_hat/p_lower/action/lo/hi。
    """
    if not (0 <= n_b_low < n_b_high):
        raise ValueError(f"需要 0 <= n_b_low < n_b_high，实际 ({n_b_low}, {n_b_high})")
    history = []
    # 下界端检查：n_b_low 已满足 → 最小可行即下界端
    ev = evaluate(n_b_low)
    satisfied = ev["p_lower"] >= TARGET_P
    history.append({
        "iter": 1, "n_b": n_b_low, "connected": ev["connected"], "trials": ev["trials"],
        "p_hat": ev["p_hat"], "p_lower": ev["p_lower"],
        "action": "high" if satisfied else "low",
        "lo": n_b_low, "hi": n_b_low if satisfied else n_b_high,
    })
    if satisfied:
        return {"n_b_min": n_b_low, "iterations": 1, "history": history}
    # 上界端检查：n_b_high 不满足 → 该 N_A 在 [0, n_b_max] 内不可行
    ev = evaluate(n_b_high)
    satisfied = ev["p_lower"] >= TARGET_P
    history.append({
        "iter": 2, "n_b": n_b_high, "connected": ev["connected"], "trials": ev["trials"],
        "p_hat": ev["p_hat"], "p_lower": ev["p_lower"],
        "action": "high" if satisfied else "low",
        "lo": n_b_low, "hi": n_b_high,
    })
    if not satisfied:
        return {"n_b_min": None, "iterations": 2, "history": history}
    lo, hi, it = n_b_low, n_b_high, 2
    while hi - lo > 1:
        if it >= max_iter:
            raise RuntimeError(
                f"二分 {max_iter} 次评估内未收敛到整数步长 1，当前区间 [{lo}, {hi}]（请增大 max_iter）")
        mid = (lo + hi) // 2
        ev = evaluate(mid)
        it += 1
        if ev["p_lower"] >= TARGET_P:
            hi, action = mid, "high"
        else:
            lo, action = mid, "low"
        history.append({
            "iter": it, "n_b": mid, "connected": ev["connected"], "trials": ev["trials"],
            "p_hat": ev["p_hat"], "p_lower": ev["p_lower"],
            "action": action, "lo": lo, "hi": hi,
        })
    return {"n_b_min": hi, "iterations": it, "history": history}


# 构造 Level 3 确认点：候选点及其 8 邻域（da∈{-step,0,step}、db∈{-50,0,50}，去重且坐标在界内）
def make_confirm_points(n_a_star, n_b_star, n_a_step, n_a_max, db=CONFIRM_DB):
    """返回候选点及其 8 邻域确认点列表（固定顺序：da 外层、db 内层，≤9 个点）。

    约束：0 <= n_a+da <= n_a_max 且 n_b+db >= 0；候选点本身（da=0, db=0）必在其中。
    """
    points = []
    for da in (-n_a_step, 0, n_a_step):
        for db_ in (-db, 0, db):
            na, nb = n_a_star + da, n_b_star + db_
            if 0 <= na <= n_a_max and nb >= 0:
                pt = (int(na), int(nb))
                if pt not in points:
                    points.append(pt)
    return points


# 写入问题 4 边界 CSV：每 (mode, N_A 网格点) 一行（n_b_min 不可行时 cost/p_hat/p_lower 记空）
def write_csv(rows, csv_path):
    """按列 boundary_mode,n_a,n_b_min,cost,p_hat,p_lower,iterations 写问题 4 边界 CSV。"""
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["boundary_mode", "n_a", "n_b_min", "cost", "p_hat", "p_lower",
                         "iterations"])
        for row in rows:
            writer.writerow([
                row["boundary_mode"], f"{row['n_a']:d}",
                "" if row["n_b_min"] is None else f"{row['n_b_min']:d}",
                "" if row["cost"] is None else f"{row['cost']:.4f}",
                "" if row["p_hat"] is None else f"{row['p_hat']:.6f}",
                "" if row["p_lower"] is None else f"{row['p_lower']:.6f}",
                f"{row['iterations']:d}",
            ])


# 统一 matplotlib 设置：Agg 后端 + 中文字体（Microsoft YaHei/SimHei）+ 负号显示
def setup_matplotlib():
    """设置 Agg 后端与中文字体（须在任何 pyplot 使用前调用一次）。"""
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False


# 绘制 Figure 6：λ_max 理论热图（N_B×N_A 颜色映射 + λ_max=1 理论可行区等高线）
def plot_lambda_max_heatmap(n_a_grid, n_b_grid, lam_grid, out_dir):
    """画 λ_max 理论筛选热图，保存 problem4_lambda_max_heatmap.png，返回路径（失败抛异常）。"""
    import matplotlib.pyplot as plt
    setup_matplotlib()
    X, Y = np.meshgrid(n_b_grid, n_a_grid)
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    pc = ax.pcolormesh(X, Y, lam_grid, shading="auto", cmap="viridis")
    cs = ax.contour(X, Y, lam_grid, levels=[1.0], colors="red", linewidths=1.6)
    ax.clabel(cs, fmt=r"$\lambda_{\max}=1$", fontsize=9)
    fig.colorbar(pc, ax=ax, label=r"$\lambda_{\max}$")
    ax.set_xlabel("N_B")
    ax.set_ylabel("N_A")
    ax.set_title(r"$\lambda_{\max}$ 理论筛选（红色等高线为理论可行区边界 $\lambda_{\max}=1$）")
    fig.tight_layout()
    png = os.path.join(out_dir, "problem4_lambda_max_heatmap.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 绘制 Figure 7：MC P_conn 热图（双模式各一子图，叠加 0.90 等高线与 Level2 边界阶梯线）
def plot_mc_pconn_heatmap(mc_grid, boundary_lines, out_dir):
    """按 mc_grid（mode->n_a_grid/n_b_grid/p_hat_grid）画 P_conn 热图，保存 problem4_mc_pconn_heatmap.png。"""
    import matplotlib.pyplot as plt
    setup_matplotlib()
    mode_keys = list(mc_grid.keys())
    fig, axes = plt.subplots(1, len(mode_keys), figsize=(13.0, 5.2), sharey=True)
    if len(mode_keys) == 1:
        axes = [axes]
    for ax, mkey in zip(axes, mode_keys):
        g = mc_grid[mkey]
        X = np.meshgrid(np.asarray(g["n_b_grid"]), np.asarray(g["n_a_grid"]))
        p = np.asarray(g["p_hat_grid"])
        pc = ax.pcolormesh(X[0], X[1], p, shading="auto", cmap="viridis", vmin=0.0, vmax=1.0)
        ax.contour(X[0], X[1], p, levels=[0.90], colors="white", linewidths=1.5)
        if mkey in boundary_lines:
            na_pts, nb_mins = boundary_lines[mkey]
            ax.step(na_pts, nb_mins, where="post", color="red", linewidth=1.6,
                    label="Level2 边界 N_B_min(N_A)")
        ax.set_xlabel("N_B")
        ax.set_title(mkey)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
        fig.colorbar(pc, ax=ax, label="P_conn")
    axes[0].set_ylabel("N_A")
    fig.suptitle("MC P_conn 热图（白线为 P_conn=0.90 等高线）")
    fig.tight_layout()
    png = os.path.join(out_dir, "problem4_mc_pconn_heatmap.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 绘制 Figure 8：成本等高线 + P=0.9 可行边界（N_A×N_B 网格 cost 等高线 + Level2 边界阶梯线 + 最优解标记）
def plot_cost_boundary(n_a_grid, n_b_grid, boundary_lines, optimal_points, out_dir):
    """画成本等高线与可行边界图，保存 problem4_cost_boundary.png，返回路径（失败抛异常）。"""
    import matplotlib.pyplot as plt
    setup_matplotlib()
    cost_grid = np.array([[total_cost(int(na), int(nb)) for nb in n_b_grid] for na in n_a_grid])
    X, Y = np.meshgrid(n_b_grid, n_a_grid)
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    cf = ax.contourf(X, Y, cost_grid, levels=20, cmap="viridis")
    ax.contour(X, Y, cost_grid, levels=10, colors="white", linewidths=0.6, alpha=0.6)
    fig.colorbar(cf, ax=ax, label="总成本 C（元）")
    colors = {"periodic_connected": "#ff7f0e", "wrapped_geometry_only": "#d62728"}
    for mkey, (na_pts, nb_mins) in boundary_lines.items():
        ax.step(na_pts, nb_mins, where="post", color=colors.get(mkey, "#333333"),
                linewidth=1.8, label=f"{mkey} 可行边界")
    for mkey, (na, nb) in optimal_points.items():
        ax.plot(nb, na, marker="*", markersize=14, color=colors.get(mkey, "#333333"),
                label=f"{mkey} 最优 ({na},{nb})")
    ax.set_xlabel("N_B")
    ax.set_ylabel("N_A")
    ax.set_title("成本等高线与 P_conn≥0.90 可行边界")
    ax.legend(fontsize=8)
    fig.tight_layout()
    png = os.path.join(out_dir, "problem4_cost_boundary.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 主流程：加载 V_AA/V_AB/V_BB → Level 1 理论筛选 → Level 2 MC 边界 + 单调性校验 →
# Level 2.5 MC 粗网格 → Level 3 成本优化与确认 → 写 CSV/JSON → 可选绘图 → 打印汇总
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
               else os.path.join(project_root, "results", "problem4"))
    os.makedirs(out_dir, exist_ok=True)
    total_t0 = time.perf_counter()

    # 数据准备：V_BB 解析值、V_AB 数值积分（固定种子可复现）、q_AB(r) 曲线（独立派生种子供论文画图）
    v_bb = v_bb_analytic()
    r_max = ab_max_contact_center_distance()
    v_ab = estimate_v_ab(r_max=r_max, dr=args.r_dr,
                         rng=np.random.default_rng(args.seed + 1000),
                         samples_per_r=args.ab_samples)
    v_ab_capsule = v_ab_capsule_approx()
    q_rng = np.random.default_rng(args.seed + 2000)
    q_r_grid = np.arange(0.0, r_max + args.r_dr * 0.5, args.r_dr)
    q_q_grid = [estimate_q_ab(float(r), q_rng, args.ab_samples) for r in q_r_grid]
    print(f"V_AA={v_aa:.4e} nm³ | V_BB={v_bb:.4e} nm³ | V_AB={v_ab:.4e} nm³ "
          f"| V_AB_capsule={v_ab_capsule:.4e} nm³ "
          f"(相对偏差 {abs(v_ab - v_ab_capsule) / v_ab_capsule * 100:.1f}%)", flush=True)

    # 双模式先 periodic 后 wrapped，与 solve_problem2/3 的枚举顺序一致
    modes = ([BoundaryMode.PERIODIC_CONNECTED, BoundaryMode.WRAPPED_GEOMETRY_ONLY]
             if args.mode == "all" else [BoundaryMode(args.mode)])

    # Level 1 理论筛选：λ_max 网格（N_A 步长 10、N_B 步长 100）与理论可行区（λ_max≥1）
    n_a_grid = np.arange(0, args.n_a_max + 1, L1_NA_STEP)
    n_b_grid = np.arange(0, args.n_b_max + 1, L1_NB_STEP)
    n_a_grid, n_b_grid, lam_grid = lambda_max_grid(n_a_grid, n_b_grid, v_aa, v_ab, v_bb)
    feas = lam_grid >= 1.0
    if feas.any():
        na_where, nb_where = np.nonzero(feas)
        feasible_region = {
            "n_a_min": int(n_a_grid[na_where.min()]), "n_a_max": int(n_a_grid[na_where.max()]),
            "n_b_min": int(n_b_grid[nb_where.min()]), "n_b_max": int(n_b_grid[nb_where.max()]),
        }
        print(f"λ_max 理论可行区（λ_max≥1）: N_A∈[{feasible_region['n_a_min']},"
              f"{feasible_region['n_a_max']}], N_B∈[{feasible_region['n_b_min']},"
              f"{feasible_region['n_b_max']}]（Level1 网格）", flush=True)
    else:
        feasible_region = None
        print("λ_max 理论可行区：空（Level1 网格上无 λ_max≥1 点）", flush=True)

    boundary_json, mc_json, confirm_json, optimal_json = None, None, None, None
    boundary_lines, optimal_points = {}, {}

    if not args.theory_only:
        # 共享进程池：Level 2/2.5/3 全部评估复用同一进程池（避免每次二分评估重建 8 个进程的开销）
        executor = (ProcessPoolExecutor(max_workers=args.workers)
                    if args.workers > 1 else None)
        # 可复现机制：主 rng 一次性派生 bisect/confirm/mc_grid 三个子种子池，评估顺序固定
        na_pts_l2 = list(range(0, args.n_a_max + 1, args.n_a_step))
        main_rng = np.random.default_rng(args.seed)
        bisect_pool = main_rng.integers(0, 2 ** 31,
                                        size=(len(modes), len(na_pts_l2),
                                              MAX_BISECT_ITER, args.trials))
        confirm_pool = main_rng.integers(0, 2 ** 31,
                                         size=(len(modes), CONFIRM_MAX_POINTS,
                                               args.confirm_trials))
        mc_na_pts = list(range(0, args.n_a_max + 1, args.mc_grid_na_step))
        mc_nb_pts = list(range(0, args.n_b_max + 1, args.mc_grid_nb_step))
        n_mc_pts = len(mc_na_pts) * len(mc_nb_pts)
        mc_grid_pool = main_rng.integers(0, 2 ** 31,
                                         size=(len(modes), n_mc_pts, args.mc_grid_trials))

        boundary_json, confirm_json, optimal_json = {}, {}, {}
        csv_rows = []
        for m_idx, mode in enumerate(modes):
            confirm_json[mode.value], optimal_json[mode.value] = None, None
            # Level 2：逐 N_A 网格点整数二分最小 N_B（判据 Wilson 单侧下界 ≥ 0.90）
            boundary_rows = []
            for k, na in enumerate(na_pts_l2):
                evaluator = make_nb_evaluator(bisect_pool[m_idx, k], int(na), mode,
                                              args.trials, args.workers, executor)
                bs = binary_search_min_nb(0, args.n_b_max, evaluator)
                last_high = next((h for h in reversed(bs["history"])
                                  if h["action"] == "high"), None)
                row = {
                    "n_a": int(na),
                    "n_b_min": bs["n_b_min"],
                    "cost": (float(total_cost(int(na), bs["n_b_min"]))
                             if bs["n_b_min"] is not None else None),
                    "p_hat": (float(last_high["p_hat"]) if last_high else None),
                    "p_lower": (float(last_high["p_lower"]) if last_high else None),
                    "iterations": bs["iterations"],
                }
                boundary_rows.append(row)
                csv_rows.append({"boundary_mode": mode.value, **row})

            # 边界单调性校验：N_A 增大 → N_B_min 不增（仅统计可行行）
            violations, prev = [], None
            for r in boundary_rows:
                if r["n_b_min"] is None:
                    continue
                if prev is not None and r["n_b_min"] > prev["n_b_min"]:
                    violations.append({"n_a_prev": prev["n_a"], "n_b_prev": prev["n_b_min"],
                                       "n_a_cur": r["n_a"], "n_b_cur": r["n_b_min"]})
                prev = r
            monotonic = not violations
            print(f"mode={mode.value}: 边界单调性 {'OK' if monotonic else '违反 ' + str(violations)}",
                  flush=True)
            boundary_json[mode.value] = {"rows": boundary_rows,
                                         "monotonic": monotonic,
                                         "monotonic_violations": violations}

            # Level 3：沿边界取成本最小候选，对候选及 8 邻域用 confirm_trials 复核
            feasible_rows = [r for r in boundary_rows if r["n_b_min"] is not None]
            candidate = min(feasible_rows, key=lambda r: r["cost"]) if feasible_rows else None
            boundary_json[mode.value]["candidate"] = candidate
            if candidate is None:
                print(f"mode={mode.value}: 无任何可行 N_A 网格点（Level 2 全部不可行），跳过 Level 3",
                      flush=True)
                continue
            print(f"mode={mode.value}: N_A=0 处 N_B_min={boundary_rows[0]['n_b_min']}, "
                  f"成本最小候选 N_A={candidate['n_a']} N_B={candidate['n_b_min']} "
                  f"cost={candidate['cost']:.4f}", flush=True)
            confirm_points = make_confirm_points(candidate["n_a"], candidate["n_b_min"],
                                                 args.n_a_step, args.n_a_max, db=CONFIRM_DB)
            confirms = []
            for pt_idx, (c_na, c_nb) in enumerate(confirm_points):
                seeds = confirm_pool[m_idx, pt_idx, :args.confirm_trials]
                connected = run_trials_mixed(seeds, c_na, c_nb, mode, args.workers, executor)
                confirms.append({
                    "n_a": c_na, "n_b": c_nb, "boundary_mode": mode.value,
                    "trials": args.confirm_trials, "connected_count": connected,
                    "p_hat": float(connected) / args.confirm_trials,
                    "p_lower": float(wilson_one_sided_lower(connected, args.confirm_trials)),
                    "cost": float(total_cost(c_na, c_nb)),
                })
            feasible_confirms = [c for c in confirms if c["p_lower"] >= TARGET_P]
            best = min(feasible_confirms, key=lambda c: c["cost"]) if feasible_confirms else None
            confirm_json[mode.value] = confirms
            optimal_json[mode.value] = ({
                "n_a": best["n_a"], "n_b": best["n_b"], "cost": best["cost"],
                "p_hat": best["p_hat"], "p_lower": best["p_lower"],
                "boundary_mode": mode.value,
            } if best else None)
            if best:
                print(f"mode={mode.value}: 最终最优 (N_A*, N_B*)=({best['n_a']},{best['n_b']}) "
                      f"cost={best['cost']:.4f} 元 p_hat={best['p_hat']:.4f} "
                      f"p_lower={best['p_lower']:.4f}", flush=True)
            else:
                print(f"mode={mode.value}: 警告！候选及邻域均未达到 p_lower>=0.90，确认点无解",
                      flush=True)
            for c in confirms:
                print(f"  confirm N_A={c['n_a']} N_B={c['n_b']} p_hat={c['p_hat']:.4f} "
                      f"p_lower={c['p_lower']:.4f} cost={c['cost']:.4f}", flush=True)
            boundary_lines[mode.value] = (
                [r["n_a"] for r in boundary_rows if r["n_b_min"] is not None],
                [r["n_b_min"] for r in boundary_rows if r["n_b_min"] is not None])
            if best:
                optimal_points[mode.value] = (best["n_a"], best["n_b"])

        # Level 2.5：MC 粗网格热图数据（逐模式、逐网格点固定顺序取 mc_grid_pool 行）
        mc_json = {}
        for m_idx, mode in enumerate(modes):
            p_grid = np.zeros((len(mc_na_pts), len(mc_nb_pts)))
            pt = 0
            for i, na in enumerate(mc_na_pts):
                for j, nb in enumerate(mc_nb_pts):
                    seeds = mc_grid_pool[m_idx, pt, :args.mc_grid_trials]
                    pt += 1
                    connected = run_trials_mixed(seeds, na, nb, mode, args.workers, executor)
                    p_grid[i, j] = float(connected) / args.mc_grid_trials
            mc_json[mode.value] = {"n_a_grid": mc_na_pts, "n_b_grid": mc_nb_pts,
                                   "p_hat_grid": p_grid.tolist()}
            print(f"mode={mode.value}: MC 粗网格 {len(mc_na_pts)}×{len(mc_nb_pts)} 点完成",
                  flush=True)

        # 写边界 CSV（每 N_A 网格边界一行）
        write_csv(csv_rows, os.path.join(out_dir, "problem4_result.csv"))
        if executor is not None:
            executor.shutdown()

    # 汇总结果 JSON（--theory-only 时 boundary/mc_grid/confirm/optimal 为 None）
    result = {
        "seed": args.seed,
        "config": vars(args),
        "V_AA": float(v_aa), "V_AB": float(v_ab), "V_BB": float(v_bb),
        "V_AB_capsule": float(v_ab_capsule),
        "q_ab": {"r_grid": q_r_grid.tolist(), "q_grid": q_q_grid},
        "theory_grid": {
            "n_a_grid": n_a_grid.tolist(), "n_b_grid": n_b_grid.tolist(),
            "lam_grid": lam_grid.tolist(), "feasible": feas.tolist(),
            "feasible_region": feasible_region,
        },
        "boundary": boundary_json,
        "mc_grid": mc_json,
        "confirm": confirm_json,
        "optimal": optimal_json,
    }
    json_path = os.path.join(out_dir, "problem4_result.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    total_time = time.perf_counter() - total_t0
    print(f"\n总耗时: {total_time:.1f}s ({total_time / 60.0:.1f} min)")
    if args.theory_only:
        print("theory-only 模式：跳过 Level 2/2.5/3（未写边界 CSV）")
    print(f"JSON 已保存: {json_path}")

    # 绘图（失败不阻断结果文件）
    if args.plot:
        try:
            png6 = plot_lambda_max_heatmap(n_a_grid, n_b_grid, lam_grid, out_dir)
            print(f"图已保存: {png6}")
        except Exception as exc:
            print(f"绘图失败（结果文件已保存）: {exc}")
        if not args.theory_only:
            try:
                png7 = plot_mc_pconn_heatmap(mc_json, boundary_lines, out_dir)
                print(f"图已保存: {png7}")
            except Exception as exc:
                print(f"绘图失败（结果文件已保存）: {exc}")
            try:
                png8 = plot_cost_boundary(n_a_grid, n_b_grid, boundary_lines,
                                          optimal_points, out_dir)
                print(f"图已保存: {png8}")
            except Exception as exc:
                print(f"绘图失败（结果文件已保存）: {exc}")


if __name__ == "__main__":
    main()