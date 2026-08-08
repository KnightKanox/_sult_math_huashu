# 生成"仿真 vs 概率云理论"对比图（问题 2）
# 左图：双 y 轴对照——仿真导通概率 P_conn（含 95% CI 误差棒，两种边界模式）
#       与概率云理论平均连接度 k̄（红色虚线），叠加 k̄=1 理论临界线，
#       直观展示"理论临界点与仿真跃迁区对齐"；
# 右图：P̂ 随 k̄ 的相图（标注各 φ 点），k̄=1 附近 P̂ 快速上升。
# 输出 results/figures/problem2_theory_vs_simulation.png（Agg 后端，中文字体，dpi=150）。
import csv
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")  # 无显示环境时使用 Agg 后端（须在导入 pyplot 前设置）
import matplotlib.pyplot as plt

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

# 理论临界点：平均连接度 k̄ = 1（平均场渗流阈值）
KBAR_CRITICAL = 1.0


# 读取 solve_problem2 结果 CSV，返回行字典列表
def read_problem2_csv(csv_path):
    """读取问题 2 结果 CSV，返回行字典列表。"""
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


# 绘制对比图（左：双轴 P̂/k̄ vs φ；右：P̂ vs k̄ 相图），保存 PNG 并返回路径
def plot_theory_vs_simulation(csv_path, out_dir):
    """生成对比图：左侧双轴对照仿真 P̂ 与理论 k̄ 并标出 k̄=1 临界线，右侧 P̂-k̄ 相图。"""
    os.makedirs(out_dir, exist_ok=True)
    rows = read_problem2_csv(csv_path)
    modes = sorted({r["boundary_mode"] for r in rows})

    # 按 phi 升序整理各边界模式数据
    by_mode = {}
    for m in modes:
        by_mode[m] = sorted([r for r in rows if r["boundary_mode"] == m],
                            key=lambda r: float(r["phi"]))

    phis = [float(r["phi"]) * 100.0 for r in by_mode[modes[0]]]  # 转为百分数
    kbar = [float(r["mean_degree_theory"]) for r in by_mode[modes[0]]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    # ---- 左图：双 y 轴对照（仿真 P̂ vs 理论 k̄）----
    for m in modes:
        g = by_mode[m]
        y = np.array([float(r["p_hat"]) for r in g])
        lo = np.array([float(r["ci_low"]) for r in g])
        hi = np.array([float(r["ci_high"]) for r in g])
        ax1.errorbar(phis, y, yerr=[y - lo, hi - y], marker="o", capsize=4,
                     label="仿真 P_conn（" + m + "）")
    ax1.set_xlabel("体积分数 φ (%)")
    ax1.set_ylabel("导通概率 P_conn（仿真）")
    ax1.set_xticks(phis)
    ax1.set_xticklabels(["{:.2f}%".format(v) for v in phis])
    ax1.set_ylim(-0.03, 1.03)
    ax1.grid(True, alpha=0.3)

    # 右轴：概率云理论平均连接度 k̄
    ax1b = ax1.twinx()
    ax1b.plot(phis, kbar, marker="s", color="#d62728", linestyle="--",
              label=r"理论 $\bar{k}$（概率云）")
    ax1b.axhline(KBAR_CRITICAL, color="gray", linestyle=":", linewidth=1.2)
    ax1b.text(phis[0] - 0.05, KBAR_CRITICAL + 0.03, r"$\bar{k}=1$（理论临界）",
              color="gray", fontsize=9)
    ax1b.set_ylabel(r"平均连接度 $\bar{k}$（概率云理论）")
    ax1b.set_ylim(0, 2.2)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=9)
    ax1.set_title(r"仿真 P_conn 与理论 $\bar{k}$ 对照（临界区对齐）")

    # ---- 右图：P̂ vs k̄ 相图（用第一种边界模式的数据）----
    g = by_mode[modes[0]]
    x = np.array([float(r["mean_degree_theory"]) for r in g])
    y = np.array([float(r["p_hat"]) for r in g])
    ax2.axvline(KBAR_CRITICAL, color="gray", linestyle=":", linewidth=1.2)
    ax2.scatter(x, y, s=70, color="#1f77b4", zorder=3)
    for xi, yi, phi in zip(x, y, phis):
        txt = "φ={:.2f}%\n$\\bar{{k}}$={:.2f}\n$\\hat{{P}}$={:.3f}".format(phi, xi, yi)
        ax2.annotate(txt, (xi, yi), textcoords="offset points",
                     xytext=(12, -10), fontsize=8)
    ax2.text(KBAR_CRITICAL + 0.03, 0.9,
             r"理论临界区($\bar{k}$=1 附近 $\hat{P}$ 快速上升)",
             color="gray", fontsize=9, va="top")
    ax2.set_xlabel(r"平均连接度 $\bar{k}$（概率云理论）")
    ax2.set_ylabel("导通概率 P_conn（仿真）")
    ax2.set_ylim(0.85, 1.01)
    ax2.grid(True, alpha=0.3)
    ax2.set_title(r"相图：仿真 $\hat{P}$ 随理论 $\bar{k}$ 的演化")

    fig.suptitle("问题2：Monte Carlo 仿真 vs 概率云理论", fontsize=13, y=1.0)
    fig.tight_layout()
    png = os.path.join(out_dir, "problem2_theory_vs_simulation.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 命令行入口：python scripts/plot_theory_vs_simulation.py
def main():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    csv_path = os.path.join(project_root, "results", "problem2", "problem2_result.csv")
    out_dir = os.path.join(project_root, "results", "figures")
    png = plot_theory_vs_simulation(csv_path, out_dir)
    print("对比图已保存: " + png)


if __name__ == "__main__":
    main()
