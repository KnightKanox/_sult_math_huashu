# 问题 2 可视化：导通概率 P_conn 与平均连接度 k̄ 随体积分数变化的图
# 输入 solve_problem2.py 输出的 CSV，输出两张 PNG：problem2_p_conn.png（含 95% CI
# 误差棒，两种边界模式各一条线）与 problem2_kbar.png（k̄ 与边界模式无关仅画一条）。
# 使用 Agg 后端（无显示环境可运行），中文字体 Microsoft YaHei/SimHei，dpi=150。
import argparse
import csv
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")  # 无显示环境时使用 Agg 后端（须在导入 pyplot 前设置）
import matplotlib.pyplot as plt

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False


# 读取 solve_problem2 结果 CSV，返回行字典列表（phi 为十进制小数）
def read_problem2_csv(csv_path):
    """读取问题 2 结果 CSV，返回行字典列表。"""
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


# 绘制图 1：导通概率 P_conn 随体积分数变化（两种边界模式各一条线，带 95% CI 误差棒）
def plot_p_conn(csv_path, out_dir):
    """读结果 CSV 画 P_conn vs φ（误差棒+双模式），保存 problem2_p_conn.png，返回路径。"""
    os.makedirs(out_dir, exist_ok=True)
    rows = read_problem2_csv(csv_path)
    phis = sorted({float(r["phi"]) for r in rows})
    x = np.array([p * 100.0 for p in phis])  # 转为百分数显示
    modes = sorted({r["boundary_mode"] for r in rows})

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for m in modes:
        g = sorted([r for r in rows if r["boundary_mode"] == m],
                   key=lambda r: float(r["phi"]))
        y = np.array([float(r["p_hat"]) for r in g])
        lo = np.array([float(r["ci_low"]) for r in g])
        hi = np.array([float(r["ci_high"]) for r in g])
        ax.errorbar(x, y, yerr=[y - lo, hi - y], marker="o", capsize=4, label=m)
    ax.set_xlabel("体积分数 φ (%)")
    ax.set_ylabel("导通概率 P_conn")
    ax.set_title("P_conn 随体积分数变化（含 95% CI）")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:.2f}%" for v in x])
    ax.set_ylim(-0.03, 1.03)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    png = os.path.join(out_dir, "problem2_p_conn.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 绘制图 2：平均连接度 k̄ 随体积分数变化（k̄ 与边界模式无关，取第一种模式的数据画单条线）
def plot_kbar(csv_path, out_dir):
    """读结果 CSV 画 k̄ vs φ（单条线），保存 problem2_kbar.png，返回路径。"""
    os.makedirs(out_dir, exist_ok=True)
    rows = read_problem2_csv(csv_path)
    phis = sorted({float(r["phi"]) for r in rows})
    x = np.array([p * 100.0 for p in phis])  # 转为百分数显示
    modes = sorted({r["boundary_mode"] for r in rows})
    g = sorted([r for r in rows if r["boundary_mode"] == modes[0]],
               key=lambda r: float(r["phi"]))

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    y = np.array([float(r["mean_degree_theory"]) for r in g])
    ax.plot(x, y, marker="s", color="#d62728", label=r"$\bar{k}$（两种边界模式相同）")
    ax.set_xlabel("体积分数 φ (%)")
    ax.set_ylabel(r"平均连接度 $\bar{k}$")
    ax.set_title(r"平均连接度 $\bar{k}$ 随体积分数变化")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:.2f}%" for v in x])
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    png = os.path.join(out_dir, "problem2_kbar.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 绘制图 3：最大连通分量占比 S_max/N 随体积分数变化（两种边界模式各一条线）
def plot_smax(csv_path, out_dir):
    """读结果 CSV 画 S_max/N vs φ（双模式），保存 problem2_smax_ratio.png，返回路径。"""
    os.makedirs(out_dir, exist_ok=True)
    rows = read_problem2_csv(csv_path)
    phis = sorted({float(r["phi"]) for r in rows})
    x = np.array([p * 100.0 for p in phis])
    modes = sorted({r["boundary_mode"] for r in rows})

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for m in modes:
        g = sorted([r for r in rows if r["boundary_mode"] == m],
                   key=lambda r: float(r["phi"]))
        y = np.array([float(r["mean_max_component_ratio"]) for r in g])
        ax.plot(x, y, marker="o", label=m)
    ax.set_xlabel("体积分数 φ (%)")
    ax.set_ylabel(r"最大连通分量占比 $S_{\max}/N$")
    ax.set_title(r"$S_{\max}/N$ 随体积分数变化")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:.2f}%" for v in x])
    ax.set_ylim(-0.03, 1.03)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    png = os.path.join(out_dir, "problem2_smax_ratio.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 命令行入口：python -m src.visualization.probability_plot --csv ... --out-dir ...
def main():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    parser = argparse.ArgumentParser(description="问题2结果可视化（读 solve_problem2 输出 CSV 绘图）")
    parser.add_argument("--csv", default=os.path.join(
        project_root, "results", "problem2", "problem2_result.csv"),
        help="solve_problem2 输出的 CSV 路径")
    parser.add_argument("--out-dir", default=os.path.join(project_root, "results", "figures"),
                        help="图片输出目录")
    args = parser.parse_args()
    png1 = plot_p_conn(args.csv, args.out_dir)
    png2 = plot_kbar(args.csv, args.out_dir)
    png3 = plot_smax(args.csv, args.out_dir)
    print(f"图已保存: {png1}")
    print(f"图已保存: {png2}")
    print(f"图已保存: {png3}")


if __name__ == "__main__":
    main()
