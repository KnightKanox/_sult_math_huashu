# 生成问题 2 的概率云图：q_AA(r) 与 4πr²q(r)
import csv
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False


# 读取 cloud_AA.csv，返回 r、q、置信区间下上界数组
def read_cloud_csv(csv_path):
    """读取 cloud_AA.csv，返回 (r, q, ci_low, ci_high) 四个 numpy 数组。"""
    r, q, lo, hi = [], [], [], []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            r.append(float(row["r_nm"]))
            q.append(float(row["q_hat"]))
            lo.append(float(row["ci_low"]))
            hi.append(float(row["ci_high"]))
    return np.array(r), np.array(q), np.array(lo), np.array(hi)


# 绘制局部连接概率云 q_AA(r)
def plot_q_aa(cloud_csv, out_dir):
    """绘制 q_AA(r) 图并保存为 problem2_q_aa.png。"""
    r, q, lo, hi = read_cloud_csv(cloud_csv)
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.plot(r, q, color="#1f77b4", linewidth=1.8, label=r"$q_{AA}(r)$")
    ax.fill_between(r, lo, hi, color="#1f77b4", alpha=0.18, label="95% CI")
    ax.set_xlabel(r"中心距 $r$ (nm)")
    ax.set_ylabel(r"局部连接概率 $q_{AA}(r)$")
    ax.set_title(r"介质 A 的局部连接概率云 $q_{AA}(r)$")
    ax.set_xlim(0, float(r.max()))
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, "problem2_q_aa.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 绘制概率云积分核 4πr²q(r)
def plot_4pir2q(cloud_csv, out_dir):
    """绘制 4πr²q_AA(r) 图并保存为 problem2_4pir2q.png。"""
    r, q, _, _ = read_cloud_csv(cloud_csv)
    kernel = 4.0 * np.pi * r * r * q
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.plot(r, kernel, color="#d62728", linewidth=1.8)
    ax.fill_between(r, 0.0, kernel, color="#d62728", alpha=0.18)
    ax.set_xlabel(r"中心距 $r$ (nm)")
    ax.set_ylabel(r"$4\pi r^2 q_{AA}(r)$")
    ax.set_title(r"概率云积分核 $4\pi r^2 q_{AA}(r)$")
    ax.set_xlim(0, float(r.max()))
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, "problem2_4pir2q.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# 命令行入口：生成两张概率云图
def main():
    """从 results/cloud/cloud_AA.csv 生成问题 2 的两张概率云图。"""
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    cloud_csv = os.path.join(project_root, "results", "cloud", "cloud_AA.csv")
    out_dir = os.path.join(project_root, "results", "figures")
    if not os.path.isfile(cloud_csv):
        print(f"缺少概率云数据: {cloud_csv}")
        sys.exit(1)
    png1 = plot_q_aa(cloud_csv, out_dir)
    png2 = plot_4pir2q(cloud_csv, out_dir)
    print("图已保存: " + png1)
    print("图已保存: " + png2)


if __name__ == "__main__":
    main()
