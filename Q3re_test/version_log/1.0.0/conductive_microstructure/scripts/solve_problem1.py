# 求解问题1：对附件三个分表，在两种边界模式下判定左右电极是否导通并输出统计与贯通路径
import argparse
import json
import os
import sys

import numpy as np

# 允许从 scripts/ 直接运行（把项目根目录加入 sys.path）
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.config import BoundaryMode
from src.graph.connectivity import analyze_group
from src.io.attachment_reader import read_attachment


# 解析命令行参数：附件路径与要运行的边界模式
def parse_args():
    parser = argparse.ArgumentParser(description="问题1：微构体导通判定（附件三组数据）")
    parser.add_argument("--xlsx", default=os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "附件.xlsx")),
        help="附件Excel路径（默认工作区根目录 附件.xlsx）")
    parser.add_argument("--mode", default="all",
                        choices=["all", "periodic_connected", "wrapped_geometry_only"],
                        help="边界模式：all=两种对比输出")
    return parser.parse_args()


# 把一组端点数据按指定边界模式做一次完整分析并格式化输出文本
def run_one_group(name, endpoints, mode):
    """分析单个分表并返回 (分析结果dict, 可打印文本)。"""
    res = analyze_group(endpoints, mode)
    lines = []
    lines.append(f"[{name}] 边界模式 {mode.value}: "
                 f"{'导通' if res['connected'] else '不导通'}")
    lines.append(f"  原始行数(片段): {res['segment_count']}, "
                 f"导体节点数: {res['node_count']}, 导体-导体连接边数: {res['edge_count']}")
    lines.append(f"  左电极直连节点数: {len(res['left_node_ids'])}, "
                 f"右电极直连节点数: {len(res['right_node_ids'])}")
    if res["merged_pairs"]:
        merged_txt = "; ".join("+".join(str(r) for r in grp)
                               for grp in res["merged_pairs"])
        lines.append(f"  跨边界合并的行组: {merged_txt}")
    if res["connected"]:
        node_path = " -> ".join(str(nd) for nd in res["path_node_ids"])
        lines.append(f"  贯通路径(节点编号): {node_path}")
        lines.append(f"  贯通路径(原始行): {res['path_row_ids']}")
    return res, "\n".join(lines)


# 主流程：读取附件，逐组逐模式分析，输出文本并保存 JSON
def main():
    args = parse_args()
    data = read_attachment(args.xlsx)
    modes = ([BoundaryMode.PERIODIC_CONNECTED, BoundaryMode.WRAPPED_GEOMETRY_ONLY]
             if args.mode == "all" else [BoundaryMode(args.mode)])
    out_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results", "problem1"))
    os.makedirs(out_dir, exist_ok=True)
    summary = {"xlsx": args.xlsx, "groups": {}}
    for name in data.keys():
        group_out = {}
        for mode in modes:
            res, text = run_one_group(name, data[name], mode)
            print(text)
            print()
            group_out[mode.value] = {
                "connected": res["connected"],
                "node_count": res["node_count"],
                "segment_count": res["segment_count"],
                "edge_count": res["edge_count"],
                "left_node_ids": res["left_node_ids"],
                "right_node_ids": res["right_node_ids"],
                "path_node_ids": res["path_node_ids"],
                "path_row_ids": res["path_row_ids"],
                "merged_pairs": [list(map(int, g)) for g in res["merged_pairs"]],
            }
        summary["groups"][name] = group_out
    json_path = os.path.join(out_dir, "problem1_result.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(f"结果已保存: {json_path}")


if __name__ == "__main__":
    main()
