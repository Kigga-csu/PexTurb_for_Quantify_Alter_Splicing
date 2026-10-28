#!/usr/bin/env python3
"""build_tpm_matrix.py — 把每个样本的 salmon quant.sf 拼成 PexTurb 需要的 TPM 矩阵。

列来源: quant.sf 的 TPM 列（第 4 列）。
输出格式（gtf_index.py / calc_psi.py 直接可用）:
    第一列: transcript_id（保留版本号，gtf_index.py 会自行去掉 .版本）
    其余列: 每个样本一列，列名必须与 anno 文件里的样本名一致

过滤规则与生产脚本 get_tpm.py 一致:
    某个转录本在 >= max_low_frac 比例的样本中 TPM < min_tpm 时，整行丢弃。
    默认: TPM < 1 的样本占比 >= 50% 则丢弃。

用法 A（推荐，目录结构 output/ours/<sample>/quant.sf）:
    python build_tpm_matrix.py --quant-root output/ours --output isomantis_tpm.txt

用法 B（显式列表，两列: sample<TAB>quant.sf）:
    python build_tpm_matrix.py --quant-list quant_files.tsv --output isomantis_tpm.txt

关闭过滤:
    --keep-all
"""
import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--quant-root", help="目录，其下每个子目录是一个样本，内含 quant.sf")
    p.add_argument("--quant-list", help="两列 TSV: sample_id<TAB>quant.sf 路径")
    p.add_argument("--output", required=True, help="输出 TPM 矩阵路径")
    p.add_argument("--min-tpm", type=float, default=1.0,
                   help="低表达判定阈值（默认 1.0）")
    p.add_argument("--max-low-frac", type=float, default=0.5,
                   help="低表达样本占比达到该值则丢弃该转录本（默认 0.5）")
    p.add_argument("--keep-all", action="store_true", help="不做低表达过滤")
    return p.parse_args()


def collect_from_root(root):
    items = []
    for name in sorted(os.listdir(root)):
        q = os.path.join(root, name, "quant.sf")
        if os.path.isfile(q):
            items.append((name, q))
    return items


def collect_from_list(path):
    items = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"列表行格式错误: {line!r}")
            items.append((parts[0], parts[1]))
    return items


def read_quant_sf(path):
    """返回 {transcript: TPM}，保持文件顺序。"""
    tpm = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if "TPM" in header:
            col = header.index("TPM")
        else:                     # 兼容无表头: Name Length EffectiveLength TPM NumReads
            col = 3
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= col:
                continue
            tpm[parts[0]] = float(parts[col])
    return tpm


def main():
    args = parse_args()
    if bool(args.quant_root) == bool(args.quant_list):
        sys.exit("请且仅请提供 --quant-root 或 --quant-list 之一")

    items = collect_from_root(args.quant_root) if args.quant_root else collect_from_list(args.quant_list)
    if not items:
        sys.exit("没有找到任何 quant.sf")
    samples = [s for s, _ in items]
    if len(set(samples)) != len(samples):
        sys.exit("样本名重复，请检查输入")

    print(f"[*] 样本数: {len(samples)}")
    mats = []
    for sample, qpath in items:
        tpm = read_quant_sf(qpath)
        mats.append(tpm)
        print(f"    {sample}: {len(tpm)} transcripts")

    ref = list(mats[0].keys())
    ref_set = set(ref)
    for sample, m in zip(samples[1:], mats[1:]):
        if set(m.keys()) != ref_set:
            missing = len(ref_set - set(m))
            extra = len(set(m) - ref_set)
            print(f"[!] {sample} 转录本集合与第一个样本不同 (missing={missing}, extra={extra})；"
                  f"按第一个样本的顺序取交集", file=sys.stderr)
            ref = [t for t in ref if t in m]
            ref_set = set(ref)

    n = len(samples)
    threshold = args.max_low_frac * n
    kept = 0
    with open(args.output, "w") as out:
        out.write("transcript_id\t" + "\t".join(samples) + "\n")
        for t in ref:
            vals = [m[t] for m in mats]
            if not args.keep_all:
                low = sum(1 for v in vals if v < args.min_tpm)
                if low >= threshold:
                    continue
            out.write(t + "\t" + "\t".join(f"{v:.6g}" for v in vals) + "\n")
            kept += 1

    print(f"[ok] {args.output}: {kept} transcripts x {n} samples"
          f"{' (未过滤)' if args.keep_all else ''}")


if __name__ == "__main__":
    main()
