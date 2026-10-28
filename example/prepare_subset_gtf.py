#!/usr/bin/env python3
"""prepare_subset_gtf.py — 从完整 GTF 中抽出示例数据涉及基因的子集 GTF。

保留规则:
  1. 第一遍：找出 TPM 矩阵中出现过的转录本，记录它们所属的 gene_id；
  2. 第二遍：写出这些 gene 的 gene / transcript / exon 记录（同基因的其它转录本也保留，
     与生产流程在完整 GTF 上建索引的行为一致）。

用法:
    python prepare_subset_gtf.py \
        --tpm example/data/simulated_tpm.txt \
        --gtf Homo_sapiens.GRCh38.77.sorted.gtf \
        --output example/data/simulated_annotation.gtf
"""
import argparse


def parse_attrs(field):
    attrs = {}
    for item in field.split("; "):
        item = item.strip()
        if not item:
            continue
        parts = item.split(" ")
        if len(parts) >= 2:
            attrs[parts[0]] = parts[1].strip('"')
    return attrs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tpm", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    tx_ids = set()
    with open(args.tpm) as f:
        f.readline()
        for line in f:
            tx = line.split("\t", 1)[0].strip()
            if tx:
                tx_ids.add(tx.split(".")[0])
    print(f"[*] TPM 中转录本数: {len(tx_ids)}")

    gene_ids = set()
    with open(args.gtf) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] not in ("transcript", "exon"):
                continue
            attrs = parse_attrs(parts[8])
            tx = attrs.get("transcript_id", "").split(".")[0]
            if tx in tx_ids:
                gene_ids.add(attrs.get("gene_id", ""))
    gene_ids.discard("")
    print(f"[*] 涉及基因数: {len(gene_ids)}")

    n_out = 0
    with open(args.gtf) as f, open(args.output, "w") as out:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            feature = parts[2]
            if feature not in ("gene", "transcript", "exon"):
                continue
            attrs = parse_attrs(parts[8])
            if attrs.get("gene_id", "") not in gene_ids:
                continue
            out.write(line if line.endswith("\n") else line + "\n")
            n_out += 1

    print(f"[ok] {args.output}: {n_out} lines")


if __name__ == "__main__":
    main()
