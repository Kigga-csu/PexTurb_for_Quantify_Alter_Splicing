#!/usr/bin/env bash
# 批量 salmon quant —— 与服务器生产脚本保持一致：
#   /data/lab/wangshixian/IsoMantis_AD_all/M_EV/quant_ours.sh
#
# 用法:
#   bash run_salmon.sh --samples samples.tsv --index /path/to/salmon_index --outdir output/ours
#
# samples.tsv 每行 2 或 3 列（制表符分隔，无表头）:
#   双端: sample_id <TAB> R1.fastq.gz <TAB> R2.fastq.gz
#   单端: sample_id <TAB> reads.fastq.gz
#
# 输出: <outdir>/<sample_id>/quant.sf
set -euo pipefail

SAMPLES=""
INDEX=""
OUTDIR="output/ours"
THREADS="${SALMON_THREADS:-8}"

while [ $# -gt 0 ]; do
    case "$1" in
        --samples) SAMPLES="$2"; shift 2 ;;
        --index)   INDEX="$2";   shift 2 ;;
        --outdir)  OUTDIR="$2";  shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
        *) echo "未知参数: $1" >&2; exit 1 ;;
    esac
done

[ -n "$SAMPLES" ] && [ -n "$INDEX" ] || { echo "必须提供 --samples 和 --index" >&2; exit 1; }
command -v salmon >/dev/null 2>&1 || { echo "找不到 salmon，请先安装/激活环境" >&2; exit 1; }

mkdir -p "$OUTDIR"

while IFS=$'\t' read -r SAMPLE R1 R2; do
    [ -z "${SAMPLE:-}" ] && continue
    case "$SAMPLE" in \#*) continue ;; esac

    OUT="$OUTDIR/$SAMPLE"
    mkdir -p "$OUT"

    if [ -n "${R2:-}" ]; then
        echo "[quant] $SAMPLE (paired-end)"
        salmon quant \
            -i "$INDEX" \
            -l A \
            -1 "$R1" -2 "$R2" \
            --gcBias --seqBias --validateMappings \
            -p "$THREADS" \
            -o "$OUT"
    else
        echo "[quant] $SAMPLE (single-end)"
        salmon quant \
            -i "$INDEX" \
            -l A \
            -r "$R1" \
            --gcBias --seqBias --validateMappings \
            -p "$THREADS" \
            -o "$OUT"
    fi
done < "$SAMPLES"

echo "[ok] quant 完成 -> $OUTDIR/<sample>/quant.sf"
