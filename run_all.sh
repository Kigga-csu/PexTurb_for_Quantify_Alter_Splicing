#!/usr/bin/env bash
# PexTurb 核心流程：GTF 建索引 -> 算 EPP/PSI -> 差异检验
#
# 从 FASTQ 开始时，先使用 01_salmon/ 量化转录本表达。
#
# 用法:
#   bash run_all.sh \
#       --tpm   isomantis_tpm.txt \
#       --gtf   Homo_sapiens.GRCh38.77.sorted.gtf \
#       --anno  samples.tsv \
#       --outdir result \
#       [--min-group-cnt 1]
#
# anno 格式: 两列 sample<TAB>group，且必须恰好两组。
#   group1 是文件里先出现的那组，delta = mean(group1) - mean(group2)。
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

TPM=""
GTF=""
ANNO=""
OUTDIR="result"
MIN_GROUP_CNT=1

while [ $# -gt 0 ]; do
    case "$1" in
        --tpm)    TPM="$2"; shift 2 ;;
        --gtf)    GTF="$2"; shift 2 ;;
        --anno)   ANNO="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --min-group-cnt) MIN_GROUP_CNT="$2"; shift 2 ;;
        -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
        *) echo "未知参数: $1" >&2; exit 1 ;;
    esac
done

[ -n "$TPM" ] && [ -n "$GTF" ] && [ -n "$ANNO" ] || {
    echo "必须提供 --tpm / --gtf / --anno" >&2; exit 1;
}

mkdir -p "$OUTDIR"
PREFIX=$(basename "$TPM" .txt)
INDEX="$OUTDIR/${PREFIX}_index.txt"
PSI="$OUTDIR/${PREFIX}_psi.txt"
OUT="$OUTDIR/${PREFIX}_pexturb.txt"

echo "==========================================================="
echo " [1/3] 根据 GTF 建事件索引 (gtf_index.py)"
echo "==========================================================="
python "$HERE/code_pexturb/gtf_index.py" \
    --tpm "$TPM" --gtf "$GTF" --index "$INDEX"

echo "==========================================================="
echo " [2/3] 计算 EPP/PSI (calc_psi.py)"
echo "==========================================================="
python "$HERE/code_pexturb/calc_psi.py" \
    --tpm "$TPM" --index "$INDEX" --psi "$PSI"

echo "==========================================================="
echo " [3/3] 差异检验 (pexturb_test.py)"
echo "==========================================================="
python "$HERE/code_pexturb/pexturb_test.py" \
    --psi "$PSI" --anno "$ANNO" --output "$OUT" \
    --min_group_cnt "$MIN_GROUP_CNT"

echo "==========================================================="
echo " 完成: $OUT"
echo " 索引: $INDEX"
echo " PSI : $PSI"
echo "==========================================================="
