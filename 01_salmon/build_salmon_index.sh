#!/usr/bin/env bash
# 构建 salmon 转录组索引（PexTurb 生产环境使用的索引参数：k=31，GENCODE 格式，含 decoy）
#
# 用法:
#   bash build_salmon_index.sh <transcripts.fa.gz> <salmon_index_dir> [decoys.txt]
#
# 示例（有 decoy，与服务器 /data/lab/wangshixian/genome/salmon_index 的配置一致）:
#   bash build_salmon_index.sh gencode.vXX.transcripts.fa.gz salmon_index decoys.txt
#
# 服务器上现成索引（无需重建）:
#   /data/lab/wangshixian/genome/salmon_index     # info.json: k=31, num_decoys=194
set -euo pipefail

TRANSCRIPTOME=${1:?需要转录组 FASTA}
INDEX_DIR=${2:?需要输出索引目录}
DECOYS=${3:-}

mkdir -p "$INDEX_DIR"

if [ -n "$DECOYS" ]; then
    salmon index \
        -t "$TRANSCRIPTOME" \
        -i "$INDEX_DIR" \
        -k 31 \
        --gencode \
        -d "$DECOYS" \
        -p "${SALMON_THREADS:-8}"
else
    salmon index \
        -t "$TRANSCRIPTOME" \
        -i "$INDEX_DIR" \
        -k 31 \
        --gencode \
        -p "${SALMON_THREADS:-8}"
fi

echo "[ok] salmon index -> $INDEX_DIR"
