#!/bin/bash
set -e
set -o pipefail

# 强制使用 UTF-8
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONIOENCODING=UTF-8

# ==============================================================================
# 通用测试脚本：比较 原始输出(Minimap2)、MAPQ过滤 和 OvFilter
# ==============================================================================

# --- 1. 核心输入文件配置（可根据需要修改） ---
# 默认使用测试目录下的数据，也可以通过命令行参数传入：
# ./run.sh reads.fastq ref.fasta truth.maf [ont|hifi]
READS_FASTQ=${1:-"path/sim_data/chr8/hifi/chr8_hifi_0001.fastq"}
REF_FASTA=${2:-"path/sim_data/chr8/hifi/chr8_hifi_0001.ref"}
TRUTH_MAF=${3:-"path/sim_data/chr8/hifi/chr8_hifi_0001.maf"}
READS_TYPE=${4:-"hifi"}  # hifi 或 ont

# --- 2. 输出目录与程序路径 ---
OUT_DIR="path/test_result_${READS_TYPE}"
OVFILTER_CPP_BIN="path/ovfilter"
EVAL_PY="path/evaluate_alignments.py"

# --- 3. 通用运行参数 ---
THREADS=32
ONLY_PRIMARY="false" # true 表示 minimap2 使用 --secondary=no

# --- 4. K-mer 参数 ---
OV_K=21

echo "========================================================="
echo "  开始执行通用评估流程"
echo "  数据类型：${READS_TYPE}"
echo "  输出目录：${OUT_DIR}"
echo "  Reads:   ${READS_FASTQ}"
echo "  Ref:     ${REF_FASTA}"
echo "  Truth:   ${TRUTH_MAF}"
echo "========================================================="

if [ ! -f "${READS_FASTQ}" ] || [ ! -f "${REF_FASTA}" ] || [ ! -f "${TRUTH_MAF}" ]; then
    echo "错误: 输入文件不存在，请检查路径！"
    echo "用法: ./run.sh [reads.fastq] [ref.fasta] [truth.maf] [ont|hifi]"
    exit 1
fi

mkdir -p "${OUT_DIR}"

MINIMAP2_PAF="${OUT_DIR}/minimap2_raw.paf"
JF_COUNT="${OUT_DIR}/mer_counts.jf"
KMER_FILE="${OUT_DIR}/target_u1.kmer"
OVFILTER_PAF="${OUT_DIR}/ovfilter.paf"
OVFILTER_TIME_LOG="${OUT_DIR}/ovfilter_time.log"

if [ "${READS_TYPE}" = "ont" ]; then
    MINIMAP2_PRESET="map-ont"
else
    MINIMAP2_PRESET="map-hifi"
fi

MINIMAP2_OPTS="-x ${MINIMAP2_PRESET}"
if [ "${ONLY_PRIMARY}" = "true" ]; then
    MINIMAP2_OPTS="${MINIMAP2_OPTS} --secondary=no"
fi

# --- Step 1: 运行 minimap2 ---
if [ ! -s "${MINIMAP2_PAF}" ]; then
    echo "[Step 1] 运行 minimap2..."
    minimap2 -t ${THREADS} ${MINIMAP2_OPTS} "${REF_FASTA}" "${READS_FASTQ}" > "${MINIMAP2_PAF}"
else
    echo "[Step 1] minimap2 输出已存在且不为空，跳过。"
fi

# --- Step 2: Jellyfish 提取频率为 1 的参考 k-mer ---
# 注意：对于 ont 数据和 hifi 数据一样，都取参考基因组频率为 1 的 kmer 作为 unikmer
if [ ! -s "${KMER_FILE}" ]; then
    echo "[Step 2] 运行 Jellyfish 提取频率为 1 的 k-mer (k=${OV_K})..."
    jellyfish count -m ${OV_K} -s 1G -t ${THREADS} -C "${REF_FASTA}" -o "${JF_COUNT}"
    jellyfish dump -c -U 1 "${JF_COUNT}" > "${KMER_FILE}"
else
    echo "[Step 2] Jellyfish 结果已存在且不为空，跳过。"
fi

# 准备 k-mer 文件格式供 OvFilter 使用（只保留序列）
awk '{print $1}' "${KMER_FILE}" > "${OUT_DIR}/ovfilter.kmer"

# 准备 FASTA（将 READS_FASTQ 转为 FASTA，并与 REF_FASTA 合并）
if [ ! -s "${OUT_DIR}/combined.fa" ]; then
    echo "[Step 3] 准备组合 FASTA 文件..."
    sed -n '1~4s/^@/>/p;2~4p' "${READS_FASTQ}" > "${OUT_DIR}/query.fa"
    awk '/^>/{print n $1; n = "\n"} !/^>/{printf "%s",$0}' "${REF_FASTA}" > "${OUT_DIR}/reference.fa"
    sed -i '/^$/d' "${OUT_DIR}/reference.fa"
    cat "${OUT_DIR}/query.fa" "${OUT_DIR}/reference.fa" > "${OUT_DIR}/combined.fa"
else
    echo "[Step 3] 组合 FASTA 文件已存在且不为空，跳过。"
fi

# --- Step 4: 运行 OvFilter ---
if [ ! -s "${OVFILTER_PAF}" ]; then
    echo "[Step 4] 运行 OvFilter..."
    /usr/bin/time -v "${OVFILTER_CPP_BIN}" \
        --read-type "${READS_TYPE}" \
        --paf "${MINIMAP2_PAF}" \
        --unikmers "${OUT_DIR}/ovfilter.kmer" \
        --reads "${OUT_DIR}/combined.fa" \
        --output "${OVFILTER_PAF}" \
        --threads ${THREADS} \
        -k ${OV_K} 2> "${OVFILTER_TIME_LOG}"
else
    echo "[Step 4] OvFilter 输出已存在且不为空，跳过。"
fi

# --- Step 5: 评估 ---
echo "[Step 5] 运行评估..."
python "${EVAL_PY}" \
    --maf "${TRUTH_MAF}" \
    --raw-paf "${MINIMAP2_PAF}" \
    --ovfilter-paf "${OVFILTER_PAF}" \
    --ovfilter-label "OvFilter" \
    --ovfilter-time-log "${OVFILTER_TIME_LOG}" \
    --read-type "${READS_TYPE}" \
    --output-report "${OUT_DIR}/evaluation_report.txt"

echo "========================================================="
echo "  测试完毕！"
echo "  评估报告已保存到: ${OUT_DIR}/evaluation_report.txt"
cat "${OUT_DIR}/evaluation_report.txt"
echo "========================================================="
