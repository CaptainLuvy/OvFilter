#!/bin/bash
set -eu -o pipefail

# ==============================================================================
# 脚本名称: run_asm_wgs.sh
# 功能描述:
#   全基因组长读装配流程，默认复用独立预处理脚本 prep_reads_wgs.sh 的输出。
#   你可以：
#     1. 先单独运行 prep_reads_wgs.sh 检查处理后的 reads
#     2. 再运行本脚本继续 unikmer -> ava overlap -> OvFilter -> Miniasm
#   也可以直接运行本脚本，让它自动调用 prep_reads_wgs.sh。
#
# 用法:
#   bash run_asm_wgs.sh <work_dir> <reads.fastq|reads.fastq.gz> <reference.fasta> <hifi|ont>
#
# 可覆盖的环境变量:
#   THREADS=32
#   KMER_LEN=21
#   MIN_READ_LEN=0
#   TARGET_BASES_GB=""      # 留空时按模式使用默认值
#   SKIP_PREP=0             # 设为 1 时跳过 prep，直接使用 work_dir/raw/reads.used.fastq
#   READS_USED_OVERRIDE=""  # 显式指定已处理好的 FASTQ
#   MINIMAP2_EXTRA_OPTS=""
#   MIN_SHARED=""
#   MIN_IDENTITY=""
#   TOLERANCE=15
#   MAX_SD=10.0
#   MIN_ALIGN_LEN=500
#   OVERHANG_TOL=5000
#
# 依赖工具:
#   - ovfilter_cpp
#   - prep_reads_wgs.sh
#   - minimap2
#   - miniasm
#   - jellyfish
#   - python3
#   - awk / grep / du
# ==============================================================================

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <work_dir> <reads.fastq(.gz)> <reference.fasta> <hifi|ont>"
    exit 1
fi

WORK_DIR="$1"
READS_IN="$2"
REF_FASTA="$3"
READ_TYPE="$4"

if [ "${READ_TYPE}" != "hifi" ] && [ "${READ_TYPE}" != "ont" ]; then
    echo "Error: read_type must be 'hifi' or 'ont'."
    exit 1
fi

if [ ! -f "${READS_IN}" ]; then
    echo "Error: reads file not found: ${READS_IN}"
    exit 1
fi

if [ ! -f "${REF_FASTA}" ]; then
    echo "Error: reference file not found: ${REF_FASTA}"
    exit 1
fi

mkdir -p "${WORK_DIR}"
WORK_DIR=$(cd "${WORK_DIR}" && pwd)
READS_IN=$(cd "$(dirname "${READS_IN}")" && pwd)/$(basename "${READS_IN}")
REF_FASTA=$(cd "$(dirname "${REF_FASTA}")" && pwd)/$(basename "${REF_FASTA}")

THREADS=${THREADS:-32}
KMER_LEN=${KMER_LEN:-21}
MIN_READ_LEN=${MIN_READ_LEN:-0}
TARGET_BASES_GB=${TARGET_BASES_GB:-}
SKIP_PREP=${SKIP_PREP:-0}
READS_USED_OVERRIDE=${READS_USED_OVERRIDE:-}
MINIMAP2_EXTRA_OPTS=${MINIMAP2_EXTRA_OPTS:-}
MIN_SHARED=${MIN_SHARED:-}
MIN_IDENTITY=${MIN_IDENTITY:-}
TOLERANCE=${TOLERANCE:-15}
MAX_SD=${MAX_SD:-10.0}
MIN_ALIGN_LEN=${MIN_ALIGN_LEN:-500}
OVERHANG_TOL=${OVERHANG_TOL:-5000}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OVFILTER_CPP_BIN="${OVFILTER_CPP_BIN:-/home/dingyc/tools/ovfilter/ovfilter_cpp}"
PREP_SCRIPT="${PREP_SCRIPT:-$SCRIPT_DIR/prep_reads_wgs.sh}"
MINIMAP2="${MINIMAP2:-minimap2}"
MINIASM="${MINIASM:-miniasm}"
JELLYFISH="${JELLYFISH:-jellyfish}"
PYTHON="${PYTHON:-python3}"

if [ -z "${TARGET_BASES_GB}" ]; then
    if [ "${READ_TYPE}" = "hifi" ]; then
        TARGET_BASES_GB="6"
    else
        TARGET_BASES_GB="12"
    fi
fi

if [ "${READ_TYPE}" = "ont" ] && [ "${MIN_READ_LEN}" = "0" ]; then
    MIN_READ_LEN=10000
fi

if [ "${READ_TYPE}" = "hifi" ]; then
    MM2_PRESET="ava-pb"
else
    MM2_PRESET="ava-ont"
fi

RAW_DIR="${WORK_DIR}/raw"
UNIKMER_DIR="${WORK_DIR}/unikmer"
ASM_DIR="${WORK_DIR}/assembly"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${RAW_DIR}" "${UNIKMER_DIR}" "${ASM_DIR}" "${LOG_DIR}"

READS_PREP="${RAW_DIR}/reads.prepared.fastq"
READS_USED="${RAW_DIR}/reads.used.fastq"
PREP_SUMMARY="${RAW_DIR}/prep_summary.txt"
PAF_RAW="${WORK_DIR}/overlaps.raw.paf"
PAF_FILTERED="${WORK_DIR}/overlaps.ovfilter.paf"
PAF_NUMERIC="${WORK_DIR}/overlaps.numeric.paf"
UNIKMER_JF="${UNIKMER_DIR}/k${KMER_LEN}.jf"
UNIKMER_HISTO="${UNIKMER_DIR}/k${KMER_LEN}.histo"
UNIKMER_FILE="${UNIKMER_DIR}/k${KMER_LEN}_unikmers.kmers"
ASM_GFA="${ASM_DIR}/assembly.ovfilter.gfa"
ASM_FASTA="${ASM_DIR}/assembly.ovfilter.fasta"
ASM_NUMERIC_GFA="${ASM_DIR}/assembly.numeric.gfa"
ASM_NUMERIC_FASTA="${ASM_DIR}/assembly.numeric.fasta"
ASM_RAW_GFA="${ASM_DIR}/assembly.raw.gfa"
ASM_RAW_FASTA="${ASM_DIR}/assembly.raw.fasta"
LOG_FILE="${LOG_DIR}/run_asm_wgs.log"

echo "=========================================================="
echo "Whole-genome assembly workflow"
echo "Work dir:        ${WORK_DIR}"
echo "Reads:           ${READS_IN}"
echo "Reference:       ${REF_FASTA}"
echo "Read type:       ${READ_TYPE}"
echo "Threads:         ${THREADS}"
echo "K-mer length:    ${KMER_LEN}"
echo "Min read len:    ${MIN_READ_LEN}"
echo "Target Gbases:   ${TARGET_BASES_GB}"
echo "Skip prep:       ${SKIP_PREP}"
echo "Minimap2 preset: ${MM2_PRESET}"
echo "=========================================================="

if [ ! -x "${OVFILTER_CPP_BIN}" ]; then
    echo "Error: ovfilter_cpp not found or not executable: ${OVFILTER_CPP_BIN}"
    exit 1
fi

if [ ! -x "${PREP_SCRIPT}" ] && [ "${SKIP_PREP}" != "1" ] && [ -z "${READS_USED_OVERRIDE}" ]; then
    echo "Error: prep_reads_wgs.sh not found or not executable: ${PREP_SCRIPT}"
    exit 1
fi

if ! command -v "${MINIMAP2}" >/dev/null 2>&1; then
    echo "Error: minimap2 not found in PATH and MINIMAP2 not set."
    exit 1
fi

if ! command -v "${MINIASM}" >/dev/null 2>&1; then
    echo "Error: miniasm not found in PATH and MINIASM not set."
    exit 1
fi

if ! command -v "${JELLYFISH}" >/dev/null 2>&1; then
    echo "Error: jellyfish not found in PATH and JELLYFISH not set."
    exit 1
fi

if ! command -v "${PYTHON}" >/dev/null 2>&1; then
    echo "Error: python3 not found in PATH and PYTHON not set."
    exit 1
fi

# --- Step 0: 预处理 reads ---
if [ -n "${READS_USED_OVERRIDE}" ]; then
    READS_USED="${READS_USED_OVERRIDE}"
    echo "[0/5] Use external processed reads: ${READS_USED}"
elif [ "${SKIP_PREP}" = "1" ]; then
    if [ ! -s "${READS_USED}" ]; then
        echo "Error: SKIP_PREP=1 but ${READS_USED} not found."
        exit 1
    fi
    echo "[0/5] Skip preprocessing and reuse ${READS_USED}"
else
    echo "[0/5] Run standalone preprocessing script..."
    THREADS="${THREADS}" \
    MIN_READ_LEN="${MIN_READ_LEN}" \
    TARGET_BASES_GB="${TARGET_BASES_GB}" \
    PYTHON="${PYTHON}" \
    "${PREP_SCRIPT}" "${WORK_DIR}" "${READS_IN}" "${READ_TYPE}"
fi

if [ ! -s "${READS_USED}" ]; then
    echo "Error: processed reads not found: ${READS_USED}"
    exit 1
fi

# --- Step 1: 生成 unikmer ---
if [ ! -s "${UNIKMER_FILE}" ]; then
    echo "[1/5] Run Jellyfish and build unikmers from processed reads..."
    "${JELLYFISH}" count -m "${KMER_LEN}" -C -s 500M -t "${THREADS}" -o "${UNIKMER_JF}" "${READS_USED}" >> "${LOG_FILE}" 2>&1
    "${JELLYFISH}" histo -o "${UNIKMER_HISTO}" -v "${UNIKMER_JF}" >> "${LOG_FILE}" 2>&1

    read LOWER UPPER < <("${PYTHON}" - "${UNIKMER_HISTO}" <<'PY'
import math
import sys

histo = sys.argv[1]
freqs = []
counts = []
with open(histo, "r", encoding="utf-8") as f:
    for line in f:
        a, b = line.strip().split()
        if a == "0":
            continue
        freqs.append(int(a))
        counts.append(int(b))

start = 0
while start + 1 < len(counts) and counts[start] / max(counts[start + 1], 1) > 1.25:
    start += 1

trunc = counts[start:]
max_idx = counts.index(max(trunc))
end = min(start + 2 * (max_idx - start), len(counts) - 1)

F = counts[start:end+1]
X = freqs[start:end+1]
total = sum(F)
mean = sum(x * f for x, f in zip(X, F)) / total
var = sum(f * ((x - mean) ** 2) for x, f in zip(X, F)) / total
sd = math.sqrt(var)

lower = max(int(math.floor(mean - 3 * sd)), 5)
upper = int(math.ceil(mean + 3 * sd))
print(lower, upper)
PY
)

    echo "[1/5] unikmer bounds: lower=${LOWER}, upper=${UPPER}"
    "${JELLYFISH}" dump -c -L "${LOWER}" -U "${UPPER}" "${UNIKMER_JF}" | awk '{print $1}' > "${UNIKMER_FILE}"
else
    echo "[1/5] Reuse existing unikmer file."
fi

# --- Step 2: all-vs-all overlaps ---
if [ ! -s "${PAF_RAW}" ]; then
    echo "[2/5] Run minimap2 ${MM2_PRESET}..."
    "${MINIMAP2}" -t "${THREADS}" ${MINIMAP2_EXTRA_OPTS} -x "${MM2_PRESET}" "${READS_USED}" "${READS_USED}" > "${PAF_RAW}" 2>> "${LOG_FILE}"
else
    echo "[2/5] Reuse existing raw PAF."
fi

# --- Step 3: OvFilter ---
if [ ! -s "${PAF_FILTERED}" ] || [ ! -s "${PAF_NUMERIC}" ]; then
    echo "[3/5] Run OvFilter..."
    OVFILTER_ARGS=(
        --read-type "${READ_TYPE}"
        --paf "${PAF_RAW}"
        --unikmers "${UNIKMER_FILE}"
        --reads "${READS_USED}"
        --output "${PAF_FILTERED}"
        --output-numeric "${PAF_NUMERIC}"
        --threads "${THREADS}"
        -k "${KMER_LEN}"
        -t "${TOLERANCE}"
        --max-sd "${MAX_SD}"
        -l "${MIN_ALIGN_LEN}"
        -e "${OVERHANG_TOL}"
        --keep-contained
        --ignore-overhang
    )

    if [ -n "${MIN_SHARED}" ]; then
        OVFILTER_ARGS+=( -s "${MIN_SHARED}" )
    fi
    if [ -n "${MIN_IDENTITY}" ]; then
        OVFILTER_ARGS+=( -i "${MIN_IDENTITY}" )
    fi

    "${OVFILTER_CPP_BIN}" "${OVFILTER_ARGS[@]}" >> "${LOG_FILE}" 2>&1
else
    echo "[3/5] Reuse existing filtered and numeric PAFs."
fi

# --- Step 4: Miniasm ---
if [ ! -s "${ASM_FASTA}" ]; then
    echo "[4/5] Run miniasm on filtered PAF..."
    "${MINIASM}" -f "${READS_USED}" "${PAF_FILTERED}" > "${ASM_GFA}" 2>> "${LOG_FILE}"
    awk '/^S/{print ">"$2"\n"$3}' "${ASM_GFA}" > "${ASM_FASTA}"
else
    echo "[4/5] Reuse existing OvFilter assembly."
fi

if [ ! -s "${ASM_NUMERIC_FASTA}" ]; then
    echo "[4b/5] Run miniasm on numeric PAF..."
    "${MINIASM}" -f "${READS_USED}" "${PAF_NUMERIC}" > "${ASM_NUMERIC_GFA}" 2>> "${LOG_FILE}"
    awk '/^S/{print ">"$2"\n"$3}' "${ASM_NUMERIC_GFA}" > "${ASM_NUMERIC_FASTA}"
else
    echo "[4b/5] Reuse existing numeric assembly."
fi

if [ ! -s "${ASM_RAW_FASTA}" ]; then
    echo "[4c/5] Run miniasm on raw PAF..."
    "${MINIASM}" -f "${READS_USED}" "${PAF_RAW}" > "${ASM_RAW_GFA}" 2>> "${LOG_FILE}"
    awk '/^S/{print ">"$2"\n"$3}' "${ASM_RAW_GFA}" > "${ASM_RAW_FASTA}"
else
    echo "[4c/5] Reuse existing raw assembly."
fi

# --- Step 5: 汇总 ---
echo "[5/5] Summarize outputs..."
calc_fasta_stats() {
    local fasta="$1"
    "${PYTHON}" - "$fasta" <<'PY'
import sys

path = sys.argv[1]
lengths = []
cur = 0
with open(path, "r", encoding="utf-8") as fin:
    for line in fin:
        if line.startswith(">"):
            if cur > 0:
                lengths.append(cur)
            cur = 0
        else:
            cur += len(line.strip())
if cur > 0:
    lengths.append(cur)

if not lengths:
    print("0\t0\t0\t0")
    sys.exit(0)

lengths.sort(reverse=True)
contigs = len(lengths)
total = sum(lengths)
longest = lengths[0]
half = total / 2
acc = 0
n50 = 0
for x in lengths:
    acc += x
    if acc >= half:
        n50 = x
        break

print(f"{contigs}\t{total}\t{longest}\t{n50}")
PY
}
read FILTERED_CONTIGS FILTERED_TOTAL FILTERED_LONGEST FILTERED_N50 < <(calc_fasta_stats "${ASM_FASTA}")
read NUMERIC_CONTIGS NUMERIC_TOTAL NUMERIC_LONGEST NUMERIC_N50 < <(calc_fasta_stats "${ASM_NUMERIC_FASTA}")
read RAW_CONTIGS RAW_TOTAL RAW_LONGEST RAW_N50 < <(calc_fasta_stats "${ASM_RAW_FASTA}")
READS_SIZE=$(du -h "${READS_USED}" | cut -f1)
PAF_SIZE=$(du -h "${PAF_RAW}" | cut -f1)
PAF_FILTERED_SIZE=$(du -h "${PAF_FILTERED}" | cut -f1)
PAF_NUMERIC_SIZE=$(du -h "${PAF_NUMERIC}" | cut -f1)

echo "=========================================================="
echo "Done."
echo "Processed reads:    ${READS_USED} (${READS_SIZE})"
if [ -f "${PREP_SUMMARY}" ]; then
    echo "Prep summary:       ${PREP_SUMMARY}"
fi
echo "Reference:          ${REF_FASTA}"
echo "Unikmer file:       ${UNIKMER_FILE}"
echo "Raw overlaps:       ${PAF_RAW} (${PAF_SIZE})"
echo "Filtered overlaps:  ${PAF_FILTERED} (${PAF_FILTERED_SIZE})"
echo "Numeric overlaps:   ${PAF_NUMERIC} (${PAF_NUMERIC_SIZE})"
echo "OvFilter Assembly:  ${ASM_FASTA}"
echo "  contigs=${FILTERED_CONTIGS}, total=${FILTERED_TOTAL}, longest=${FILTERED_LONGEST}, N50=${FILTERED_N50}"
echo "Numeric Assembly:   ${ASM_NUMERIC_FASTA}"
echo "  contigs=${NUMERIC_CONTIGS}, total=${NUMERIC_TOTAL}, longest=${NUMERIC_LONGEST}, N50=${NUMERIC_N50}"
echo "Raw assembly:       ${ASM_RAW_FASTA}"
echo "  contigs=${RAW_CONTIGS}, total=${RAW_TOTAL}, longest=${RAW_LONGEST}, N50=${RAW_N50}"
echo "Contig comparison:  Raw=${RAW_CONTIGS}, Numeric=${NUMERIC_CONTIGS}, OvFilter=${FILTERED_CONTIGS}"
echo "Log file:           ${LOG_FILE}"
echo "=========================================================="
