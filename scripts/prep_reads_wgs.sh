#!/bin/bash
set -eu -o pipefail

# ==============================================================================
# 脚本名称: prep_reads_wgs.sh
# 功能描述:
#   为全基因组长读实验做独立预处理，便于单独检查预处理后的数据。
#
# 处理内容:
#   1. 解压 FASTQ（若输入为 .gz）
#   2. 按最短读长过滤
#   3. 按目标总碱基数下采样
#
# 说明:
#   - 当前默认实现是轻量级预处理，不依赖额外专门工具。
#   - 若环境中提供 rasusa / seqtk，可切换为对应下采样工具。
#   - 若后续需要更激进的 ONT 质量筛选，可再接入 filtlong。
#
# 用法:
#   bash prep_reads_wgs.sh <work_dir> <reads.fastq|reads.fastq.gz> <hifi|ont>
# ==============================================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <work_dir> <reads.fastq(.gz)> <hifi|ont>"
    exit 1
fi

WORK_DIR="$1"
READS_IN="$2"
READ_TYPE="$3"

if [ "${READ_TYPE}" != "hifi" ] && [ "${READ_TYPE}" != "ont" ]; then
    echo "Error: read_type must be 'hifi' or 'ont'."
    exit 1
fi

if [ ! -f "${READS_IN}" ]; then
    echo "Error: reads file not found: ${READS_IN}"
    exit 1
fi

mkdir -p "${WORK_DIR}"
WORK_DIR=$(cd "${WORK_DIR}" && pwd)
READS_IN=$(cd "$(dirname "${READS_IN}")" && pwd)/$(basename "${READS_IN}")

THREADS=${THREADS:-32}
MIN_READ_LEN=${MIN_READ_LEN:-0}
TARGET_BASES_GB=${TARGET_BASES_GB:-}
DOWNSAMPLE_TOOL=${DOWNSAMPLE_TOOL:-builtin}
PYTHON=${PYTHON:-python3}
RASUSA=${RASUSA:-rasusa}
SEQTK=${SEQTK:-seqtk}
FORCE_PREP_SUMMARY=${FORCE_PREP_SUMMARY:-0}

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

RAW_DIR="${WORK_DIR}/raw"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${RAW_DIR}" "${LOG_DIR}"

READS_PREP="${RAW_DIR}/reads.prepared.fastq"
READS_USED="${RAW_DIR}/reads.used.fastq"
SUMMARY_FILE="${RAW_DIR}/prep_summary.txt"
LOG_FILE="${LOG_DIR}/prep_reads_wgs.log"

echo "==========================================================" | tee -a "${LOG_FILE}"
echo "Preprocessing whole-genome reads" | tee -a "${LOG_FILE}"
echo "Work dir:         ${WORK_DIR}" | tee -a "${LOG_FILE}"
echo "Reads:            ${READS_IN}" | tee -a "${LOG_FILE}"
echo "Read type:        ${READ_TYPE}" | tee -a "${LOG_FILE}"
echo "Threads:          ${THREADS}" | tee -a "${LOG_FILE}"
echo "Min read len:     ${MIN_READ_LEN}" | tee -a "${LOG_FILE}"
echo "Target Gbases:    ${TARGET_BASES_GB}" | tee -a "${LOG_FILE}"
echo "Downsample tool:  ${DOWNSAMPLE_TOOL}" | tee -a "${LOG_FILE}"
echo "==========================================================" | tee -a "${LOG_FILE}"

if [ ! -s "${READS_PREP}" ]; then
    echo "[1/3] 解压并按长度过滤 reads..." | tee -a "${LOG_FILE}"
    "${PYTHON}" - "${READS_IN}" "${READS_PREP}" "${MIN_READ_LEN}" <<'PY' | tee -a "${LOG_FILE}"
import gzip
import sys

in_path, out_path, min_len = sys.argv[1], sys.argv[2], int(sys.argv[3])

def opener(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

kept = 0
total = 0
kept_bases = 0
total_bases = 0
with opener(in_path, "rt") as fin, open(out_path, "w", encoding="utf-8") as fout:
    while True:
        h = fin.readline()
        if not h:
            break
        s = fin.readline()
        p = fin.readline()
        q = fin.readline()
        if not s or not p or not q:
            raise RuntimeError("FASTQ truncated")
        seq = s.rstrip("\r\n")
        total += 1
        total_bases += len(seq)
        if len(seq) < min_len:
            continue
        fout.write(h)
        fout.write(s)
        fout.write(p)
        fout.write(q)
        kept += 1
        kept_bases += len(seq)

print(f"[prep] total_reads={total} total_bases={total_bases} kept_reads={kept} kept_bases={kept_bases} min_len={min_len}")
PY
else
    echo "[1/3] 已存在 prepared reads，跳过。" | tee -a "${LOG_FILE}"
fi

if [ ! -s "${READS_USED}" ]; then
    echo "[2/3] 下采样到目标总碱基数..." | tee -a "${LOG_FILE}"
    TARGET_BASES=$(python3 - "${TARGET_BASES_GB}" <<'PY'
import sys
print(int(float(sys.argv[1]) * 1_000_000_000))
PY
)

    if [ "${DOWNSAMPLE_TOOL}" = "rasusa" ] && command -v "${RASUSA}" >/dev/null 2>&1; then
        "${RASUSA}" reads -i "${READS_PREP}" -o "${READS_USED}" -g "${TARGET_BASES}" -c 1 >/dev/null 2>>"${LOG_FILE}"
        echo "[downsample] used rasusa" | tee -a "${LOG_FILE}"
    elif [ "${DOWNSAMPLE_TOOL}" = "seqtk" ] && command -v "${SEQTK}" >/dev/null 2>&1; then
        TOTAL_BASES=$("${PYTHON}" - "${READS_PREP}" <<'PY'
import sys
total = 0
with open(sys.argv[1], "r", encoding="utf-8") as fin:
    while True:
        h = fin.readline()
        if not h:
            break
        s = fin.readline()
        p = fin.readline()
        q = fin.readline()
        total += len(s.rstrip("\r\n"))
print(total)
PY
)
        if [ "${TOTAL_BASES}" -le "${TARGET_BASES}" ]; then
            cp "${READS_PREP}" "${READS_USED}"
            echo "[downsample] total_bases<=target_bases, copied all reads" | tee -a "${LOG_FILE}"
        else
            FRACTION=$("${PYTHON}" - "${TOTAL_BASES}" "${TARGET_BASES}" <<'PY'
import sys
total = int(sys.argv[1])
target = int(sys.argv[2])
print(target / total)
PY
)
            "${SEQTK}" sample -s42 "${READS_PREP}" "${FRACTION}" > "${READS_USED}"
            echo "[downsample] used seqtk fraction=${FRACTION}" | tee -a "${LOG_FILE}"
        fi
    else
        "${PYTHON}" - "${READS_PREP}" "${READS_USED}" "${TARGET_BASES_GB}" <<'PY' | tee -a "${LOG_FILE}"
import random
import sys

in_path, out_path, target_gb = sys.argv[1], sys.argv[2], float(sys.argv[3])
target_bases = int(target_gb * 1_000_000_000)

total_bases = 0
with open(in_path, "r", encoding="utf-8") as fin:
    while True:
        h = fin.readline()
        if not h:
            break
        s = fin.readline()
        p = fin.readline()
        q = fin.readline()
        if not s or not p or not q:
            raise RuntimeError("FASTQ truncated while counting")
        total_bases += len(s.rstrip("\r\n"))

if total_bases <= target_bases:
    import shutil
    shutil.copyfile(in_path, out_path)
    print(f"[downsample] total_bases={total_bases} <= target_bases={target_bases}, copied all reads")
    sys.exit(0)

prob = target_bases / total_bases
random.seed(42)
kept_reads = 0
kept_bases = 0
with open(in_path, "r", encoding="utf-8") as fin, open(out_path, "w", encoding="utf-8") as fout:
    while True:
        h = fin.readline()
        if not h:
            break
        s = fin.readline()
        p = fin.readline()
        q = fin.readline()
        if not s or not p or not q:
            raise RuntimeError("FASTQ truncated while sampling")
        if random.random() <= prob:
            fout.write(h)
            fout.write(s)
            fout.write(p)
            fout.write(q)
            kept_reads += 1
            kept_bases += len(s.rstrip("\r\n"))

print(f"[downsample] total_bases={total_bases} target_bases={target_bases} prob={prob:.6f} kept_reads={kept_reads} kept_bases={kept_bases}")
PY
        echo "[downsample] used builtin sampler" | tee -a "${LOG_FILE}"
    fi
else
    echo "[2/3] 已存在 sampled reads，跳过。" | tee -a "${LOG_FILE}"
fi

echo "[3/3] 汇总预处理结果..." | tee -a "${LOG_FILE}"
if [ "${FORCE_PREP_SUMMARY}" = "1" ] || [ ! -s "${SUMMARY_FILE}" ]; then
"${PYTHON}" - "${READS_PREP}" "${READS_USED}" "${SUMMARY_FILE}" <<'PY' | tee -a "${LOG_FILE}"
import sys

prep_path, used_path, summary_path = sys.argv[1], sys.argv[2], sys.argv[3]

def summarize(path):
    reads = 0
    bases = 0
    min_len = None
    max_len = 0
    with open(path, "r", encoding="utf-8") as fin:
        while True:
            h = fin.readline()
            if not h:
                break
            s = fin.readline()
            p = fin.readline()
            q = fin.readline()
            if not s or not p or not q:
                raise RuntimeError(f"FASTQ truncated: {path}")
            l = len(s.rstrip("\r\n"))
            reads += 1
            bases += l
            min_len = l if min_len is None else min(min_len, l)
            max_len = max(max_len, l)
    mean_len = bases / reads if reads else 0.0
    return reads, bases, min_len or 0, max_len, mean_len

prep = summarize(prep_path)
used = summarize(used_path)

with open(summary_path, "w", encoding="utf-8") as fout:
    fout.write("dataset\treads\tbases\tmin_len\tmax_len\tmean_len\n")
    fout.write(f"prepared\t{prep[0]}\t{prep[1]}\t{prep[2]}\t{prep[3]}\t{prep[4]:.2f}\n")
    fout.write(f"used\t{used[0]}\t{used[1]}\t{used[2]}\t{used[3]}\t{used[4]:.2f}\n")

print(f"[summary] prepared_reads={prep[0]} prepared_bases={prep[1]} prepared_mean_len={prep[4]:.2f}")
print(f"[summary] used_reads={used[0]} used_bases={used[1]} used_mean_len={used[4]:.2f}")
print(f"[summary] summary_file={summary_path}")
PY
else
echo "[3/3] Reuse existing prep summary." | tee -a "${LOG_FILE}"
fi

echo "Done."
echo "Prepared reads: ${READS_PREP}"
echo "Used reads:     ${READS_USED}"
echo "Summary:        ${SUMMARY_FILE}"
