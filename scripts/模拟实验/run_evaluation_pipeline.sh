#!/bin/bash
set -e
set -o pipefail

# 强制使用 UTF-8 终端与 Python 输出编码，避免中文乱码。
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONIOENCODING=UTF-8

# ==============================================================================
# 完整评估流程：
#   Minimap2 -> Jellyfish -> RAfilter -> OvFilter -> Evaluation
# 使用 run 目录下的 .done 文件支持断点续跑。
# ==============================================================================

# --- 1. 核心输入文件：按当前实验修改这些路径 ---
READS_FASTQ="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp/sim_data/chr8/ont/chr8_ont_0001_complete.fastq"
REF_FASTA="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp/censat_ref_data/chr8_cen_wide.fasta"
TRUTH_MAF="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp/sim_data/chr8/ont/chr8_ont_0001.maf"

# --- 2. 输出目录与程序路径 ---
OUT_DIR="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp/result/chr8/ont30x_secondary"
RAFILTER_BIN="/home/dingyc/tools/RAfilter/src/rafilter"
OVFILTER_CPP_BIN="/home/dingyc/tools/ovfilter/ovfilter_cpp"
EVAL_PY="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp/result/evaluate_alignments.py"

# --- 3. 通用运行参数 ---
THREADS=32
READS_TYPE="ont"     # hifi 或 ont
ONLY_PRIMARY="false" # true 表示 minimap2 使用 --secondary=no

# --- 4. OvFilter 参数 ---
OV_K=21
# 若留空，则交给 ovfilter.cpp 根据模式自动取默认值：
#   hifi -> min_shared = 15
#   ont  -> min_shared = 6
OV_MIN_SHARED=""
# 若留空，则交给 ovfilter.cpp 根据模式自动取默认值：
#   hifi -> min_identity = 0
#   ont  -> 从原始 PAF 计算简单均值并四舍五入到两位小数
OV_MIN_IDENTITY=""
OV_TOLERANCE=15
OV_MAX_SD=10.0
OV_MIN_ALIGN_LEN=0
OV_OVERHANG=5000
OV_EXTRA_ARGS="--keep-contained --ignore-overhang"

# 可选：额外跑一组 ONT 特调版本（改变参考 k-mer 频次上限）
OVFILTER_USE_ONT_TUNED_KMERS="true"
OVFILTER_TUNED_MAX_COUNT=1
OVFILTER_TUNED_LABEL="OvFilter（ONT特调）"

# 固定输出目录结构，便于断点续跑。
RUN_OUT_DIR="${OUT_DIR}/run"
RAFILTER_OUT_DIR="${RUN_OUT_DIR}/rafilter"
OVFILTER_OUT_DIR="${RUN_OUT_DIR}/ovfilter"
MINIMAP2_OUT_DIR="${RUN_OUT_DIR}/minimap2"
TEMP_DIR="${RUN_OUT_DIR}/tmp"
LOG_DIR="${RUN_OUT_DIR}/logs"

mkdir -p "${RAFILTER_OUT_DIR}" "${OVFILTER_OUT_DIR}" "${MINIMAP2_OUT_DIR}" "${TEMP_DIR}" "${LOG_DIR}"

MINIMAP2_PAF="${MINIMAP2_OUT_DIR}/minimap2_raw.paf"
QUERY_FA="${TEMP_DIR}/query.fa"
REF_FA="${TEMP_DIR}/reference.fa"
JF_COUNT="${TEMP_DIR}/mer_counts.jf"
KMER_FILE="${TEMP_DIR}/target.kmer"
RAFILTER_BUILD_TIME_LOG="${LOG_DIR}/rafilter_build_time.log"
TIME_LOG="${LOG_DIR}/rafilter_time.log"
MINIMAP2_LOG="${LOG_DIR}/minimap2.stderr.log"

OVFILTER_DEFAULT_LABEL="OvFilter"
OVFILTER_DEFAULT_TAG="${READS_TYPE}_default"
OVFILTER_DEFAULT_PAF="${OVFILTER_OUT_DIR}/ovfilter_${OVFILTER_DEFAULT_TAG}.paf"
OVFILTER_DEFAULT_TIME_LOG="${LOG_DIR}/ovfilter_time_${OVFILTER_DEFAULT_TAG}.log"
OVFILTER_DEFAULT_DONE_FILE="${LOG_DIR}/step6_ovfilter_${OVFILTER_DEFAULT_TAG}.done"
OVFILTER_DEFAULT_KMER_FILE="${TEMP_DIR}/ovfilter_target_${OVFILTER_DEFAULT_TAG}.kmer"

if [ "${READS_TYPE}" = "ont" ] && [ "${OVFILTER_USE_ONT_TUNED_KMERS}" = "true" ]; then
    OVFILTER_TUNED_TAG="ont_u${OVFILTER_TUNED_MAX_COUNT}"
    OVFILTER_TUNED_PAF="${OVFILTER_OUT_DIR}/ovfilter_${OVFILTER_TUNED_TAG}.paf"
    OVFILTER_TUNED_TIME_LOG="${LOG_DIR}/ovfilter_time_${OVFILTER_TUNED_TAG}.log"
    OVFILTER_TUNED_DONE_FILE="${LOG_DIR}/step6_ovfilter_${OVFILTER_TUNED_TAG}.done"
    OVFILTER_TUNED_KMER_FILE="${TEMP_DIR}/ovfilter_target_${OVFILTER_TUNED_TAG}.kmer"
else
    OVFILTER_TUNED_TAG=""
    OVFILTER_TUNED_PAF=""
    OVFILTER_TUNED_TIME_LOG=""
    OVFILTER_TUNED_DONE_FILE=""
    OVFILTER_TUNED_KMER_FILE=""
fi

if [ "${READS_TYPE}" = "ont" ]; then
    MINIMAP2_PRESET="map-ont"
else
    MINIMAP2_PRESET="map-hifi"
fi

MINIMAP2_OPTS="-x ${MINIMAP2_PRESET}"
if [ "${ONLY_PRIMARY}" = "true" ]; then
    MINIMAP2_OPTS="${MINIMAP2_OPTS} --secondary=no"
    echo "提示：ONLY_PRIMARY=true，minimap2 只输出 primary alignments。"
fi

echo "========================================================="
echo "  开始执行 RAfilter / OvFilter 完整评估流程"
echo "  输出目录：${RUN_OUT_DIR}"
echo "========================================================="

cd "${RUN_OUT_DIR}"

# --- Step 1: 运行 minimap2 生成原始 PAF ---
if [ ! -f "${LOG_DIR}/step1_minimap2.done" ]; then
    echo "[Step 1] 运行 minimap2..."
    if [ ! -f "${REF_FASTA}" ]; then
        echo "[Step 1] ERROR: REF_FASTA 不存在：${REF_FASTA}" >&2
        exit 1
    fi
    if [ ! -f "${READS_FASTQ}" ]; then
        echo "[Step 1] ERROR: READS_FASTQ 不存在：${READS_FASTQ}" >&2
        exit 1
    fi

    python - "${READS_FASTQ}" <<'PY'
import sys
if len(sys.argv) < 2:
    sys.exit(1)
path = sys.argv[1]
max_records = 2000
max_line_len = 2000000

def bad(msg: str) -> None:
    sys.stderr.write(f"[Step 1] ERROR: FASTQ 格式有问题：{msg}\n")
    sys.stderr.write("[Step 1] 建议：先把输入修成严格的 4 行 FASTQ。\n")
    sys.exit(1)

def read_line(f):
    line = f.readline()
    if not line:
        return None
    if len(line) > max_line_len:
        bad(f"单行长度超过 {max_line_len} 字节，可能缺少换行或文件损坏。")
    return line

with open(path, "rb") as f:
    for i in range(max_records):
        h = read_line(f)
        if h is None:
            if i == 0:
                bad("文件为空。")
            break
        s = read_line(f)
        p = read_line(f)
        q = read_line(f)
        if s is None or p is None or q is None:
            bad(f"在第 {i+1} 条记录处提前 EOF（总行数不是 4 的倍数）。")
        if not h.startswith(b"@"):
            bad(f"第 {i+1} 条记录头行不是 @ 开头。")
        if not p.startswith(b"+"):
            bad(f"第 {i+1} 条记录第三行不是 + 开头。")
        s2 = s.rstrip(b"\r\n")
        q2 = q.rstrip(b"\r\n")
        if len(s2) == 0:
            bad(f"第 {i+1} 条记录序列为空。")
        if len(s2) != len(q2):
            bad(f"第 {i+1} 条记录序列长度 {len(s2)} 与质量值长度 {len(q2)} 不一致。")
sys.stderr.write(f"[Step 1] FASTQ 预检查通过（抽查前 {max_records} 条记录）。\n")
PY

    echo "[Step 1] 输入文件大小："
    ls -lh "${REF_FASTA}" "${READS_FASTQ}"
    echo "[Step 1] minimap2 日志：${MINIMAP2_LOG}"

    : > "${MINIMAP2_PAF}"
    : > "${MINIMAP2_LOG}"

    minimap2 -t ${THREADS} ${MINIMAP2_OPTS} "${REF_FASTA}" "${READS_FASTQ}" \
        > "${MINIMAP2_PAF}" \
        2> >(tee -a "${MINIMAP2_LOG}" >&2) &
    MINIMAP2_PID=$!

    while kill -0 ${MINIMAP2_PID} 2>/dev/null; do
        if [ -f "${MINIMAP2_PAF}" ]; then
            PAF_BYTES=$(wc -c < "${MINIMAP2_PAF}" 2>/dev/null || echo 0)
            echo "[Step 1] minimap2 运行中... 当前 PAF 字节数=${PAF_BYTES}"
        else
            echo "[Step 1] minimap2 运行中... PAF 还未生成"
        fi
        sleep 20
    done

    wait ${MINIMAP2_PID}

    if [ ! -s "${MINIMAP2_PAF}" ]; then
        echo "[Step 1] ERROR: minimap2 输出为空：${MINIMAP2_PAF}" >&2
        echo "[Step 1] minimap2 stderr 末尾 80 行：" >&2
        tail -n 80 "${MINIMAP2_LOG}" >&2 || true
        exit 1
    fi

    touch "${LOG_DIR}/step1_minimap2.done"
    echo "[Step 1] 完成"
else
    echo "[Step 1] minimap2 已完成，跳过。"
fi

# --- Step 2: 转成 RAfilter 需要的单行 FASTA ---
if [ ! -f "${LOG_DIR}/step2_format_fasta.done" ]; then
    echo "[Step 2] 转换 reads/reference 为单行 FASTA..."
    echo "  转换 reads..."
    sed -n '1~4s/^@/>/p;2~4p' "${READS_FASTQ}" > "${QUERY_FA}"
    echo "  转换 reference..."
    awk '/^>/{print n $1; n = "\n"} !/^>/{printf "%s",$0}' "${REF_FASTA}" > "${REF_FA}"
    sed -i '/^$/d' "${REF_FA}"
    touch "${LOG_DIR}/step2_format_fasta.done"
    echo "[Step 2] 完成"
else
    echo "[Step 2] FASTA 转换已完成，跳过。"
fi

# --- Step 3: Jellyfish 统计参考 k-mer ---
if [ ! -f "${LOG_DIR}/step3_jellyfish.done" ]; then
    echo "[Step 3] 运行 Jellyfish（k=${OV_K}）..."
    jellyfish count -m ${OV_K} -s 1G -t ${THREADS} -C "${REF_FASTA}" -o "${JF_COUNT}"
    if [ "${READS_TYPE}" = "ont" ]; then
        echo "  提取 ONT 用参考 rare k-mer（频次 <= 3）..."
        jellyfish dump -c -U 3 "${JF_COUNT}" > "${KMER_FILE}"
    else
        echo "  提取 HiFi 用参考 unique k-mer（频次 == 1）..."
        jellyfish dump -c -U 1 "${JF_COUNT}" > "${KMER_FILE}"
    fi
    touch "${LOG_DIR}/step3_jellyfish.done"
    echo "[Step 3] 完成"
else
    echo "[Step 3] Jellyfish 已完成，跳过。"
fi

# --- Step 4: RAfilter build ---
if [ ! -f "${LOG_DIR}/step4_rafilter_build.done" ]; then
    echo "[Step 4] 运行 RAfilter build..."
    /usr/bin/time -v "${RAFILTER_BIN}" build -t ${THREADS} -q "${QUERY_FA}" -r "${REF_FA}" -o "${TEMP_DIR}/" "${KMER_FILE}" 2> "${RAFILTER_BUILD_TIME_LOG}"
    touch "${LOG_DIR}/step4_rafilter_build.done"
    echo "[Step 4] 完成"
else
    echo "[Step 4] RAfilter build 已完成，跳过。"
fi

# --- Step 5: RAfilter filter ---
if [ ! -f "${LOG_DIR}/step5_rafilter_filter.done" ]; then
    echo "[Step 5] 运行 RAfilter filter，并记录时间/内存..."
    /usr/bin/time -v "${RAFILTER_BIN}" filter -o "${RAFILTER_OUT_DIR}/" -p --threshold 12 "${TEMP_DIR}/ref.pos" "${TEMP_DIR}/query.pos" "${MINIMAP2_PAF}" 2> "${TIME_LOG}"
    touch "${LOG_DIR}/step5_rafilter_filter.done"
    echo "[Step 5] 完成"
else
    echo "[Step 5] RAfilter filter 已完成，跳过。"
fi

run_ovfilter_variant() {
    local label="$1"
    local paf_out="$2"
    local time_log="$3"
    local done_file="$4"
    local kmer_source="$5"
    local params_file="${done_file}.params"
    local min_shared_stamp="${OV_MIN_SHARED:-auto_by_mode}"
    local min_identity_stamp="${OV_MIN_IDENTITY:-auto_by_mode}"
    local params_stamp="k=${OV_K};read_type=${READS_TYPE};s=${min_shared_stamp};t=${OV_TOLERANCE};sd=${OV_MAX_SD};l=${OV_MIN_ALIGN_LEN};i=${min_identity_stamp};e=${OV_OVERHANG};kmer=${kmer_source};extra=${OV_EXTRA_ARGS}"
    local min_shared_arg=""
    local min_identity_arg=""

    if [ -n "${OV_MIN_SHARED}" ]; then
        min_shared_arg="-s ${OV_MIN_SHARED}"
    fi
    if [ -n "${OV_MIN_IDENTITY}" ]; then
        min_identity_arg="-i ${OV_MIN_IDENTITY}"
    fi

    if [ -f "${done_file}" ] && [ -f "${params_file}" ] && [ "$(cat "${params_file}")" = "${params_stamp}" ]; then
        echo "[Step 6] ${label} 已完成，跳过。"
        return
    fi

    if [ -f "${done_file}" ] || [ -f "${params_file}" ]; then
        echo "[Step 6] ${label} 检测到参数变化，重新运行。"
        rm -f "${done_file}" "${params_file}"
    fi

    echo "[Step 6] 运行 ${label}..."
    echo "[Step 6] read_type=${READS_TYPE}, min_shared=${min_shared_stamp}, min_identity=${min_identity_stamp}"
    awk '{print $1}' "${kmer_source}" > "${TEMP_DIR}/ovfilter.kmer"
    cat "${QUERY_FA}" "${REF_FA}" > "${TEMP_DIR}/combined_for_ovfilter.fa"

    /usr/bin/time -v "${OVFILTER_CPP_BIN}" \
        --read-type "${READS_TYPE}" \
        --paf "${MINIMAP2_PAF}" \
        --unikmers "${TEMP_DIR}/ovfilter.kmer" \
        --reads "${TEMP_DIR}/combined_for_ovfilter.fa" \
        --output "${paf_out}" \
        --threads ${THREADS} \
        -k ${OV_K} \
        ${min_shared_arg} \
        -t ${OV_TOLERANCE} \
        --max-sd ${OV_MAX_SD} \
        -l ${OV_MIN_ALIGN_LEN} \
        ${min_identity_arg} \
        -e ${OV_OVERHANG} \
        ${OV_EXTRA_ARGS} 2> "${time_log}"

    printf '%s\n' "${params_stamp}" > "${params_file}"
    touch "${done_file}"
    echo "[Step 6] ${label} 完成"
}

# --- Step 6: 默认 OvFilter ---
cp "${KMER_FILE}" "${OVFILTER_DEFAULT_KMER_FILE}"
run_ovfilter_variant "${OVFILTER_DEFAULT_LABEL}" "${OVFILTER_DEFAULT_PAF}" "${OVFILTER_DEFAULT_TIME_LOG}" "${OVFILTER_DEFAULT_DONE_FILE}" "${OVFILTER_DEFAULT_KMER_FILE}"

# --- Step 6b: 可选的 ONT 特调 OvFilter ---
if [ -n "${OVFILTER_TUNED_TAG}" ]; then
    if [ ! -f "${OVFILTER_TUNED_KMER_FILE}" ]; then
        echo "[Step 6b] 为 ${OVFILTER_TUNED_LABEL} 生成 ONT 特调参考 k-mer（jellyfish dump -c -U ${OVFILTER_TUNED_MAX_COUNT}）..."
        jellyfish dump -c -U ${OVFILTER_TUNED_MAX_COUNT} "${JF_COUNT}" > "${OVFILTER_TUNED_KMER_FILE}"
    fi
    run_ovfilter_variant "${OVFILTER_TUNED_LABEL}" "${OVFILTER_TUNED_PAF}" "${OVFILTER_TUNED_TIME_LOG}" "${OVFILTER_TUNED_DONE_FILE}" "${OVFILTER_TUNED_KMER_FILE}"
fi

# --- Step 7: 评估 ---
echo "[Step 7] 运行评估脚本..."
RAFILTER_FINAL_PAF="${RAFILTER_OUT_DIR}/rafiltered.paf"
if [ ! -f "${RAFILTER_FINAL_PAF}" ]; then
    RAFILTER_FINAL_PAF=$(ls "${RAFILTER_OUT_DIR}"/*rafiltered*.paf 2>/dev/null | head -n 1 || echo "")
fi

REPORT_FILE="${RUN_OUT_DIR}/evaluation_report.txt"
OV_MIN_SHARED_DESC="${OV_MIN_SHARED:-auto_by_mode}"
OV_MIN_IDENTITY_DESC="${OV_MIN_IDENTITY:-auto_by_mode}"
{
    echo "========================================================="
    echo "  运行记录（时间：$(date "+%Y-%m-%d %H:%M:%S")）"
    echo "========================================================="
    echo "1. 输入文件："
    echo "   - 模拟 reads: ${READS_FASTQ}"
    echo "   - 参考 fasta: ${REF_FASTA}"
    echo "   - 真实 MAF:   ${TRUTH_MAF}"
    echo ""
    echo "2. 程序路径："
    echo "   - RAfilter: ${RAFILTER_BIN}"
    echo "   - OvFilter(C++): ${OVFILTER_CPP_BIN}"
    echo ""
    echo "3. OvFilter 模式与覆盖参数："
    echo "   - read_type: ${READS_TYPE}"
    echo "   - min_shared override: ${OV_MIN_SHARED_DESC}"
    echo "   - min_identity override: ${OV_MIN_IDENTITY_DESC}"
    echo ""
    echo "4. 关键命令："
    echo "   - Minimap2:"
    echo "     minimap2 -t ${THREADS} ${MINIMAP2_OPTS} ${REF_FASTA} ${READS_FASTQ} > ${MINIMAP2_PAF}"
    echo "   - Jellyfish:"
    echo "     jellyfish count -m ${OV_K} -s 1G -t ${THREADS} -C ${REF_FASTA} -o ${JF_COUNT}"
    echo "     jellyfish dump -c ${JF_COUNT} > ${KMER_FILE}"
    echo "   - RAfilter build:"
    echo "     /usr/bin/time -v ${RAFILTER_BIN} build -t ${THREADS} -q ${QUERY_FA} -r ${REF_FA} -o ${TEMP_DIR}/ ${KMER_FILE}"
    echo "   - RAfilter filter:"
    echo "     ${RAFILTER_BIN} filter -o ${RAFILTER_OUT_DIR}/ -p --threshold 12 ${TEMP_DIR}/ref.pos ${TEMP_DIR}/query.pos ${MINIMAP2_PAF}"
    echo "   - ${OVFILTER_DEFAULT_LABEL}:"
    echo "     ${OVFILTER_CPP_BIN} --read-type ${READS_TYPE} --paf ${MINIMAP2_PAF} --unikmers ${OVFILTER_DEFAULT_KMER_FILE} --reads ${TEMP_DIR}/combined_for_ovfilter.fa --output ${OVFILTER_DEFAULT_PAF} --threads ${THREADS} -k ${OV_K} -t ${OV_TOLERANCE} --max-sd ${OV_MAX_SD} -l ${OV_MIN_ALIGN_LEN} -e ${OV_OVERHANG} ${OV_EXTRA_ARGS}"
    if [ -n "${OV_MIN_SHARED}" ]; then
        echo "     显式覆盖：-s ${OV_MIN_SHARED}"
    fi
    if [ -n "${OV_MIN_IDENTITY}" ]; then
        echo "     显式覆盖：-i ${OV_MIN_IDENTITY}"
    fi
    if [ -n "${OVFILTER_TUNED_TAG}" ]; then
        echo "   - ${OVFILTER_TUNED_LABEL}:"
        echo "     jellyfish dump -c -U ${OVFILTER_TUNED_MAX_COUNT} ${JF_COUNT} > ${OVFILTER_TUNED_KMER_FILE}"
        echo "     ${OVFILTER_CPP_BIN} --read-type ${READS_TYPE} --paf ${MINIMAP2_PAF} --unikmers ${OVFILTER_TUNED_KMER_FILE} --reads ${TEMP_DIR}/combined_for_ovfilter.fa --output ${OVFILTER_TUNED_PAF} --threads ${THREADS} -k ${OV_K} -t ${OV_TOLERANCE} --max-sd ${OV_MAX_SD} -l ${OV_MIN_ALIGN_LEN} -e ${OV_OVERHANG} ${OV_EXTRA_ARGS}"
        if [ -n "${OV_MIN_SHARED}" ]; then
            echo "     显式覆盖：-s ${OV_MIN_SHARED}"
        fi
        if [ -n "${OV_MIN_IDENTITY}" ]; then
            echo "     显式覆盖：-i ${OV_MIN_IDENTITY}"
        fi
    fi
    echo "========================================================="
    echo ""
} > "${REPORT_FILE}"

if [ -n "${OVFILTER_TUNED_TAG}" ]; then
    python "${EVAL_PY}" \
        --maf "${TRUTH_MAF}" \
        --raw-paf "${MINIMAP2_PAF}" \
        --rafilter-paf "${RAFILTER_FINAL_PAF}" \
        --ovfilter-paf "${OVFILTER_DEFAULT_PAF}" \
        --ovfilter-label "${OVFILTER_DEFAULT_LABEL}" \
        --ovfilter-time-log "${OVFILTER_DEFAULT_TIME_LOG}" \
        --ovfilter-extra-paf "${OVFILTER_TUNED_PAF}" \
        --ovfilter-extra-label "${OVFILTER_TUNED_LABEL}" \
        --ovfilter-extra-time-log "${OVFILTER_TUNED_TIME_LOG}" \
        --rafilter-build-time-log "${RAFILTER_BUILD_TIME_LOG}" \
        --time-log "${TIME_LOG}" \
        --read-type "${READS_TYPE}" \
        --output-report "${RUN_OUT_DIR}/temp_eval_result.txt"
else
    python "${EVAL_PY}" \
        --maf "${TRUTH_MAF}" \
        --raw-paf "${MINIMAP2_PAF}" \
        --rafilter-paf "${RAFILTER_FINAL_PAF}" \
        --ovfilter-paf "${OVFILTER_DEFAULT_PAF}" \
        --ovfilter-label "${OVFILTER_DEFAULT_LABEL}" \
        --ovfilter-time-log "${OVFILTER_DEFAULT_TIME_LOG}" \
        --rafilter-build-time-log "${RAFILTER_BUILD_TIME_LOG}" \
        --time-log "${TIME_LOG}" \
        --read-type "${READS_TYPE}" \
        --output-report "${RUN_OUT_DIR}/temp_eval_result.txt"
fi

cat "${RUN_OUT_DIR}/temp_eval_result.txt" >> "${REPORT_FILE}"
rm "${RUN_OUT_DIR}/temp_eval_result.txt"

echo "========================================================="
echo "  全部流程执行完毕"
echo "  结果目录：${RUN_OUT_DIR}"
echo "  最终评估报告：${REPORT_FILE}"
echo "========================================================="
