#!/bin/bash
# ==============================================================================
# 脚本名称: run_ablation_map.sh
# 功能: 针对 Read-to-Ref (map) 模式下的 PAF 进行消融实验对比。
# 实验设计:
#   1. Raw: 无任何过滤 (Minimap2 原始输出)
#   2. Basic: 只有基础过滤 (自比对、overhang=1000、默认identity)
#   3. Unikmer: 只有 Unikmer 相关共线性过滤 (包含 RPS 特赦，关闭 overhang 和 identity 过滤)
#   4. Full: 全功能过滤 (基础过滤 + Unikmer + RPS特赦，关闭 overhang 过滤防止影响召回)# ==============================================================================
#
# 用法: bash run_ablation_map.sh [数据类型(hifi/ont)] [染色体(chr8/chr10...)]
# 示例: bash run_ablation_map.sh hifi chr8
# ==============================================================================

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export PYTHONIOENCODING=UTF-8

# 接收数据类型参数 (hifi 或 ont)
READS_TYPE=${1:-"hifi"}
# 接收染色体参数 (默认为 chr8)
CHR_NAME=${2:-"chr8"}

# 基础路径配置
BASE_DIR="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp"
OVFILTER_CPP_BIN="/home/dingyc/tools/ovfilter/ovfilter_cpp"
EVAL_PY="${BASE_DIR}/result/evaluate_alignments.py"

if [ "$READS_TYPE" == "ont" ]; then
    echo "当前模式: ONT 数据, 染色体: ${CHR_NAME}"
    TRUTH_MAF="${BASE_DIR}/sim_data/${CHR_NAME}/ont/${CHR_NAME}_ont_0001.maf"
    RAW_PAF="${BASE_DIR}/result/${CHR_NAME}/ont30x_secondary/run/minimap2/minimap2_raw.paf"
    UNIKMER_FILE="${BASE_DIR}/result/${CHR_NAME}/ont30x_secondary/run/tmp/ovfilter.kmer"
    COMBINED_FA="${BASE_DIR}/result/${CHR_NAME}/ont30x_secondary/run/tmp/combined_for_ovfilter.fa"
    OUT_DIR="${BASE_DIR}/result/${CHR_NAME}/ont30x_secondary/ablation_map"
    MIN_SHARED=6
else
    echo "当前模式: HiFi 数据, 染色体: ${CHR_NAME}"
    TRUTH_MAF="${BASE_DIR}/sim_data/${CHR_NAME}/hifi/${CHR_NAME}_hifi_0001.maf"
    RAW_PAF="${BASE_DIR}/result/${CHR_NAME}/hifi30x_secondary/run/minimap2/minimap2_raw.paf"
    UNIKMER_FILE="${BASE_DIR}/result/${CHR_NAME}/hifi30x_secondary/run/tmp/ovfilter.kmer"
    COMBINED_FA="${BASE_DIR}/result/${CHR_NAME}/hifi30x_secondary/run/tmp/combined_for_ovfilter.fa"
    OUT_DIR="${BASE_DIR}/result/${CHR_NAME}/hifi30x_secondary/ablation_map"
    MIN_SHARED=1
fi

mkdir -p "${OUT_DIR}"
LOG_FILE="${OUT_DIR}/ablation_results.log"
echo "Group Precision Recall F1_Score" > "${LOG_FILE}"

# 固定参数
THREADS=32
OV_K=21
OV_TOLERANCE=15
OV_OVERHANG=1000
RPS_FULL=0.08
# 实验设计中提到基础过滤包括“自比对、overhang、identity”，
# 这里直接留空，因为 ovfilter.cpp 已经默认保留 contained 了。
OV_EXTRA_ARGS=""

echo "========================================================="
echo " 开始消融实验 (Map 模式: ${READS_TYPE})"
echo " 固定参数: Min_Shared=${MIN_SHARED}, Overhang=${OV_OVERHANG}, RPS=${RPS_FULL}"
echo "========================================================="

# 辅助函数: 运行评估并提取指标
run_eval() {
    local group_name=$1
    local paf_file=$2
    local eval_out="${OUT_DIR}/eval_${group_name}.txt"
    
    echo "   正在评估 ${group_name} ..."
    python "${EVAL_PY}" \
        --maf "${TRUTH_MAF}" \
        --raw-paf "${RAW_PAF}" \
        --rafilter-paf "${RAW_PAF}" \
        --ovfilter-paf "${paf_file}" \
        --ovfilter-label "${group_name}" \
        --read-type "${READS_TYPE}" \
        --output-report "${eval_out}" > /dev/null 2>&1
        
    prec=$(awk -v prefix="${group_name}" '$0 ~ prefix && $0 ~ /\|/ {print $(NF-4); exit}' "${eval_out}")
    rec=$(awk -v prefix="${group_name}" '$0 ~ prefix && $0 ~ /\|/ {print $(NF-2); exit}' "${eval_out}")
    f1=$(awk -v prefix="${group_name}" '$0 ~ prefix && $0 ~ /\|/ {print $NF; exit}' "${eval_out}")
    
    echo "   结果 [${group_name}]: Precision=${prec}, Recall=${rec}, F1=${f1}"
    echo "${group_name} ${prec} ${rec} ${f1}" >> "${LOG_FILE}"
}

# ---------------------------------------------------------
# Group 1: Raw (无任何过滤)
# ---------------------------------------------------------
echo "组别 1: Raw (无任何过滤)"
run_eval "Raw" "${RAW_PAF}"

# ---------------------------------------------------------
# Group 2: Basic (只有基础过滤: Self, Overhang, Identity)
# ---------------------------------------------------------
echo "组别 2: Basic (只有基础过滤)"
PAF_BASIC="${OUT_DIR}/overlaps_basic.paf"
${OVFILTER_CPP_BIN} \
    --read-type ${READS_TYPE} \
    --paf "${RAW_PAF}" \
    --unikmers "${UNIKMER_FILE}" \
    --reads "${COMBINED_FA}" \
    --output "${PAF_BASIC}" \
    --threads ${THREADS} \
    -k ${OV_K} \
    -e ${OV_OVERHANG} \
    --numeric-only \
    ${OV_EXTRA_ARGS} > /dev/null 2>&1

run_eval "Basic" "${PAF_BASIC}"

# ---------------------------------------------------------
# Group 3: Unikmer (只有 Unikmer共线性及 RPS特赦, 无基础数值过滤)
# ---------------------------------------------------------
echo "组别 3: Unikmer (只有 Unikmer共线性 + RPS特赦, 关闭基础过滤)"
PAF_UNIKMER="${OUT_DIR}/overlaps_unikmer.paf"
${OVFILTER_CPP_BIN} \
    --read-type ${READS_TYPE} \
    --paf "${RAW_PAF}" \
    --unikmers "${UNIKMER_FILE}" \
    --reads "${COMBINED_FA}" \
    --output "${PAF_UNIKMER}" \
    --threads ${THREADS} \
    -k ${OV_K} \
    -s ${MIN_SHARED} \
    -t ${OV_TOLERANCE} \
    --ignore-overhang \
    -i 0.0 \
    -R ${RPS_FULL} > /dev/null 2>&1

run_eval "Unikmer" "${PAF_UNIKMER}"

# ---------------------------------------------------------
# Group 4: Full (全功能过滤: 基础 + Unikmer + RPS特赦)
# ---------------------------------------------------------
echo "组别 4: Full (全功能过滤: 基础 + Unikmer + RPS=${RPS_FULL}, 关闭 overhang)"
PAF_FULL="${OUT_DIR}/overlaps_full.paf"
${OVFILTER_CPP_BIN} \
    --read-type ${READS_TYPE} \
    --paf "${RAW_PAF}" \
    --unikmers "${UNIKMER_FILE}" \
    --reads "${COMBINED_FA}" \
    --output "${PAF_FULL}" \
    --threads ${THREADS} \
    -k ${OV_K} \
    -s ${MIN_SHARED} \
    -t ${OV_TOLERANCE} \
    --ignore-overhang \
    -R ${RPS_FULL} \
    ${OV_EXTRA_ARGS} > /dev/null 2>&1

run_eval "Full" "${PAF_FULL}"

echo "========================================================="
echo " ${CHR_NAME} ${READS_TYPE}消融实验完成！结果如下："
echo "========================================================="
column -t "${LOG_FILE}"
echo "========================================================="
