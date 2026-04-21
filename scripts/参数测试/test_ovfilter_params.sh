#!/bin/bash
# 移除 set -e 以免中间命令报错直接退出脚本
# set -e
# set -o pipefail

# 尝试使用通用的 UTF-8，避免 C.UTF-8 找不到的 warning
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export PYTHONIOENCODING=UTF-8

# ==============================================================================
# 脚本名称: test_ovfilter_params.sh
# 功能: 针对给定的原始 PAF 和真实 MAF，测试 ovfilter 不同参数组合 (Identity, Min_Shared, Max_SD) 
#       对 Precision, Recall, F1 的影响，并找出最优组合。
# ==============================================================================

# --- 1. 核心输入文件 (根据你的实际情况修改) ---
# 注意：为了节省时间，这里假设 minimap2 raw paf, kmer 文件, 和 combined fasta 已经生成好了。
# 如果没有，请先运行一次 run_evaluation_pipeline.sh
BASE_DIR="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp"
TRUTH_MAF="${BASE_DIR}/sim_data/chr8/hifi/chr8_hifi_0001.maf"
RAW_PAF="${BASE_DIR}/result/chr8/hifi30x_secondary/run/minimap2/minimap2_raw.paf"
UNIKMER_FILE="${BASE_DIR}/result/chr8/hifi30x_secondary/run/tmp/ovfilter.kmer"
COMBINED_FA="${BASE_DIR}/result/chr8/hifi30x_secondary/run/tmp/combined_for_ovfilter.fa"

# 程序和脚本路径
OVFILTER_CPP_BIN="/home/dingyc/tools/ovfilter/ovfilter_cpp"
EVAL_PY="${BASE_DIR}/result/evaluate_alignments.py"

# 输出目录
OUT_DIR="${BASE_DIR}/result/chr8/hifi30x_secondary/param_test_ovfilter"
mkdir -p "${OUT_DIR}"
LOG_FILE="${OUT_DIR}/ovfilter_param_search.log"

# --- 2. 固定参数 ---
THREADS=32
READS_TYPE="hifi"     # hifi 或 ont
OV_K=21
OV_TOLERANCE=15
OV_MIN_ALIGN_LEN=0
OV_OVERHANG=5000
OV_EXTRA_ARGS="--keep-contained --ignore-overhang"

# --- 3. 网格搜索的参数空间 (可根据需要修改范围) ---
MIN_SHARED_VALS=(1 3 6 9)
# MIN_SHARED_VALS=(6)
# 测试 RPS (Repeat Penalty Score) 阈值
RPS_VALS=("0.0" "0.04" "0.08" "0.12" "0.16")

echo "Min_Shared RPS_Threshold Precision Recall F1_Score" > "${LOG_FILE}"

echo "========================================================="
echo " 开始 ovfilter 参数网格搜索 (使用 RPS 特赦)"
echo " 输出目录：${OUT_DIR}"
echo "========================================================="

# --- 4. 核心网格搜索循环 ---
for shared in "${MIN_SHARED_VALS[@]}"; do
    for rps in "${RPS_VALS[@]}"; do
        echo "---------------------------------------------------------"
        echo "测试组合: Min_Shared=${shared}, RPS_Threshold=${rps}"
        
        # 定义当前组合的输出文件
        PREFIX="sh${shared}_rps${rps}"
        OUT_PAF="${OUT_DIR}/ovfilter_${PREFIX}.paf"
        EVAL_OUT="${OUT_DIR}/eval_${PREFIX}.txt"
        
        # 运行 ovfilter
        ${OVFILTER_CPP_BIN} \
            --read-type ${READS_TYPE} \
            --paf "${RAW_PAF}" \
            --unikmers "${UNIKMER_FILE}" \
            --reads "${COMBINED_FA}" \
            --output "${OUT_PAF}" \
            --threads ${THREADS} \
            -k ${OV_K} \
            -s ${shared} \
            -t ${OV_TOLERANCE} \
            -R ${rps} \
            -l ${OV_MIN_ALIGN_LEN} \
            -e ${OV_OVERHANG} \
            ${OV_EXTRA_ARGS}
                        
                    if [ $? -ne 0 ]; then
                        echo "❌ 警告: ovfilter 运行失败，跳过该参数组合。"
                        continue
                    fi
                        
                    # 运行评估脚本 (去掉隐藏输出以便排错)
                    python "${EVAL_PY}" \
                        --maf "${TRUTH_MAF}" \
                        --raw-paf "${RAW_PAF}" \
                        --rafilter-paf "${RAW_PAF}" \
                        --ovfilter-paf "${OUT_PAF}" \
                        --ovfilter-label "Test_${PREFIX}" \
                        --read-type "${READS_TYPE}" \
                        --output-report "${EVAL_OUT}"
                        
                    if [ ! -f "${EVAL_OUT}" ]; then
                        echo "❌ 警告: evaluate_alignments.py 没有生成报告文件，跳过该参数组合。"
                        continue
                    fi
                        
                    # 从 evaluate_alignments.py 的输出中提取指标
                    prec=$(awk -v prefix="Test_${PREFIX}" '$0 ~ prefix && $0 ~ /\|/ {print $(NF-4)}' "${EVAL_OUT}")
                    rec=$(awk -v prefix="Test_${PREFIX}" '$0 ~ prefix && $0 ~ /\|/ {print $(NF-2)}' "${EVAL_OUT}")
                    f1=$(awk -v prefix="Test_${PREFIX}" '$0 ~ prefix && $0 ~ /\|/ {print $NF}' "${EVAL_OUT}")
                    
                    if [ -z "$prec" ] || [ -z "$rec" ] || [ -z "$f1" ]; then
                        echo "警告: 无法解析评估结果。请检查 evaluate_alignments.py 的输出格式。"
                        prec="0"
                        rec="0"
                        f1="0"
                    fi

                    echo "结果: Precision=${prec}, Recall=${rec}, F1=${f1}"
                    echo "${shared} ${rps} ${prec} ${rec} ${f1}" >> "${LOG_FILE}"
                    
                    # 清理临时 PAF 以节省空间
                    rm -f "${OUT_PAF}" "${EVAL_OUT}"
    done
done

echo "========================================================="
echo "所有测试完成！结果如下："
echo "========================================================="
column -t "${LOG_FILE}"
echo "========================================================="

# --- 5. 自动找出指标最高的参数组合 ---
echo "✨ 最优参数组合汇总 ✨"

# 1. Highest F1 (列索引改变：shared=1, rps=2, prec=3, rec=4, f1=5)
BEST_F1_RESULT=$(tail -n +2 "${LOG_FILE}" | sort -k5 -n -r | head -n 1)
F1_SH=$(echo "$BEST_F1_RESULT" | awk '{print $1}')
F1_RPS=$(echo "$BEST_F1_RESULT" | awk '{print $2}')
F1_VAL=$(echo "$BEST_F1_RESULT" | awk '{print $5}')

# 2. Highest Precision
BEST_PREC_RESULT=$(tail -n +2 "${LOG_FILE}" | sort -k3 -n -r | head -n 1)
PREC_SH=$(echo "$BEST_PREC_RESULT" | awk '{print $1}')
PREC_RPS=$(echo "$BEST_PREC_RESULT" | awk '{print $2}')
PREC_VAL=$(echo "$BEST_PREC_RESULT" | awk '{print $3}')

# 3. Highest Recall
BEST_REC_RESULT=$(tail -n +2 "${LOG_FILE}" | sort -k4 -n -r | head -n 1)
REC_SH=$(echo "$BEST_REC_RESULT" | awk '{print $1}')
REC_RPS=$(echo "$BEST_REC_RESULT" | awk '{print $2}')
REC_VAL=$(echo "$BEST_REC_RESULT" | awk '{print $4}')

echo "[最高 F1 Score]: $F1_VAL (Min_Shared=$F1_SH, RPS_Threshold=$F1_RPS)"
echo "[最高 Precision]: $PREC_VAL (Min_Shared=$PREC_SH, RPS_Threshold=$PREC_RPS)"
echo "[最高 Recall]: $REC_VAL (Min_Shared=$REC_SH, RPS_Threshold=$REC_RPS)"
echo "========================================================="

BEST_FILE="${OUT_DIR}/best_params_sh${F1_SH}_rps${F1_RPS}.txt"
echo "=== Best Parameters Summary ===" > "$BEST_FILE"
echo "" >> "$BEST_FILE"
echo "[Highest F1 Score]" >> "$BEST_FILE"
echo "Min_Shared: $F1_SH, RPS_Threshold: $F1_RPS, F1: $F1_VAL" >> "$BEST_FILE"
echo "" >> "$BEST_FILE"
echo "[Highest Precision]" >> "$BEST_FILE"
echo "Min_Shared: $PREC_SH, RPS_Threshold: $PREC_RPS, Precision: $PREC_VAL" >> "$BEST_FILE"
echo "" >> "$BEST_FILE"
echo "[Highest Recall]" >> "$BEST_FILE"
echo "Min_Shared: $REC_SH, RPS_Threshold: $REC_RPS, Recall: $REC_VAL" >> "$BEST_FILE"
