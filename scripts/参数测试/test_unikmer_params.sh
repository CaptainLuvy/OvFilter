#!/bin/bash
# ==============================================================================
# 脚本名称: test_unikmer_params.sh
# 功能: 评估不同 k 值和阈值提取 Unikmer 的精确率、召回率和 F1 分数。
# 逻辑: 
#   1. 使用 Jellyfish 从【参考基因组 (Reference)】提取单拷贝 (count=1) 的 k-mer 作为【真值集合 (Ground Truth)】。
#   2. 使用 Jellyfish 从【Reads 数据】提取 k-mer 频率分布。
#   3. 根据指定的 k 值和不同的阈值窗口 (如 mu ± 1*sd, mu ± 2*sd, mu ± 3*sd) 提取【测试集合 (Test Set)】。
#   4. 使用 Python 脚本比对真值集合和测试集合，计算 Precision, Recall, F1。
# ==============================================================================

set -e

# --- 配置区 (根据实际情况修改) ---
REF_FILE="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/chm13/regions/chrX_57319763_61427195/ref.fa"  # 替换为你的参考基因组文件
READS_FILE="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/chm13/regions/chrX_57319763_61427195/reads.fq"    # 替换为你的 reads 文件
OUT_DIR="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/chm13/regions/chrX_57319763_61427195/param_test"
THREADS=16

# 测试的 K 值列表
K_VALUES=(17 19 21 23 25)

# 测试的阈值窗口乘数 (如 1 表示 mu ± 1*sd)
SD_MULTIPLIERS=(1 2 3 4)

mkdir -p "$OUT_DIR"
LOG_FILE="$OUT_DIR/evaluation_results.txt"
echo "K_Value SD_Multiplier True_Unikmers Extracted_Unikmers Correct_Unikmers Precision Recall F1_Score" > "$LOG_FILE"

# --- 核心循环 ---
for k in "${K_VALUES[@]}"; do
    echo "========================================="
    echo "正在测试 K = $k"
    echo "========================================="
    
    K_DIR="$OUT_DIR/k${k}"
    mkdir -p "$K_DIR"
    
    # 1. 提取真值集合 (Ground Truth)
    # 从参考基因组中严格提取 count = 1 的 k-mer
    TRUTH_COUNT="$K_DIR/ref_k${k}.jf"
    TRUTH_KMERS="$K_DIR/ref_k${k}_truth.kmers"
    
    echo "[1/4] 提取参考基因组单拷贝真值 (k=$k)..."
    if [ ! -s "$TRUTH_KMERS" ]; then
        jellyfish count -m "$k" -C -s 100M -t "$THREADS" -o "$TRUTH_COUNT" "$REF_FILE"
        # 仅 dump 出 count 为 1 的 k-mer (单拷贝)
        jellyfish dump -c -L 1 -U 1 "$TRUTH_COUNT" | awk '{print $1}' > "$TRUTH_KMERS"
    fi
    num_truth=$(wc -l < "$TRUTH_KMERS")
    echo "      真值 Unikmer 数量: $num_truth"
    
    # 2. 从 Reads 中计算频率分布
    READS_COUNT="$K_DIR/reads_k${k}.jf"
    READS_HISTO="$K_DIR/reads_k${k}.histo"
    READS_XCEL="$K_DIR/reads_k${k}.xlsx"
    
    echo "[2/4] 计算 Reads 的 K-mer 分布 (k=$k)..."
    if [ ! -s "$READS_HISTO" ]; then
        jellyfish count -m "$k" -C -o "$READS_COUNT" -c 3 -s 10M --disk -t "$THREADS" "$READS_FILE"
        jellyfish histo -o "$READS_HISTO" -v "$READS_COUNT"
    fi
    
    # 3. 使用 Python 计算统计量 (mu 和 sd)
    # 假设你之前用的 freq_distribution_kmers.py 支持传入文件并输出 mean 和 sd
    # 我们这里调用它来获取 mu 和 sd
    echo "[3/4] 估算频率峰值和标准差..."
    OUTPUT=$(python3 /homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/freq_distribution_kmers.py "$READS_HISTO" "$READS_XCEL")
    
    # 从 python 输出中提取 mean 和 sd (适配格式: "mean: 30.3 std: 7.29")
    mu=$(echo "$OUTPUT" | awk '/mean:/ {print $2}')
    sd=$(echo "$OUTPUT" | awk '/std:/ {print $4}')
    
    if [ -z "$mu" ] || [ -z "$sd" ]; then
        echo "      [Debug] Python 输出为: $OUTPUT"
        echo "警告: 无法从 python 脚本输出中解析 mu 或 sd。跳过此 k 值。"
        continue
    fi
    echo "      mu = $mu, sd = $sd"
    
    # 4. 循环测试不同的阈值乘数
    for mult in "${SD_MULTIPLIERS[@]}"; do
        echo "  --- 测试阈值: mu ± ${mult}*sd ---"
        
        # 计算上下界
        lower=$(echo "$mu - $mult * $sd" | bc | awk '{printf "%d", $1}')
        upper=$(echo "$mu + $mult * $sd" | bc | awk '{printf "%d", $1+0.5}')
        
        if [ "$lower" -lt 1 ]; then lower=1; fi
        
        echo "      提取窗口: [$lower, $upper]"
        
        TEST_KMERS="$K_DIR/test_k${k}_m${mult}.kmers"
        jellyfish dump -c -t -L "$lower" -U "$upper" "$READS_COUNT" | awk '{print $1}' > "$TEST_KMERS"
        
        # 5. 使用 Python 脚本计算评估指标
        # 编写一个临时的内联 python 脚本来快速计算集合交集
        EVAL_OUT=$(python3 -c "
import sys

truth_file = '$TRUTH_KMERS'
test_file = '$TEST_KMERS'

# 读取真值集合
with open(truth_file, 'r') as f:
    truth_set = set(line.strip() for line in f)

# 读取测试集合
with open(test_file, 'r') as f:
    test_set = set(line.strip() for line in f)

# 计算交集
correct = len(truth_set.intersection(test_set))
num_truth = len(truth_set)
num_test = len(test_set)

precision = correct / num_test if num_test > 0 else 0
recall = correct / num_truth if num_truth > 0 else 0
f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

print(f'{num_truth} {num_test} {correct} {precision:.4f} {recall:.4f} {f1:.4f}')
")
        
        # 将结果写入日志
        echo "$k $mult $EVAL_OUT" >> "$LOG_FILE"
        echo "      结果: Precision=$(echo $EVAL_OUT | awk '{print $4}'), Recall=$(echo $EVAL_OUT | awk '{print $5}'), F1=$(echo $EVAL_OUT | awk '{print $6}')"
    done
done

echo "========================================="
echo "所有测试完成！结果如下："
echo "========================================="
# 使用 column 命令格式化输出日志文件内容
column -t "$LOG_FILE"
echo "========================================="
echo "结果已同时保存在 $LOG_FILE"

# 自动找出指标最高的参数组合并打印和保存
echo "========================================="
echo "⭐ 最优参数组合汇总 ⭐"

# 1. Highest F1
BEST_F1_RESULT=$(tail -n +2 "$LOG_FILE" | sort -k8 -n -r | head -n 1)
F1_K=$(echo "$BEST_F1_RESULT" | awk '{print $1}')
F1_MULT=$(echo "$BEST_F1_RESULT" | awk '{print $2}')
F1_VAL=$(echo "$BEST_F1_RESULT" | awk '{print $8}')

# 2. Highest Precision
BEST_PREC_RESULT=$(tail -n +2 "$LOG_FILE" | sort -k6 -n -r | head -n 1)
PREC_K=$(echo "$BEST_PREC_RESULT" | awk '{print $1}')
PREC_MULT=$(echo "$BEST_PREC_RESULT" | awk '{print $2}')
PREC_VAL=$(echo "$BEST_PREC_RESULT" | awk '{print $6}')

# 3. Highest Recall
BEST_REC_RESULT=$(tail -n +2 "$LOG_FILE" | sort -k7 -n -r | head -n 1)
REC_K=$(echo "$BEST_REC_RESULT" | awk '{print $1}')
REC_MULT=$(echo "$BEST_REC_RESULT" | awk '{print $2}')
REC_VAL=$(echo "$BEST_REC_RESULT" | awk '{print $7}')

echo "[最高 F1 Score]: $F1_VAL (K=$F1_K, 阈值=mu±${F1_MULT}*sd)"
echo "[最高 Precision]: $PREC_VAL (K=$PREC_K, 阈值=mu±${PREC_MULT}*sd)"
echo "[最高 Recall]: $REC_VAL (K=$REC_K, 阈值=mu±${REC_MULT}*sd)"
echo "========================================="

# 将最优结果单独保存到一个文件中
# 使用 F1 最高的参数来命名文件，让结果一目了然
BEST_FILE="$OUT_DIR/best_params_k${F1_K}_sd${F1_MULT}.txt"
echo "=== Best Parameters Summary ===" > "$BEST_FILE"
echo "" >> "$BEST_FILE"
echo "[Highest F1 Score]" >> "$BEST_FILE"
echo "K_Value: $F1_K, SD_Multiplier: $F1_MULT, F1: $F1_VAL" >> "$BEST_FILE"
echo "Full_Stats: $BEST_F1_RESULT" >> "$BEST_FILE"
echo "" >> "$BEST_FILE"
echo "[Highest Precision]" >> "$BEST_FILE"
echo "K_Value: $PREC_K, SD_Multiplier: $PREC_MULT, Precision: $PREC_VAL" >> "$BEST_FILE"
echo "Full_Stats: $BEST_PREC_RESULT" >> "$BEST_FILE"
echo "" >> "$BEST_FILE"
echo "[Highest Recall]" >> "$BEST_FILE"
echo "K_Value: $REC_K, SD_Multiplier: $REC_MULT, Recall: $REC_VAL" >> "$BEST_FILE"
echo "Full_Stats: $BEST_REC_RESULT" >> "$BEST_FILE"

echo "你可以使用这些数据使用 matplotlib 绘制曲线图了。"
echo "========================================="
