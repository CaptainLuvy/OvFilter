#!/bin/bash

# ==============================================================================
# 脚本名称: run_ovfilter_benchmark.sh
# 功能: 验证 PAF 过滤效果。
#       1. 自动生成 PAF (如果不存在): Minimap2 -> PAF
#       2. 对比组装结果: "Minimap2 -> Miniasm" vs "Minimap2 -> 过滤 -> Miniasm"
# ==============================================================================

# --- 1. 参数配置 (请根据实际情况修改这里) ---

# 输入文件
READS_FILE="/homeb/dingyc/fsa/hifi/chm13/ovfilter/chrX/repeat_data/chrX.fasta"          # 你的 reads 文件
UNIKMER_FILE="/homeb/dingyc/fsa/hifi/chm13/ovfilter/chrX/unikmer/k21_unikmers.kmers"       # 你的独特 kmer 文件

# PAF 文件路径设置
# 脚本会优先查找此文件。如果不存在，则检查 $OUT_DIR 下是否有 overlaps.paf。
# 如果都不存在，则调用 minimap2 生成。
RAW_PAF_FILE="overlaps.paf"

# 输出目录
OUT_DIR="/homeb/dingyc/fsa/hifi/chm13/ovfilter/chrX/output"

# 工具路径
FILTER_SCRIPT="/home/dingyc/tools/ovfilter/ovfilter.py"       # 你的 Python 过滤脚本路径
MINIASM_BIN="miniasm"             # miniasm 可执行文件
MINIMAP2_BIN="minimap2"           # minimap2 可执行文件

# Minimap2 参数
THREADS=8
PRESET="ava-pb"                   # ava-pb: PacBio HiFi/CLR, ava-ont: Nanopore

# 过滤参数 (传给 ovfilter.py)
K_SIZE=21
MIN_SHARED=15
TOLERANCE=15
MAX_SD=10.0
# KEEP_CONTAINED_FLAG="--keep-contained" # 如果需要保留被包含的重叠，取消注释此行

# ==============================================================================

# 创建输出目录
mkdir -p "$OUT_DIR"

# 0. 检查依赖
echo "Checking dependencies..."
for tool in "$MINIASM_BIN" "$MINIMAP2_BIN"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "[Error] 未找到 $tool，请确保已安装并在 PATH 中。"
        # exit 1 # 暂时不强制退出，因为可能是在 Windows Git Bash 下，某些环境配置不同
    fi
done

if [ ! -f "$FILTER_SCRIPT" ]; then
    echo "[Error] 未找到 Python 过滤脚本: $FILTER_SCRIPT"
    exit 1
fi

if [ ! -f "$READS_FILE" ]; then
    echo "[Error] 未找到 Reads 文件: $READS_FILE"
    exit 1
fi

echo "=========================================================="
echo "   开始对比实验: Raw PAF vs Filtered PAF"
echo "=========================================================="

# --- 步骤 1: 准备 PAF 文件 ---
TARGET_PAF=""

# 优先级 1: 用户指定的 RAW_PAF_FILE
if [ -f "$RAW_PAF_FILE" ]; then
    echo "[Step 1] 使用指定的现有 PAF 文件: $RAW_PAF_FILE"
    TARGET_PAF="$RAW_PAF_FILE"

# 优先级 2: 输出目录下的 overlaps.paf
elif [ -f "$OUT_DIR/overlaps.paf" ]; then
    echo "[Step 1] 使用输出目录中现有的 PAF 文件: $OUT_DIR/overlaps.paf"
    TARGET_PAF="$OUT_DIR/overlaps.paf"

# 优先级 3: 运行 Minimap2 生成
else
    echo "[Step 1] 未找到 PAF 文件，运行 Minimap2 生成..."
    echo "   Command: $MINIMAP2_BIN -x $PRESET -t $THREADS $READS_FILE $READS_FILE > $OUT_DIR/overlaps.paf"
    
    $MINIMAP2_BIN -x "$PRESET" -t "$THREADS" "$READS_FILE" "$READS_FILE" > "$OUT_DIR/overlaps.paf"
    
    if [ $? -ne 0 ]; then 
        echo "[Error] Minimap2 运行失败！"
        exit 1
    fi
    TARGET_PAF="$OUT_DIR/overlaps.paf"
    echo "   -> PAF 生成完成: $TARGET_PAF"
fi


# --- 步骤 2: 运行你的过滤脚本 ---
echo "[Step 2] 正在运行 Python 脚本过滤 PAF..."
start_time=$(date +%s)

# 自动检测 python 命令
PYTHON_BIN="python3"
if ! command -v python3 &> /dev/null; then
    PYTHON_BIN="python"
fi

$PYTHON_BIN "$FILTER_SCRIPT" \
    --paf "$TARGET_PAF" \
    --unikmers "$UNIKMER_FILE" \
    --reads "$READS_FILE" \
    --output "$OUT_DIR/filtered.paf" \
    --discarded "$OUT_DIR/discarded.paf" \
    --k $K_SIZE \
    --min-shared $MIN_SHARED \
    --tolerance $TOLERANCE \
    --max-sd $MAX_SD \
    $KEEP_CONTAINED_FLAG

if [ $? -ne 0 ]; then echo "过滤脚本运行失败！"; exit 1; fi
end_time=$(date +%s)
echo "   -> 过滤完成，耗时 $((end_time - start_time)) 秒。"


# --- 步骤 3: 运行 Miniasm 组装 (对照组: 原始 PAF) ---
echo "[Step 3] 正在组装原始数据 (Raw PAF)..."
# miniasm -f reads.fasta overlaps.paf > output.gfa
$MINIASM_BIN -f "$READS_FILE" "$TARGET_PAF" > "$OUT_DIR/raw.gfa" 2> "$OUT_DIR/miniasm_raw.log"

# 将 GFA 转为 FASTA (提取 S 行)
awk '/^S/{print ">"$2"\n"$3}' "$OUT_DIR/raw.gfa" > "$OUT_DIR/raw.asm.fasta"


# --- 步骤 4: 运行 Miniasm 组装 (实验组: 过滤后 PAF) ---
echo "[Step 4] 正在组装过滤数据 (Filtered PAF)..."
$MINIASM_BIN -f "$READS_FILE" "$OUT_DIR/filtered.paf" > "$OUT_DIR/filtered.gfa" 2> "$OUT_DIR/miniasm_filtered.log"

# 将 GFA 转为 FASTA
awk '/^S/{print ">"$2"\n"$3}' "$OUT_DIR/filtered.gfa" > "$OUT_DIR/filtered.asm.fasta"


# --- 步骤 5: 统计对比结果 ---
echo "=========================================================="
echo "                  实验结果对比表"
echo "=========================================================="
printf "%-15s %-10s %-15s %-15s %-15s\n" "Group" "Contigs" "Total_Len(bp)" "Max_Len(bp)" "N50(bp)"
echo "----------------------------------------------------------"

# 定义统计函数 (使用 awk)
calc_stats() {
    local file=$1
    local label=$2
    if [ ! -s "$file" ]; then
        printf "%-15s %-10s %-15s %-15s %-15s\n" "$label" "0" "0" "0" "0"
        return
    fi
    awk -v lbl="$label" '
    BEGIN {len=0; n=0}
    /^>/ {if(len>0){L[n++]=len}; len=0; next}
    {len+=length($0)}
    END {
        if(len>0){L[n++]=len}
        asort(L)
        total=0;
        for(i=1;i<=n;i++){total+=L[i]}
        half=total/2
        sum=0; n50=0
        for(i=n;i>=1;i--){
            sum+=L[i]
            if(sum>=half){n50=L[i]; break}
        }
        printf "%-15s %-10d %-15d %-15d %-15d\n", lbl, n, total, L[n], n50
    }
    ' "$file"
}

# 计算并打印
calc_stats "$OUT_DIR/raw.asm.fasta" "Raw(原始)"
calc_stats "$OUT_DIR/filtered.asm.fasta" "Filtered(过滤后)"

echo "----------------------------------------------------------"
echo "结果文件已保存在: $OUT_DIR"
