#!/bin/bash
set -eu -o pipefail

# ==============================================================================
# 脚本名称: run_asm_single.sh
# 功能描述: 
#   针对单个基因组区域（如着丝粒）执行完整的组装流程，并对比 OvFilter 过滤前后的效果。
#   该脚本自动化执行从 Reads 比对到最终组装的全过程。
#
# 核心流程:
#   1. 环境检查: 验证输入目录、reads 文件及 unikmer 文件是否存在。
#   2. 原始重叠生成 (Minimap2): 使用 ava-pb 模式生成全对全比对结果 (overlaps.raw.paf)。
#   3. 重叠过滤 (OvFilter): 调用 ovfilter.py 基于 unikmer 和几何距离标准差清洗重叠 (overlaps.ovfilter.paf)。
#   4. 过滤后组装 (Miniasm): 使用清洗后的重叠图进行组装，生成 GFA 和 FASTA (assembly.fasta)。
#   5. 原始组装对比 (Miniasm): 使用原始重叠图进行组装，生成对照组结果 (assembly.raw.fasta)。
#   6. 结果统计: 输出过滤前后组装结果的 Contig 数量对比。
#
# 输入依赖:
#   - 区域目录需包含: reads.fq
#   - 区域目录需包含: unikmer/k21_unikmers.kmers (需预先生成)
#
# 用法: bash run_asm_single.sh <region_dir>
# 示例: bash run_asm_single.sh centromere_data/regions/chr10_39133793_42426237
# ==============================================================================

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <region_directory>"
    echo "Example: $0 hg002/regions/chr10_39133793_42426237"
    exit 1
fi

REGION_DIR="$1"

# 转为绝对路径，防止目录切换导致找不到
if [ -d "$REGION_DIR" ]; then
    REGION_DIR=$(cd "$REGION_DIR" && pwd)
else
    echo "Error: Directory '$REGION_DIR' not found."
    exit 1
fi

REGION_NAME=$(basename "$REGION_DIR")

# --- 配置参数 ---
THREADS=${THREADS:-32}

# Minimap2 额外参数
# -f 0.001: 过滤掉最高频 1% 的 k-mer (平衡速度与着丝粒区域的连通性)
MINIMAP2_OPTS=${MINIMAP2_OPTS:-"-f 0.001"}

MINIMAP2="minimap2"
MINIASM="miniasm"
PYTHON="python3"
AWK="awk"

# 自动定位 ovfilter.py (假设在 scripts/ 上一级)
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OVFILTER_PY="$SCRIPT_DIR/../ovfilter.py"

# Ovfilter 参数 (可根据需要调整)
KMER_LEN=21
MAX_SD=10.0
MIN_SHARED=15
TOLERANCE=15
OVERHANG_TOL=5000 # 保持放宽 (原5000) 以减少误杀
MIN_IDENTITY=0    # 新增: 最小 Identity 过滤，剔除低质量比对

# --- 环境检查 ---
if [ ! -f "$OVFILTER_PY" ]; then
    echo "Error: ovfilter.py not found at $OVFILTER_PY"
    exit 1
fi

echo "=========================================================="
echo "   Processing Single Region: $REGION_NAME"
echo "   Dir: $REGION_DIR"
echo "=========================================================="

# 定义文件路径
READS_FILE="$REGION_DIR/reads.fq"
UNIKMER_FILE="$REGION_DIR/unikmer/k${KMER_LEN}_unikmers.kmers"
PAF_RAW="$REGION_DIR/overlaps.raw.paf"
PAF_FILTERED="$REGION_DIR/overlaps.ovfilter.paf"
PAF_NUMERIC="$REGION_DIR/overlaps.numeric.paf"
PAF_DISCARDED="$REGION_DIR/overlaps.discarded.paf"
ASM_GFA="$REGION_DIR/assembly.gfa"
ASM_FASTA="$REGION_DIR/assembly.fasta"
ASM_NUMERIC_GFA="$REGION_DIR/assembly.numeric.gfa"
ASM_NUMERIC_FASTA="$REGION_DIR/assembly.numeric.fasta"
LOG_FILE="$REGION_DIR/assembly_single.log"

echo "   Logs will be saved to: $LOG_FILE"

# 检查输入
if [ ! -f "$READS_FILE" ]; then
    echo "Error: Reads file not found: $READS_FILE"
    exit 1
fi

if [ ! -f "$UNIKMER_FILE" ]; then
    echo "Error: Unikmer file not found: $UNIKMER_FILE"
    echo "Please run unikmer extraction first."
    exit 1
fi

# --- Step 1: Minimap2 (智能跳过) ---
if [ -f "$PAF_RAW" ]; then
    size=$(du -h "$PAF_RAW" | cut -f1)
    echo "   [1/4] Found existing raw PAF ($size). SKIPPING minimap2."
else
    echo "   [1/4] Running minimap2 (ava-pb)..."
    $MINIMAP2 -t "$THREADS" $MINIMAP2_OPTS -x ava-pb "$READS_FILE" "$READS_FILE" > "$PAF_RAW" 2>> "$LOG_FILE"
fi

# --- Step 2: Ovfilter (Run Once, Generate Both Outputs) ---
if [ -f "$PAF_FILTERED" ] && [ -s "$PAF_FILTERED" ] && [ -f "$PAF_NUMERIC" ] && [ -s "$PAF_NUMERIC" ]; then
     echo "   [2/4] Found existing filtered and numeric PAFs. SKIPPING ovfilter."
else
    echo "   [2/4] Running ovfilter (generating both filtered and numeric outputs)..."
    # Added --keep-contained to preserve connectivity
    if $PYTHON "$OVFILTER_PY" \
        -p "$PAF_RAW" \
        -u "$UNIKMER_FILE" \
        -r "$READS_FILE" \
        -o "$PAF_FILTERED" \
        --output-numeric "$PAF_NUMERIC" \
        -k "$KMER_LEN" \
        -s "$MIN_SHARED" \
        -t "$TOLERANCE" \
        --max-sd "$MAX_SD" \
        --overhang-tolerance "$OVERHANG_TOL" \
        --min-identity "$MIN_IDENTITY" \
        --keep-contained \
        -T "$THREADS" \
        >> "$LOG_FILE" 2>&1; then
        
        echo "         Ovfilter completed successfully."
    else
        echo "   Error: ovfilter failed. Last 10 lines of log:"
        echo "   ----------------------------------------"
        tail -n 10 "$LOG_FILE"
        echo "   ----------------------------------------"
        exit 1
    fi
fi

# --- Step 3: Miniasm (Filtered) ---
echo "   [3/4] Running miniasm (Filtered)..."
$MINIASM -f "$READS_FILE" "$PAF_FILTERED" > "$ASM_GFA" 2>> "$LOG_FILE"

# --- Step 4: GFA to FASTA (Filtered) ---
echo "   [4/4] Converting GFA to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' "$ASM_GFA" > "$ASM_FASTA"

# --- Step 4.5: Miniasm & GFA to FASTA (Numeric Only) ---
echo "   [4.5/4] Running miniasm (Numeric Only)..."
$MINIASM -f "$READS_FILE" "$PAF_NUMERIC" > "$ASM_NUMERIC_GFA" 2>> "$LOG_FILE"
awk '/^S/{print ">"$2"\n"$3}' "$ASM_NUMERIC_GFA" > "$ASM_NUMERIC_FASTA"
num_numeric_contigs=$(grep -c "^>" "$ASM_NUMERIC_FASTA")

# --- Step 5: Miniasm (Raw - for comparison) ---
ASM_RAW_GFA="$REGION_DIR/assembly.raw.gfa"
ASM_RAW_FASTA="$REGION_DIR/assembly.raw.fasta"

echo "   [5/5] Running miniasm on RAW PAF for comparison..."
if [ -f "$PAF_RAW" ]; then
    $MINIASM -f "$READS_FILE" "$PAF_RAW" > "$ASM_RAW_GFA" 2>> "$LOG_FILE"
    awk '/^S/{print ">"$2"\n"$3}' "$ASM_RAW_GFA" > "$ASM_RAW_FASTA"
    num_raw_contigs=$(grep -c "^>" "$ASM_RAW_FASTA")
else
    echo "   Warning: Raw PAF not found, skipping raw assembly."
    num_raw_contigs="N/A"
fi

num_contigs=$(grep -c "^>" "$ASM_FASTA")
echo "=========================================================="
echo "   Done!"
echo "   Filtered Assembly: $ASM_FASTA ($num_contigs contigs)"
echo "   Numeric Assembly:  $ASM_NUMERIC_FASTA ($num_numeric_contigs contigs)"
echo "   Raw Assembly:      $ASM_RAW_FASTA ($num_raw_contigs contigs)"
echo "=========================================================="
