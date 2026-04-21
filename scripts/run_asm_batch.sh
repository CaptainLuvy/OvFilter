#!/bin/bash
set -eu -o pipefail

# ==============================================================================
# 脚本名称: run_asm_batch.sh
# 功能描述: 
#   批量组装脚本。遍历指定目录下的所有着丝粒区域数据，针对每个区域执行：
#   Minimap2 (全对全比对) -> OvFilter (重叠过滤) -> Miniasm (组装)
#   同时会生成 Raw Assembly (不过滤) 作为对照。
#
# 输入:
#   - centromere_data/regions/ 下的各个区域子目录
#   - 每个子目录需包含: reads.fq 和 unikmer/k21_unikmers.kmers
#
# 输出:
#   - assembly.fasta (OvFilter 过滤后的组装结果)
#   - assembly.raw.fasta (原始组装结果，用于对比)
#   - 中间文件: overlaps.raw.paf, overlaps.ovfilter.paf, assembly.gfa
# ==============================================================================

# --- 1. 参数配置 ---

# 数据根目录 (默认 chm13/regions，可通过第一个参数传入)
DATA_DIR="${1:-chm13/regions}"

# 线程数
THREADS=${THREADS:-32}

# Minimap2 额外参数
# -f 0.01: 过滤掉最高频 1% 的 k-mer (平衡速度与着丝粒区域的连通性)
MINIMAP2_OPTS=${MINIMAP2_OPTS:-"-f 0.01"}

# 工具路径
MINIMAP2="minimap2"
MINIASM="miniasm"
PYTHON="python3"
AWK="awk"

# OvFilter 脚本路径
# 假设当前脚本在 scripts/ 目录下，ovfilter.py 在上一级目录
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OVFILTER_PY="$SCRIPT_DIR/../ovfilter.py"

# OvFilter 核心参数
KMER_LEN=21
MAX_SD=10.0      # 几何距离标准差阈值
MIN_SHARED=15    # 最小共享 unikmer 数量
TOLERANCE=15     # 几何距离容忍度
OVERHANG_TOL=5000
MIN_IDENTITY=0.8 # 新增: 最小 Identity 过滤

# ==============================================================================

# --- 2. 环境检查 ---
echo "[Check] 正在检查依赖工具..."
for tool in "$MINIMAP2" "$MINIASM" "$PYTHON" "$AWK"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: 未找到工具 '$tool'，请确保其在 PATH 环境变量中。"
        exit 1
    fi
done

if [ ! -f "$OVFILTER_PY" ]; then
    echo "Error: 未找到 ovfilter.py 脚本，路径: $OVFILTER_PY"
    # 尝试备用路径 (可选)
    # OVFILTER_PY="/path/to/ovfilter.py"
    exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
    echo "Error: 数据目录 '$DATA_DIR' 不存在。"
    exit 1
fi

echo "=========================================================="
echo "   开始批量组装流程 (Batch Assembly Pipeline)"
echo "   数据目录: $DATA_DIR"
echo "   流程: Minimap2 (ava-pb) -> OvFilter -> Miniasm"
echo "   (+ Raw Assembly 对照)"
echo "=========================================================="

# --- 3. 遍历处理每个区域 ---

count=0
# 遍历 DATA_DIR 下的所有条目
for region_path in "$DATA_DIR"/*; do
    if [ ! -d "$region_path" ]; then continue; fi
    
    region_name=$(basename "$region_path")
    echo "-> 正在处理区域: $region_name"
    
    # 定义关键文件路径
    READS_FILE="$region_path/reads.fq"
    UNIKMER_FILE="$region_path/unikmer/k${KMER_LEN}_unikmers.kmers"
    
    # 检查输入文件是否存在
    if [ ! -f "$READS_FILE" ]; then
        echo "   [Skip] 未找到 reads 文件: $READS_FILE"
        continue
    fi
    
    if [ ! -f "$UNIKMER_FILE" ]; then
        echo "   [Skip] 未找到 unikmer 文件: $UNIKMER_FILE"
        echo "   (请先运行 unikmer 提取步骤)"
        continue
    fi
    
    # 定义输出文件
    PAF_RAW="$region_path/overlaps.raw.paf"
    PAF_FILTERED="$region_path/overlaps.ovfilter.paf"
    ASM_GFA="$region_path/assembly.gfa"
    ASM_FASTA="$region_path/assembly.fasta"
    LOG_FILE="$region_path/assembly.log"
    
    # Raw 对照组文件
    ASM_RAW_GFA="$region_path/assembly.raw.gfa"
    ASM_RAW_FASTA="$region_path/assembly.raw.fasta"
    
    echo "   日志将保存至: $LOG_FILE"
    
    # --- Step 1: Minimap2 All-vs-All 比对 ---
    if [ ! -f "$PAF_RAW" ]; then
        echo "   [1/5] 运行 Minimap2 (ava-pb)..."
        # 使用 ava-pb 模式针对 PacBio HiFi/CCS 数据
        $MINIMAP2 -t "$THREADS" $MINIMAP2_OPTS -x ava-pb "$READS_FILE" "$READS_FILE" > "$PAF_RAW" 2>> "$LOG_FILE"
    else
        echo "   [1/5] 原始比对文件已存在，跳过。"
    fi
    
    # --- Step 2: OvFilter 过滤 ---
    if [ ! -f "$PAF_FILTERED" ]; then
        echo "   [2/5] 运行 OvFilter 过滤..."
        if $PYTHON "$OVFILTER_PY" \
            -p "$PAF_RAW" \
            -u "$UNIKMER_FILE" \
            -r "$READS_FILE" \
            -o "$PAF_FILTERED" \
            -k "$KMER_LEN" \
            -s "$MIN_SHARED" \
            -t "$TOLERANCE" \
            --max-sd "$MAX_SD" \
            --overhang-tolerance "$OVERHANG_TOL" \
            --min-identity "$MIN_IDENTITY" \
            -T "$THREADS" \
            >> "$LOG_FILE" 2>&1; then
            
            echo "         OvFilter 完成。"
        else
            echo "   Error: OvFilter 运行失败，请检查日志。"
            continue
        fi
    else
        echo "   [2/5] 过滤后比对文件已存在，跳过。"
    fi
    
    # --- Step 3: Miniasm 组装 (Filtered) ---
    if [ ! -f "$ASM_GFA" ]; then
        echo "   [3/5] 运行 Miniasm 组装 (Filtered)..."
        # Miniasm 需要原始 reads 和 paf 文件
        $MINIASM -f "$READS_FILE" "$PAF_FILTERED" > "$ASM_GFA" 2>> "$LOG_FILE"
    else
        echo "   [3/5] 组装图 (Filtered GFA) 已存在，跳过。"
    fi
    
    # --- Step 4: GFA 转 FASTA (Filtered) ---
    if [ ! -f "$ASM_FASTA" ]; then
        echo "   [4/5] 将 GFA 转换为 FASTA..."
        # 提取 GFA 中的 S 行 (Sequence)
        awk '/^S/{print ">"$2"\n"$3}' "$ASM_GFA" > "$ASM_FASTA"
        
        # 简单统计 Contig 数量
        num_contigs=$(grep -c "^>" "$ASM_FASTA")
        echo "   完成 (Filtered)。生成了 $num_contigs 个 Contigs。"
    else
        echo "   [4/5] 组装结果 (Filtered FASTA) 已存在，跳过。"
    fi
    
    # --- Step 5: Miniasm 组装 (Raw / Unfiltered) ---
    if [ ! -f "$ASM_RAW_FASTA" ]; then
        echo "   [5/5] 运行 Miniasm 组装 (Raw/Unfiltered) 作为对照..."
        
        if [ -f "$PAF_RAW" ]; then
            $MINIASM -f "$READS_FILE" "$PAF_RAW" > "$ASM_RAW_GFA" 2>> "$LOG_FILE"
            awk '/^S/{print ">"$2"\n"$3}' "$ASM_RAW_GFA" > "$ASM_RAW_FASTA"
            
            num_raw_contigs=$(grep -c "^>" "$ASM_RAW_FASTA")
            echo "   完成 (Raw)。生成了 $num_raw_contigs 个 Contigs。"
        else
            echo "   Warning: 原始 PAF 文件丢失，无法生成 Raw 对照。"
        fi
    else
        echo "   [5/5] Raw 对照组装结果已存在，跳过。"
    fi
    
    count=$((count + 1))
    echo "   区域 $region_name 处理完毕。"
    echo "----------------------------------------------------------"
    
done

echo "全部完成！共处理了 $count 个区域。"
