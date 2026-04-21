#!/bin/bash
set -eu -o pipefail

# ==============================================================================
# 脚本名称: run_hifiasm_assembly.sh
# 功能描述: 
#   批量或针对特定区域运行 Hifiasm 组装。
#
# 用法:
#   1. 针对特定区域 (推荐):
#      ./run_hifiasm_assembly.sh /path/to/regions/chr8_xxxx
#
#   2. 批量运行整个目录 (智能识别):
#      ./run_hifiasm_assembly.sh /path/to/regions_root
#
#   3. 批量运行默认目录 (需修改脚本中的 INPUT_ROOT):
#      ./run_hifiasm_assembly.sh [filter_pattern]
#
# 输出:
#   - 在目标目录下生成 hifiasm/ 子目录
# ==============================================================================

# --- 默认配置 ---
INPUT_ROOT="hg002/regions" 
THREADS=32
HIFIASM="hifiasm"
AWK="awk"

# ==============================================================================

# 获取参数
TARGET="${1:-}"

# 定义一个函数来处理单个区域
process_region() {
    local region_dir="$1"
    local region_name=$(basename "$region_dir")
    local reads_file="$region_dir/reads.fq"
    
    # 简单的检查，如果不是目录则跳过
    if [ ! -d "$region_dir" ]; then return; fi
    
    # 检查是否包含 reads.fq，如果不包含说明不是有效的区域目录
    if [ ! -f "$reads_file" ]; then
        # 只有在 verbose 模式或者非批量遍历时才报错，避免刷屏
        # echo "   [Skip] reads.fq not found in $region_name"
        return
    fi

    echo "-> Processing Region: $region_name"
    
    # 创建输出目录
    local asm_dir="$region_dir/hifiasm"
    mkdir -p "$asm_dir"
    
    # 检查结果是否存在 (断点续传核心逻辑)
    local final_fasta="$asm_dir/hifiasm.p_ctg.fa"
    if [ -f "$final_fasta" ] && [ -s "$final_fasta" ]; then
        echo "   [Skip] Assembly already exists: $final_fasta"
        return
    fi
    
    # 运行 Hifiasm
    local prefix="$asm_dir/hifiasm"
    echo "   Running hifiasm..."
    $HIFIASM -o "$prefix" -t "$THREADS" -f0 "$reads_file" > "$asm_dir/hifiasm.log" 2>&1
    
    # GFA 转 FASTA
    local gfa_file=""
    if [ -f "${prefix}.bp.p_ctg.gfa" ]; then
        gfa_file="${prefix}.bp.p_ctg.gfa"
    elif [ -f "${prefix}.p_ctg.gfa" ]; then
        gfa_file="${prefix}.p_ctg.gfa"
    else
        echo "   [Error] Hifiasm GFA output not found!"
        return
    fi
    
    echo "   Converting GFA to FASTA..."
    $AWK '/^S/{print ">"$2;print $3}' "$gfa_file" > "$final_fasta"
    
    local num_contigs=$(grep -c "^>" "$final_fasta" || true)
    echo "   [Done] Assembled $num_contigs contigs."
    echo "   Output: $final_fasta"
}

# ==============================================================================

# 检查依赖
if ! command -v "$HIFIASM" &> /dev/null; then
    echo "Error: Tool '$HIFIASM' not found."
    exit 1
fi

echo "=========================================================="
echo "   Start Hifiasm Assembly Task"
echo "=========================================================="

if [ -n "$TARGET" ] && [ -d "$TARGET" ]; then
    # 传入了一个存在的目录路径
    
    # 判断该目录是 "单个区域" 还是 "包含多个区域的父目录"
    # 判断依据: 是否存在 reads.fq
    
    if [ -f "$TARGET/reads.fq" ]; then
        # --- 情况 A: 单个区域 ---
        echo "Mode: Single Region (Path provided)"
        process_region "$TARGET"
    else
        # --- 情况 B: 父目录 (批量) ---
        echo "Mode: Batch Processing (Parent Directory provided)"
        echo "Input Root: $TARGET"
        
        count=0
        for region_dir in "$TARGET"/*; do
            [ -d "$region_dir" ] || continue
            process_region "$region_dir"
            count=$((count+1))
        done
        
        if [ "$count" -eq 0 ]; then
             echo "Warning: No subdirectories found in $TARGET"
        fi
    fi

else
    # --- 情况 C: 默认配置 / 关键词筛选 ---
    echo "Mode: Batch Processing (Default Config)"
    echo "Input Root: $INPUT_ROOT"
    
    if [ ! -d "$INPUT_ROOT" ]; then
        echo "Error: Input root '$INPUT_ROOT' does not exist."
        echo "Tip: Provide full path as argument: ./run_hifiasm_assembly.sh /path/to/regions"
        exit 1
    fi

    if [ -n "$TARGET" ]; then
        echo "Filter: Only regions matching '*$TARGET*'"
    else
        echo "Filter: ALL regions"
    fi

    for region_dir in "$INPUT_ROOT"/*; do
        [ -d "$region_dir" ] || continue
        region_name=$(basename "$region_dir")
        
        # 筛选
        if [ -n "$TARGET" ] && [[ "$region_name" != *"$TARGET"* ]]; then
            continue
        fi
        
        process_region "$region_dir"
    done
fi

echo "=========================================================="
echo "   All Tasks Completed!"
echo "=========================================================="
