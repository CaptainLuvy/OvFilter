#!/bin/bash
set -eu -o pipefail

# ==============================================================================
# 脚本名称: extract_centromere_reads.sh
# 功能描述: 
#   复刻 mm2-ivh 流程，专门用于提取 CHM13 着丝粒 (HOR) 区域的 HiFi Reads 及其对应的 Unikmer 集合。
#   该脚本是数据准备阶段的核心工具，为后续的组装实验提供标准化的输入数据。
#
# 核心流程:
#   1. 数据准备: 自动下载 CHM13v2.0 参考基因组及 CenSat (着丝粒/卫星序列) 注释文件。
#   2. 目标区域定义: 
#      - 筛选注释中的 HOR (高阶重复序列) 区域。
#      - 合并间距 < 5Mb 的相邻区域。
#      - 对合并后的区域两侧各外扩 500kb，确保包含侧翼序列。
#   3. 全基因组比对: 使用 Minimap2 (map-pb) 将 HiFi Reads 比对到参考基因组。
#   4. Reads 提取: 遍历每个目标区域，提取落在该区域内的 Reads，保存为 FASTQ。
#   5. Unikmer 提取: 针对提取出的 Reads，计算 k-mer 频率分布，利用 Python 脚本 (freq_distribution_kmers.py)
#      自动确定频率阈值 (mean ± std)，并提取高置信度的 Unikmers。
#
# 输入:
#   - HiFi Reads 文件 (支持通配符)
#   - (可选) CHM13 参考基因组 (若无则自动下载)
#
# 输出:
#   - centromere_data/regions/<区域名>/reads.fq
#   - centromere_data/regions/<区域名>/unikmer/k21_unikmers.kmers
#   - 中间文件: mapping.paf, targets.bed 等
# ==============================================================================

# --- 1. 用户参数配置 (请修改这里) ---
# 自动下载官方参考基因组以确保染色体名称匹配 (chr1 vs chr1)
REF_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
REF_FILE="chm13v2.0.fa"

# 指定要提取的染色体 (空表示所有染色体)
# 格式: 竖线分隔的正则模式，例如: "chr1|chr8|chr21"
# 如果留空 "" 则提取所有染色体
TARGET_CHROMS="chr1|chr8|chr21"

# 输入文件路径 (绝对路径或相对路径)
READS_FILE="/homeb/dingyc/data/biodata/hg002/hifi/66/*.fastq.gz"      # 您的 HiFi Reads 文件

# 输出目录
OUT_DIR="hg002"

# 工具路径 (如果不在 PATH 中，请修改为绝对路径)
MINIMAP2="minimap2"
BEDTOOLS="bedtools"
SAMTOOLS="samtools"
SEQKIT="seqkit"
WGET="wget"
JELLYFISH="jellyfish"
PYTHON="python3"

# 脚本所在目录 (用于定位 python 脚本)
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
FREQ_DIST_PY="$SCRIPT_DIR/freq_distribution_kmers.py"


# 线程数
THREADS=32

# 原作者流程中的参数
JOIN_DIST=5000000    # 合并距离: 5Mb
MARGIN_SIZE=500000   # 外扩距离: 500kb
BED_PATTERN="^hor_"  # 筛选模式: HOR (Higher Order Repeats)

# ==============================================================================

# 0. 环境检查
echo "[Check] Checking dependencies..."
for tool in "$MINIMAP2" "$BEDTOOLS" "$SAMTOOLS" "$SEQKIT" "$WGET"; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: Tool '$tool' not found in PATH."
        exit 1
    fi
done

mkdir -p "$OUT_DIR"
# 转为绝对路径方便后续处理
OUT_DIR=$(cd "$OUT_DIR" && pwd)

echo "=========================================================="
echo "   Start extracting Centromere (HOR) reads"
echo "=========================================================="

# 0.1 展开 Reads 文件路径 (处理通配符)
# 使用数组来存储展开后的文件名，兼容包含空格等情况
READS_FILES_ARRAY=()
# 这里的 trick 是利用 shell 的 glob expansion
if compgen -G "$READS_FILE" > /dev/null; then
    for f in $READS_FILE; do
        READS_FILES_ARRAY+=("$f")
    done
else
    echo "Error: No files matching pattern '$READS_FILE' found."
    exit 1
fi

echo "Found ${#READS_FILES_ARRAY[@]} read files:"
printf "  - %s\n" "${READS_FILES_ARRAY[@]}"

# 1. Check/Download Reference Genome
if [ ! -f "$REF_FILE" ]; then
    echo "Reference file '$REF_FILE' not found."
    REF_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
    echo "Downloading standard T2T-CHM13v2.0 reference to match annotation..."
    if [ ! -f "chm13v2.0.fa.gz" ]; then
        $WGET -q "$REF_URL" -O "chm13v2.0.fa.gz"
    fi
    echo "Decompressing reference..."
    gzip -d "chm13v2.0.fa.gz"
    REF_FILE="$(pwd)/chm13v2.0.fa"
fi

if [ ! -f "$REF_FILE" ]; then echo "Error: Reference file '$REF_FILE' not found."; exit 1; fi

# --- Step 1: 准备目标区域 (Target Regions) ---
echo "[Step 1] Preparing target regions (HOR)..."

ANNOT_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed"
ANNOT_FILE="$OUT_DIR/chm13_censat.bed"
TARGETS_BED="$OUT_DIR/targets.bed"

# 1.1 下载注释
if [ ! -f "$ANNOT_FILE" ]; then
    echo "  -> Downloading annotation..."
    $WGET -q "$ANNOT_URL" -O "$ANNOT_FILE"
else
    echo "  -> Annotation file exists."
fi

# 1.2 生成参考基因组索引 (如果不存在)
REF_IDX="${REF_FILE}.fai"
if [ ! -f "$REF_IDX" ]; then
    echo "  -> Indexing reference genome..."
    $SAMTOOLS faidx "$REF_FILE"
fi

# 1.3 执行筛选、合并、外扩 (完全复刻 run_local_assembly.sh 逻辑)
# 逻辑: Filter (^hor_) -> Merge (5Mb) -> Slop (500kb)
echo "  -> Generating target BED (Filter -> Merge -> Slop)..."

# 临时文件
TMP_FLT="$OUT_DIR/annot.flt.bed"
TMP_MRG="$OUT_DIR/annot.merged.bed"

# Filter
if [ -n "$TARGET_CHROMS" ]; then
    echo "  -> Filtering for chromosomes: $TARGET_CHROMS"
    # 使用 awk 匹配染色体名称
    awk -v p="$BED_PATTERN" -v c="$TARGET_CHROMS" '$1 ~ "^("c")$" && $4 ~ p' "$ANNOT_FILE" > "$TMP_FLT"
else
    echo "  -> Extracting all chromosomes..."
    awk -v p="$BED_PATTERN" '$4 ~ p' "$ANNOT_FILE" > "$TMP_FLT"
fi

# Merge (按染色体合并，距离<5Mb的合并)
# -c 4 -o distinct 是为了保留合并后的区域名称列表
$BEDTOOLS merge -i "$TMP_FLT" -d "$JOIN_DIST" -c 4 -o distinct > "$TMP_MRG"

# Slop (外扩)
# -g 指定基因组大小文件(.fai)以防越界
$BEDTOOLS slop -i "$TMP_MRG" -g "$REF_IDX" -b "$MARGIN_SIZE" > "$TARGETS_BED"

REGION_COUNT=$(wc -l < "$TARGETS_BED")
echo "  -> Generated $REGION_COUNT target regions in: $TARGETS_BED"


# --- Step 2: 全基因组比对 (Mapping) ---
echo "[Step 2] Mapping reads to reference..."
# 注意：原脚本针对 ONT 用 map-ont，这里改为 map-pb
# 输出 PAF 格式
PAF_FILE="$OUT_DIR/mapping.paf"

if [ ! -f "$PAF_FILE" ]; then
    echo "  -> Running minimap2 (map-pb)..."
    # Pass the array of files
    $MINIMAP2 -t "$THREADS" -x map-pb "$REF_FILE" "${READS_FILES_ARRAY[@]}" > "$PAF_FILE"
else
    echo "  -> PAF file exists, skipping mapping."
fi


# --- Step 3: 提取 Reads (Extraction) ---
echo "[Step 3] Extracting reads for each region..."

# 创建输出子目录
DATA_DIR="$OUT_DIR/regions"
mkdir -p "$DATA_DIR"

# 遍历 BED 文件
# BED 格式: chr start end names
while read -r chrom start end names; do
    # 构造一个安全的文件名 (将逗号等特殊字符替换为下划线)
    region_name=$(echo "$names" | sed 's/,/_/g' | cut -c 1-50) # 截断一下避免文件名过长
    # 或者简单点，用 coordinates 命名
    safe_name="${chrom}_${start}_${end}"
    
    echo "  -> Processing: $chrom:$start-$end (HORs: $region_name)"
    
    REGION_SUBDIR="$DATA_DIR/$safe_name"
    mkdir -p "$REGION_SUBDIR"
    
    READ_LIST="$REGION_SUBDIR/read_ids.txt"
    READS_OUT="$REGION_SUBDIR/reads.fq"
    
    # 筛选逻辑 (原脚本逻辑简化版):
    # 找出比对到该染色体，且比对区间与目标区间有交集的 reads
    # 原脚本逻辑：$6==chr && $8 < end && $9 >= start (Reads mapped within the bin)
    # 这里我们沿用这个逻辑，确保 reads 主要位于该区域内
    
    awk -v c="$chrom" -v s="$start" -v e="$end" \
        '$6 == c && $8 <= e && $9 >= s { print $1 }' "$PAF_FILE" | sort | uniq > "$READ_LIST"
    
    count=$(wc -l < "$READ_LIST")
    echo "     Found $count reads."
    
    if [ "$count" -gt 0 ]; then
        # 使用 seqkit 提取
        # 如果 reads 文件很大，这一步可能会花点时间
        $SEQKIT grep -f "$READ_LIST" "${READS_FILES_ARRAY[@]}" > "$READS_OUT"
        
        # 顺便生成该区域的参考序列 (用于后续组装)
        $SAMTOOLS faidx "$REF_FILE" "${chrom}:${start}-${end}" > "$REGION_SUBDIR/ref.fa"
        
        echo "     Saved to: $REGION_SUBDIR/reads.fq"
        
        # --- Step 3.1: 提取 Unique K-mers ---
        # 仿照 extractUnikmers.sh
        echo "     -> Extracting unique k-mers..."
        
        KMER_DIR="$REGION_SUBDIR/unikmer"
        mkdir -p "$KMER_DIR"
        
        K=21
        OUT_COUNT="$KMER_DIR/k${K}_count.jf"
        OUT_HISTO="$KMER_DIR/k${K}_counts.histo"
        OUT_XCEL="$KMER_DIR/k${K}_distribution.xlsx"
        UNIKMERS="$KMER_DIR/k${K}_unikmers.kmers"
        
        # 1. Count k-mers
        # -s 10M 对于单个区域应该够了，如果不够可以调大
        $JELLYFISH count -m $K -C -o "$OUT_COUNT" -c 3 -s 10000000 --disk -t "$THREADS" "$READS_OUT"
        
        # 2. Generate Histogram
        $JELLYFISH histo -o "$OUT_HISTO" -v "$OUT_COUNT"
        
        # 3. Calculate Bounds (using python script)
        # 注意: 这里的 python 脚本需要 xlsxwriter 和 numpy
        if OUTPUT=$($PYTHON "$FREQ_DIST_PY" "$OUT_HISTO" "$OUT_XCEL"); then
             # extract the lower and upper bound values from the output
             # Python output format:
             # mean: ... std: ...
             # lower: ...
             # upper: ...
             
             lower=$(echo "$OUTPUT" | grep "lower:" | cut -d':' -f2 | xargs)
             upper=$(echo "$OUTPUT" | grep "upper:" | cut -d':' -f2 | xargs)
             
             echo "        Bounds: [$lower, $upper]"
             
             # 4. Dump Unikmers
             $JELLYFISH dump -c -t -L "$lower" -U "$upper" "$OUT_COUNT" | awk '{print $1}' > "$UNIKMERS"
             
             echo "        Unikmers saved to: $UNIKMERS"
        else
             echo "        Warning: Python script failed. Skipping unikmer extraction."
        fi
        
    else
        echo "     No reads found, skipping."
        rm -rf "$REGION_SUBDIR"
    fi

done < "$TARGETS_BED"

echo "=========================================================="
echo "   All Done!"
echo "   Extracted data is in: $DATA_DIR"
echo "=========================================================="
