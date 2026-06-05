#!/bin/bash
set -e

# ==============================================================================
# OvFilter Automated Test & Evaluation Pipeline
# ==============================================================================
# This script runs the complete pipeline: Minimap2 -> Jellyfish -> OvFilter.
# It also supports optional comparison with RAfilter if installed.
# Finally, it evaluates precision/recall against a ground-truth MAF file.
# ==============================================================================

# --- Configuration ---
# Assuming the user extracted sim_data.tar.gz into test/
DATA_DIR="sim_data/chr1/hifi"
REF_FA="${DATA_DIR}/chr1_hifi_0001.ref"
READS_FQ="${DATA_DIR}/chr1_hifi_0001.fastq"
TRUTH_MAF="${DATA_DIR}/chr1_hifi_0001.maf"

# Note: The test dataset is not included in the GitHub repository due to its large size.
# You can download the test dataset from Zenodo: https://doi.org/10.5281/zenodo.20554435

OVFILTER_BIN="../ovfilter"
EVAL_PY="./evaluate_alignments.py"

# Set this to the RAfilter executable path to enable comparison (e.g., "/path/to/rafilter")
RAFILTER_BIN=""

THREADS=8
KMER_SIZE=21
READ_TYPE="hifi" # "hifi" or "ont"

# --- Intermediate & Output Files ---
OUT_DIR="results"
mkdir -p "${OUT_DIR}"

MINIMAP2_PAF="${OUT_DIR}/minimap2_raw.paf"
JF_COUNT="${OUT_DIR}/mer_counts.jf"
TARGET_KMER="${OUT_DIR}/target.kmer"
COMBINED_FA="${OUT_DIR}/combined.fa"

OVFILTER_PAF="${OUT_DIR}/ovfilter.paf"
OVFILTER_TIME="${OUT_DIR}/ovfilter_time.log"

RAFILTER_OUT_DIR="${OUT_DIR}/rafilter_out"
RAFILTER_BUILD_TIME="${OUT_DIR}/rafilter_build.log"
RAFILTER_TIME="${OUT_DIR}/rafilter_filter.log"

# --- 1. Check Input Files and Dependencies ---
if [ ! -f "${REF_FA}" ] || [ ! -f "${READS_FQ}" ]; then
    echo "Error: Missing ${REF_FA} or ${READS_FQ}."
    echo "Please download and extract the dataset into test/sim_data/ first."
    exit 1
fi

if ! command -v minimap2 &> /dev/null; then
    echo "Error: minimap2 could not be found. Please install it or add it to your PATH."
    exit 1
fi

if ! command -v jellyfish &> /dev/null; then
    echo "Error: jellyfish could not be found. Please install it or add it to your PATH."
    exit 1
fi

if [ ! -f "${OVFILTER_BIN}" ]; then
    echo "Error: ${OVFILTER_BIN} not found. Please compile OvFilter first (run 'make' in the parent directory)."
    exit 1
fi

if [ ! -f "${EVAL_PY}" ]; then
    echo "Warning: Evaluation script ${EVAL_PY} not found. Evaluation will be skipped."
fi

echo "=========================================="
echo " Starting OvFilter Evaluation Pipeline"
echo "=========================================="

# --- 2. Run Minimap2 ---
echo "[1/5] Running Minimap2..."
minimap2 -t ${THREADS} -x map-${READ_TYPE} "${REF_FA}" "${READS_FQ}" > "${MINIMAP2_PAF}"

# --- 3. Run Jellyfish ---
echo "[2/5] Extracting unique k-mers using Jellyfish..."
jellyfish count -m ${KMER_SIZE} -s 10M -t ${THREADS} -C "${REF_FA}" -o "${JF_COUNT}"
if [ "${READ_TYPE}" = "ont" ]; then
    jellyfish dump -c -U 3 "${JF_COUNT}" > "${TARGET_KMER}"
else
    jellyfish dump -c -U 1 "${JF_COUNT}" > "${TARGET_KMER}"
fi

# --- 4. Run RAfilter (Optional) ---
USE_RAFILTER=false
if [ -n "${RAFILTER_BIN}" ] && [ -x "${RAFILTER_BIN}" ]; then
    echo "[3/5] Running RAfilter for comparison..."
    USE_RAFILTER=true
    mkdir -p "${RAFILTER_OUT_DIR}"
    
    # RAfilter requires single-line fasta
    QUERY_FA="${OUT_DIR}/query.fa"
    REF_SINGLE_FA="${OUT_DIR}/ref_single.fa"
    sed -n '1~4s/^@/>/p;2~4p' "${READS_FQ}" > "${QUERY_FA}"
    awk '/^>/{print n $1; n = "\n"} !/^>/{printf "%s",$0}' "${REF_FA}" | sed '/^$/d' > "${REF_SINGLE_FA}"

    /usr/bin/time -v "${RAFILTER_BIN}" build -t ${THREADS} -q "${QUERY_FA}" -r "${REF_SINGLE_FA}" -o "${RAFILTER_OUT_DIR}/" "${TARGET_KMER}" 2> "${RAFILTER_BUILD_TIME}"
    /usr/bin/time -v "${RAFILTER_BIN}" filter -o "${RAFILTER_OUT_DIR}/" -p --threshold 12 "${RAFILTER_OUT_DIR}/ref.pos" "${RAFILTER_OUT_DIR}/query.pos" "${MINIMAP2_PAF}" 2> "${RAFILTER_TIME}"
else
    echo "[3/5] Skipping RAfilter (RAFILTER_BIN not set or not executable)."
fi

# --- 5. Run OvFilter ---
echo "[4/5] Running OvFilter..."
cat "${READS_FQ}" "${REF_FA}" > "${COMBINED_FA}"
/usr/bin/time -v "${OVFILTER_BIN}" --read-type ${READ_TYPE} -p "${MINIMAP2_PAF}" -u "${TARGET_KMER}" -r "${COMBINED_FA}" -o "${OVFILTER_PAF}" -T ${THREADS} -k ${KMER_SIZE} 2> "${OVFILTER_TIME}"

# --- 6. Evaluation ---
echo "[5/5] Evaluating results..."
if [ -f "${TRUTH_MAF}" ] && [ -f "${EVAL_PY}" ]; then
    EVAL_CMD=(python "${EVAL_PY}" 
        --maf "${TRUTH_MAF}" 
        --raw-paf "${MINIMAP2_PAF}" 
        --ovfilter-paf "${OVFILTER_PAF}" 
        --ovfilter-label "OvFilter" 
        --ovfilter-time-log "${OVFILTER_TIME}" 
        --read-type "${READ_TYPE}" 
        --output-report "${OUT_DIR}/evaluation_report.txt")

    if [ "${USE_RAFILTER}" = true ]; then
        RAFILTER_FINAL_PAF=$(ls "${RAFILTER_OUT_DIR}"/*rafiltered*.paf | head -n 1)
        EVAL_CMD+=(--rafilter-paf "${RAFILTER_FINAL_PAF}" 
                   --rafilter-build-time-log "${RAFILTER_BUILD_TIME}" 
                   --time-log "${RAFILTER_TIME}")
    fi

    "${EVAL_CMD[@]}"
    
    echo "=========================================="
    echo " Evaluation finished. Report saved to: ${OUT_DIR}/evaluation_report.txt"
    cat "${OUT_DIR}/evaluation_report.txt"
else
    echo "Warning: Evaluation skipped. ${TRUTH_MAF} or ${EVAL_PY} not found."
fi
