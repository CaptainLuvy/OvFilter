#!/bin/bash

# ==============================================================================
# Script: run_analyze_fn.sh
# Purpose: Run the analyze_fn.py error analysis script with the specific file paths
#          from the ovfilter parameter search tests.
# ==============================================================================

# Common Paths (based on test_ovfilter_params.sh)
BASE_DIR="/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp"
TRUTH_MAF="${BASE_DIR}/sim_data/chr8/hifi/chr8_hifi_0001.maf"
RAW_PAF="${BASE_DIR}/result/chr8/hifi30x_secondary/run/minimap2/minimap2_raw.paf"
OUT_DIR="${BASE_DIR}/result/chr8/hifi30x_secondary/param_test_ovfilter"

# We will analyze the result of the best parameter configuration (Ext_Ratio=0.0)
FILTERED_PAF="${OUT_DIR}/ovfilter_sh1_ext0.0.paf"
# If the PAF was deleted by the grid search script, the user must re-run ovfilter to generate it, or provide the stats PAF if it's identical
# Let's check if it exists, if not we will run ovfilter to generate it
if [ ! -f "$FILTERED_PAF" ]; then
    echo "Filtered PAF not found. Regenerating it with optimal parameters (Min_Shared=1, Ext_Ratio=0.0)..."
    OVFILTER_CPP_BIN="/home/dingyc/tools/ovfilter/ovfilter_cpp"
    UNIKMER_FILE="${BASE_DIR}/result/chr8/hifi30x_secondary/run/tmp/ovfilter.kmer"
    COMBINED_FA="${BASE_DIR}/result/chr8/hifi30x_secondary/run/tmp/combined_for_ovfilter.fa"
    
    ${OVFILTER_CPP_BIN} \
        --read-type hifi \
        --paf "${RAW_PAF}" \
        --unikmers "${UNIKMER_FILE}" \
        --reads "${COMBINED_FA}" \
        --output "${FILTERED_PAF}" \
        --threads 32 \
        -k 21 \
        -s 1 \
        -t 15 \
        -x 0.0 \
        -l 0 \
        -e 5000 \
        --keep-contained --ignore-overhang
fi

# Output prefix for the detailed TSV files
OUT_PREFIX="${OUT_DIR}/error_analysis_ext0.0"

echo "========================================="
echo "Running Error Analysis for ovfilter..."
echo "Truth MAF:    $TRUTH_MAF"
echo "Raw PAF:      $RAW_PAF"
echo "Filtered PAF: $FILTERED_PAF"
echo "Output Prefix:$OUT_PREFIX"
echo "========================================="

# Run the python script
python /home/dingyc/tools/ovfilter/analyze_fn.py \
    --maf "$TRUTH_MAF" \
    --raw-paf "$RAW_PAF" \
    --filtered-paf "$FILTERED_PAF" \
    --out-prefix "$OUT_PREFIX"

echo ""
echo "Analysis complete! Check the output directory for the generated TSV files:"
echo "  ${OUT_PREFIX}_FN.tsv  (False Negatives - True alignments we missed)"
echo "  ${OUT_PREFIX}_FP.tsv  (False Positives - Bad alignments we kept)"
echo "  ${OUT_PREFIX}_TP_sample.tsv (Sample of True Positives)"
