# -*- coding: utf-8 -*-
import gzip
import os

# ================= Configuration =================
# Input filenames (renamed by user)
INPUT_HIFI_GZ = "hifi_raw.fastq.gz"
INPUT_ONT_GZ = "ont_raw.fastq.gz"

# Output filenames
OUTPUT_HIFI = "sample_profile_hifi.fastq"
OUTPUT_ONT = "sample_profile_ont.fastq"

# Extract 40,000 lines (approx 10,000 reads)
LINE_COUNT = 40000
# =============================================

def extract_fastq(input_gz, output_fastq, description):
    print(f"\n[{description}] Extracting...")
    print(f"  Input: {input_gz}")
    print(f"  Output: {output_fastq}")
    
    if not os.path.exists(input_gz):
        print(f"  Error: File {input_gz} not found. Please check the filename.")
        return

    try:
        # Use errors='ignore' to handle potentially truncated gzip files safely
        with gzip.open(input_gz, 'rt', encoding='utf-8', errors='ignore') as f_in:
            with open(output_fastq, 'w', encoding='utf-8') as f_out:
                count = 0
                for line in f_in:
                    f_out.write(line)
                    count += 1
                    if count >= LINE_COUNT:
                        break
        print(f"  Success! Extracted {count} lines.")
        
    except Exception as e:
        print(f"  Extraction failed: {e}")
        print("  (If 'CRC check failed' or 'Not a gzipped file' appears, it's expected for incomplete downloads. Check if output file has content.)")

print("================ Start Local Extraction ================")

# 1. Extract HiFi
extract_fastq(INPUT_HIFI_GZ, OUTPUT_HIFI, "PacBio HiFi")

# 2. Extract ONT
extract_fastq(INPUT_ONT_GZ, OUTPUT_ONT, "ONT Ultra-Long")

print("\n================ Task Finished ================")
print(f"Please check {OUTPUT_HIFI} and {OUTPUT_ONT}")
