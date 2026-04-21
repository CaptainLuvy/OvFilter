#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: stats.py
Description: 
    Calculate and compare assembly statistics (N50, NG50, Length, etc.) 
    for OvFilter, Raw Assembly, and Hifiasm.

Usage:
    python stats.py -d <region_directory>

Output:
    - Prints a comparison table to stdout.
    - Saves "stats.txt" in the region directory.
"""

import argparse
import os
import sys

def get_fasta_lengths(fasta_file):
    """
    Parse a FASTA file and return a list of sequence lengths.
    """
    if not os.path.exists(fasta_file):
        return []

    lengths = []
    current_len = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line)
        if current_len > 0:
            lengths.append(current_len)
    return sorted(lengths, reverse=True)

def calculate_n50(lengths, genome_size=None):
    """
    Calculate N50 and L50. 
    If genome_size is provided, calculates NG50 and LG50.
    """
    if not lengths:
        return 0, 0

    total_len = sum(lengths)
    target_size = genome_size if genome_size else total_len
    
    half_size = target_size / 2.0
    current_sum = 0
    n50 = 0
    l50 = 0
    
    for i, l in enumerate(lengths):
        current_sum += l
        if current_sum >= half_size:
            n50 = l
            l50 = i + 1
            break
            
    return n50, l50

def get_reference_length(region_dir):
    """
    Try to find reference file and calculate its length.
    Checks ref.fa and ref_region.fasta.
    """
    candidates = [
        os.path.join(region_dir, "ref.fa"),
        os.path.join(region_dir, "ref_region.fasta")
    ]
    
    for ref_file in candidates:
        if os.path.exists(ref_file) and os.path.getsize(ref_file) > 0:
            lengths = get_fasta_lengths(ref_file)
            return sum(lengths)
    return None

def format_row(name, count, total_len, max_len, n50, ng50):
    return f"{name:<15} | {count:>8} | {total_len:>15,} | {max_len:>12,} | {n50:>12,} | {ng50:>12}"

def main():
    parser = argparse.ArgumentParser(description="Calculate Assembly Statistics (N50, NG50).")
    parser.add_argument("-d", "--region_dir", required=True, help="Region data directory")
    args = parser.parse_args()
    
    region_dir = os.path.abspath(args.region_dir)
    region_name = os.path.basename(region_dir)
    
    # Define files
    files = {
        "OvFilter": os.path.join(region_dir, "assembly.fasta"),
        "NumericOnly": os.path.join(region_dir, "assembly.numeric.fasta"),
        "Raw": os.path.join(region_dir, "assembly.raw.fasta"),
        "Hifiasm": os.path.join(region_dir, "hifiasm", "hifiasm.p_ctg.fa")
    }
    
    # Get Reference Length for NG50
    ref_len = get_reference_length(region_dir)
    ref_str = f"{ref_len:,}" if ref_len else "N/A"
    
    # Header
    output_lines = []
    output_lines.append(f"=== Assembly Statistics: {region_name} ===")
    output_lines.append(f"Reference Length: {ref_str} bp")
    output_lines.append("-" * 95)
    output_lines.append(f"{'Tool':<15} | {'Contigs':>8} | {'Total bp':>15} | {'Max bp':>12} | {'N50':>12} | {'NG50':>12}")
    output_lines.append("-" * 95)
    
    # Process each tool
    for name, path in files.items():
        lengths = get_fasta_lengths(path)
        
        if not lengths:
            output_lines.append(f"{name:<15} | {'MISSING':>8} | {'-':>15} | {'-':>12} | {'-':>12} | {'-':>12}")
            continue
            
        count = len(lengths)
        total_len = sum(lengths)
        max_len = lengths[0]
        
        n50, _ = calculate_n50(lengths) # Standard N50 (based on assembly size)
        
        ng50_val = "-"
        if ref_len:
            ng50, _ = calculate_n50(lengths, genome_size=ref_len)
            ng50_val = f"{ng50:,}"
            
        output_lines.append(format_row(name, count, total_len, max_len, n50, ng50_val))

    output_lines.append("-" * 95)
    
    # Print and Save
    result_text = "\n".join(output_lines)
    print(result_text)
    
    out_file = os.path.join(region_dir, "stats.txt")
    with open(out_file, 'w') as f:
        f.write(result_text)
    print(f"\n[Saved] Statistics saved to: {out_file}")

if __name__ == "__main__":
    main()
