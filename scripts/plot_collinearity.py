#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: plot_collinearity.py
Description: 
    Plot collinearity dotplot between assembly contigs and reference genome for a specific region.
    
    Key Features:
    1. Auto Alignment: Uses minimap2 to align assembly.fasta to ref.fa.
    2. Multi-Comparison: Plots "Filtered Assembly", "Raw Assembly", and "Hifiasm Assembly" side-by-side.
    3. Smart Plotting:
       - Automatically extracts or uses existing reference sequence.
       - Parses PAF alignment format.
       - Intelligently sorts Contigs along Y-axis for better visualization.
       - Blue lines = Forward alignment, Red lines = Reverse complement.

Usage:
    python plot_collinearity.py -d hg002/regions/chr1_123_456

Input Requirements:
    - region_dir (-d): 
        - Must contain ref.fa (Reference).
        - Optional: assembly.fasta (OvFilter result).
        - Optional: assembly.raw.fasta (Raw result).
        - Optional: hifiasm/hifiasm.p_ctg.fa (Hifiasm result).

Output:
    - Generates collinearity.png in the region_dir.
"""

import argparse
import os
import subprocess
import sys
import matplotlib
matplotlib.use('Agg') # Backend setting, no GUI required
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --- Helper Functions ---

def check_dependency(cmd):
    """Check if command exists in PATH"""
    try:
        subprocess.check_call(f"which {cmd}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except:
        print(f"Error: Dependency '{cmd}' not found. Please ensure it is in your PATH.")
        sys.exit(1)

def run_command(cmd):
    """Execute shell command and handle errors"""
    print(f"[Exec] {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed: {e}")
        sys.exit(1)

def extract_reference(ref_file, region_str, output_file):
    """
    Extract specific region sequence from full genome FASTA.
    """
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        print(f"[Skip] Reference file exists: {output_file}")
        return

    parts = region_str.split("_")
    if len(parts) < 3:
        print(f"Error: Invalid region name format '{region_str}'. Expected format: chr_start_end (e.g., chr1_1000_2000)")
        sys.exit(1)
        
    start = parts[-2]
    end = parts[-1]
    chrom = "_".join(parts[:-2])
    
    region_query = f"{chrom}:{start}-{end}"
    print(f"Extracting region from reference: {region_query}")
    
    cmd = f"samtools faidx \"{ref_file}\" \"{region_query}\" > \"{output_file}\""
    run_command(cmd)

def run_minimap2(target, query, output_paf, threads=4):
    """
    Run minimap2 for sequence alignment.
    """
    if os.path.exists(output_paf) and os.path.getsize(output_paf) > 0:
        print(f"[Skip] Alignment file exists: {output_paf}")
        return

    # Use asm5 mode for genome-to-genome alignment
    cmd = f"minimap2 -x asm5 -t {threads} \"{target}\" \"{query}\" > \"{output_paf}\" 2> /dev/null"
    run_command(cmd)

def parse_paf(paf_file):
    """
    Parse PAF format alignment file.
    """
    alignments = []
    if not os.path.exists(paf_file):
        return alignments
        
    with open(paf_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12: continue
            
            q_name = parts[0]
            q_len = int(parts[1])
            q_start = int(parts[2])
            q_end = int(parts[3])
            strand = parts[4]
            t_name = parts[5]
            t_len = int(parts[6])
            t_start = int(parts[7])
            t_end = int(parts[8])
            
            matches = int(parts[9])
            block_len = int(parts[10])
            
            # Filter: Ignore fragments shorter than 1kb
            if block_len < 1000: continue
            
            alignments.append({
                'q_name': q_name, 'q_start': q_start, 'q_end': q_end, 'q_len': q_len,
                't_name': t_name, 't_start': t_start, 't_end': t_end, 't_len': t_len,
                'strand': strand, 'matches': matches, 'block_len': block_len
            })
    return alignments

def plot_dotplot(alignments, title, ax, ref_len):
    """
    Plot collinearity dotplot.
    """
    if not alignments:
        ax.text(0.5, 0.5, "No Alignments Found", ha='center', va='center', transform=ax.transAxes, fontsize=12, color='red')
        ax.set_title(title)
        ax.set_xlim(0, ref_len)
        return

    # --- 1. Smart Contig Layout (Y-axis) ---
    contig_center_map = {}
    contig_len_map = {}
    
    for aln in alignments:
        q = aln['q_name']
        t_center = (aln['t_start'] + aln['t_end']) / 2
        
        if q not in contig_center_map:
            contig_center_map[q] = []
            contig_len_map[q] = aln['q_len']
            
        contig_center_map[q].append(t_center)
        
    sorted_contigs = sorted(contig_center_map.keys(), 
                          key=lambda k: sum(contig_center_map[k])/len(contig_center_map[k]))
    
    y_offsets = {}
    current_y = 0
    padding = ref_len * 0.01 
    
    for q in sorted_contigs:
        y_offsets[q] = current_y
        current_y += contig_len_map[q] + padding
        
    total_y_height = max(current_y, 100) # Ensure non-zero height
    
    # --- 2. Plotting ---
    for aln in alignments:
        q = aln['q_name']
        y_base = y_offsets[q]
        
        x_start = aln['t_start']
        x_end = aln['t_end']
        
        if aln['strand'] == '+':
            y_start = y_base + aln['q_start']
            y_end = y_base + aln['q_end']
            color = '#1f77b4' # Blue
        else:
            y_start = y_base + aln['q_end']
            y_end = y_base + aln['q_start']
            color = '#d62728' # Red
            
        ax.plot([x_start, x_end], [y_start, y_end], color=color, alpha=0.8, linewidth=1.5)
        
    # Contig boundary lines
    for q in sorted_contigs:
        y_line = y_offsets[q] + contig_len_map[q] + (padding/2)
        ax.axhline(y=y_line, color='gray', linestyle=':', alpha=0.3, linewidth=0.5)

    # Axis settings
    ax.set_title(title, fontsize=12, pad=10, fontweight='bold')
    ax.set_xlabel("Reference Position (bp)")
    ax.set_ylabel("Assembly Contigs (Stacked bp)")
    ax.set_xlim(0, ref_len)
    ax.set_ylim(0, total_y_height)
    ax.grid(True, linestyle='--', alpha=0.2)
    
    # Statistics Box
    num_contigs = len(sorted_contigs)
    total_len = sum(contig_len_map.values())
    stats_text = f"Contigs: {num_contigs}\nTotal Len: {total_len:,} bp"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

def main():
    parser = argparse.ArgumentParser(description="Plot collinearity between assembly and reference.")
    parser.add_argument("-d", "--region_dir", required=True, help="Region data directory")
    parser.add_argument("-r", "--reference", help="Full reference genome FASTA (optional)")
    parser.add_argument("-o", "--output", help="Output image filename")
    
    args = parser.parse_args()
    
    check_dependency("samtools")
    check_dependency("minimap2")
    
    region_dir = os.path.abspath(args.region_dir)
    region_name = os.path.basename(region_dir)
    
    if not args.output:
        args.output = os.path.join(region_dir, "collinearity.png")
    
    # Define File Paths
    filtered_asm = os.path.join(region_dir, "assembly.fasta")
    raw_asm = os.path.join(region_dir, "assembly.raw.fasta")
    hifiasm_asm = os.path.join(region_dir, "hifiasm", "hifiasm.p_ctg.fa")
    
    # Temporary files
    ref_region_fasta = os.path.join(region_dir, "ref_region.fasta")
    paf_filtered = os.path.join(region_dir, "aln_filtered.paf")
    paf_raw = os.path.join(region_dir, "aln_raw.paf")
    paf_hifiasm = os.path.join(region_dir, "aln_hifiasm.paf")
    
    print(f"=== Processing Region: {region_name} ===")
    
    # --- 1. Prepare Reference ---
    local_ref = os.path.join(region_dir, "ref.fa")
    final_ref_fasta = ""

    if os.path.exists(local_ref) and os.path.getsize(local_ref) > 0:
        print(f"Found local reference: {local_ref}")
        final_ref_fasta = local_ref
    elif args.reference:
        if not os.path.exists(ref_region_fasta) or os.path.getsize(ref_region_fasta) == 0:
            extract_reference(args.reference, region_name, ref_region_fasta)
        final_ref_fasta = ref_region_fasta
    else:
        print("Error: Reference not found. Please provide --reference or ensure 'ref.fa' exists.")
        sys.exit(1)
        
    ref_len = 0
    with open(final_ref_fasta, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                ref_len += len(line.strip())
    print(f"Reference Length: {ref_len:,} bp")

    # --- 2. Alignments ---
    # Filtered
    if os.path.exists(filtered_asm):
        print("Aligning Filtered Assembly...")
        run_minimap2(final_ref_fasta, filtered_asm, paf_filtered)
    else:
        print(f"Warning: Filtered Assembly not found ({filtered_asm})")
    
    # Raw
    if os.path.exists(raw_asm):
        print("Aligning Raw Assembly...")
        run_minimap2(final_ref_fasta, raw_asm, paf_raw)
    
    # Hifiasm
    if os.path.exists(hifiasm_asm):
        print("Aligning Hifiasm Assembly...")
        run_minimap2(final_ref_fasta, hifiasm_asm, paf_hifiasm)
    else:
        print(f"Warning: Hifiasm Assembly not found ({hifiasm_asm})")

    # --- 3. Plotting ---
    print("Generating Collinearity Plot...")
    
    # Adjust layout to 3 columns
    fig, axes = plt.subplots(1, 3, figsize=(36, 10))
    (ax1, ax2, ax3) = axes
    
    # Plot 1: Filtered (OvFilter)
    alns_filtered = parse_paf(paf_filtered)
    plot_dotplot(alns_filtered, "OvFilter Assembly", ax1, ref_len)
    
    # Plot 2: Raw (No Filter)
    alns_raw = parse_paf(paf_raw)
    plot_dotplot(alns_raw, "Raw Assembly (No Filter)", ax2, ref_len)
    
    # Plot 3: Hifiasm
    alns_hifiasm = parse_paf(paf_hifiasm)
    plot_dotplot(alns_hifiasm, "Hifiasm Assembly", ax3, ref_len)
    
    # Overall Title
    plt.suptitle(f"Assembly Comparison: {region_name}\nReference Length: {ref_len:,} bp", fontsize=16)
    plt.tight_layout()
    
    try:
        plt.savefig(args.output, dpi=300)
        print(f"\n[Success] Plot saved to: {args.output}")
    except Exception as e:
        print(f"Error: Failed to save plot: {e}")

if __name__ == "__main__":
    main()
