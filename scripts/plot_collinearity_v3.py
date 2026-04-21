#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: plot_collinearity_v3.py
Description: 
    Plot collinearity dotplot between assembly contigs and reference genome.
    
    v3 Update:
    - Auto-Flip: Automatically detects if a contig is predominantly Reverse (-) mapped.
      If so, it visually "flips" the contig on the Y-axis so it appears diagonal (positive slope),
      making visual comparison easier.
    - Preserves Color: Flipped contigs are still colored Red to indicate they are reverse complements.
    - Layout: 3 Rows (Tools) x 2 Columns (Strands).
    
    Key Features:
    1. Auto Alignment: Uses minimap2.
    2. Smart Plotting:
       - Auto-flips reverse contigs for consistent diagonal visualization.
       - Separates Forward/Reverse strands into subplots.
       - Skips redundant steps if files exist.

Usage:
    python plot_collinearity_v3.py -d hg002/regions/chr1_123_456

Output:
    - Generates collinearity_v3.png
"""

import argparse
import os
import subprocess
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# --- Helper Functions ---

def run_command(cmd):
    print(f"[Exec] {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed: {e}")
        sys.exit(1)

def extract_reference(ref_file, region_str, output_file):
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        print(f"[Skip] Reference file exists: {output_file}")
        return

    parts = region_str.split("_")
    if len(parts) < 3:
        print(f"Error: Invalid region name format '{region_str}'")
        sys.exit(1)
        
    start = parts[-2]
    end = parts[-1]
    chrom = "_".join(parts[:-2])
    
    region_query = f"{chrom}:{start}-{end}"
    print(f"Extracting region from reference: {region_query}")
    
    cmd = f"samtools faidx \"{ref_file}\" \"{region_query}\" > \"{output_file}\""
    run_command(cmd)

def run_minimap2(target, query, output_paf, threads=4):
    if os.path.exists(output_paf) and os.path.getsize(output_paf) > 0:
        print(f"[Skip] Alignment file exists: {output_paf}")
        return

    cmd = f"minimap2 -x asm5 -t {threads} \"{target}\" \"{query}\" > \"{output_paf}\" 2> /dev/null"
    run_command(cmd)

def parse_paf(paf_file):
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
            
            if block_len < 1000: continue
            
            alignments.append({
                'q_name': q_name, 'q_start': q_start, 'q_end': q_end, 'q_len': q_len,
                't_name': t_name, 't_start': t_start, 't_end': t_end, 't_len': t_len,
                'strand': strand, 'matches': matches, 'block_len': block_len
            })
    return alignments

def calculate_layout(alignments, ref_len):
    """
    Calculate Y-axis layout and detect which contigs need flipping.
    """
    if not alignments:
        return [], {}, 100, {}, set()

    contig_center_map = {}
    contig_len_map = {}
    contig_strand_score = {} # +1 for forward bases, -1 for reverse bases
    
    for aln in alignments:
        q = aln['q_name']
        t_center = (aln['t_start'] + aln['t_end']) / 2
        
        if q not in contig_center_map:
            contig_center_map[q] = []
            contig_len_map[q] = aln['q_len']
            contig_strand_score[q] = 0
            
        contig_center_map[q].append(t_center)
        
        # Calculate strand score based on match length
        score = aln['block_len']
        if aln['strand'] == '-':
            score = -score
        contig_strand_score[q] += score
        
    # Identify contigs to flip (mostly reverse aligned)
    flipped_contigs = set()
    for q, score in contig_strand_score.items():
        if score < 0:
            flipped_contigs.add(q)

    # Sort contigs by their average mapping position on reference
    # For flipped contigs, the t_center logic still holds (center of alignment on ref)
    sorted_contigs = sorted(contig_center_map.keys(), 
                          key=lambda k: sum(contig_center_map[k])/len(contig_center_map[k]))
    
    y_offsets = {}
    current_y = 0
    padding = ref_len * 0.01 
    
    for q in sorted_contigs:
        y_offsets[q] = current_y
        current_y += contig_len_map[q] + padding
        
    total_y_height = max(current_y, 100)
    
    return sorted_contigs, y_offsets, total_y_height, contig_len_map, flipped_contigs

def plot_single_panel(ax, alignments, title, ref_len, layout_info, target_strand):
    """
    Plot a single dotplot panel for a specific strand.
    """
    sorted_contigs, y_offsets, total_y_height, contig_len_map, flipped_contigs = layout_info
    
    strand_alns = [a for a in alignments if a['strand'] == target_strand]
    
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.set_xlim(0, ref_len)
    ax.set_ylim(0, total_y_height)
    ax.grid(True, linestyle='--', alpha=0.2)
    
    if not strand_alns:
        ax.text(0.5, 0.5, "No Alignments", ha='center', va='center', transform=ax.transAxes, fontsize=10, color='gray')
    
    for aln in strand_alns:
        q = aln['q_name']
        if q not in y_offsets: continue 
        
        y_base = y_offsets[q]
        q_len = contig_len_map[q]
        is_flipped = q in flipped_contigs
        
        # Original coordinates
        q_start = aln['q_start']
        q_end = aln['q_end']
        
        # Apply Flip if needed
        if is_flipped:
            # If contig is flipped, we treat q_len as 0 and 0 as q_len visually
            # Transform: new_q = q_len - old_q
            vis_q_start = q_len - q_start
            vis_q_end = q_len - q_end
        else:
            vis_q_start = q_start
            vis_q_end = q_end
            
        x_start = aln['t_start']
        x_end = aln['t_end']
        
        # Calculate Y coordinates based on strand AND flip status
        # Logic:
        # Normal +: Slope > 0. y_start < y_end.
        # Normal -: Slope < 0. y_start > y_end.
        # Flipped +: Effectively becomes -, Slope < 0.
        # Flipped -: Effectively becomes +, Slope > 0.
        
        if aln['strand'] == '+':
            # Base logic: y follows q
            y_s = y_base + vis_q_start
            y_e = y_base + vis_q_end
            color = '#1f77b4' # Blue
        else:
            # Base logic: y is inverted relative to q (start maps to end)
            # Standard: y_start = y_base + q_end, y_end = y_base + q_start
            # With vis_q (which might be flipped):
            y_s = y_base + vis_q_end
            y_e = y_base + vis_q_start
            color = '#d62728' # Red
            
        ax.plot([x_start, x_end], [y_s, y_e], color=color, alpha=0.8, linewidth=1.5)
        
    # Plot contig boundaries only if not too many contigs (avoid gray wash for reads)
    if len(sorted_contigs) < 200:
        padding = ref_len * 0.01
        for q in sorted_contigs:
            y_line = y_offsets[q] + contig_len_map[q] + (padding/2)
            ax.axhline(y=y_line, color='gray', linestyle=':', alpha=0.3, linewidth=0.5)
        
        # Mark flipped contigs
        if q in flipped_contigs:
            # Add a small marker or text? Maybe too cluttered.
            pass

    ax.tick_params(axis='both', which='major', labelsize=8)

def main():
    parser = argparse.ArgumentParser(description="Plot collinearity v3 (Auto-Flip).")
    parser.add_argument("-d", "--region_dir", required=True, help="Region data directory")
    parser.add_argument("-r", "--reference", help="Full reference genome FASTA (optional)")
    parser.add_argument("-o", "--output", help="Output image filename")
    
    args = parser.parse_args()
    
    region_dir = os.path.abspath(args.region_dir)
    region_name = os.path.basename(region_dir)
    
    if not args.output:
        args.output = os.path.join(region_dir, "collinearity_v3.png")
    
    # Define File Paths
    filtered_asm = os.path.join(region_dir, "assembly.fasta")
    raw_asm = os.path.join(region_dir, "assembly.raw.fasta")
    hifiasm_asm = os.path.join(region_dir, "hifiasm", "hifiasm.p_ctg.fa")
    
    # Temporary files
    ref_region_fasta = os.path.join(region_dir, "ref_region.fasta")
    paf_filtered = os.path.join(region_dir, "aln_filtered.paf")
    paf_raw = os.path.join(region_dir, "aln_raw.paf")
    paf_hifiasm = os.path.join(region_dir, "aln_hifiasm.paf")
    
    print(f"=== Processing Region (v3): {region_name} ===")
    
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
    if os.path.exists(filtered_asm): run_minimap2(final_ref_fasta, filtered_asm, paf_filtered)
    if os.path.exists(raw_asm): run_minimap2(final_ref_fasta, raw_asm, paf_raw)
    if os.path.exists(hifiasm_asm): run_minimap2(final_ref_fasta, hifiasm_asm, paf_hifiasm)

    # --- 3. Plotting ---
    print("Generating Collinearity Plot v3 (Auto-Flipped)...")
    
    fig, axes = plt.subplots(3, 2, figsize=(12, 18), facecolor='white')
    
    # Row 1: OvFilter
    alns_filtered = parse_paf(paf_filtered)
    layout_filtered = calculate_layout(alns_filtered, ref_len)
    plot_single_panel(axes[0, 0], alns_filtered, "OvFilter (Forward/+)", ref_len, layout_filtered, '+')
    plot_single_panel(axes[0, 1], alns_filtered, "OvFilter (Reverse/-)", ref_len, layout_filtered, '-')
    
    # Row 2: Raw
    alns_raw = parse_paf(paf_raw)
    layout_raw = calculate_layout(alns_raw, ref_len)
    plot_single_panel(axes[1, 0], alns_raw, "Raw Assembly (Forward/+)", ref_len, layout_raw, '+')
    plot_single_panel(axes[1, 1], alns_raw, "Raw Assembly (Reverse/-)", ref_len, layout_raw, '-')
    
    # Row 3: Hifiasm
    alns_hifiasm = parse_paf(paf_hifiasm)
    layout_hifiasm = calculate_layout(alns_hifiasm, ref_len)
    plot_single_panel(axes[2, 0], alns_hifiasm, "Hifiasm (Forward/+)", ref_len, layout_hifiasm, '+')
    plot_single_panel(axes[2, 1], alns_hifiasm, "Hifiasm (Reverse/-)", ref_len, layout_hifiasm, '-')
    
    # Global Labels
    fig.text(0.5, 0.01, 'Reference Position (bp)', ha='center', fontsize=12)
    fig.text(0.01, 0.5, 'Assembly Contigs (Stacked bp)', va='center', rotation='vertical', fontsize=12)
    fig.text(0.5, 0.95, f"Assembly Comparison v3 (Auto-Flipped): {region_name}\nReference Length: {ref_len:,} bp", ha='center', fontsize=16)
    
    plt.tight_layout(rect=[0.02, 0.02, 1, 0.93])
    
    try:
        plt.savefig(args.output, dpi=300)
        print(f"\n[Success] Plot saved to: {args.output}")
    except Exception as e:
        print(f"Error: Failed to save plot: {e}")

if __name__ == "__main__":
    main()
