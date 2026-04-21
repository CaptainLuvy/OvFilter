#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
from typing import Any, Dict, List, Set, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run_command(cmd: str) -> None:
    print(f"[Exec] {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed: {e}")
        sys.exit(1)


def extract_reference(ref_file: str, region_str: str, output_file: str) -> None:
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


def run_minimap2(target: str, query: str, output_paf: str, threads: int) -> None:
    if os.path.exists(output_paf) and os.path.getsize(output_paf) > 0:
        print(f"[Skip] Alignment file exists: {output_paf}")
        return

    cmd = f"minimap2 -x asm5 -t {threads} \"{target}\" \"{query}\" > \"{output_paf}\" 2> /dev/null"
    run_command(cmd)


def parse_paf(paf_file: str) -> List[Dict[str, Any]]:
    alignments: List[Dict[str, Any]] = []
    if not os.path.exists(paf_file):
        return alignments

    with open(paf_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue

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

            if block_len < 1000:
                continue

            alignments.append(
                {
                    "q_name": q_name,
                    "q_start": q_start,
                    "q_end": q_end,
                    "q_len": q_len,
                    "t_name": t_name,
                    "t_start": t_start,
                    "t_end": t_end,
                    "t_len": t_len,
                    "strand": strand,
                    "matches": matches,
                    "block_len": block_len,
                }
            )
    return alignments


def calculate_layout(
    alignments: List[Dict[str, Any]], ref_len: int
) -> Tuple[List[str], Dict[str, float], float, Dict[str, int], Set[str]]:
    if not alignments:
        return [], {}, 100.0, {}, set()

    contig_center_map: Dict[str, List[float]] = {}
    contig_len_map: Dict[str, int] = {}
    contig_strand_score: Dict[str, int] = {}

    for aln in alignments:
        q = aln["q_name"]
        t_center = (aln["t_start"] + aln["t_end"]) / 2

        if q not in contig_center_map:
            contig_center_map[q] = []
            contig_len_map[q] = aln["q_len"]
            contig_strand_score[q] = 0

        contig_center_map[q].append(t_center)

        score = aln["block_len"]
        if aln["strand"] == "-":
            score = -score
        contig_strand_score[q] += score

    flipped_contigs = {q for q, score in contig_strand_score.items() if score < 0}

    sorted_contigs = sorted(
        contig_center_map.keys(), key=lambda k: sum(contig_center_map[k]) / len(contig_center_map[k])
    )

    y_offsets: Dict[str, float] = {}
    current_y = 0.0
    padding = ref_len * 0.01

    for q in sorted_contigs:
        y_offsets[q] = current_y
        current_y += contig_len_map[q] + padding

    total_y_height = max(current_y, 100.0)
    return sorted_contigs, y_offsets, total_y_height, contig_len_map, flipped_contigs


def plot_single_panel(
    ax,
    alignments: List[Dict[str, Any]],
    title: str,
    ref_len: int,
    layout_info: Tuple[List[str], Dict[str, float], float, Dict[str, int], Set[str]],
    target_strand: str,
) -> None:
    sorted_contigs, y_offsets, total_y_height, contig_len_map, flipped_contigs = layout_info

    strand_alns = [a for a in alignments if a["strand"] == target_strand]

    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlim(0, ref_len)
    ax.set_ylim(0, total_y_height)
    ax.grid(True, linestyle="--", alpha=0.2)

    if not strand_alns:
        ax.text(
            0.5,
            0.5,
            "No Alignments",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=10,
            color="gray",
        )

    for aln in strand_alns:
        q = aln["q_name"]
        if q not in y_offsets:
            continue

        y_base = y_offsets[q]
        q_len = contig_len_map[q]
        is_flipped = q in flipped_contigs

        q_start = aln["q_start"]
        q_end = aln["q_end"]

        if is_flipped:
            vis_q_start = q_len - q_start
            vis_q_end = q_len - q_end
        else:
            vis_q_start = q_start
            vis_q_end = q_end

        x_start = aln["t_start"]
        x_end = aln["t_end"]

        if aln["strand"] == "+":
            y_s = y_base + vis_q_start
            y_e = y_base + vis_q_end
            color = "#1f77b4"
        else:
            y_s = y_base + vis_q_end
            y_e = y_base + vis_q_start
            color = "#d62728"

        ax.plot([x_start, x_end], [y_s, y_e], color=color, alpha=0.8, linewidth=1.5)

    if len(sorted_contigs) < 200:
        padding = ref_len * 0.01
        for q in sorted_contigs:
            y_line = y_offsets[q] + contig_len_map[q] + (padding / 2)
            ax.axhline(y=y_line, color="gray", linestyle=":", alpha=0.3, linewidth=0.5)

    ax.tick_params(axis="both", which="major", labelsize=8)


def read_fasta_length(fasta_path: str) -> int:
    ref_len = 0
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                ref_len += len(line.strip())
    return ref_len


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot collinearity (Raw vs OvFilter), auto-flip reverse contigs.")
    parser.add_argument("-d", "--region_dir", required=True, help="Region data directory")
    parser.add_argument("-r", "--reference", help="Full reference genome FASTA (optional)")
    parser.add_argument("-o", "--output", help="Output image filename")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Threads for minimap2")

    args = parser.parse_args()

    region_dir = os.path.abspath(args.region_dir)
    region_name = os.path.basename(region_dir)

    if not args.output:
        args.output = os.path.join(region_dir, "collinearity_v4.png")

    filtered_asm = os.path.join(region_dir, "assembly.fasta")
    raw_asm = os.path.join(region_dir, "assembly.raw.fasta")

    ref_region_fasta = os.path.join(region_dir, "ref_region.fasta")
    paf_raw = os.path.join(region_dir, "aln_raw.paf")
    paf_filtered = os.path.join(region_dir, "aln_filtered.paf")

    print(f"=== Processing Region (v4): {region_name} ===")

    local_ref = os.path.join(region_dir, "ref.fa")
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

    ref_len = read_fasta_length(final_ref_fasta)
    print(f"Reference Length: {ref_len:,} bp")

    if os.path.exists(raw_asm):
        run_minimap2(final_ref_fasta, raw_asm, paf_raw, threads=args.threads)
    else:
        print(f"Warning: Raw assembly not found: {raw_asm}")

    if os.path.exists(filtered_asm):
        run_minimap2(final_ref_fasta, filtered_asm, paf_filtered, threads=args.threads)
    else:
        print(f"Warning: OvFilter assembly not found: {filtered_asm}")

    print("Generating Collinearity Plot v4 (Raw row 1, OvFilter row 2)...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 12), facecolor="white")

    alns_raw = parse_paf(paf_raw)
    layout_raw = calculate_layout(alns_raw, ref_len)
    plot_single_panel(axes[0, 0], alns_raw, "Raw Assembly (Forward/+)", ref_len, layout_raw, "+")
    plot_single_panel(axes[0, 1], alns_raw, "Raw Assembly (Reverse/-)", ref_len, layout_raw, "-")

    alns_filtered = parse_paf(paf_filtered)
    layout_filtered = calculate_layout(alns_filtered, ref_len)
    plot_single_panel(axes[1, 0], alns_filtered, "OvFilter (Forward/+)", ref_len, layout_filtered, "+")
    plot_single_panel(axes[1, 1], alns_filtered, "OvFilter (Reverse/-)", ref_len, layout_filtered, "-")

    fig.text(0.5, 0.01, "Reference Position (bp)", ha="center", fontsize=12)
    fig.text(0.01, 0.5, "Assembly Contigs (Stacked bp)", va="center", rotation="vertical", fontsize=12)
    fig.text(
        0.5,
        0.98,
        f"Assembly Comparison v4 (Auto-Flipped): {region_name}\nReference Length: {ref_len:,} bp",
        ha="center",
        fontsize=14,
    )

    plt.tight_layout(rect=[0.02, 0.02, 1, 0.94])

    try:
        plt.savefig(args.output, dpi=300)
        print(f"\n[Success] Plot saved to: {args.output}")
    except Exception as e:
        print(f"Error: Failed to save plot: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

