import argparse
import os

def parse_maf(maf_path):
    truth = {}
    with open(maf_path, 'r', encoding='utf-8') as f:
        block_ref = None
        block_read = None

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('a'):
                block_ref = None
                block_read = None
            elif line.startswith('s'):
                parts = line.split()
                if not block_ref:
                    block_ref = parts
                elif not block_read:
                    block_read = parts

                    ref_name = block_ref[1]
                    ref_start = int(block_ref[2])
                    ref_len = int(block_ref[3])
                    ref_strand = block_ref[4]
                    
                    ref_end = ref_start + ref_len

                    read_name = block_read[1]
                    clean_read_name = read_name.split('/')[0]
                    q_strand = block_read[4]

                    relative_strand = '+' if ref_strand == q_strand else '-'
                    truth.setdefault(clean_read_name, []).append({
                        'ref_name': ref_name,
                        'ref_start': ref_start,
                        'ref_end': ref_end,
                        'strand': relative_strand,
                    })
    return truth

def is_true_alignment(pred_aln, truth_alignments, overlap_ratio=0.99):
    for truth_aln in truth_alignments:
        name_match = (truth_aln['ref_name'] == pred_aln['ref_name']) or (truth_aln['ref_name'] == 'ref')
        if not name_match or truth_aln['strand'] != pred_aln['strand']:
            continue
        
        overlap_start = max(pred_aln['ref_start'], truth_aln['ref_start'])
        overlap_end = min(pred_aln['ref_end'], truth_aln['ref_end'])
        overlap_len = max(0, overlap_end - overlap_start)
        
        aln_len = pred_aln['ref_end'] - pred_aln['ref_start']
        if aln_len > 0 and (overlap_len / aln_len) >= overlap_ratio:
            return True
    return False

def extract_tags(parts):
    tags = {}
    for p in parts[12:]:
        if p.startswith('rl:i:'):
            tags['rl'] = int(p.split(':')[2])
        elif p.startswith('NM:i:'):
            tags['NM'] = int(p.split(':')[2])
        elif p.startswith('de:f:'):
            tags['de'] = float(p.split(':')[2])
    return tags

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--maf', required=True)
    parser.add_argument('--raw-paf', required=True)
    parser.add_argument('--filtered-paf', required=True)
    parser.add_argument('--overlap-ratio', type=float, default=0.99)
    parser.add_argument('--out-prefix', default='error_analysis', help='Prefix for output TSV files')
    args = parser.parse_args()

    print("Parsing MAF truth...")
    truth = parse_maf(args.maf)

    print("Parsing filtered PAF to get survived alignments...")
    survived_ids = set()
    with open(args.filtered_paf, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12: continue
            read_name = parts[0].split('/')[0]
            ref_name = parts[5]
            ref_start = parts[7]
            ref_end = parts[8]
            # Create a unique signature for the alignment
            sig = f"{read_name}_{ref_name}_{ref_start}_{ref_end}"
            survived_ids.add(sig)

    print("Analyzing Raw PAF...")
    tp_stats = {'rl': [], 'mapq': [], 'len': [], 'de': [], 'id': [], 'max_oh': []}
    fn_stats = {'rl': [], 'mapq': [], 'len': [], 'de': [], 'id': [], 'max_oh': []}
    fp_stats = {'rl': [], 'mapq': [], 'len': [], 'de': [], 'id': [], 'max_oh': []}

    fn_records = []
    fp_records = []
    tp_records = []

    with open(args.raw_paf, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12: continue
            
            read_name = parts[0].split('/')[0]
            read_len = int(parts[1])
            read_start = int(parts[2])
            read_end = int(parts[3])
            strand = parts[4]
            ref_name = parts[5]
            ref_len = int(parts[6])
            ref_start = int(parts[7])
            ref_end = int(parts[8])
            matches = int(parts[9])
            block_len = int(parts[10])
            mapq = int(parts[11])
            aln_len = ref_end - ref_start
            
            # Identity and Overhangs calculation
            identity = matches / block_len if block_len > 0 else 0
            q_left = read_start
            q_right = read_len - read_end
            t_left = ref_start if strand == '+' else ref_len - ref_end
            t_right = ref_len - ref_end if strand == '+' else ref_start
            
            true_left_overhang = min(q_left, t_left)
            true_right_overhang = min(q_right, t_right)
            max_oh = max(true_left_overhang, true_right_overhang)
            
            sig = f"{read_name}_{ref_name}_{ref_start}_{ref_end}"
            is_survived = sig in survived_ids
            
            pred_aln = {
                'ref_name': ref_name,
                'ref_start': ref_start,
                'ref_end': ref_end,
                'strand': strand
            }
            
            truth_alns = truth.get(read_name, [])
            is_true = is_true_alignment(pred_aln, truth_alns, args.overlap_ratio)
            
            tags = extract_tags(parts)
            rl = tags.get('rl', 0)
            de = tags.get('de', 0.0)
            
            record = f"{read_name}\t{read_len}\t{read_start}\t{read_end}\t{strand}\t{ref_name}\t{ref_len}\t{ref_start}\t{ref_end}\t{aln_len}\t{identity:.4f}\t{max_oh}\t{mapq}\t{rl}\t{de}"
            
            if is_true and is_survived:
                # True Positive
                tp_stats['rl'].append(rl)
                tp_stats['mapq'].append(mapq)
                tp_stats['len'].append(aln_len)
                tp_stats['de'].append(de)
                tp_stats['id'].append(identity)
                tp_stats['max_oh'].append(max_oh)
                if len(tp_records) < 1000: # limit to avoid huge files
                    tp_records.append(record)
            elif is_true and not is_survived:
                # False Negative (Missed true alignments)
                fn_stats['rl'].append(rl)
                fn_stats['mapq'].append(mapq)
                fn_stats['len'].append(aln_len)
                fn_stats['de'].append(de)
                fn_stats['id'].append(identity)
                fn_stats['max_oh'].append(max_oh)
                fn_records.append(record)
            elif not is_true and is_survived:
                # False Positive (Bad alignments kept)
                fp_stats['rl'].append(rl)
                fp_stats['mapq'].append(mapq)
                fp_stats['len'].append(aln_len)
                fp_stats['de'].append(de)
                fp_stats['id'].append(identity)
                fp_stats['max_oh'].append(max_oh)
                fp_records.append(record)

    def print_stats(name, stats_dict):
        count = len(stats_dict['rl'])
        print(f"--- {name} (Count: {count}) ---")
        if count == 0:
            return
        
        rl_arr = stats_dict['rl']
        len_arr = stats_dict['len']
        mapq_arr = stats_dict['mapq']
        de_arr = stats_dict['de']
        id_arr = stats_dict['id']
        max_oh_arr = stats_dict['max_oh']
        
        rl_ratio = [r / max(1, l) for r, l in zip(rl_arr, len_arr)]
        
        avg_mapq = sum(mapq_arr) / count
        avg_len = sum(len_arr) / count
        avg_de = sum(de_arr) / count
        avg_rl = sum(rl_arr) / count
        avg_id = sum(id_arr) / count
        avg_max_oh = sum(max_oh_arr) / count
        avg_rl_ratio = sum(rl_ratio) / count
        
        pct_rl_gt_0 = sum(1 for r in rl_arr if r > 0) / count * 100
        pct_rl_ratio_gt_50 = sum(1 for rr in rl_ratio if rr > 0.5) / count * 100
        pct_rl_ratio_gt_90 = sum(1 for rr in rl_ratio if rr > 0.9) / count * 100
        
        pct_id_gt_99 = sum(1 for i in id_arr if i > 0.99) / count * 100
        pct_id_gt_995 = sum(1 for i in id_arr if i > 0.995) / count * 100
        
        print(f"  Avg MAPQ: {avg_mapq:.1f}")
        print(f"  Avg Align Len: {avg_len:.1f} bp")
        print(f"  Avg Identity: {avg_id:.4f} (>99.0%: {pct_id_gt_99:.1f}%, >99.5%: {pct_id_gt_995:.1f}%)")
        print(f"  Avg Max Overhang: {avg_max_oh:.1f} bp")
        print(f"  Avg Divergence (de): {avg_de:.4f}")
        print(f"  Avg RL (Repeat Len): {avg_rl:.1f} bp")
        print(f"  Avg RL Ratio (RL / Align Len): {avg_rl_ratio:.3f}")
        print(f"  % Alignments with RL > 0: {pct_rl_gt_0:.1f}%")
        print(f"  % Alignments with RL Ratio > 50%: {pct_rl_ratio_gt_50:.1f}%")
        print(f"  % Alignments with RL Ratio > 90%: {pct_rl_ratio_gt_90:.1f}%")
        print("")

    print("\n=== Analysis Results ===")
    print_stats("True Positives (TP - Correctly Kept)", tp_stats)
    print_stats("False Negatives (FN - Wrongly Discarded)", fn_stats)
    print_stats("False Positives (FP - Wrongly Kept)", fp_stats)

    def write_records(filename, records):
        with open(filename, 'w') as f:
            f.write("read_name\tread_len\tread_start\tread_end\tstrand\tref_name\tref_len\tref_start\tref_end\taln_len\tidentity\tmax_oh\tmapq\trl\tde\n")
            for r in records:
                f.write(r + '\n')

    print(f"\nWriting detailed records to {args.out_prefix}_*.tsv ...")
    write_records(f"{args.out_prefix}_FN.tsv", fn_records)
    write_records(f"{args.out_prefix}_FP.tsv", fp_records)
    write_records(f"{args.out_prefix}_TP_sample.tsv", tp_records)
    print("Done.")

if __name__ == '__main__':
    main()