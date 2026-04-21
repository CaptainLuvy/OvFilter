import sys
from analyze_fn import parse_maf, is_true_alignment, extract_tags

truth = parse_maf("/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp/sim_data/chr8/hifi/chr8_hifi_0001.maf")

fn_mapq = []
fn_de = []
tn_mapq = []
tn_de = []

with open("/homeb/dingyc/fsa/hifi/chm13/ovfilter/scripts/sim_exp/result/chr8/hifi30x_secondary/run/minimap2/minimap2_raw.paf") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 12: continue
        
        read_name = parts[0].split('/')[0]
        ref_name = parts[5]
        ref_start = int(parts[7])
        ref_end = int(parts[8])
        strand = parts[4]
        mapq = int(parts[11])
        
        pred_aln = {'ref_name': ref_name, 'ref_start': ref_start, 'ref_end': ref_end, 'strand': strand}
        is_true = is_true_alignment(pred_aln, truth.get(read_name, []), 0.99)
        
        tags = extract_tags(parts)
        rl = tags.get('rl', 0)
        de = tags.get('de', 0.0)
        
        # We are interested in highly repetitive regions where rl > 5000
        align_len = ref_end - ref_start
        if rl > 10000 or (rl > 5000 and align_len > 0 and rl / align_len > 0.5):
            if is_true:
                fn_mapq.append(mapq)
                fn_de.append(de)
            else:
                tn_mapq.append(mapq)
                tn_de.append(de)

print(f"True Alignments in repeat regions (FN candidates): {len(fn_mapq)}")
print(f"  Avg MAPQ: {sum(fn_mapq)/max(1, len(fn_mapq)):.2f}, MAPQ>0: {sum(1 for m in fn_mapq if m>0)}, MAPQ>=5: {sum(1 for m in fn_mapq if m>=5)}, MAPQ>=10: {sum(1 for m in fn_mapq if m>=10)}")
print(f"  Avg de: {sum(fn_de)/max(1, len(fn_de)):.4f}")

print(f"False Alignments in repeat regions (TN candidates): {len(tn_mapq)}")
print(f"  Avg MAPQ: {sum(tn_mapq)/max(1, len(tn_mapq)):.2f}, MAPQ>0: {sum(1 for m in tn_mapq if m>0)}, MAPQ>=5: {sum(1 for m in tn_mapq if m>=5)}, MAPQ>=10: {sum(1 for m in tn_mapq if m>=10)}")
print(f"  Avg de: {sum(tn_de)/max(1, len(tn_de)):.4f}")