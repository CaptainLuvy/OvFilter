import argparse
from pathlib import Path


def parse_dataset(spec):
    parts = spec.split('::')
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            '--dataset format must be LABEL::PAF_PATH'
        )
    label, paf_path = parts
    return {
        'label': label,
        'paf_path': Path(paf_path),
    }


def load_paf_records(paf_path):
    records = []
    with paf_path.open('r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 12:
                continue
            try:
                matches = int(cols[9])
                block_len = int(cols[10])
            except ValueError:
                continue
            if block_len <= 0:
                continue
            records.append({
                'matches': matches,
                'block_len': block_len,
                'identity': matches / block_len,
            })
    return records


def quantile(sorted_values, q):
    if not sorted_values:
        return float('nan')
    if q <= 0:
        return sorted_values[0]
    if q >= 1:
        return sorted_values[-1]
    pos = (len(sorted_values) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = pos - lo
    return sorted_values[lo] * (1.0 - frac) + sorted_values[hi] * frac


def summarize_records(records):
    n = len(records)
    identities = sorted(r['identity'] for r in records)
    total_matches = sum(r['matches'] for r in records)
    total_block_len = sum(r['block_len'] for r in records)
    simple_mean = sum(identities) / n if n else float('nan')
    weighted_mean = (total_matches / total_block_len) if total_block_len > 0 else float('nan')
    block_len_mean = (total_block_len / n) if n else float('nan')
    return {
        'n': n,
        'block_len_mean': block_len_mean,
        'simple_mean_identity': simple_mean,
        'weighted_mean_identity': weighted_mean,
        'median_identity': quantile(identities, 0.5),
        'q25_identity': quantile(identities, 0.25),
        'q75_identity': quantile(identities, 0.75),
        'q10_identity': quantile(identities, 0.10),
        'q90_identity': quantile(identities, 0.90),
        'min_identity': identities[0] if identities else float('nan'),
        'max_identity': identities[-1] if identities else float('nan'),
    }


def main():
    parser = argparse.ArgumentParser(
        description='Summarize identity distribution statistics from one or more PAF files.'
    )
    parser.add_argument(
        '--dataset',
        action='append',
        type=parse_dataset,
        required=True,
        help='dataset spec: LABEL::PAF_PATH',
    )
    parser.add_argument('--min-block-len', type=int, default=0, help='ignore alignments shorter than this block length')
    parser.add_argument('--out-tsv', type=Path, default=None, help='optional TSV output path')
    args = parser.parse_args()

    header = [
        'label', 'n', 'block_len_mean',
        'simple_mean_identity', 'weighted_mean_identity',
        'median_identity', 'q10_identity', 'q25_identity', 'q75_identity', 'q90_identity',
        'min_identity', 'max_identity',
    ]
    rows = []

    print('Identity summary from PAF:')
    if args.min_block_len > 0:
        print('  min_block_len filter = {}'.format(args.min_block_len))
    print()

    for ds in args.dataset:
        if not ds['paf_path'].exists():
            raise SystemExit('PAF file not found: {}'.format(ds['paf_path']))
        records = load_paf_records(ds['paf_path'])
        if args.min_block_len > 0:
            records = [r for r in records if r['block_len'] >= args.min_block_len]
        if not records:
            raise SystemExit('No valid PAF records found after filtering: {}'.format(ds['paf_path']))
        stats = summarize_records(records)
        row = [
            ds['label'],
            str(stats['n']),
            '{:.1f}'.format(stats['block_len_mean']),
            '{:.4f}'.format(stats['simple_mean_identity']),
            '{:.4f}'.format(stats['weighted_mean_identity']),
            '{:.4f}'.format(stats['median_identity']),
            '{:.4f}'.format(stats['q10_identity']),
            '{:.4f}'.format(stats['q25_identity']),
            '{:.4f}'.format(stats['q75_identity']),
            '{:.4f}'.format(stats['q90_identity']),
            '{:.4f}'.format(stats['min_identity']),
            '{:.4f}'.format(stats['max_identity']),
        ]
        rows.append(row)

    widths = [max(len(str(cell)) for cell in col) for col in zip(header, *rows)]
    print('  '.join(str(h).ljust(w) for h, w in zip(header, widths)))
    for row in rows:
        print('  '.join(str(cell).ljust(w) for cell, w in zip(row, widths)))

    if args.out_tsv:
        args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.out_tsv.open('w', encoding='utf-8') as out:
            out.write('\t'.join(header) + '\n')
            for row in rows:
                out.write('\t'.join(row) + '\n')
        print('\nSaved TSV to: {}'.format(args.out_tsv))


if __name__ == '__main__':
    main()
