# -*- coding: utf-8 -*-
import argparse
import unicodedata
import os
from collections import defaultdict


def parse_maf(maf_path):
    """解析 PBSIM2 MAF 文件，获取每条 read 的真实比对位置。"""
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
                    ref_src_size = int(block_ref[5])
                    if ref_strand == '-':
                        ref_start = ref_src_size - ref_start - ref_len
                    ref_end = ref_start + ref_len

                    read_name = block_read[1]
                    clean_read_name = read_name.split('/')[0]
                    q_start = int(block_read[2])
                    q_len = int(block_read[3])
                    q_strand = block_read[4]
                    q_src_size = int(block_read[5])
                    if q_strand == '-':
                        q_start_fwd = q_src_size - q_start - q_len
                    else:
                        q_start_fwd = q_start

                    relative_strand = '+' if ref_strand == q_strand else '-'
                    truth.setdefault(clean_read_name, []).append({
                        'ref_name': ref_name,
                        'ref_start': ref_start,
                        'ref_end': ref_end,
                        'q_start': q_start_fwd,
                        'q_total': q_src_size,
                        'strand': relative_strand,
                    })
    return truth


def count_true_alignments(truth):
    return sum(len(t_list) for t_list in truth.values())


def prediction_matches_truth(pred_aln, truth_aln, overlap_ratio):
    name_match = (truth_aln['ref_name'] == pred_aln['ref_name']) or (truth_aln['ref_name'] == 'ref')
    if not name_match or truth_aln['strand'] != pred_aln['strand']:
        return False

    overlap_start = max(pred_aln['ref_start'], truth_aln['ref_start'])
    overlap_end = min(pred_aln['ref_end'], truth_aln['ref_end'])
    overlap_len = max(0, overlap_end - overlap_start)
    if overlap_len <= 0:
        return False

    aln_len = pred_aln['ref_end'] - pred_aln['ref_start']
    return aln_len > 0 and (overlap_len / aln_len) >= overlap_ratio


def max_bipartite_match_count(pred_alignments, truth_alignments, overlap_ratio):
    if not pred_alignments or not truth_alignments:
        return 0

    adjacency = []
    for pred_aln in pred_alignments:
        candidates = []
        for truth_idx, truth_aln in enumerate(truth_alignments):
            if prediction_matches_truth(pred_aln, truth_aln, overlap_ratio):
                candidates.append(truth_idx)
        adjacency.append(candidates)

    truth_to_pred = {}

    def dfs(pred_idx, seen_truth):
        for truth_idx in adjacency[pred_idx]:
            if truth_idx in seen_truth:
                continue
            seen_truth.add(truth_idx)
            if truth_idx not in truth_to_pred or dfs(truth_to_pred[truth_idx], seen_truth):
                truth_to_pred[truth_idx] = pred_idx
                return True
        return False

    matched = 0
    for pred_idx in range(len(pred_alignments)):
        if dfs(pred_idx, set()):
            matched += 1
    return matched


def eval_paf(paf_path, truth, mapq_thresh=0, overlap_ratio=0.95):
    if not paf_path or not os.path.exists(paf_path):
        return None

    total_alignments = 0
    predictions_by_read = defaultdict(list)

    with open(paf_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            mapq = int(parts[11])
            if mapq < mapq_thresh:
                continue

            read_name = parts[0].split('/')[0]
            total_alignments += 1
            predictions_by_read[read_name].append({
                'ref_name': parts[5],
                'ref_start': int(parts[7]),
                'ref_end': int(parts[8]),
                'strand': parts[4],
            })

    total_true_alignments = count_true_alignments(truth)
    matched_alignments = 0
    for read_name, pred_alignments in predictions_by_read.items():
        truth_alignments = truth.get(read_name, [])
        matched_alignments += max_bipartite_match_count(pred_alignments, truth_alignments, overlap_ratio)

    precision = matched_alignments / total_alignments if total_alignments > 0 else 0.0
    recall = matched_alignments / total_true_alignments if total_true_alignments > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        'total_alignments': total_alignments,
        'matched_alignments': matched_alignments,
        'total_true_alignments': total_true_alignments,
        'precision': precision,
        'recall': recall,
        'f1': f1,
    }


def parse_time_log(log_path):
    time_sec = 0.0
    peak_mem_kb = 0
    if not log_path or not os.path.exists(log_path):
        return None
    with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if 'Elapsed (wall clock) time' in line:
                val = line.split('):', 1)[1].strip()
                parts = val.split(':')
                if len(parts) == 3:
                    time_sec = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
                elif len(parts) == 2:
                    time_sec = int(parts[0]) * 60 + float(parts[1])
            elif 'Maximum resident set size (kbytes):' in line:
                peak_mem_kb = int(line.split(':', 1)[1].strip())
    return time_sec, peak_mem_kb


def main():
    parser = argparse.ArgumentParser(description='使用 PBSIM2 MAF 真实数据评估 Minimap2、RAfilter 与 OvFilter 的比对结果')
    parser.add_argument('--maf', required=True, help='PBSIM2 真实位置 MAF 文件')
    parser.add_argument('--raw-paf', required=True, help='原始的 Minimap2 PAF 文件')
    parser.add_argument('--rafilter-paf', required=True, help='RAfilter 过滤后的 PAF 文件')
    parser.add_argument('--ovfilter-paf', help='OvFilter 默认版本过滤后的 PAF 文件（可选）')
    parser.add_argument('--ovfilter-label', default='OvFilter', help='OvFilter 默认版本的展示标签')
    parser.add_argument('--ovfilter-time-log', help="由 '/usr/bin/time -v ovfilter ...' 生成的默认 OvFilter 日志文件")
    parser.add_argument('--ovfilter-extra-paf', help='额外 OvFilter 版本过滤后的 PAF 文件（可选）')
    parser.add_argument('--ovfilter-extra-label', default='OvFilter Extra', help='额外 OvFilter 版本的展示标签')
    parser.add_argument('--ovfilter-extra-time-log', help='额外 OvFilter 版本的时间与内存日志文件')
    parser.add_argument('--rafilter-build-time-log', help="由 '/usr/bin/time -v rafilter build ...' 生成的日志文件")
    parser.add_argument('--time-log', help="由 '/usr/bin/time -v rafilter filter ...' 生成的日志文件")
    parser.add_argument('--read-type', choices=['hifi', 'ont'], default='ont', help='模拟 reads 的类型，决定交集判定阈值（hifi: 0.99, ont: 0.95）')
    parser.add_argument('--output-report', help='[可选] 将评估结果保存到指定的文本文件中')
    args = parser.parse_args()

    overlap_ratio = 0.99 if args.read_type == 'hifi' else 0.95
    output_lines = []

    def log_and_print(text=''):
        print(text)
        output_lines.append(text)

    truth = parse_maf(args.maf)
    total_true_alignments = count_true_alignments(truth)
    log_and_print(f'正在从 {args.maf} 加载 MAF 真实数据 ...')
    log_and_print(f'成功从 MAF 文件中加载了 {len(truth)} 条真实 reads，共 {total_true_alignments} 条真实对齐记录。\n')

    log_and_print(f'评估 Minimap2（原始结果）[判断标准: 参考名称与链方向一致，且交叠占比 >= {overlap_ratio * 100:.1f}%]...')
    raw_res = eval_paf(args.raw_paf, truth, mapq_thresh=0, overlap_ratio=overlap_ratio)
    log_and_print(f'评估 Minimap2（MAPQ >= 40 过滤后）[判断标准: 参考名称与链方向一致，且交叠占比 >= {overlap_ratio * 100:.1f}%]...')
    mapq_res = eval_paf(args.raw_paf, truth, mapq_thresh=40, overlap_ratio=overlap_ratio)
    log_and_print(f'评估 RAfilter（过滤后）[判断标准: 参考名称与链方向一致，且交叠占比 >= {overlap_ratio * 100:.1f}%]...')
    ra_res = eval_paf(args.rafilter_paf, truth, mapq_thresh=0, overlap_ratio=overlap_ratio)

    evaluated_methods = []
    if args.ovfilter_paf and os.path.exists(args.ovfilter_paf):
        log_and_print(f'评估 {args.ovfilter_label}（过滤后）[判断标准: 参考名称与链方向一致，且交叠占比 >= {overlap_ratio * 100:.1f}%]...')
        evaluated_methods.append((args.ovfilter_label, eval_paf(args.ovfilter_paf, truth, mapq_thresh=0, overlap_ratio=overlap_ratio), args.ovfilter_time_log))
    if args.ovfilter_extra_paf and os.path.exists(args.ovfilter_extra_paf):
        log_and_print(f'评估 {args.ovfilter_extra_label}（过滤后）[判断标准: 参考名称与链方向一致，且交叠占比 >= {overlap_ratio * 100:.1f}%]...')
        evaluated_methods.append((args.ovfilter_extra_label, eval_paf(args.ovfilter_extra_paf, truth, mapq_thresh=0, overlap_ratio=overlap_ratio), args.ovfilter_extra_time_log))

    log_and_print('\n说明：Precision / Recall / F1 统一按 alignment-level 的一对一匹配结果统计。')

    def display_width(text):
        width = 0
        for ch in str(text):
            width += 2 if unicodedata.east_asian_width(ch) in ('W', 'F') else 1
        return width

    def pad_cell(text, width):
        text = str(text)
        return text + ' ' * max(0, width - display_width(text))

    col_method = 25
    col_total = 10
    col_matched = 12
    col_true = 10
    col_metric = 10

    header = (
        f"{pad_cell('Method', col_method)} | "
        f"{pad_cell('Total Aln', col_total)} | "
        f"{pad_cell('Matched Aln', col_matched)} | "
        f"{pad_cell('True Aln', col_true)} | "
        f"{pad_cell('Precision', col_metric)} | "
        f"{pad_cell('Recall', col_metric)} | "
        f"{pad_cell('F1 Score', col_metric)}"
    )
    separator = '-' * len(header)

    log_and_print('\n' + '=' * len(header))
    log_and_print(header)
    log_and_print(separator)

    def print_row(name, res):
        if res:
            line = (
                f"{pad_cell(name, col_method)} | "
                f"{pad_cell(res['total_alignments'], col_total)} | "
                f"{pad_cell(res['matched_alignments'], col_matched)} | "
                f"{pad_cell(res['total_true_alignments'], col_true)} | "
                f"{pad_cell('{:.4f}'.format(res['precision']), col_metric)} | "
                f"{pad_cell('{:.4f}'.format(res['recall']), col_metric)} | "
                f"{pad_cell('{:.4f}'.format(res['f1']), col_metric)}"
            )
            log_and_print(line)
        else:
            log_and_print(f"{pad_cell(name, col_method)} | 文件未找到或为空。")

    print_row('Minimap2 (Raw)', raw_res)
    print_row('Minimap2 (MAPQ >= 40)', mapq_res)
    print_row('RAfilter', ra_res)
    for name, res, _ in evaluated_methods:
        print_row(name, res)
    log_and_print('=' * 95)

    if args.time_log or args.rafilter_build_time_log:
        total_ra_time = 0.0
        max_ra_mem = 0
        if args.rafilter_build_time_log and os.path.exists(args.rafilter_build_time_log):
            res_build = parse_time_log(args.rafilter_build_time_log)
            if res_build:
                total_ra_time += res_build[0]
                max_ra_mem = max(max_ra_mem, res_build[1])
        if args.time_log and os.path.exists(args.time_log):
            res_filter = parse_time_log(args.time_log)
            if res_filter:
                total_ra_time += res_filter[0]
                max_ra_mem = max(max_ra_mem, res_filter[1])
        if total_ra_time > 0 or max_ra_mem > 0:
            log_and_print('\n--- RAfilter 性能消耗（包含 Build + Filter）---')
            log_and_print(f'总运行时间 (Wall clock): {total_ra_time:.2f} 秒')
            log_and_print(f'峰值内存 (RSS):    {max_ra_mem / 1024:.2f} MB')

    for name, _, time_log in evaluated_methods:
        res = parse_time_log(time_log)
        if res:
            t, m = res
            log_and_print(f'\n--- {name} 性能消耗 ---')
            log_and_print(f'运行时间 (Wall clock): {t:.2f} 秒')
            log_and_print(f'峰值内存 (RSS):    {m / 1024:.2f} MB')

    if args.output_report:
        with open(args.output_report, 'w', encoding='utf-8') as f:
            f.write('\n'.join(output_lines) + '\n')
        print(f'\n[INFO] 评估结果已保存至: {args.output_report}')


if __name__ == '__main__':
    main()

