import sys
import csv
import argparse
import os
import math
import gzip
import multiprocessing
import array
from functools import partial, lru_cache

# Global variable for read-only data sharing in worker processes
GLOBAL_R2U_MAP = None

# Try importing fast parsers; prompt if Biopython is missing
try:
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError:
    print("Error: Please install Biopython (pip install biopython)")
    sys.exit(1)

# --- High-performance helper functions ---

# Precompiled translation table for fast reverse complement
TRANS_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def get_canonical_hash(kmer_str):
    """
    Calculate hash (integer) for canonical K-mer.
    Storing hashes instead of strings reduces memory usage by 3-5x and speeds up comparisons.
    """
    # Fast reverse complement
    rev_str = kmer_str.translate(TRANS_TABLE)[::-1]
    # Compare lexicographically, take hash of the smaller one
    if kmer_str < rev_str:
        return hash(kmer_str)
    else:
        return hash(rev_str)

def open_file(filename, mode='rt'):
    """Automatically handle gzip file opening"""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def load_unikmers(unikmer_file_path):
    """
    Load unique k-mers, stored as a Set of Hashes.
    """
    unikmers = set()
    
    try:
        with open_file(unikmer_file_path) as f:
            for line in f:
                line = line.strip()
                if not line: continue
                # Compatible with CSV or plain text
                kmer = line.split(',')[0] if ',' in line else line
                # Store hash
                unikmers.add(get_canonical_hash(kmer))
    except Exception as e:
        print(f"Failed to read Unikmer file: {e}")
        sys.exit(1)
            
    return unikmers

def map_reads_to_unikmers(reads_file_path, unikmers, k):
    """
    Scan Reads using fast parser, return {read_id: (sorted_hits, u_to_idx)}
    """
    read_to_unikmer_map = {}
    
    # 1. Determine file format and parser
    is_fastq = False
    
    # Auto-detect format (based on file content)
    try:
        # Read first char to judge
        with open_file(reads_file_path, 'rt') as f_check:
            first_char = f_check.read(1)
            if first_char == '@':
                is_fastq = True
            elif first_char == '>':
                is_fastq = False
            else:
                # Unable to identify, fallback to suffix check
                if reads_file_path.endswith(".fastq") or reads_file_path.endswith(".fq") or ".fastq" in reads_file_path:
                    is_fastq = True
    except Exception:
         # Read failed, fallback to suffix check
        if reads_file_path.endswith(".fastq") or reads_file_path.endswith(".fq") or ".fastq" in reads_file_path:
            is_fastq = True
            
    # 2. Scan file
    count = 0
    mapped_count = 0
    
    with open_file(reads_file_path) as handle:
        # Use generator for extremely low memory usage
        iterator = FastqGeneralIterator(handle) if is_fastq else SimpleFastaParser(handle)
        
        for title, seq, *_ in iterator:
            # Extract ID (take content before first space, compatible with PAF format)
            read_id = title.split()[0]
            seq_len = len(seq)
            
            if seq_len < k: continue
            
            hits = []
            # Convert to uppercase to avoid hash differences due to case inconsistency
            seq = seq.upper()
            
            # Sliding window scan
            for i in range(seq_len - k + 1):
                kmer = seq[i:i+k]
                k_hash = get_canonical_hash(kmer)
                
                if k_hash in unikmers:
                    hits.append((k_hash, i))
            
            if hits:
                read_to_unikmer_map[read_id] = hits
                mapped_count += 1
            
            count += 1

    # 3. Preprocessing: Compressed storage (Memory Optimization)
    # Convert List of Tuples to array.array, and remove dict index to save memory
    processed_map = {}
    
    # Determine array type code
    # 'q': signed long long (8 bytes) for hash
    # 'i': signed int (4 bytes) for position
    
    for rid, hits in read_to_unikmer_map.items():
        # Sort by position
        sorted_hits = sorted(hits, key=lambda x: x[1])
        
        # Separate hash and pos
        hashes = [h for h, p in sorted_hits]
        positions = [p for h, p in sorted_hits]
        
        # Create compressed arrays
        hash_arr = array.array('q', hashes)
        pos_arr = array.array('i', positions)
        
        # Store only arrays, not dict (generated on demand in worker)
        processed_map[rid] = (hash_arr, pos_arr)
        
    return processed_map

# --- Filtering logic (Logic unchanged, adapted for Hash) ---

def check_overhang(q_len, q_start, q_end, t_len, t_start, t_end, strand, overhang_tolerance):
    """Check Overhang"""
    tol = overhang_tolerance
    if strand == '+':
        if (q_len - q_end <= tol) and (t_start <= tol): return True
        if (q_start <= tol) and (t_len - t_end <= tol): return True
    else: # strand == '-'
        if (q_start <= tol) and (t_start <= tol): return True
        if (q_len - q_end <= tol) and (t_len - t_end <= tol): return True
    return False

def calculate_sd(deltas):
    """Calculate Standard Deviation"""
    n = len(deltas)
    if n <= 1: return 0.0
    mean = sum(deltas) / n
    variance = sum((x - mean) ** 2 for x in deltas) / n
    return math.sqrt(variance)

@lru_cache(maxsize=1024)
def get_cached_index(read_id):
    """
    Cached function used inside Worker.
    Restores u_to_idx dictionary from global GLOBAL_R2U_MAP on demand.
    """
    if read_id not in GLOBAL_R2U_MAP:
        return None, None, None
        
    hash_arr, pos_arr = GLOBAL_R2U_MAP[read_id]
    
    # Fast index reconstruction
    # Note: If duplicate hashes exist, this takes the index of the last occurrence, consistent with previous logic
    u_to_idx = {h: i for i, h in enumerate(hash_arr)}
    
    return hash_arr, pos_arr, u_to_idx

def verify_overlap(q_id, t_id, q_start, q_end, t_start, t_end, strand, r2u_map_ignored, min_shared, tolerance, max_sd):
    """
    Verify overlap.
    Returns (is_valid, reason_string, shared_count, sd_value)
    """
    # Get data (using cache)
    q_hashes, q_pos, q_u_to_idx = get_cached_index(q_id)
    t_hashes, t_pos, t_u_to_idx = get_cached_index(t_id)
    
    if q_hashes is None or t_hashes is None:
        return False, "NoUnikmer", 0, 0.0
    
    # Set operations (based on keys view, very fast)
    shared_hashes = q_u_to_idx.keys() & t_u_to_idx.keys()
    if not shared_hashes:
        return False, "UnikmerFail_NoShared", 0, 0.0
        
    if len(shared_hashes) < min_shared:
        # Return raw shared count even if failed
        return False, "UnikmerFail_MinShared", len(shared_hashes), 0.0

    verified_pairs = set()
    best_fail_reason = "UnikmerFail_Tolerance"
    best_chain_len = 0
    best_chain_sd = 0.0
    
    for seed_hash in shared_hashes:
        q_idx = q_u_to_idx[seed_hash]
        t_idx = t_u_to_idx[seed_hash]
        
        q_pos_seed = q_pos[q_idx]
        t_pos_seed = t_pos[t_idx]
        
        # Region check
        if not (q_start <= q_pos_seed < q_end): continue
        if not (t_start <= t_pos_seed < t_end): continue
        if (q_idx, t_idx) in verified_pairs: continue
            
        current_chain = set()
        current_chain.add((q_idx, t_idx))
        
        # Extend forward
        curr_q, curr_t = q_idx, t_idx
        while True:
            next_q = curr_q + 1
            next_t = curr_t + 1 if strand == '+' else curr_t - 1
            
            if next_q >= len(q_hashes): break
            if strand == '+' and next_t >= len(t_hashes): break
            if strand == '-' and next_t < 0: break
            
            u_q = q_hashes[next_q]
            u_t = t_hashes[next_t]
            
            if u_q != u_t: break # Hash mismatch
            
            delta_q = q_pos[next_q] - q_pos[curr_q]
            delta_t = t_pos[next_t] - t_pos[curr_t]
            
            dist_diff = abs(delta_q - delta_t) if strand == '+' else abs(delta_q - abs(delta_t))
            if dist_diff > tolerance: break

            current_chain.add((next_q, next_t))
            curr_q, curr_t = next_q, next_t
            
        # Extend backward
        curr_q, curr_t = q_idx, t_idx
        while True:
            prev_q = curr_q - 1
            prev_t = curr_t - 1 if strand == '+' else curr_t + 1
            
            if prev_q < 0: break
            if strand == '+' and prev_t < 0: break
            if strand == '-' and prev_t >= len(t_hashes): break
            
            u_q = q_hashes[prev_q]
            u_t = t_hashes[prev_t]
            
            if u_q != u_t: break
            
            delta_q = q_pos[curr_q] - q_pos[prev_q]
            delta_t = t_pos[curr_t] - t_pos[prev_t]
            
            dist_diff = abs(delta_q - delta_t) if strand == '+' else abs(delta_q - abs(delta_t))
            if dist_diff > tolerance: break
                
            current_chain.add((prev_q, prev_t))
            curr_q, curr_t = prev_q, prev_t
            
        verified_pairs.update(current_chain)
        
        # Chain length and SD verification
        chain_len = len(current_chain)
        if chain_len > best_chain_len:
            best_chain_len = chain_len
            
        if chain_len >= min_shared:
            chain_deltas = []
            for (qi, ti) in current_chain:
                pq = q_pos[qi]
                pt = t_pos[ti]
                chain_deltas.append(pq - pt if strand == '+' else pq + pt)
            
            sd = calculate_sd(chain_deltas)
            if sd <= max_sd:
                return True, "Pass", chain_len, sd
            else:
                best_fail_reason = "UnikmerFail_MaxSD"
                best_chain_sd = sd # Store the failed SD
            
    return False, best_fail_reason, best_chain_len, best_chain_sd

# --- Parallel processing helper functions ---

def init_worker(r2u_map):
    """
    Worker process initialization function.
    Receives data from main process and stores in global variable.
    On Windows, this avoids massive overhead of serializing data for each task call.
    """
    global GLOBAL_R2U_MAP
    GLOBAL_R2U_MAP = r2u_map

def process_chunk(chunk_data):
    """
    Worker process task function.
    chunk_data: (rows, args_dict)
    """
    rows, args = chunk_data
    results = []
    
    # Get data from global variable
    r2u_map = GLOBAL_R2U_MAP
    
    # Extract parameters
    keep_contained = args['keep_contained']
    overhang_tolerance = args['overhang_tolerance']
    min_align_len = args['min_align_len']
    min_identity = args['min_identity']
    min_shared = args['min_shared']
    tolerance = args['tolerance']
    max_sd = args['max_sd']
    numeric_only = args.get('numeric_only', False)
    
    for row in rows:
        if not row or len(row) < 12: continue
        
        # Parse row
        try:
            q_name, t_name = row[0], row[5]
            q_len, q_start, q_end = int(row[1]), int(row[2]), int(row[3])
            t_len, t_start, t_end = int(row[6]), int(row[7]), int(row[8])
            # PAF format spec:
            # Col 9: Number of matching bases in the mapping
            # Col 10: Number bases, including gaps, in the mapping (block length)
            matches, block_len = int(row[9]), int(row[10])
            strand = row[4]
        except ValueError:
            continue

        reason = None
        
        # 1. Self alignment
        if q_name == t_name:
            reason = "Self"
            
        # 2. Contained relationship
        elif not keep_contained:
            slack = overhang_tolerance
            is_q_in_t = (t_len > q_len) and (q_end - q_start >= q_len - slack)
            is_t_in_q = (q_len > t_len) and (t_end - t_start >= t_len - slack)
            if is_q_in_t or is_t_in_q:
                reason = "Contained"
        
        # 3. Basic filtering
        if not reason:
            if block_len < min_align_len:
                reason = "Short"
            elif (matches / block_len) < min_identity:
                reason = "LowIdt"
        
        # 4. Overhang
        if not reason:
            if not check_overhang(q_len, q_start, q_end, t_len, t_start, t_end, strand, overhang_tolerance):
                reason = "Overhang"
        
        # Capture numeric pass state BEFORE unikmer verification
        is_numeric_pass = (reason is None)
        
        # 5. Unikmer verification
        
        # Stats placeholders
        stat_shared = 0
        stat_sd = 0.0
        
        if is_numeric_pass and not numeric_only:
            if q_name not in r2u_map or t_name not in r2u_map:
                reason = "NoUnikmer"
            else:
                verify_res, verify_detail, stat_shared, stat_sd = verify_overlap(q_name, t_name, q_start, q_end, t_start, t_end, strand, r2u_map, 
                                    min_shared, tolerance, max_sd)
                if not verify_res:
                    reason = verify_detail # e.g., "UnikmerFail_Shared", "UnikmerFail_SD", "UnikmerFail_Dist"
        
        # Record result (row content, rejection reason, numeric_pass_status, stats)
        # Calculate Identity
        idt = 0.0
        if block_len > 0:
            idt = matches / block_len
            
        results.append((row, reason, is_numeric_pass, stat_shared, stat_sd, idt))
        
    return results

def process_paf(paf_file, output_file, r2u_map, args):
    """Process PAF file (Supports parallel stream processing, reduces memory usage)"""
    global GLOBAL_R2U_MAP
    
    # Prepare argument dictionary
    args_dict = {
        'keep_contained': args.keep_contained,
        'overhang_tolerance': args.overhang_tolerance,
        'min_align_len': args.min_align_len,
        'min_identity': args.min_identity,
        'min_shared': args.min_shared,
        'tolerance': args.tolerance,
        'max_sd': args.max_sd,
        'numeric_only': args.numeric_only
    }
    
    # Prepare output files
    f_out = open_file(output_file, 'wt') if output_file.endswith('.gz') else open(output_file, 'w', newline='')
    
    # Optional: Numeric Only Output
    f_numeric = None
    numeric_writer = None
    if args.output_numeric:
        f_numeric = open_file(args.output_numeric, 'wt') if args.output_numeric.endswith('.gz') else open(args.output_numeric, 'w', newline='')
        numeric_writer = csv.writer(f_numeric, delimiter='\t')
        
    writer = csv.writer(f_out, delimiter='\t')
    
    # Define batch reader generator
    def batch_reader(file_path, size, params):
        with open_file(file_path, 'rt') as f:
            reader = csv.reader(f, delimiter='\t')
            batch = []
            for r in reader:
                batch.append(r)
                if len(batch) >= size:
                    yield (batch, params)
                    batch = []
            if batch:
                yield (batch, params)

    # Result writing helper function
    def write_results(chunk_res):
        for r_row, r_reason, is_numeric_pass, stat_shared, stat_sd, idt in chunk_res:
            
            # Write to numeric-only file if applicable
            if is_numeric_pass and numeric_writer:
                numeric_writer.writerow(r_row)

            if not r_reason:
                writer.writerow(r_row)

    # Parallel or serial processing
    if args.threads > 1:
        # Set global variable
        GLOBAL_R2U_MAP = r2u_map
        
        pool = None
        is_fork_safe = False
        
        # Environment check
        try:
            start_method = multiprocessing.get_start_method(allow_none=True)
            if start_method == 'fork' or (start_method is None and sys.platform != 'win32'):
                is_fork_safe = True
        except:
            if sys.platform != 'win32':
                is_fork_safe = True
                
        if is_fork_safe:
            pool = multiprocessing.Pool(processes=args.threads)
        else:
            pool = multiprocessing.Pool(processes=args.threads, initializer=init_worker, initargs=(r2u_map,))
            
        # Use imap for stream processing
        chunk_size = 50000  # Fixed batch size to avoid memory backlog
        batch_gen = batch_reader(paf_file, chunk_size, args_dict)
        
        try:
            with pool:
                # imap yields results in order as soon as they are ready
                # Use a larger chunksize for imap to reduce IPC overhead
                for chunk_results in pool.imap(process_chunk, batch_gen, chunksize=1):
                    write_results(chunk_results)
        except Exception as e:
            raise e
            
    else:
        GLOBAL_R2U_MAP = r2u_map
        chunk_size = 5000
        batch_gen = batch_reader(paf_file, chunk_size, args_dict)
        for args_tuple in batch_gen:
            res = process_chunk(args_tuple)
            write_results(res)

    f_out.close()
    if f_numeric: f_numeric.close()

# --- Main program entry ---
def main():
    parser = argparse.ArgumentParser(description="Filter overlap.paf files using unique k-mers (Unikmers) [Optimized]")
    parser.add_argument("-p", "--paf", required=True, help="Input PAF file (supports .gz)")
    parser.add_argument("-u", "--unikmers", required=True, help="Unique k-mer file")
    parser.add_argument("-r", "--reads", required=True, help="Reads file (FASTA/FASTQ, supports .gz)")
    parser.add_argument("-o", "--output", required=True, help="Output PAF file (supports .gz)")
    parser.add_argument("-T", "--threads", type=int, default=1, help="Number of threads for parallel processing (Default: 1)")
    
    parser.add_argument("-k", type=int, default=21, help="K-mer length (21)")
    parser.add_argument("-s", "--min-shared", type=int, default=15, help="Minimum shared count (15)")
    parser.add_argument("-t", "--tolerance", type=int, default=15, help="Position tolerance (15)")
    parser.add_argument("--max-sd", type=float, default=10.0, help="Max SD (10.0)")
    parser.add_argument("--keep-contained", action="store_true", help="Keep contained alignments")
    
    parser.add_argument("-l", "--min-align-len", type=int, default=500, help="Minimum alignment length (500)")
    parser.add_argument("-i", "--min-identity", type=float, default=0.8, help="Minimum identity (0.0)")
    parser.add_argument("-e", "--overhang-tolerance", type=int, default=5000, help="Overhang tolerance (1000)")
    parser.add_argument("--numeric-only", action="store_true", help="Skip unikmer verification, use only numeric filters")
    parser.add_argument("--output-numeric", help="[Optional] Output file for overlaps passing only numeric filters")
    
    args = parser.parse_args()
    
    if args.numeric_only:
        r2u_map = {}
    else:
        unikmers = load_unikmers(args.unikmers)
        if not unikmers: return
        r2u_map = map_reads_to_unikmers(args.reads, unikmers, args.k)
        
    process_paf(args.paf, args.output, r2u_map, args)

if __name__ == "__main__":
    main()
