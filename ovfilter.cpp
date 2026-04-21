#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include <future>
#include <cstdint>
#include <chrono>
#include <cstring>
#include <utility>
#include <zlib.h>

#include <iomanip>

using namespace std;

// Filter stats structure
struct FilterStats {
    long long self = 0;
    long long contained = 0;
    long long short_align = 0;
    long long low_identity = 0;
    long long overhang = 0;
    long long no_unikmer = 0;
    long long u_no_shared = 0;
    long long u_min_shared = 0;
    long long u_no_valid_seed = 0;
    long long u_chain_break = 0;
    long long u_tolerance = 0;
    long long u_max_sd = 0;
    long long unknown = 0;

    void add(const FilterStats& o) {
        self += o.self;
        contained += o.contained;
        short_align += o.short_align;
        low_identity += o.low_identity;
        overhang += o.overhang;
        no_unikmer += o.no_unikmer;
        u_no_shared += o.u_no_shared;
        u_min_shared += o.u_min_shared;
        u_no_valid_seed += o.u_no_valid_seed;
        u_chain_break += o.u_chain_break;
        u_tolerance += o.u_tolerance;
        u_max_sd += o.u_max_sd;
        unknown += o.unknown;
    }
};

// Wrapper for reading plain text or gzip files
class GzFile {
    gzFile file;
public:
    GzFile(const string& path) {
        file = gzopen(path.c_str(), "r");
        if (!file) {
            cerr << "Error opening file: " << path << endl;
            exit(1);
        }
        gzbuffer(file, 128 * 1024);
    }
    ~GzFile() { if (file) gzclose(file); }
    
    bool getline(string& line) {
        line.clear();
        char buf[8192];
        while (gzgets(file, buf, sizeof(buf))) {
            line += buf;
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
        }
        return !line.empty();
    }
};


// Fast k-mer encoding (2-bit)
// A=00(0), C=01(1), G=10(2), T=11(3)
inline uint8_t char_to_val(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // N or other
    }
}

inline uint8_t comp_val(uint8_t v) {
    return v ^ 3; // 0->3, 1->2, 2->1, 3->0
}

// Convert k-mer string to canonical 64-bit integer
uint64_t get_canonical_hash(const string& kmer) {
    uint64_t fwd = 0;
    uint64_t rev = 0;
    int k = kmer.length();
    for (int i = 0; i < k; ++i) {
        uint8_t v = char_to_val(kmer[i]);
        if (v > 3) v = 0; // treat N as A for simplicity, or handle better
        fwd = (fwd << 2) | v;
        rev = rev | ((uint64_t)comp_val(v) << (2 * i));
    }
    return min(fwd, rev);
}

struct KmerHit {
    uint32_t hash_idx;
    uint32_t pos;
};

int global_shift_bits = 0;

int global_prefix_bits = 0;
vector<uint32_t> global_prefix_idx;
vector<uint64_t> global_unikmers;

// Load unique k-mers into a sorted vector (faster lookup, less memory than unordered_set)
vector<uint64_t> load_unikmers(const string& filepath, int k) {
    vector<uint64_t> unikmers;
    // Load unikmers as hashes directly, fast I/O
    FILE* fp = fopen(filepath.c_str(), "r");
    if (!fp) {
        cerr << "Failed to open unikmers file: " << filepath << endl;
        return unikmers;
    }
    char line_buf[256];
    while (fgets(line_buf, sizeof(line_buf), fp)) {
        string kmer;
        for (int i = 0; line_buf[i] && line_buf[i] != ',' && line_buf[i] != '\n'; ++i) {
            kmer += line_buf[i];
        }
        if (kmer.length() >= 10) {
            unikmers.push_back(get_canonical_hash(kmer));
        }
    }
    fclose(fp);
    sort(unikmers.begin(), unikmers.end());
    auto last = unique(unikmers.begin(), unikmers.end());
    unikmers.erase(last, unikmers.end());
    
    // Build prefix index for O(1) fast lookup with minimal memory
    global_shift_bits = (2 * k > 24) ? (2 * k - 24) : 0;
    global_prefix_bits = (2 * k > 24) ? 24 : (2 * k);
    global_prefix_idx.assign((1 << global_prefix_bits) + 1, 0);
    
    uint32_t current_prefix = 0;
    for (uint32_t i = 0; i < unikmers.size(); ++i) {
        uint32_t prefix = unikmers[i] >> global_shift_bits;
        while (current_prefix < prefix) {
            current_prefix++;
            global_prefix_idx[current_prefix] = i;
        }
    }
    while (current_prefix < (1 << global_prefix_bits)) {
        current_prefix++;
        global_prefix_idx[current_prefix] = unikmers.size();
    }
    
    cout << "Loaded " << unikmers.size() << " unique k-mer hashes." << endl;
    return unikmers;
}

// Global fast mapping optimization
struct ReadHits {
    vector<KmerHit> hits;
    vector<pair<uint32_t, uint32_t>> hits_by_hash;
};
vector<ReadHits> global_read_hits;
unordered_map<string, int> read_id_to_idx;

// Read FASTA/FASTQ and map k-mers
void map_reads_to_unikmers(const string& filepath, const vector<uint64_t>& unikmers, int k, const unordered_set<string>& relevant_reads) {
    GzFile f(filepath);
    string line, id, seq;
    int count = 0;
    
    auto process_read = [&](const string& read_id, const string& sequence) {
        if (sequence.length() < k) return;
        vector<KmerHit> hits;
        
        uint64_t fwd = 0, rev = 0;
        uint64_t mask = (1ULL << (2 * k)) - 1;
        
        int valid_bases = 0;
        int seq_len = sequence.length();
        for (int i = 0; i < seq_len; ++i) {
            uint8_t v = char_to_val(sequence[i]);
            if (v > 3) {
                valid_bases = 0;
                fwd = 0; rev = 0;
                continue;
            }
            fwd = ((fwd << 2) | v) & mask;
            rev = (rev >> 2) | ((uint64_t)comp_val(v) << (2 * (k - 1)));
            valid_bases++;
            
            if (valid_bases >= k) {
                uint64_t can_hash = min(fwd, rev);
                uint32_t prefix = can_hash >> global_shift_bits;
                uint32_t start_idx = global_prefix_idx[prefix];
                uint32_t end_idx = global_prefix_idx[prefix + 1];
                int hash_idx = -1;
                // Linear scan in the extremely small range
                for (uint32_t idx = start_idx; idx < end_idx; ++idx) {
                    if (unikmers[idx] == can_hash) {
                        hash_idx = idx;
                        break;
                    }
                }
                if (hash_idx != -1) {
                    hits.push_back({(uint32_t)hash_idx, (uint32_t)(i - k + 1)});
                }
            }
        }
        if (!hits.empty()) {
            vector<pair<uint32_t, uint32_t>> hits_by_hash;
            hits_by_hash.reserve(hits.size());
            for (uint32_t i = 0; i < hits.size(); ++i) {
                hits_by_hash.push_back({hits[i].hash_idx, i});
            }
            sort(hits_by_hash.begin(), hits_by_hash.end(), [](const pair<uint32_t, uint32_t>& a, const pair<uint32_t, uint32_t>& b) {
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second; // largest original_idx first
            });
            auto it_unique = unique(hits_by_hash.begin(), hits_by_hash.end(), [](const pair<uint32_t, uint32_t>& a, const pair<uint32_t, uint32_t>& b) {
                return a.first == b.first;
            });
            hits_by_hash.erase(it_unique, hits_by_hash.end());
            hits_by_hash.shrink_to_fit();

            hits.shrink_to_fit();
            int idx = -1;
            auto it = read_id_to_idx.find(read_id);
            if (it != read_id_to_idx.end()) {
                idx = it->second;
                global_read_hits[idx].hits = move(hits);
                global_read_hits[idx].hits_by_hash = move(hits_by_hash);
            } else {
                idx = global_read_hits.size();
                global_read_hits.push_back({move(hits), move(hits_by_hash)});
                read_id_to_idx[read_id] = idx;
            }
        }
        count++;
        if (count % 1000 == 0) cout << "Scanned " << count << " relevant reads...\r" << flush;
    };

    bool is_fastq = (filepath.find(".fastq") != string::npos || filepath.find(".fq") != string::npos);
    bool current_read_relevant = false;
    
    while (f.getline(line)) {
        if (line.empty()) continue;
        if ((is_fastq && line[0] == '@') || (!is_fastq && line[0] == '>')) {
            if (!id.empty() && current_read_relevant) {
                process_read(id, seq);
            }
            seq.clear();
            
            size_t space_pos = line.find(' ');
            if (space_pos != string::npos) {
                id = line.substr(1, space_pos - 1);
            } else {
                id = line.substr(1);
            }
            current_read_relevant = (relevant_reads.find(id) != relevant_reads.end());
        } else if (is_fastq && line[0] == '+') {
            f.getline(line); // skip quality
        } else {
            if (current_read_relevant) {
                seq += line;
            }
        }
    }
    if (!id.empty() && current_read_relevant) process_read(id, seq);
    cout << "\nProcessing complete. Reads with Unikmers: " << global_read_hits.size() << endl;
}

// Check overhang using Minimap2's definition:
// A true overhang is the minimum of the unaligned lengths at the same side of the alignment.
// If one sequence aligns to its end, the overhang on that side is 0 (or small), which is valid.
bool check_overhang(int q_len, int q_start, int q_end, int t_len, int t_start, int t_end, char strand, int tol) {
    int q_left_len = q_start;
    int q_right_len = q_len - q_end;
    int t_left_len = (strand == '+') ? t_start : (t_len - t_end);
    int t_right_len = (strand == '+') ? (t_len - t_end) : t_start;
    
    int true_left_overhang = std::min(q_left_len, t_left_len);
    int true_right_overhang = std::min(q_right_len, t_right_len);
    
    return (true_left_overhang <= tol && true_right_overhang <= tol);
}

double calculate_sd(const vector<int>& deltas) {
    if (deltas.size() <= 1) return 0.0;
    double sum = 0;
    for (int d : deltas) sum += d;
    double mean = sum / deltas.size();
    double var = 0;
    for (int d : deltas) var += (d - mean) * (d - mean);
    return sqrt(var / deltas.size());
}

// Verify overlap logic
bool verify_overlap(const string& q_id, int q_len, const string& t_id, int t_len, int q_start, int q_end, int t_start, int t_end, char strand, 
                    int min_shared, int tolerance, double max_sd, string& reason, int rl_tag = 0, int align_len = 0, int mapq = 0, double rps_threshold = 0.0) {
    
    auto q_it = read_id_to_idx.find(q_id);
    auto t_it = read_id_to_idx.find(t_id);
    if (q_it == read_id_to_idx.end() || t_it == read_id_to_idx.end()) {
        reason = "NoUnikmer";
        return false;
    }

    const auto& q_hits = global_read_hits[q_it->second].hits;
    const auto& t_hits = global_read_hits[t_it->second].hits;
    const auto& q_by_hash = global_read_hits[q_it->second].hits_by_hash;
    const auto& t_by_hash = global_read_hits[t_it->second].hits_by_hash;

    thread_local vector<pair<uint32_t, uint32_t>> shared_pairs;
    shared_pairs.clear();

    int i = 0, j = 0;
    while (i < q_by_hash.size() && j < t_by_hash.size()) {
        if (q_by_hash[i].first < t_by_hash[j].first) i++;
        else if (q_by_hash[i].first > t_by_hash[j].first) j++;
        else {
            shared_pairs.push_back({q_by_hash[i].second, t_by_hash[j].second});
            i++; j++;
        }
    }

    if (shared_pairs.empty()) {
        reason = "UnikmerFail_NoShared";
        return false;
    }
    if (shared_pairs.size() < min_shared) {
        reason = "UnikmerFail_MinShared";
        return false;
    }

    // Sort shared_pairs by q_idx to find contiguous chains in a single pass
    sort(shared_pairs.begin(), shared_pairs.end(), [](const pair<uint32_t, uint32_t>& a, const pair<uint32_t, uint32_t>& b) {
        return a.first < b.first;
    });

    int best_chain_len = 0;
    bool saw_any_valid_seed = false;
    bool saw_chain_break = false;
    bool saw_tolerance_violation = false;
    bool saw_max_sd_violation = false;
    reason = "UnikmerFail_ChainBreak";

    int n_shared = shared_pairs.size();
    int idx = 0;
    while (idx < n_shared) {
        int start_idx = idx;
        bool has_valid_seed = false;
        
        while (idx < n_shared) {
            int q_idx = shared_pairs[idx].first;
            int t_idx = shared_pairs[idx].second;
            int q_pos = q_hits[q_idx].pos;
            int t_pos = t_hits[t_idx].pos;
            
            if (q_pos >= q_start && q_pos < q_end && 
                t_pos >= t_start && t_pos < t_end) {
                has_valid_seed = true;
                saw_any_valid_seed = true;
            }
            
            if (idx + 1 < n_shared) {
                int next_q = shared_pairs[idx+1].first;
                int next_t = shared_pairs[idx+1].second;
                
                if (next_q == q_idx + 1 && next_t == (strand == '+' ? t_idx + 1 : t_idx - 1)) {
                    int next_q_pos = q_hits[next_q].pos;
                    int next_t_pos = t_hits[next_t].pos;
                    
                    int delta_q = next_q_pos - q_pos;
                    int delta_t = next_t_pos - t_pos;
                    int dist_diff = (strand == '+') ? abs(delta_q - delta_t) : abs(delta_q - abs(delta_t));
                    
                    if (dist_diff <= tolerance) {
                        idx++;
                        continue;
                    }
                    saw_tolerance_violation = true;
                } else {
                    saw_chain_break = true;
                }
            }
            break; // Chain ends
        }
        
        int chain_len = idx - start_idx + 1;
        if (has_valid_seed) {
            if (chain_len > best_chain_len) best_chain_len = chain_len;
            
            if (chain_len >= min_shared) {
                thread_local vector<int> chain_deltas;
                chain_deltas.clear();
                chain_deltas.reserve(chain_len);
                for (int k = start_idx; k <= idx; ++k) {
                    int pq = q_hits[shared_pairs[k].first].pos;
                    int pt = t_hits[shared_pairs[k].second].pos;
                    chain_deltas.push_back((strand == '+') ? (pq - pt) : (pq + pt));
                }
                double sd = calculate_sd(chain_deltas);
                if (sd <= max_sd) {
                    reason = "Pass";
                    return true;
                } else {
                    saw_max_sd_violation = true;
                }
            }
        }
        
        idx++;
    }

    string fail_reason = "UnikmerFail_Unknown";
    if (saw_max_sd_violation) fail_reason = "UnikmerFail_MaxSD";
    else if (saw_tolerance_violation) fail_reason = "UnikmerFail_Tolerance";
    else if (!saw_any_valid_seed) fail_reason = "UnikmerFail_NoValidSeed";
    else if (saw_chain_break || best_chain_len < min_shared) fail_reason = "UnikmerFail_ChainBreak";

    // Only apply Repeat Rescue to cases where valid shared seeds were found but fell outside 
    // the boundary (NoValidSeed). This targets collinear signals truncated by repetitive regions.
    if (fail_reason == "UnikmerFail_NoValidSeed" && rps_threshold > 0.0) {
        // mapq > 0 filters cross-copy FP in reference mapping (always true/ineffective in ava mode)
        if (align_len > 0 && rl_tag > 0 && mapq > 0) {
            int q_left_len = q_start;
            int q_right_len = q_len - q_end;
            int t_left_len = (strand == '+') ? t_start : (t_len - t_end);
            int t_right_len = (strand == '+') ? (t_len - t_end) : t_start;
            
            int true_left_overhang = std::min(q_left_len, t_left_len);
            int true_right_overhang = std::min(q_right_len, t_right_len);
            int max_overhang = std::max(true_left_overhang, true_right_overhang);
            
            // 计算悬垂比例和重复区比例，两者相乘得到一个“重复区惩罚分数” (Repeat Penalty Score, RPS)。
            // RPS 越低越好（代表更符合真实比对特征），只有低于阈值 (rps_threshold) 才能被特赦。
            // 真比对 (FN)：悬垂小，能覆盖大部分重复区 -> RPS 低，被捞回。
            // 假比对 (FP)：悬垂大，深陷巨大的重复区 -> RPS 高，被拒绝。
            double overhang_ratio = static_cast<double>(max_overhang) / align_len;
            double repeat_ratio = static_cast<double>(rl_tag) / align_len;
            
            double rps = overhang_ratio * repeat_ratio;
            
            if (rps < rps_threshold) {
                reason = "Pass_RepeatRescue";
                return true;
            }
        }
    }

    reason = fail_reason;
    return false;
}

// Fast line splitting that avoids allocating strings where possible
void split_paf_line_fast(const string& line, vector<pair<int, int>>& cols) {
    cols.clear();
    int start = 0;
    int len = line.length();
    for (int i = 0; i < len; ++i) {
        if (line[i] == '\t') {
            cols.push_back({start, i - start});
            start = i + 1;
        }
    }
    cols.push_back({start, len - start});
}

// Convert string to int fast
inline int fast_stoi_sub(const string& s, int start, int len) {
    int res = 0;
    for (int i = 0; i < len; ++i) {
        char c = s[start + i];
        if (c >= '0' && c <= '9') {
            res = res * 10 + (c - '0');
        }
    }
    return res;
}

double round_to_two_decimals(double x) {
    return std::round(x * 100.0) / 100.0;
}

double compute_simple_mean_identity_from_paf(const string& paf_file) {
    GzFile paf(paf_file);
    string line;
    vector<pair<int, int>> cols;
    double total_identity = 0.0;
    long long count = 0;

    while (paf.getline(line)) {
        if (line.empty()) continue;
        split_paf_line_fast(line, cols);
        if (cols.size() < 11) continue;

        int matches = fast_stoi_sub(line, cols[9].first, cols[9].second);
        int block_len = fast_stoi_sub(line, cols[10].first, cols[10].second);
        if (block_len <= 0) continue;

        total_identity += static_cast<double>(matches) / block_len;
        count++;
    }

    if (count == 0) {
        cerr << "Error: no valid PAF records found for identity mean calculation." << endl;
        exit(1);
    }

    return round_to_two_decimals(total_identity / count);
}

int main(int argc, char* argv[]) {
    // Parse args manually for simplicity
    string paf_file, unikmers_file, reads_file, output_file, output_numeric_file;
    string read_type = "hifi";
    int k = 21, min_shared = 15, tolerance = 15, min_align_len = 0, overhang_tolerance = 1000;
    double min_identity = 0.0, max_sd = 10.0, rps_threshold = 0.08;
    bool min_shared_explicit = false;
    bool min_identity_explicit = false;
    int configured_min_shared = min_shared;
    string min_shared_source = "default_hifi";
    double configured_min_identity = min_identity;
    string min_identity_source = "default_hifi";
    bool keep_contained = true, ignore_overhang = true, numeric_only = false;
    int threads = 1;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--paf" || arg == "-p") paf_file = argv[++i];
        else if (arg == "--unikmers" || arg == "-u") unikmers_file = argv[++i];
        else if (arg == "--reads" || arg == "-r") reads_file = argv[++i];
        else if (arg == "--output" || arg == "-o") output_file = argv[++i];
        else if (arg == "--output-numeric") output_numeric_file = argv[++i];
        else if (arg == "-k") k = stoi(argv[++i]);
        else if (arg == "-s" || arg == "--min-shared") { min_shared = stoi(argv[++i]); min_shared_explicit = true; }
        else if (arg == "-t" || arg == "--tolerance") tolerance = stoi(argv[++i]);
        else if (arg == "--max-sd") max_sd = stod(argv[++i]);
        else if (arg == "-l" || arg == "--min-align-len") min_align_len = stoi(argv[++i]);
        else if (arg == "-i" || arg == "--min-identity") { min_identity = stod(argv[++i]); min_identity_explicit = true; }
        else if (arg == "-e" || arg == "--overhang-tolerance") overhang_tolerance = stoi(argv[++i]);
        else if (arg == "--read-type" || arg == "--mode") read_type = argv[++i];
        else if (arg == "--threads" || arg == "-T") threads = stoi(argv[++i]);
        else if (arg == "--drop-contained") keep_contained = false;
        else if (arg == "--keep-overhang") ignore_overhang = false;
        else if (arg == "--numeric-only") numeric_only = true;
        else if (arg == "-R" || arg == "--rps-threshold") rps_threshold = stod(argv[++i]);
    }

    if (paf_file.empty() || output_file.empty() || reads_file.empty() || unikmers_file.empty()) {
        cerr << "Usage: " << argv[0] << " --paf in.paf --unikmers in.kmer --reads in.fa --output out.paf [--output-numeric numeric.paf] [--read-type hifi|ont] [options]\n";
        return 1;
    }

    if (numeric_only && !output_numeric_file.empty()) {
        cerr << "Error: --numeric-only cannot be combined with --output-numeric. "
             << "Use --output to specify the numeric PAF path in numeric-only mode." << endl;
        return 1;
    }

    if (read_type != "hifi" && read_type != "ont") {
        cerr << "Error: --read-type must be either 'hifi' or 'ont'." << endl;
        return 1;
    }

    configured_min_shared = min_shared;
    if (read_type == "ont") {
        if (min_shared_explicit) {
            min_shared_source = "explicit_argument";
        } else {
            min_shared = 6;
            min_shared_source = "default_ont";
        }
    } else {
        if (min_shared_explicit) {
            min_shared_source = "explicit_argument";
        } else {
            min_shared = 15;
            min_shared_source = "default_hifi";
        }
    }

    configured_min_identity = min_identity;
    if (read_type == "ont") {
        if (min_identity_explicit) {
            min_identity_source = "explicit_argument";
        } else {
            min_identity = compute_simple_mean_identity_from_paf(paf_file);
            min_identity_source = "simple_mean_from_raw_paf";
        }
    } else {
        min_identity_source = min_identity_explicit ? "explicit_argument" : "default_hifi";
    }

    cout << fixed << setprecision(2);
    cout << "Read type: " << read_type
         << ", configured min_shared: " << configured_min_shared
         << ", effective min_shared: " << min_shared
         << " (" << min_shared_source << ")"
         << ", configured min_identity: " << configured_min_identity
         << ", effective min_identity: " << min_identity
         << " (" << min_identity_source << ")" << endl;

    auto start_time = chrono::high_resolution_clock::now();

    // Pass 1: Scan PAF for relevant reads
    unordered_set<string> relevant_reads;
    cout << "Pass 1: Scanning PAF for relevant reads..." << endl;
    {
        GzFile paf(paf_file);
        string line;
        vector<pair<int, int>> cols;
        while (paf.getline(line)) {
            if (line.empty()) continue;
            split_paf_line_fast(line, cols);
            if (cols.size() < 12) continue;
            relevant_reads.insert(line.substr(cols[0].first, cols[0].second));
            relevant_reads.insert(line.substr(cols[5].first, cols[5].second));
        }
    }
    cout << "Found " << relevant_reads.size() << " unique reads in PAF." << endl;

    global_unikmers = load_unikmers(unikmers_file, k);
    map_reads_to_unikmers(reads_file, global_unikmers, k, relevant_reads);

    // Free memory not needed anymore for the next phases
    relevant_reads.clear();
    
    // We only need the hash_idx stored in KmerHit, not the actual 64-bit kmers
    global_unikmers.clear();
    global_unikmers.shrink_to_fit();
    global_prefix_idx.clear();
    global_prefix_idx.shrink_to_fit();

    // Stream PAF in chunks
    cout << "Reading and filtering PAF file..." << endl;
    GzFile paf2(paf_file);
    string line2;
    int kept_count = 0;
    int numeric_kept_count = 0;
    int total_processed = 0;
    ofstream out(output_file);
    ofstream out_numeric;
    bool write_numeric_output = !output_numeric_file.empty();
    if (write_numeric_output) {
        out_numeric.open(output_numeric_file);
        if (!out_numeric.is_open()) {
            cerr << "Error: cannot open numeric output file: " << output_numeric_file << endl;
            return 1;
        }
    }
    FilterStats global_stats;
    
    vector<string> chunk_lines;
    int chunk_capacity = 200000;
    chunk_lines.reserve(chunk_capacity);
    
    auto process_and_write_chunk = [&]() {
        if (chunk_lines.empty()) return;
        vector<char> keep(chunk_lines.size(), 0);
        vector<char> keep_numeric;
        if (write_numeric_output) {
            keep_numeric.assign(chunk_lines.size(), 0);
        }
        
        auto process_chunk = [&](int start, int end, FilterStats* local_stats) {
            vector<pair<int, int>> cols;
            cols.reserve(32);
            for (int i = start; i < end; ++i) {
                const string& paf_line = chunk_lines[i];
                split_paf_line_fast(paf_line, cols);
                if (cols.size() < 12) continue;
                
                string q_name = paf_line.substr(cols[0].first, cols[0].second);
                string t_name = paf_line.substr(cols[5].first, cols[5].second);
            int q_len = fast_stoi_sub(paf_line, cols[1].first, cols[1].second);
            int q_start = fast_stoi_sub(paf_line, cols[2].first, cols[2].second);
            int q_end = fast_stoi_sub(paf_line, cols[3].first, cols[3].second);
            char strand = paf_line[cols[4].first];
            int t_len = fast_stoi_sub(paf_line, cols[6].first, cols[6].second);
            int t_start = fast_stoi_sub(paf_line, cols[7].first, cols[7].second);
            int t_end = fast_stoi_sub(paf_line, cols[8].first, cols[8].second);
            int matches = fast_stoi_sub(paf_line, cols[9].first, cols[9].second);
            int block_len = fast_stoi_sub(paf_line, cols[10].first, cols[10].second);
            int mapq = fast_stoi_sub(paf_line, cols[11].first, cols[11].second);
            
            // Extract rl tag and alignment length for potential rescue
            int rl_tag = 0;
            int align_len = t_end - t_start;
            for (size_t c = 12; c < cols.size(); ++c) {
                string tag = paf_line.substr(cols[c].first, cols[c].second);
                if (tag.rfind("rl:i:", 0) == 0) {
                    rl_tag = std::stoi(tag.substr(5));
                    break;
                }
            }
            
            string reason = "";
            if (q_name == t_name) reason = "Self";
            else if (!keep_contained) {
                bool is_contained = false;
                if (strand == '+') {
                    if ((q_start <= t_start && q_len - q_end <= t_len - t_end) || 
                        (t_start <= q_start && t_len - t_end <= q_len - q_end)) {
                        is_contained = true;
                    }
                } else {
                    if ((q_start <= t_len - t_end && q_len - q_end <= t_start) || 
                        (t_len - t_end <= q_start && t_start <= q_len - q_end)) {
                        is_contained = true;
                    }
                }
                if (is_contained) reason = "Contained";
            }
            
            if (reason.empty() && !ignore_overhang) {
                if (!check_overhang(q_len, q_start, q_end, t_len, t_start, t_end, strand, overhang_tolerance)) {
                    reason = "Overhang";
                }
            }
            
            if (reason.empty() && block_len < min_align_len) reason = "ShortAlign";
            double curr_identity = block_len > 0 ? ((double)matches / block_len) : 0.0;
            if (reason.empty() && curr_identity < min_identity) reason = "LowIdentity";
            
            bool numeric_pass = reason.empty();

            if (reason.empty() && !numeric_only) {
                verify_overlap(q_name, q_len, t_name, t_len, q_start, q_end, t_start, t_end, strand, 
                               min_shared, tolerance, max_sd, reason, rl_tag, align_len, mapq, rps_threshold);
            }
            
            if (reason == "Self") local_stats->self++;
            else if (reason == "Contained") local_stats->contained++;
            else if (reason == "ShortAlign") local_stats->short_align++;
            else if (reason == "LowIdentity") local_stats->low_identity++;
            else if (reason == "Overhang") local_stats->overhang++;
            else if (reason == "NoUnikmer") local_stats->no_unikmer++;
            else if (reason == "UnikmerFail_NoShared") local_stats->u_no_shared++;
            else if (reason == "UnikmerFail_MinShared") local_stats->u_min_shared++;
            else if (reason == "UnikmerFail_NoValidSeed") local_stats->u_no_valid_seed++;
            else if (reason == "UnikmerFail_ChainBreak") local_stats->u_chain_break++;
            else if (reason == "UnikmerFail_Tolerance") local_stats->u_tolerance++;
            else if (reason == "UnikmerFail_MaxSD") local_stats->u_max_sd++;
            else if (reason == "Pass_RepeatRescue") {
                // Treated as Pass, but we could add a counter if we want to track rescued reads
            }
            else if (!reason.empty() && reason != "Pass") local_stats->unknown++;

            if (write_numeric_output && numeric_pass) {
                keep_numeric[i] = 1;
            }

            if (reason.empty() || reason == "Pass" || reason == "Pass_RepeatRescue") {
                keep[i] = 1;
            }
        }
    };

        vector<thread> pool;
        vector<FilterStats> thread_stats(threads);
        int chunk_size_per_thread = chunk_lines.size() / threads;
        if (chunk_size_per_thread == 0) chunk_size_per_thread = 1;
        for (int i = 0; i < threads; ++i) {
            int start = i * chunk_size_per_thread;
            int end = (i == threads - 1) ? chunk_lines.size() : start + chunk_size_per_thread;
            if (start >= chunk_lines.size()) break;
            pool.emplace_back(process_chunk, start, end, &thread_stats[i]);
        }
        for (auto& t : pool) t.join();

        for (int i = 0; i < threads; ++i) {
            global_stats.add(thread_stats[i]);
        }

        for (size_t i = 0; i < chunk_lines.size(); ++i) {
            if (write_numeric_output && keep_numeric[i]) {
                out_numeric << chunk_lines[i] << "\n";
                numeric_kept_count++;
            }
            if (keep[i]) {
                out << chunk_lines[i] << "\n";
                kept_count++;
            }
        }
        chunk_lines.clear();
    };

    cout << "Filtering alignments using " << threads << " threads..." << endl;
    while (paf2.getline(line2)) {
        if (line2.empty()) continue;
        chunk_lines.push_back(line2);
        total_processed++;
        if (chunk_lines.size() >= chunk_capacity) {
            process_and_write_chunk();
        }
    }
    if (!chunk_lines.empty()) {
        process_and_write_chunk();
    }
    
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end_time - start_time;
    
    // Write statistics file
    string stats_file = output_file + ".stats.txt";
    ofstream f_stats(stats_file);
    if (f_stats.is_open()) {
        f_stats << "=== Ovfilter Statistics ===\n";
        f_stats << "Read Type:        " << read_type << "\n";
        f_stats << "Min Shared Source: " << min_shared_source << "\n";
        f_stats << "Configured Min Shared: " << configured_min_shared << "\n";
        f_stats << "Effective Min Shared:  " << min_shared << "\n";
        f_stats << "Identity Source:  " << min_identity_source << "\n";
        f_stats << "Configured Min Identity: " << configured_min_identity << "\n";
        f_stats << "Effective Min Identity:  " << min_identity << "\n";
        f_stats << "Total Processed: " << total_processed << "\n";
        double kept_pct = total_processed > 0 ? (double)kept_count / total_processed * 100.0 : 0.0;
        f_stats << "Kept Overlaps:   " << kept_count << " (" << fixed << setprecision(2) << kept_pct << "%)\n";
        if (write_numeric_output) {
            double numeric_kept_pct = total_processed > 0 ? (double)numeric_kept_count / total_processed * 100.0 : 0.0;
            f_stats << "Numeric Kept:    " << numeric_kept_count << " (" << fixed << setprecision(2) << numeric_kept_pct << "%)\n";
        }
        f_stats << "Discarded:       " << (total_processed - kept_count) << "\n\n";
        f_stats << "--- Discard Reasons ---\n";

        auto print_stat = [&](const string& name, long long count) {
            double pct = total_processed > 0 ? (double)count / total_processed * 100.0 : 0.0;
            f_stats << left << setw(20) << name << ": " << left << setw(12) << count 
                    << " (" << fixed << setprecision(2) << pct << "%)\n";
        };

        print_stat("self", global_stats.self);
        print_stat("contained", global_stats.contained);
        print_stat("short", global_stats.short_align);
        print_stat("identity", global_stats.low_identity);
        print_stat("overhang", global_stats.overhang);
        print_stat("no_unikmer", global_stats.no_unikmer);

        long long total_unikmer_fail = global_stats.u_no_shared + global_stats.u_min_shared + 
                                       global_stats.u_no_valid_seed + global_stats.u_chain_break +
                                       global_stats.u_tolerance + global_stats.u_max_sd + global_stats.unknown;
        print_stat("unikmer_total", total_unikmer_fail);

        auto print_sub_stat = [&](const string& name, long long count) {
            double pct = total_processed > 0 ? (double)count / total_processed * 100.0 : 0.0;
            f_stats << "  - " << left << setw(16) << name << ": " << left << setw(12) << count 
                    << " (" << fixed << setprecision(2) << pct << "%)\n";
        };

        print_sub_stat("no_shared", global_stats.u_no_shared);
        print_sub_stat("min_shared", global_stats.u_min_shared);
        print_sub_stat("no_valid_seed", global_stats.u_no_valid_seed);
        print_sub_stat("chain_break", global_stats.u_chain_break);
        print_sub_stat("tolerance", global_stats.u_tolerance);
        print_sub_stat("max_sd", global_stats.u_max_sd);
        print_sub_stat("unknown", global_stats.unknown);

        f_stats.close();
        cout << "Statistics saved to: " << stats_file << endl;
    } else {
        cerr << "Failed to write statistics file: " << stats_file << endl;
    }

    cout << "\n\n=== Filtering Complete ===" << endl;
    cout << "Total Processed: " << total_processed << ", Kept: " << kept_count << endl;
    cout << "Discard details: {'self': " << global_stats.self 
         << ", 'contained': " << global_stats.contained 
         << ", 'short': " << global_stats.short_align 
         << ", 'identity': " << global_stats.low_identity 
         << ", 'overhang': " << global_stats.overhang 
         << ", 'no_unikmer': " << global_stats.no_unikmer 
         << ", 'unikmer_no_shared': " << global_stats.u_no_shared 
         << ", 'unikmer_min_shared': " << global_stats.u_min_shared 
         << ", 'unikmer_no_valid_seed': " << global_stats.u_no_valid_seed 
         << ", 'unikmer_chain_break': " << global_stats.u_chain_break 
         << ", 'unikmer_tolerance': " << global_stats.u_tolerance 
         << ", 'unikmer_max_sd': " << global_stats.u_max_sd 
         << ", 'unknown': " << global_stats.unknown << "}" << endl;
    cout << "Time elapsed: " << diff.count() << " s" << endl;

    return 0;
}



