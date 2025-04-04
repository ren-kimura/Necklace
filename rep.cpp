#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <cstdint>
#include <unistd.h>
#include <memory_resource>
#include <string>
#include <vector>
#include <queue>
#include <array> // sada
#include <tuple> // sada
#include <unordered_map>
#include <unordered_set>
#include <algorithm>


using namespace std;
namespace fs = std::filesystem;
namespace chrono = std::chrono;

using INT = uint64_t;
using VSTR = vector<string>;
using Vint = vector<uint8_t>;
using VINT = vector<INT>;
using PINT = pair<INT, INT>;
using VTINT = vector<tuple<INT, INT, INT, INT>>;
using VVINT = vector<VINT>;
using REP = pair<string, VINT>;
using USET = unordered_set<INT>;

unordered_set<char> base = {'A', 'C', 'G', 'T'};
const INT INF = UINT64_MAX;

size_t get_available_memory() {
    return sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGE_SIZE);
}

void progress(INT now, INT total, const string& task) {
    const int disp_interval = total / 100;
    if (total >= 100 && now % disp_interval != 0) return;
    double progress = (double)(now) / (double)(total);
    int width = 30;
    int pos = (int)(width * progress);
    cout << "\r" << task << " [";
    for (int i = 0; i < width; ++i) {
        if (i < pos)        cout << "=";
        else if (i == pos)  cout << ">";
        else                cout << " ";
    }
    cout << "] " << (int)(progress * 100.0) << "%";
    cout.flush();
}

void finished(const string& task) {
    const int width = 30;
    cout << "\r" << task << " [";
    for (int i = 0; i < width - 1; ++i)
        cout << "=";
    cout << ">] 100%";
}

string remove_extension(string str, const string& ex) {
    if (str.size() >= ex.size() && str.rfind(ex) == str.size() - ex.size()) {
        str.erase(str.size() - ex.size());
    }
    return str;
}

string path_to_filename(string path, const char& dlm) {
    size_t pos = path.find_last_of(dlm);
    if (pos != string::npos)  path = path.substr(pos + 1);
    return path;
}

string get_current_timestamp() {
    auto now = chrono::system_clock::now();
    auto in_time_t = chrono::system_clock::to_time_t(now);
    stringstream ss;
    ss << put_time(localtime(&in_time_t), "%Y-%m-%d %X");
    return ss.str();
}

size_t get_memory_usage() {
    size_t memory_used = sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGE_SIZE);
    return memory_used / (1024 * 1024);
}

void log_time_and_memory(const string& task_name, const chrono::duration<double>& elapsed_time, size_t memory_usage, ofstream& logfile) {
    string task_name_trimmed = task_name.substr(0, 20);

    logfile << left << setw(20) << task_name_trimmed
            << setw(20) << fixed << setprecision(4) << elapsed_time.count()
            << setw(20) << memory_usage << "\n";

    cout << left << setw(20) << task_name_trimmed
         << setw(20) << fixed << setprecision(4) << elapsed_time.count()
         << setw(20) << memory_usage << "\n";
}

void print_table_header(ofstream& logfile) {
    const int col_width = 20;

    logfile << string(col_width * 3, '-') << "\n";
    logfile << left << setw(col_width) << "Task name"
            << setw(col_width) << "Time (s)"
            << setw(col_width) << "Memory (MB)" << endl;
    logfile << string(col_width * 3, '-') << "\n"; 

    cout << string(col_width * 3, '-') << "\n";
    cout << left << setw(col_width) << "Task name"
         << setw(col_width) << "Time (s)"
         << setw(col_width) << "Memory (MB)" << endl;
    cout << string(col_width * 3, '-') << "\n"; 
}

void to_uppercase(VSTR& strs) {
    for (auto& str: strs)
        transform(str.begin(), str.end(), str.begin(),
                [](unsigned char c) {return toupper(c);});
}

uint64_t encode_kmer(const string& kmer) {
    uint64_t val = 0;
    for (char c : kmer) {
        if (c == 'A')      val = (val << 2) | 0; // A -> 00
        else if (c == 'C') val = (val << 2) | 1; // C -> 01
        else if (c == 'G') val = (val << 2) | 2; // G -> 10
        else if (c == 'T') val = (val << 2) | 3; // T -> 11
        else return UINT64_MAX;
    }
    return val;
}

string decode_kmer(uint64_t hash, INT K) {
    string kmer;
    for (INT i = 0; i < K; ++i) {
        int base = hash & 3;  // take rightmost 2 bits
        if (base == 0)        kmer = 'A' + kmer;  // 00 -> A
        else if (base == 1)   kmer = 'C' + kmer;  // 01 -> C
        else if (base == 2)   kmer = 'G' + kmer;  // 10 -> G
        else if (base == 3)   kmer = 'T' + kmer;  // 11 -> T
        hash >>= 2;
    }
    return kmer;
}

char decode_base(uint64_t val) {
    return "ACGT"[val];
}

INT rolling_valid_kmer_start(const string& read, INT start, INT K) {
    INT len = read.size();
    INT count = 0;
    for (INT j = start; j < len; ++j) {
        if (read[j] == 'A' || read[j] == 'C' ||
            read[j] == 'G' || read[j] == 'T') {
            ++count;
            if (count == K) return j - K + 1;
        } else {
            count = 0; // reset if read[j] is non-ACGT
        }
    }
    return INF;
}

class EDBG{
public:
    string filename;
    INT K;

    size_t total_memory = get_available_memory();
    size_t pool_size = (total_memory > (64L * 1024 * 1024 * 1024)) ? (16L * 1024 * 1024 * 1024) : // 16GB if more than 64GB
                       (total_memory > (16L * 1024 * 1024 * 1024)) ? (4L * 1024 * 1024 * 1024) :  // 4GB if more than 16GB
                                                                      (1L * 1024 * 1024 * 1024);  // else 1GB
    size_t reserve_size = (total_memory > (64L * 1024 * 1024 * 1024)) ? 500'000'000 :
                          (total_memory > (16L * 1024 * 1024 * 1024)) ? 100'000'000 :
                                                                         10'000'000;
    
    pmr::monotonic_buffer_resource pool{pool_size};
    using UMAP = pmr::unordered_map<INT, INT>;
    pmr::unordered_map<INT, int64_t> k_minus_1mers;
    USET kmers;

    EDBG(string _filename, INT _K)
        : filename(_filename), K(_K), k_minus_1mers(&pool)
    {
        k_minus_1mers.reserve(reserve_size);
        cout << "\nUsing pool size: " << pool_size / (1024 * 1024) << " MB\n";
        cout << "K-mers reserve size: " << reserve_size << "\n";
    };

    INT process() {
        string timestamp = get_current_timestamp();
        cout << "Timestamp: " << timestamp << "\n";
        size_t initial_memory = get_memory_usage();
        cout << "Initial Memory: " << initial_memory << " MB\n\n";
        
        INT S = get_kmers(k_minus_1mers, kmers);
        return S;
    }

    INT get_kmers(pmr::unordered_map<INT, int64_t>& k_minus_1mers, USET& kmers) {
        ifstream inputFile(filename, ios::in | ios::binary);
        if (!inputFile) {
            cerr << "Error opening input file." << "\n";
            exit(1);
        }

        const size_t BUF_SIZE = 64 * 1024 * 1024; // 64MB
        char *buffer = new char[BUF_SIZE];
        inputFile.rdbuf()->pubsetbuf(buffer, BUF_SIZE);

        VSTR reads;
        string line, read;
        bool is_new_read = true;
        while (getline(inputFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!read.empty()) {
                    reads.emplace_back(read);
                    read.clear();
                }
                is_new_read = true;
                continue;
            }
            if (is_new_read) {
                read = line;
                is_new_read = false;
            } else read += line;
        }
        if (!read.empty()) reads.emplace_back(read);
        inputFile.close();
        delete[] buffer;

        INT tlen = 0;
        to_uppercase(reads);
        for (string& read: reads)
            tlen += read.size();
        
        INT id = 0, footing = 0;
        uint64_t mask = (1ULL << (2 * (K - 1))) - 1; // 2K-bits mask

        for (INT i = 0; i < (INT)(reads.size()); ++i) {
            const string& read = reads[i];
            INT len = read.size();
            if (len < K) continue;
            uint64_t hash;

            INT j = rolling_valid_kmer_start(read, 0, K);
            if (j == INF) continue; // if not any valid k-mer, go to next read
            hash = encode_kmer(read.substr(j, K));
            auto it = kmers.find(hash);
            if (it == kmers.end()) {
                kmers.insert(hash);
                --k_minus_1mers[hash / 4];
                ++k_minus_1mers[hash & mask];
                ++id;
                if (id == INF) {
                    cerr << "Too much k-mers!!\n\n";
                    exit(1);
                }
            }

            // update kmer by rolling hash
            for (++j; j <= len - K; ++j) {
                char c_in = read[j + K - 1];

                if (c_in == 'A')        hash = ((hash & mask) << 2) | 0;
                else if (c_in == 'C')   hash = ((hash & mask) << 2) | 1;
                else if (c_in == 'G')   hash = ((hash & mask) << 2) | 2;
                else if (c_in == 'T')   hash = ((hash & mask) << 2) | 3;
                else {
                    j = rolling_valid_kmer_start(read, j + 1, K);
                    if (j == INF) break;
                    hash = encode_kmer(read.substr(j, K));
                }

                // register k-mer to kmers
                auto it = kmers.find(hash);
                if (it == kmers.end()) {
                    kmers.insert(hash);
                    --k_minus_1mers[hash / 4];
                    ++k_minus_1mers[hash & mask];
                    ++id;
                    if (id == INF) {
                        cerr << "Too much k-mers!!\n\n";
                        exit(1);
                    }
                }

                progress(footing + j, tlen, "Detecting k-mers");
            }
            footing += len;
        }
        finished("Detecting k-mers");

        INT N = k_minus_1mers.size();
        INT M = kmers.size();
        cout << "\nTotal number of (k-1)-mers: " << N << "\n"
            << "Total number of k-mers: " << M << "\n";

        // calc number of breaking edges
        INT S = 0;
        for (const auto& entry: k_minus_1mers) {
            int64_t val = entry.second;
            if (val > 0) S += val;
        }
        cout << "Number of breaking edges: " << S << "\n\n";
        return S * (K - 1) + M;
    }
};

class NDBG{
public:
    string filename, logfilename;
    INT K;
    INT option;

    size_t total_memory = get_available_memory();
    size_t pool_size = (total_memory > (64L * 1024 * 1024 * 1024)) ? (16L * 1024 * 1024 * 1024) : // 16GB if more than 64GB
                       (total_memory > (16L * 1024 * 1024 * 1024)) ? (4L * 1024 * 1024 * 1024) :  // 4GB if more than 16GB
                                                                      (1L * 1024 * 1024 * 1024);  // else 1GB
    size_t reserve_size = (total_memory > (64L * 1024 * 1024 * 1024)) ? 500'000'000 :
                          (total_memory > (16L * 1024 * 1024 * 1024)) ? 100'000'000 :
                                                                         10'000'000;
    
    pmr::monotonic_buffer_resource pool{pool_size};
    using UMAP = pmr::unordered_map<INT, INT>;
    UMAP kmers;
    UMAP heads;

    NDBG(string _filename, INT _K, INT _option)
        : filename(_filename), K(_K), option(_option), kmers(&pool), heads(&pool)
    {
        if (option != 0 && option != 1 && option != 2) {
            cerr << "\nInvalid option value\n";
            exit(1);
        }
        logfilename = remove_extension(filename, ".fa") + ".ours.k" + to_string(K) + ".opt" + to_string(option) + ".log";
        kmers.reserve(reserve_size);
        cout << "\nUsing pool size: " << pool_size / (1024 * 1024) << " MB\n";
        cout << "K-mers reserve size: " << reserve_size << "\n";
    };

    REP process() {
        ofstream logfile(logfilename, ios::trunc);
        if (!logfile) {
            cerr << "Error: Could not open " << logfilename << " for writing.\n";
            exit(1);
        }
        string timestamp = get_current_timestamp();
        logfile << "Timestamp: " << timestamp << "\n";
        cout << "Timestamp: " << timestamp << "\n";

        size_t initial_memory = get_memory_usage();
        logfile << "Initial Memory: " << initial_memory << " MB\n\n";
        cout << "Initial Memory: " << initial_memory << " MB\n\n";
        VINT kmerv; 
        
        auto start_time = chrono::high_resolution_clock::now();
        INT N = get_kmers(kmers, kmerv);
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> get_kmers_time = end_time - start_time;
        size_t get_kmers_memory = get_memory_usage();

        VVINT adj(N), inv_adj(N);

        start_time = chrono::high_resolution_clock::now();
        INT E = add_edges(kmers, kmerv, N, adj, inv_adj);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> add_edges_time = end_time - start_time;
        size_t add_edges_memory = get_memory_usage();

        VINT match_u, match_v;

        start_time = chrono::high_resolution_clock::now();
        INT M = hopcroft_karp(N, adj, match_u, match_v);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> hopcroft_karp_time = end_time - start_time;
        size_t hopcroft_karp_memory = get_memory_usage();

        VVINT cycles, paths;

        start_time = chrono::high_resolution_clock::now();
        PINT CnP = decompose(N, match_u, match_v, cycles, paths);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> decompose_time = end_time - start_time;
        size_t decompose_memory = get_memory_usage();

        cout << "(#k-mers, #edges, #matching, #cycles, #paths) = (" 
            << N << ", " << E << ", " << M << ", "
            << CnP.first << ", " << CnP.second << ")\n";

        REP rep;
        chrono::duration<double> align_time;
        size_t align_memory = 0;
        if (option == 0) {
            // plain
            start_time = chrono::high_resolution_clock::now();
            rep = plain(kmerv, cycles, paths);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
            align_memory = get_memory_usage();
        } else if (option == 1) {
            // unsorted
            start_time = chrono::high_resolution_clock::now();
            rep = unsorted(kmers, N, kmerv, inv_adj, cycles, paths, CnP);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
            align_memory = get_memory_usage();
        } else {
            // sorted
            start_time = chrono::high_resolution_clock::now();
            rep = sorted(kmers, N, kmerv, inv_adj, cycles, paths, CnP);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
            align_memory = get_memory_usage();
        }

        // display benchmarks
        print_table_header(logfile);
        log_time_and_memory("get_kmers", get_kmers_time, get_kmers_memory, logfile);
        log_time_and_memory("add_edges", add_edges_time, add_edges_memory, logfile);
        log_time_and_memory("hopcroft_karp", hopcroft_karp_time, hopcroft_karp_memory, logfile);
        log_time_and_memory("decompose", decompose_time, decompose_memory, logfile);
        log_time_and_memory((option == 0 ? "plain" : option == 1 ? "unsorted" : "sorted"),
                            align_time, align_memory, logfile);
        logfile << string(60, '-');
        cout << string(60, '-');

        size_t final_memory = get_memory_usage();
        logfile << "\n\nFinal Memory: " << final_memory << " MB\n";
        cout << "\nFinal Memory: " << final_memory << " MB\n";

        auto total_time = get_kmers_time.count() + add_edges_time.count()
                          + hopcroft_karp_time.count() + decompose_time.count()
                          + align_time.count();
        logfile << "Total time: " << total_time << " s\n";
        cout << "Total time: " << total_time << " s\n";
        logfile.close();

        return rep;
    }

    INT get_kmers(UMAP& kmers, VINT& kmerv) {
        ifstream inputFile(filename, ios::in | ios::binary);
        if (!inputFile) {
            cerr << "Error opening input file." << "\n";
            exit(1);
        }

        const size_t BUF_SIZE = 64 * 1024 * 1024; // 64MB
        char *buffer = new char[BUF_SIZE];
        inputFile.rdbuf()->pubsetbuf(buffer, BUF_SIZE);

        VSTR reads;
        string line, read;
        bool is_new_read = true;
        while (getline(inputFile, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!read.empty()) {
                    reads.emplace_back(read);
                    read.clear();
                }
                is_new_read = true;
                continue;
            }
            if (is_new_read) {
                read = line;
                is_new_read = false;
            } else read += line;
        }
        if (!read.empty()) reads.emplace_back(read);
        inputFile.close();
        delete[] buffer;

        INT tlen = 0;
        to_uppercase(reads);
        for (string& read: reads)
            tlen += read.size();
        
        INT id = 0, footing = 0;
        uint64_t mask = (1ULL << (2 * (K - 1))) - 1; // 2K-bits mask

        for (INT i = 0; i < (INT)(reads.size()); ++i) {
            const string& read = reads[i];
            INT len = read.size();
            if (len < K) continue;
            uint64_t hash;

            INT j = rolling_valid_kmer_start(read, 0, K);
            if (j == INF) continue; // if not any valid k-mer, go to next read
            hash = encode_kmer(read.substr(j, K));
            auto it = kmers.find(hash);
            if (it == kmers.end()) {
                kmers[hash] = id;
                kmerv.emplace_back(hash);
                ++id;
                if (id == INF) {
                    cerr << "Too much k-mers!!\n\n";
                    exit(1);
                }
            }

            // update kmer by rolling hash
            for (++j; j <= len - K; ++j) {
                char c_in = read[j + K - 1];

                if (c_in == 'A')        hash = ((hash & mask) << 2) | 0;
                else if (c_in == 'C')   hash = ((hash & mask) << 2) | 1;
                else if (c_in == 'G')   hash = ((hash & mask) << 2) | 2;
                else if (c_in == 'T')   hash = ((hash & mask) << 2) | 3;
                else {
                    j = rolling_valid_kmer_start(read, j + 1, K);
                    if (j == INF) break;
                    hash = encode_kmer(read.substr(j, K));
                }

                // register k-mer to kmers
                auto it = kmers.find(hash);
                if (it == kmers.end()) {
                    kmers[hash] = id;
                    kmerv.emplace_back(hash);
                    ++id;
                    if (id == INF) {
                        cerr << "Too much k-mers!!\n\n";
                        exit(1);
                    }
                }

                progress(footing + j, tlen, "Detecting k-mers");
            }
            footing += len;
        }
        finished("Detecting k-mers");

        INT N = kmers.size();
        cout << "\nTotal number of k-mers: " << N << "\n\n";
        return N;
    }

    INT forward(UMAP& kmers, const VINT& kmerv, INT id, char c) {
        INT hash = kmerv[id]; 
        if (hash == INF) return INF;

        INT mask = (1ULL << (2 * (K - 1))) - 1;
        hash = (hash & mask) << 2;

        if (c == 'A')       hash |= 0;
        else if (c == 'C')  hash |= 1;
        else if (c == 'G')  hash |= 2;
        else if (c == 'T')  hash |= 3;
        else return INF; // invalid for non-ACGT c

        if (kmers.find(hash) != kmers.end() && kmers[hash] != id)
            return kmers[hash];
        
        return INF; // no branch to c
    }

    INT add_edges(UMAP& kmers, const VINT& kmerv, const INT& N, VVINT& adj, VVINT& inv_adj) {
        INT cnt = 0, E = 0;
        adj.resize(N);
        inv_adj.resize(N);
        for (const auto& entry : kmers) {
            auto id = entry.second;
            for (auto const c : base) {
                auto next_id = forward(kmers, kmerv, id, c);
                if (next_id != INF && next_id != id) {
                    adj[id].emplace_back(next_id);
                    inv_adj[next_id].emplace_back(id);
                    ++E;
                }
            }
            progress(++cnt, N, "Constructing de Bruijn graph");
        }
        finished("Constructing de Bruijn graph");

        cout << "\nTotal number of edges: " << E << "\n\n";
        return E;
    }

    bool bfs (const INT& N, const VVINT& adj, 
            const VINT& match_u, const VINT& match_v, VINT& dist) {
        queue<INT> q;
        bool found_augpath = false;
        for (INT u = 0; u < N; ++u) {
            if (match_u[u] == INF) {
                dist[u] = 0;
                q.push(u);
            } else dist[u] = INF;
        }
        while (!q.empty()) {
            INT u = q.front();
            q.pop();
            for (auto v : adj[u]) {
                if (match_v[v] == INF) {
                    found_augpath = true;
                } else if (dist[match_v[v]] == INF) {
                    dist[match_v[v]] = dist[u] + 1;
                    q.push(match_v[v]);
                }
            }
        }
        return found_augpath;
    }

    bool dfs (const VVINT& adj, VINT& match_u, VINT& match_v, 
            VINT& dist, INT u) {
        for (auto v : adj[u]) {
            if (match_v[v] == INF || (dist[match_v[v]] == dist[u] + 1 && dfs(adj, match_u, match_v, dist, match_v[v]))) {
                match_u[u] = v;
                match_v[v] = u;
                return true;
            }
        }
        dist[u] = INF; // not found any augpath
        return false;
    }

    INT hopcroft_karp(const INT& N, const VVINT& adj, VINT& match_u,
                VINT& match_v) {
        VINT dist(N); // distance of bfs
        INT M = 0;
        match_u.assign(N, INF);
        match_v.assign(N, INF);
        while (bfs(N, adj, match_u, match_v, dist)) {
            for (INT u = 0; u < N; ++u) {
                if (match_u[u] == INF && dfs(adj, match_u, match_v, dist, u))
                    M++; // found an augpath
            }
            progress(M, N, "Maximum matching");
        }
        finished("Maximum matching");
        cout << "\nFound maximum matching of size " << M << "\n\n";

        return M;
    }

    VINT dfs_forward(const VINT& match_u, const INT& start, Vint& seen_u, Vint& seen_v) {
        VINT trl;
        INT crnt = start;

        while(1) {
            if (seen_u[crnt]) break;
            seen_u[crnt] = 1;
            trl.emplace_back(crnt);

            auto next = match_u[crnt];
            if (next == INF) break;
            seen_v[next] = 1;
            crnt = next;
        }
        return trl;
    }

    VINT dfs_backward(const VINT& match_v, const INT& start, Vint& seen_u, Vint& seen_v) {
        VINT lrt;
        INT crnt = start;

        while(1) {
            if (seen_v[crnt]) break;
            seen_v[crnt] = 1;
            lrt.emplace_back(crnt);

            auto prev = match_v[crnt];
            if (prev == INF) break;
            seen_u[prev] = 1;
            crnt = prev;
        }

        return lrt;
    }

    PINT decompose(const INT& N, const VINT& match_u, const VINT& match_v,
                VVINT& cycles, VVINT& paths) {
        Vint seen_u(N, 0);
        Vint seen_v(N, 0);

        for (INT u = 0; u < N; ++u) {
            if (seen_u[u]) continue;
            VINT trl = dfs_forward(match_u, u, seen_u, seen_v);
            if (trl.size() > 1 && match_u[trl.back()] == trl.front())
                cycles.emplace_back(trl);
            else {
                VINT lrt = dfs_backward(match_v, trl.front(), seen_u, seen_v);
                if (!lrt.empty()) {
                    reverse(lrt.begin(), lrt.end());
                    lrt.pop_back();
                    trl.insert(trl.begin(), lrt.begin(), lrt.end());
                }
                paths.emplace_back(trl);
            }
            progress(u, N, "Decomposing");
        }
        finished("Decomposing");
        if (has_dup(cycles, paths))
            cout << "\nDuplicates found.\n";
        else
            cout << "\nNo duplicates found.\n";
        INT C = cycles.size(), P = paths.size();
        cout << "Found " << C << " cycles and " << P << " paths.\n\n";

        return {C, P};
    }

    bool has_dup(const VVINT& cycles, const VVINT& paths) {
        USET seen;
        auto check = [&](const VVINT& v) {
            for (const auto& vec: v)
                for (const auto& x: vec)
                    if (!seen.insert(x).second)
                        return true;
            return false;
        };
        return check(cycles) || check(paths);
    }

    REP plain(const VINT& kmerv, const VVINT& cycles, const VVINT& paths) {
        string txt;
        VINT pnt;
        cout << "Generating SPSS...";
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                auto c = decode_base(kmerv[node] % 4);
                txt += c;
            } txt += "$";
        } txt += "$";
        for (const auto& path: paths) {
            for (const auto& node: path) {
                auto hash = kmerv[node];
                if (node == path.front())
                    txt += decode_kmer(hash, K).substr(0, K - 1);
                txt += decode_base(hash % 4);
            } txt += "$";
        } if (!txt.empty()) txt.pop_back();
        cout << "\rCompleted generating SPSS.\n";
        return {txt, {}};
    }

    REP unsorted(UMAP& kmers, const INT& N, const VINT& kmerv,
                const VVINT& inv_adj, const VVINT& cycles, const VVINT& paths, const PINT& CnP) {
        INT P = CnP.second;  
        INT S = CnP.first + P + N;
        USET self_paths;
        for (INT i = 0; i < P; ++i)
            if (inv_adj[paths[i][0]].empty())
                self_paths.insert(i);
            else heads[paths[i][0]] = i; // record paths' heads
        INT pos = 0;
        VINT pnt(P, -1);
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                for (const auto& c: base) {
                    auto next = forward(kmers, kmerv, node, c);
                    if (next == INF) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pnt[heads[next]] = pos;
                    heads.erase(it);
                    if (heads.empty()) goto end;
                } pos++;
                progress(pos, S, "Pointing");
            } pos++;
        }
        for (const auto& path: paths) {
            for (const auto& node: path) {
                for (const auto& c: base) {
                    auto next = forward(kmers, kmerv, node, c);
                    if (next == INF) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pnt[heads[next]] = pos;
                    heads.erase(it);
                    if (heads.empty()) goto end;
                } pos++;
                progress(pos, S, "Pointing");
            } pos++;
        }
        for (const auto& pid: self_paths) {
            pnt[pid] = pos;
            pos += paths[pid].size() + 1;
        }
        finished("Pointing");
        end:
        if (heads.empty()) 
            cout << "\nAll paths are pointed.\n";
        else cout << "\n" << heads.size() << " paths are not pointed.\n";
        // generate representation
        string txt;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                auto c = decode_base(kmerv[node] % 4);
                txt += c;
            } txt += "$";
        } for (INT i = 0; i < P; ++i) {
            auto path = paths[i];
            if (self_paths.find(i) != self_paths.end())
                txt += "$" + decode_kmer(kmerv[path[0]], K).substr(0, K - 1);
            for (const auto& node: path) {
                auto c = decode_base(kmerv[node] % 4);
                txt += c;
            } txt += "$";
        } if(!txt.empty()) txt.pop_back();

        return {txt, pnt};
    }

    bool has_pointer_cycle(
        INT current, VVINT& paths, UMAP& kmers, const VINT& kmerv, VINT& new_cycle,
        Vint& visited
    ) {
        if (visited[current]) return true;
        visited[current] = 1;
        auto& path = paths[current];

        for (INT i = 0; i < (INT)(path.size()); ++i) {
            auto node = path[i];
            for (const auto& c: base) {
                auto next = forward(kmers, kmerv, node, c);
                if (next == INF) continue;

                auto it = heads.find(next);
                if (it != heads.end()) {
                    INT next_path = it->second;
                    if (next_path == current) {
                        new_cycle.assign(path.begin(), path.begin() + i + 1);
                        return true;
                    }
                    if (has_pointer_cycle(next_path, paths, kmers, kmerv, new_cycle, visited)) {
                        new_cycle.insert(new_cycle.begin(), path.begin(), path.begin() + i + 1);
                        return true;
                    }
                }
            }
        }
        return false; // no pointer cycle found
    }

    REP sorted(UMAP& kmers, const INT& N, const VINT& kmerv,
                VVINT& inv_adj, VVINT& cycles, VVINT& paths, const PINT& CnP) {
        INT P = CnP.second;
        INT S = CnP.first + P + N, ofst = 0, from = 0;
        USET self_paths;
        for (INT i = 0; i < P; ++i) {
            if (inv_adj[paths[i][0]].empty()){
                self_paths.insert(i);
            }
            else heads[paths[i][0]] = i; // record paths' heads
        }

        cout << "\n# non-self paths: " << heads.size()
             << "\n# self paths: " << self_paths.size() << "\n";

        VINT pord; // sorted path ID order
        VINT pntc, pntp; // record the relative positions of pointees

        // add pointers from paths to cycles
        INT pos = 0, dist = 0, exit_code = -1;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                for (const auto& c: base) {
                    auto next = forward(kmers, kmerv, node, c);
                    if (next == INF) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pord.emplace_back(it->second);
                    pntc.emplace_back(dist);
                    dist = 0;
                    heads.erase(it);
                    if ((INT)pord.size() == P)
                        {exit_code = 0; goto end;}
                } ++pos; ++dist;
                progress(pos, S, "Pointing");
            } ++pos; ++dist;
        }
        ofst = dist; // distance from pos of the last pointee in cycles and the end
        cout << "ofst(right after finishing scanning cycles): " << ofst << "\n";
        dist = 0;

        if (!self_paths.empty()) {
            for (const auto& pid: self_paths) {
                pord.emplace_back(pid);
                pntp.emplace_back(dist);
                dist = K;
                for (const auto& node: paths[pid]) {
                    for (const auto& c: base) {
                        auto next = forward(kmers, kmerv, node, c);
                        if (next == INF) continue;
                        auto it = heads.find(next);
                        if (it == heads.end()) continue;
                        pord.emplace_back(it->second);
                        pntp.emplace_back(dist);
                        dist = 0;
                        heads.erase(it);
                        if ((INT)pord.size() == P) 
                            {exit_code = 3; goto end;}
                    } ++pos; ++dist;
                    progress(pos, S, "Pointing");
                } ++pos; ++dist;
            }
        }
        
        while ((INT)(pord.size()) < P) {         
            INT bound = (INT)pord.size();

            for (INT i = from; i < bound; ++i) {
                for (const auto& node: paths[pord[i]]) {
                    for (const auto& c: base) {
                        auto next = forward(kmers, kmerv, node, c);
                        if (next == INF) continue;
                        auto it = heads.find(next);
                        if (it == heads.end()) continue;
                        pord.emplace_back(it->second);
                        pntp.emplace_back(dist);
                        dist = 0;
                        ++bound;
                        heads.erase(it);
                        if ((INT)pord.size() == P) 
                            {exit_code = 1; goto end;}
                    } ++pos; ++dist;
                    progress(pos, S, "Pointing");
                } ++pos; ++dist;
            } 
            from = bound;

            cout << "\n# remaining paths: " << heads.size() << "\n"
                 << "Resolving pointer cycles...\n";
            int flg = 1;
            for (auto it1 = heads.begin(); it1 != heads.end(); ++it1) {
                VINT new_cycle;
                auto start = it1->second;
                Vint visited(P, 0);

                if (has_pointer_cycle(start, paths, kmers, kmerv, new_cycle, visited)
                    && !new_cycle.empty()) {
                    cycles.emplace_back(new_cycle); // append new_cycle to cycles if found

                    INT R = (INT)(new_cycle.size()), i = 0, distb = 0;
                    while (i < R) {
                        auto it = heads.find(new_cycle[i]);
                        if (it != heads.end()) {
                            distb = 0;
                            auto pid = it->second;
                            heads.erase(it);
                            auto& path = paths[pid];
                            INT j = 0;
                            while (i < R && j < (INT)path.size() && new_cycle[i] == path[j]) {
                                for (const auto& c: base) {
                                    auto next = forward(kmers, kmerv, new_cycle[i], c);
                                    if (next == INF || next == new_cycle[(i + 1) % R]) continue;
                                    auto itit = heads.find(next);
                                    if (itit == heads.end()) continue;
                                    pord.emplace_back(itit->second);
                                    if (flg) {
                                        pntc.emplace_back(ofst + i);
                                        --flg; distb = 0;
                                    } else {
                                        pntc.emplace_back(distb);
                                        distb = 0;
                                    }
                                    heads.erase(itit);
                                }
                                ++i; ++j; ++distb;
                            }
                            pord.emplace_back(pid);
                            if (flg) {
                                pntc.emplace_back(ofst + i - 1);
                                --flg;
                            } else {
                                pntc.emplace_back(distb);
                            }
                            path.erase(path.begin(), path.begin() + j);
                        }
                    }
                    ofst = distb;

                    if ((INT)pord.size() == P) {
                        exit_code = 2;
                        goto end;
                    }
                }
            }
            cout << "\n# remaining paths: " << heads.size() << "\n"
                 << "Resolved\n";
        }
        end:
        finished("Pointing");
        if (heads.empty()) 
            cout << "\nAll paths are pointed.\n";
        else cout << "\n" << heads.size() << " paths are not pointed.\n";
        cout << "\nResult code: " << exit_code << "\n";
        cout << "0: All pointers to cycles\n"
             << "1: At least one pointer to a path. (No pointer cycle)\n"
             << "2: At least one pointer cycle in the initial decomposition.\n"
             << "3: At least one self-path which has no pointer to it.\n"
             << "-1: No paths || Self-pointing paths only || Unexpected behavior (Exception)\n\n";

        INT N_count = 0; // verify the resulting cumulative length of cycles & paths is N
        if (pntp.size()) pntp[0] += ofst; // because there is a delimiter in between
        
        // representation
        string txt;
        VINT pnt;
        pnt.insert(pnt.end(), pntc.begin(), pntc.end());
        pnt.insert(pnt.end(), pntp.begin(), pntp.end());

        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                ++N_count;
                auto c = decode_base(kmerv[node] % 4);
                txt += c;
            } txt += "$";
        } for (const auto& pid: pord) {
            if (self_paths.find(pid) != self_paths.end())
                txt += "$" + decode_kmer(kmerv[paths[pid][0]], K).substr(0, K - 1);
            for (const auto& node: paths[pid]) {
                ++N_count;
                auto c = decode_base(kmerv[node] % 4);
                txt += c;
            } txt += "$";
        } if(!txt.empty()) txt.pop_back();
        
        return {txt, pnt};
    }

    string cp_spell_out(VINT& kmerv, const VINT& walk, bool is_cycle) {
        string s;
        for (auto id: walk) {
            if (id == walk.front() && !is_cycle) {
                s += decode_kmer(kmerv[id], K).substr(0, K - 1);
            }
            s += decode_base(kmerv[id] % 4);
        }
        return s;
    }

    REP forest(UMAP& kmers, const INT& N, const VINT& kmerv,
        VVINT& inv_adj, VVINT& cycles, VVINT& paths, const PINT& CnP) {
        string txt;
        INT P = CnP.second;
        USET roots;
        for (INT i = 0; i < P; ++i) {
            if (inv_adj[paths[i][0]].empty()){
                roots.insert(i);
            }
            else heads[paths[i][0]] = i; // record paths' heads
        }
        Vint embedded(P, 0); // embedded[pid] = 1 iff path pid is already embedded
        // embed_from_cycle(kmers, kmerv, paths, embedded, cycles);
        // embed_from_path(kmerv, paths, roots, embedded);

        return {txt, {}};
    }

    void write(const REP& rep) {
        if (rep.first.empty()) {
            cerr << "Error: txt is empty.\n";
            return;
        }
        string txtfilename = remove_extension(filename, ".fa") + ".ours.k" + to_string(K) + ".opt" + to_string(option) + ".txt";
        ofstream txtfile(txtfilename);
        if (!txtfile) {
            cerr << "Error: Could not open file " << txtfilename << " for writing.\n";
            return;
        }
        txtfile << rep.first;
        txtfile.close();
        cout << "\nFile created: " << path_to_filename(txtfilename, '/');
        fs::path txtfilepath = txtfilename;
        auto txtfilesize = fs::file_size(txtfilepath);
        cout << " (File Size: " << txtfilesize << " B)\n";

        if (rep.second.empty()) {
            cerr << "Note: pnt is empty.\n";
            return;
        }
        string pntfilename = remove_extension(filename, ".fa") + ".ours.k" + to_string(K) + ".opt" + to_string(option) + ".pnt.txt";
        ofstream pntfile(pntfilename);
        if (!pntfile) {
            cerr << "Error: Could not open file " << pntfilename << " for writing.\n";
            return;
        }
        for (const auto& p: rep.second)
            pntfile << p << ",";
        pntfile.close();
        cout << "File created: " << path_to_filename(pntfilename, '/');
        fs::path pntfilepath = pntfilename;
        auto pntfilesize = fs::file_size(pntfilepath);
        cout << " (File Size: " << pntfilesize << " B)\n";

        // record output size on logfile
        ofstream logfile(logfilename, ios::app);
        if (!logfile) {
            cerr << "Error: Could not open file " << logfilename << " for writing.\n";
            return;
        }
        logfile << "\nFile created: " << path_to_filename(txtfilename, '/')
                << " (File Size: " << txtfilesize << " B)\n"
                << "File created: " << path_to_filename(pntfilename, '/')
                << " (File Size: " << pntfilesize << " B)\n";
        logfile.close();
    }
};

void print_usage(string cmdname) {
    cerr << "Usage: " << cmdname << "[0(node-centric) / 1(edge-centric)] [****.fa] [k] ([0 / 1 / 2])\n"
        << "Output options(node-centric) are...\n"
        << "0: explicit representation without pointers (SPSS)\n"
        << "1: representation with unsorted pointers\n"
        << "2: representation with sorted pointers\n"
        << "\nFor edge-centric, set 1st arg = 0\n";
}

int main(int argc, char *argv[]) {
    if (argc < 4 || 5 < argc) {
        print_usage(argv[0]);
        return 1;
    }
    string _filename = argv[2];
    INT _K = atoi(argv[3]);

    // calc size of undirected eulertigs
    if (atoi(argv[1]) == 1) {
        if (argc != 4) {
            print_usage(argv[0]);
            return 1;
        }
        EDBG edbg = EDBG(_filename, _K);
        INT S = edbg.process();
        cout << "Total base count: " << S << "\n";
    }
    // ours
    else if (atoi(argv[1]) == 0) {
        if (argc != 5) {
            print_usage(argv[0]);
            return 1;
        }
        INT _option = atoi(argv[4]);

        NDBG ndbg = NDBG(_filename, _K, _option);
        REP rep = ndbg.process();
        ndbg.write(rep);
    }
    else {
        print_usage(argv[0]);
        return 1;
    }    

    return 0;
}