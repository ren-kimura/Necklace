#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <cstdint>
#include <unistd.h>
#include <memory_resource>
#include <sys/resource.h>
#include <string>
#include <vector>
#include <stack>
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

vector<char> base = {'A', 'C', 'G', 'T'};
const INT INF = UINT64_MAX;

size_t get_available_memory() {
    return sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGE_SIZE);
}

size_t get_peak_memory_kb() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // in kilobytes
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

void log_time(const string& task_name, const chrono::duration<double>& elapsed_time, ofstream& logfile) {
    string task_name_trimmed = task_name.substr(0, 20);

    logfile << left << setw(20) << task_name_trimmed
            << setw(20) << fixed << setprecision(4) << elapsed_time.count()
            << "\n";

    cout << left << setw(20) << task_name_trimmed
         << setw(20) << fixed << setprecision(4) << elapsed_time.count()
         << "\n";
}

void print_table_header(ofstream& logfile) {
    const int col_width = 20;

    logfile << string(30, '-') << "\n";
    logfile << left << setw(col_width) << "Task name"
            << setw(col_width) << "Time (s)"
            << endl;
    logfile << string(30, '-') << "\n"; 

    cout << string(30, '-') << "\n";
    cout << left << setw(col_width) << "Task name"
         << setw(col_width) << "Time (s)"
         << endl;
    cout << string(30, '-') << "\n"; 
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
        return S * (K - 1) + M + (S - 1); // prefix + non-prefix + delim
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
        if (option != 0 && option != 1 && option != 2 && option != 3) {
            cerr << "\nInvalid option value\n";
            exit(1);
        }
        logfilename = remove_extension(filename, ".fa") + ".o" + to_string(K) + "." + to_string(option) + ".log";
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

        VINT kmerv; 
        
        auto start_time = chrono::high_resolution_clock::now();
        INT N = get_kmers(kmers, kmerv);
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> get_kmers_time = end_time - start_time;

        VVINT adj(N), inv_adj(N);

        start_time = chrono::high_resolution_clock::now();
        INT E = add_edges(kmers, kmerv, N, adj, inv_adj);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> add_edges_time = end_time - start_time;

        VINT match_u, match_v;

        start_time = chrono::high_resolution_clock::now();
        INT M = hopcroft_karp(N, adj, match_u, match_v);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> hopcroft_karp_time = end_time - start_time;

        VVINT cycles, paths;

        start_time = chrono::high_resolution_clock::now();
        PINT CnP = decompose(N, match_u, match_v, cycles, paths);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> decompose_time = end_time - start_time;

        cout << "\n(#k-mers, #edges, #matching, #cycles, #paths)\n= (" 
            << N << ", " << E << ", " << M << ", "
            << CnP.first << ", " << CnP.second << ")\n";

        REP rep;
        chrono::duration<double> align_time;
        if (option == 0) {
            // plain
            start_time = chrono::high_resolution_clock::now();
            rep = plain(kmerv, cycles, paths);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else if (option == 1) {
            // unsorted
            start_time = chrono::high_resolution_clock::now();
            rep = unsorted(kmers, N, kmerv, inv_adj, cycles, paths, CnP);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else if (option == 2) {
            // sorted
            start_time = chrono::high_resolution_clock::now();
            rep = sorted(kmers, N, kmerv, inv_adj, cycles, paths, CnP);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else if (option == 3) {
            // tree (BP)
            start_time = chrono::high_resolution_clock::now();
            rep = forest(kmers, kmerv, inv_adj, cycles, paths, CnP);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else {
            cerr << "Invalid option !\n";
            exit(1);
        }

        // display benchmarks
        print_table_header(logfile);
        log_time("get_kmers", get_kmers_time, logfile);
        log_time("add_edges", add_edges_time, logfile);
        log_time("hopcroft_karp", hopcroft_karp_time, logfile);
        log_time("decompose", decompose_time, logfile);
        log_time((option == 0 ? "plain" : option == 1 ? "unsorted" : option == 2 ? "sorted" : "tree (BP)"),
                            align_time, logfile);
        logfile << string(30, '-');
        cout << string(30, '-');

        auto total_time = get_kmers_time.count() + add_edges_time.count()
                          + hopcroft_karp_time.count() + decompose_time.count()
                          + align_time.count();
        logfile << "\nTotal time: " << total_time << " s\n";
        cout << "\nTotal time: " << total_time << " s\n";

        size_t peak_memory_kb = get_peak_memory_kb();
        logfile << "Peak memory: " << peak_memory_kb << " KB\n";
        cout << "Peak memory: " << peak_memory_kb << " KB\n";
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

        // Error handling
        Vint seen(N, 0); INT base_count = 0;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                if (seen[node]) {
                    cerr << "\nNODE" << node << " IS DUPLICATED\n";
                    exit(1);
                }
                seen[node] = 1;
                ++base_count;
            }
        }
        for (const auto& path: paths) {
            for (const auto& node: path) {
                if (seen[node]) {
                    cerr << "\nNODE" << node << " IS DUPLICATED\n";
                    exit(1);
                }
                seen[node] = 1;
                ++base_count;
            }
        }
        if (base_count != N) {
            cerr << "\nWRONG NUMBER OF BASES IN CYCLES AND PATHS\n";
            exit(1);
        }
        cout << "\nNo duplicates found. Total number of bases is correct.\n";

        INT C = cycles.size(), P = paths.size();
        cout << "Found " << C << " cycles and " << P << " paths.\n";

        return {C, P};
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
        } 
        if (!txt.empty()) txt.pop_back();
        if (paths.empty()) txt.pop_back();
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

    VINT find_new_cycle(INT start, VVINT& paths, Vint& in_pord, UMAP& kmers, const VINT& kmerv, Vint& visited) {
        struct Frame {
            INT current;
            INT idx;
            INT b_idx;
        };
        vector<Frame> stack;
        VINT cycle;

        visited[start] = 1;
        stack.push_back({start, 0, 0});

        while (!stack.empty()) {
            Frame &top = stack.back();
            const VINT &cur_path = paths[top.current];

            // rm from stack if explored all the nodes in cur_path
            if (top.idx >= (INT)cur_path.size()) {
                stack.pop_back();
                if (!cycle.empty()) cycle.pop_back();
                continue;
            }

            // when entering new node (on start of this frame)
            if (top.idx == 0 && top.b_idx == 0) cycle.push_back(top.current);

            // try next base
            if (top.b_idx < 4) {
                char c = base[top.b_idx];
                ++top.b_idx;
                
                INT cur_node = cur_path[top.idx];
                INT next_node = forward(kmers, kmerv, cur_node, c);
                if (next_node == INF) continue;
                auto it = heads.find(next_node);
                if (it == heads.end()) continue;
                INT next_pid = it->second;
                if (in_pord[next_pid]) continue;
                if (next_pid == start) return cycle;
                // skip visited path
                if (visited[next_pid]) continue;

                visited[next_pid] = 1;
                stack.push_back({next_pid, 0, 0});
            } else {
                // already tried all bases. proceed to next node
                ++top.idx;
                top.b_idx = 0;
            }
        }
        return {};
    }

    REP sorted(UMAP& kmers, const INT& N, const VINT& kmerv,
                VVINT& inv_adj, VVINT& cycles, VVINT& paths, const PINT& CnP) {
        INT C = CnP.first, P = CnP.second, S = C + P + N, ncycs = 0;
        USET roots;
        for (INT i = 0; i < P; ++i) {
            if (inv_adj[paths[i][0]].empty()) roots.insert(i);
            else heads[paths[i][0]] = i;
        }
        cout << "\n(#non-root paths, #root paths) = (" << heads.size() << ", "
             << roots.size() << ")\n";

        queue<INT> q;
        vector<pair<uint8_t, INT>> ord; // sorted pairs (0 or 1, id), where 0 is cycle, 1 is path
        Vint in_ord(P, 0); // 0: path is in ord, 1: not in ord yet
        VINT pnt;

        INT pos = 0;
        // bfs from root paths
        for (const auto& pid: roots) {
            q.push(pid);
            ord.emplace_back(1, pid);
            in_ord[pid] = 1;
            pos += K; // extra "$" and K-1 letters before a root path

            while (!q.empty()) {
                auto cpid = q.front();
                q.pop();
                for (const auto& node: paths[cpid]) {
                    for (const auto& c: base) {
                        auto next = forward(kmers, kmerv, node, c);
                        if (next == INF) continue;
                        auto it = heads.find(next);
                        if (it == heads.end()) continue;
                        auto npid = it->second;
                        if (in_ord[npid]) continue;
                        q.push(npid);
                        ord.emplace_back(1, npid);
                        in_ord[npid] = 1;
                        pnt.emplace_back(pos);
                    }
                    ++pos;
                    progress(pos, S, "Pointing");
                }
                ++pos;
            }
        }
        // bfs from cycles
        for (INT i = 0; i < C; ++i) {
            auto cycle = cycles[i];
            ord.emplace_back(0, i);
            ++pos; // need flag "*" for a cycle
            for (const auto& node: cycle) {
                for (const auto& c: base) {
                    auto next = forward(kmers, kmerv, node, c);
                    if (next == INF) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    auto npid = it->second;
                    if (in_ord[npid]) continue;
                    q.push(npid);
                    ord.emplace_back(1, npid);
                    in_ord[npid] = 1;
                    pnt.emplace_back(pos);                    
                }
                ++pos;
                progress(pos, S, "Pointing");
            }
            ++pos;
            while (!q.empty()) {
                auto cpid = q.front();
                q.pop();
                for (const auto& node: paths[cpid]) {
                    for (const auto& c: base) {
                        auto next = forward(kmers, kmerv, node, c);
                        if (next == INF) continue;
                        auto it = heads.find(next);
                        if (it == heads.end()) continue;
                        auto npid = it->second;
                        if (in_ord[npid]) continue;
                        q.push(npid);
                        ord.emplace_back(1, npid);
                        in_ord[npid] = 1;
                        pnt.emplace_back(pos);
                    }
                    ++pos;
                    progress(pos, S, "Pointing");
                }
                ++pos;
            }
        }

        // find a new cycle and bfs from it
        while (INT res = count(in_ord.begin(), in_ord.end(), 0)) {
            cout << "\n\n#residue paths = " << res << "\n"
                 << "Searching for a new cycle...";
            for (INT i = 0; i < P; ++i) {
                if (in_ord[i]) continue;
                Vint vis(P, 0);
                VINT pids = find_new_cycle(i, paths, in_ord, kmers, kmerv, vis);
                if (!pids.empty()) {
                    cout << "\rFound a cycle of " << pids.size() << " paths          \n";
                    ++ncycs;
                    VINT cycle;
                    INT p = (INT)pids.size();
                    for (INT j = 0; j < p; ++j) {
                        auto& path = paths[pids[j]];
                        auto h = (INT)path.size();
                        for (INT k = 0; k < h; ++k) {
                            auto node = path[k];
                            cycle.emplace_back(node);
                            for (auto& c: base) {
                                auto next = forward(kmers, kmerv, node, c);
                                if (next == INF) continue;
                                if (next == paths[pids[(j + 1) % p]][0]) {
                                    // trim cycle portion from path, and renew heads
                                    heads.erase(path[0]);
                                    path.erase(path.begin(), path.begin() + k + 1);
                                    heads[path[0]] = pids[j];
                                    goto endj;
                                }
                            }
                        }
                        endj:;
                    }
                    cycles.emplace_back(cycle);
                    ord.emplace_back(0, (INT)cycles.size() - 1);
                    ++pos; // need flag "*" for a cycle
                    for (const auto& node: cycle) {
                        for (const auto& c: base) {
                            auto next = forward(kmers, kmerv, node, c);
                            if (next == INF) continue;
                            auto it = heads.find(next);
                            if (it == heads.end()) continue;
                            auto npid = it->second;
                            if (in_ord[npid]) continue;
                            q.push(npid);
                            ord.emplace_back(1, npid);
                            in_ord[npid] = 1;
                            pnt.emplace_back(pos);                    
                        }
                        ++pos;
                        progress(pos, S, "Pointing");
                    }
                    ++pos;
                    while (!q.empty()) {
                        auto cpid = q.front();
                        q.pop();
                        for (const auto& node: paths[cpid]) {
                            for (const auto& c: base) {
                                auto next = forward(kmers, kmerv, node, c);
                                if (next == INF) continue;
                                auto it = heads.find(next);
                                if (it == heads.end()) continue;
                                auto npid = it->second;
                                if (in_ord[npid]) continue;
                                q.push(npid);
                                ord.emplace_back(1, npid);
                                in_ord[npid] = 1;
                                pnt.emplace_back(pos);
                            }
                            ++pos;
                            progress(pos, S, "Pointing");
                        }
                        ++pos;
                    }
                    break;
                }
            }
        }
        finished("Pointing");
        
        // for check
        INT res = count(in_ord.begin(), in_ord.end(), 0);
        if (res == 0) cout << "\n\nAll paths are pointed.\n\n";
        else {
            cout << "\n\n" << res << " paths are not pointed.\n";
            for (INT i = 0; i < P; ++i) if (!in_ord[i])
                cout << "path" << i << " not added\n";
            cout << "\n\n";
        }
        cout << "#new cycles = " << ncycs << "\n";

        // take difference of pnt
        for (INT i = (INT)pnt.size() - 1; i > 0; --i) {
            pnt[i] -= pnt[i - 1];
        }

        // generate text rep
        string txt;
        for (const auto& e: ord) {
            auto is_path = e.first;
            auto id = e.second;
            if (is_path) {
                auto path = paths[id];
                if (roots.find(id) != roots.end())
                    txt += "$" + decode_kmer(kmerv[path[0]], K).substr(0, K - 1);
                for (const auto& node: path)
                    txt += decode_base(kmerv[node] % 4);
            } else {
                auto cycle = cycles[id];
                txt += "*";
                for (const auto& node: cycle)
                    txt += decode_base(kmerv[node] % 4);
            }
            txt += "$";
        }
        if (!txt.empty()) txt.pop_back();

        return {txt, pnt};
    }

    string write_subtree(INT pid, VVINT& paths, UMAP& kmers, const VINT& kmerv, 
                         UMAP& heads, Vint& visited) {
        if (visited[pid]) return "";
        visited[pid] = 1;
        auto path = paths[pid];

        string s;        
        for (const auto& node: path) {
            s += decode_base(kmerv[node] % 4);
            for (const auto& c: base) {
                auto next = forward(kmers, kmerv, node, c);
                if (next == INF) continue;
                auto it = heads.find(next);
                if (it == heads.end()) continue;
                auto npid = it->second;
                if (visited[npid]) continue;
                s += "(" + write_subtree(npid, paths, kmers, kmerv, heads, visited) + ")";
            }
        }
        return s;
    }

    REP forest(UMAP& kmers, const VINT& kmerv,
        VVINT& inv_adj, VVINT& cycles, VVINT& paths, const PINT& CnP) {
        string txt;
        INT P = CnP.second, S = CnP.first + P, n = 0, ncycs = 0;
        USET roots;        
        for (INT i = 0; i < P; ++i) {
            if (inv_adj[paths[i][0]].empty()) roots.insert(i);
            else heads[paths[i][0]] = i;
        }
        cout << "\n(#non-root paths, #root paths) = (" << heads.size() << ", "
             << roots.size() << ")\n";

        Vint visited(P, 0); // 1 iff path already embedded

        for (const auto& pid: roots) {
            txt += decode_kmer(kmerv[paths[pid][0]], K).substr(0, K - 1)
                   + write_subtree(pid, paths, kmers, kmerv, heads, visited)
                   + ",";
            ++n; progress(n, S, "Finding pseudoforest");
        }

        for (const auto& cycle: cycles) {
            txt += "*";
            for (const auto& node: cycle) {
                txt += decode_base(kmerv[node] % 4);
                for (const auto& c: base) {
                    auto next = forward(kmers, kmerv, node, c);
                    if (next == INF) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    auto npid = it->second;
                    if (visited[npid]) continue;
                    txt += "(" + write_subtree(npid, paths, kmers, kmerv, heads, visited) + ")";
                }
            }
            txt += ","; ++n; progress(n, S, "Finding pseudoforest");
        }
        cout << "\n";

        while (INT res = count(visited.begin(), visited.end(), 0)) {
            cout << "\r#residue paths = " << res << ", searching for a new cycle..." << flush;
            for (INT i = 0; i < P; ++i) {
                if (visited[i]) continue;
                Vint vis(P, 0);
                VINT pids = find_new_cycle(i, paths, visited, kmers, kmerv, vis);
                if (!pids.empty()) {
                    cout << "\r#residue paths = " << count(visited.begin(), visited.end(), 0) << ", found a cycle of " << pids.size() << " paths          " << flush;
                    ++ncycs;
                    VINT cycle;
                    INT p = (INT)pids.size();
                    for (INT j = 0; j < p; ++j) {
                        auto& path = paths[pids[j]];
                        auto h = (INT)path.size();
                        for (INT k = 0; k < h; ++k) {
                            auto node = path[k];
                            cycle.emplace_back(node);
                            for (auto& c: base) {
                                auto next = forward(kmers, kmerv, node, c);
                                if (next == INF) continue;
                                if (next == paths[pids[(j + 1) % p]][0]) {
                                    // trim cycle portion from path, and renew heads
                                    heads.erase(path[0]);
                                    path.erase(path.begin(), path.begin() + k + 1);
                                    heads[path[0]] = pids[j];
                                    goto endj;
                                }
                            }
                        }
                        endj:;
                    }
                    cycles.emplace_back(cycle);
                    txt += "*";
                    for (const auto& node: cycle) {
                        for (const auto& c: base) {
                            auto next = forward(kmers, kmerv, node, c);
                            if (next == INF) continue;
                            auto it = heads.find(next);
                            if (it == heads.end()) continue;
                            auto npid = it->second;
                            if (visited[npid]) continue;
                            txt += "(" + write_subtree(npid, paths, kmers, kmerv, heads, visited) + ")";
                        }
                    }
                    txt += ","; ++n; progress(n, S, "Finding pseudoforest");
                    break;
                }
            }
            // debug
            if ((INT)count(visited.begin(), visited.end(), 0) == res) {
                cerr << "ERROR: no new cycle found in remaining paths\n";
                exit(1);
            }
        }
        finished("Finding pseudoforest");
        cout << "\n";

        if (!txt.empty()) txt.pop_back();
        return {txt, {}};
    }

    void write(const REP& rep) {
        // --- Write .txt file ---
        if (rep.first.empty()) {
            cerr << "Note: txt is empty.\n";
        }
    
        string txtfilename = remove_extension(filename, ".fa") + ".o" + to_string(K) + "." + to_string(option) + ".txt";
        ofstream txtfile(txtfilename);
        if (!txtfile) {
            cerr << "Error: Could not open file " << txtfilename << " for writing.\n";
            return;
        }
        txtfile << rep.first;
        txtfile.close();
    
        fs::path txtfilepath = txtfilename;
        auto txtfilesize = fs::file_size(txtfilepath);
        cout << "\nFile created: " << path_to_filename(txtfilename, '/')
             << " (File Size: " << txtfilesize << " B)\n";
    
        // --- Write .pnt.txt only if rep.second is non-empty ---
        string pntfilename;
        size_t pntfilesize = 0;
        bool wrote_pnt = false;
    
        if (rep.second.empty()) {
            cerr << "Note: pnt is empty.\n";
        } else {
            pntfilename = remove_extension(filename, ".fa") + ".o" + to_string(K) + "." + to_string(option) + "p.txt";
            ofstream pntfile(pntfilename);
            if (!pntfile) {
                cerr << "Error: Could not open file " << pntfilename << " for writing.\n";
                return;
            }
            for (const auto& p : rep.second)
                pntfile << p << ",";
            pntfile.close();
    
            fs::path pntfilepath = pntfilename;
            pntfilesize = fs::file_size(pntfilepath);
            cout << "File created: " << path_to_filename(pntfilename, '/')
                 << " (File Size: " << pntfilesize << " B)\n";
    
            wrote_pnt = true;
        }
    
        // --- Log file ---
        ofstream logfile(logfilename, ios::app);
        if (!logfile) {
            cerr << "Error: Could not open file " << logfilename << " for writing.\n";
            return;
        }
        logfile << "\nFile created: " << path_to_filename(txtfilename, '/')
                << " (File Size: " << txtfilesize << " B)\n";
        if (wrote_pnt) {
            logfile << "File created: " << path_to_filename(pntfilename, '/')
                    << " (File Size: " << pntfilesize << " B)\n";
        }
        logfile.close();
    }    
};

void print_usage(string cmdname) {
    cerr << "Usage: " << cmdname << "[0(node-centric) / 1(edge-centric)] [****.fa] [k] ([0 / 1 / 2 / 3])\n"
        << "Output options(node-centric) are...\n"
        << "0: explicit representation without pointers (SPSS)\n"
        << "1: representation with unsorted pointers\n"
        << "2: representation with sorted pointers\n"
        << "3: tree representation (BP)\n"
        << "\nFor edge-centric, set 1st arg = 0\n";
}

int main(int argc, char *argv[]) {
    if (argc < 4 || 5 < argc) {
        print_usage(argv[0]);
        return 1;
    }
    
    string _filename = argv[2];
    filesystem::path filepath(_filename);
    string ext = filepath.extension().string();
    if (ext.empty()) {
        _filename += ".fa";
        cout << "No extension detected. Appending '.fa': " << _filename << "\n";
    } else if (ext != ".fa" && ext != ".fna" && ext != ".fasta") {
        cerr << "Invalid extension (must be \".fa\", \".fna\" or \".fasta\")\n\n";
        print_usage(argv[0]);
        return 1;
    }

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