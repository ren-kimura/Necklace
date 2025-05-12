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
using VVINT = vector<VINT>;
using REP = pair<string, VINT>;
using USET = unordered_set<INT>;

vector<char> base = {'A', 'C', 'G', 'T'};
const INT INF = UINT64_MAX;

size_t get_available_memory() {
    std::ifstream meminfo("/proc/meminfo");
    std::string line;
    size_t available_kb = 0;

    while (std::getline(meminfo, line)) {
        if (line.find("MemAvailable:") == 0) {
            std::istringstream iss(line);
            std::string key, unit;
            iss >> key >> available_kb >> unit;
            break;
        }
    }

    return available_kb * 1024; // bytes
}

size_t estimate_pool_size_from_file(const string& filename, INT K) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file) {
        cerr << "Error: Cannot open file " << filename << "\n";
        exit(1);
    }

    const size_t BUF_SIZE = 64 * 1024 * 1024;
    vector<char> buffer(BUF_SIZE);
    file.rdbuf()->pubsetbuf(buffer.data(), BUF_SIZE);

    string line, read;
    size_t base_count = 0;
    bool is_new_read = true;

    while (getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!read.empty()) {
                base_count += read.size();
                read.clear();
            }
            is_new_read = true;
        } else {
            if (is_new_read) {
                read = line;
                is_new_read = false;
            } else {
                read += line;
            }
        }
    }
    if (!read.empty()) base_count += read.size();
    file.close();

    size_t kmer_upper_bound = static_cast<size_t>(1ULL << (2 * K)); // 4^K
    base_count = std::min(base_count, kmer_upper_bound);

    const double RECORD_SIZE = 48.0;
    const double SAFETY_MARGIN = 2.0;
    const size_t MIN_POOL = 64L * 1024 * 1024;

    size_t estimate = static_cast<size_t>(base_count * RECORD_SIZE * SAFETY_MARGIN);

    size_t available_memory = get_available_memory();
    size_t dynamic_max_pool = static_cast<size_t>(available_memory * 0.7);

    estimate = std::clamp(estimate, MIN_POOL, dynamic_max_pool);

    return estimate;
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
    cout << ">] 100%\n";
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
    string filename, logfilename;
    INT K;    
    size_t pool_size;
    size_t reserve_size;
    
    pmr::monotonic_buffer_resource pool;
    pmr::unordered_set<INT> kmers;
    pmr::unordered_map<INT, short> adj;
    pmr::unordered_map<INT, int> outcnt;
    pmr::unordered_map<INT, VINT> dummy_map;

    struct pair_hash {
        size_t operator()(pair<INT, INT> const &p) const {
            return hash<INT>()(p.first) ^ (hash<INT>()(p.second) << 1);
        }
    };
    pmr::unordered_multiset<pair<INT, INT>, pair_hash> used_dummy_set;

    EDBG(string _filename, INT _K)
        : filename(_filename), K(_K),
          pool_size(estimate_pool_size_from_file(_filename, _K)),
          reserve_size(pool_size / 32),
          pool(pool_size),
          adj(&pool), outcnt(&pool), dummy_map(&pool), used_dummy_set(&pool)
    {
        logfilename = remove_extension(filename, ".fa") + ".ue" + to_string(K) + ".log";
        kmers = pmr::unordered_set<INT>(&pool);
        kmers.reserve(reserve_size);
        cout << "\nUsing pool size: " << pool_size / (1024 * 1024) << " MB\n";
        cout << "K-mers reserve size: " << reserve_size << "\n";
    };

    string process() {
        ofstream logfile(logfilename, ios::trunc);
        if (!logfile) {
            cerr << "Error: Could not open " << logfilename << " for writing\n";
            exit(1);
        }
        string timestamp = get_current_timestamp();
        logfile << "Timestamp: " << timestamp << "\n";
        cout << "Timestamp: " << timestamp << "\n";
        
        auto start_time = chrono::high_resolution_clock::now();
        get_kmers();
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> get_kmers_time = end_time - start_time;
        cout << "Collected " << kmers.size() << " unique k-mers\n\n";

        adj.reserve(kmers.size());
        outcnt.reserve(kmers.size());

        start_time = chrono::high_resolution_clock::now();
        build_adj();
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> build_adj_time = end_time - start_time;

        dummy_map.reserve(kmers.size());

        start_time = chrono::high_resolution_clock::now();
        balance_graph();
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> balance_graph_time = end_time - start_time;

        start_time = chrono::high_resolution_clock::now();
        VINT tour = find_eulerian_cycle();
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> find_eulerian_cycle_time = end_time - start_time;

        start_time = chrono::high_resolution_clock::now();
        string txt = spell_out(tour);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> spell_out_time = end_time - start_time;

        print_table_header(logfile);
        log_time("get_kmers", get_kmers_time, logfile);
        log_time("build_adj", build_adj_time, logfile);
        log_time("balance_graph", balance_graph_time, logfile);
        log_time("find_eulerian_cycle", find_eulerian_cycle_time, logfile);
        log_time("spell_out", spell_out_time, logfile);
        logfile << string(30, '-');
        cout << string(30, '-');
        auto total_time = get_kmers_time.count() + build_adj_time.count()
                        + balance_graph_time.count() + find_eulerian_cycle_time.count()
                        + spell_out_time.count();
        logfile << "\nTotal time: " << total_time << " s\n";
        cout << "\nTotal time: " << total_time << " s\n";

        size_t peak_memory_kb = get_peak_memory_kb();
        logfile << "Peak memory: " << peak_memory_kb << " KB\n";
        cout << "Peak memory: " << peak_memory_kb << " KB\n";
        logfile.close();

        return txt;
    }

    void get_kmers() {
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
    }

    void build_adj() {
        INT mask = (1ULL << (2 * (K - 1))) - 1;
        INT i = 0, N = kmers.size();

        for (const auto& kmer: kmers) {
            INT pre = kmer >> 2;        // prefix
            INT suf = kmer & mask;      // suffix
            int last = kmer & 3;
            int first = (kmer >> (2 * (K - 1))) & 3;

            adj[pre] |= (1 << (last + 4));  // out: ms4-bits
            adj[suf] |= (1 << first);       // in:  ls4-bits
            ++outcnt[pre];

            progress(++i, N, "Building EDBG");
        }
        finished("Building EDBG"); cout << "\n";
    }

    void balance_graph() {
        VINT sources, sinks;
        INT i = 0, N = adj.size();

        for (const auto& [node, bits]: adj) {
            int odeg = 0, ideg = 0;
            for (int i = 0; i < 4; ++i) {
                if (bits & (1 << (i + 4)))  ++odeg;
                if (bits & (1 << i))        ++ideg;
            }
            if (odeg > ideg) {
                int diff = odeg - ideg;
                for (int j = 0; j < diff; ++j) sources.emplace_back(node);
            } else if (ideg > odeg) {
                int diff = ideg - odeg;
                for (int j = 0; j < diff; ++j) sinks.emplace_back(node);
            }
            progress(++i, N, "Counting imbalances");
        }
        finished("Counting imbalances"); cout << "\n";

        if (sources.size() != sinks.size()) {
            cerr << "ERROR: mismatch between sources and sinks!\n";
            exit(1);
        }
        auto I = sources.size();
        for (size_t i = 0; i < I; ++i) {
            INT u = sinks[i];
            INT v = sources[i];
            ++outcnt[u];
            dummy_map[u].push_back(v);
            progress(i, I, "Adding dummy edges");
        }
        finished("Adding dummy edges");
        cout << "Added " << sources.size() << " dummy edges to balance the graph\n\n";
    }

    VINT find_eulerian_cycle() {
        size_t m = 0, M = 0; // M: total moves
        for (auto& [node, cnt]: outcnt) M += cnt;
        auto dmap = dummy_map;

        INT start = INF;
        if (!dmap.empty()) {
            start = dmap.begin()->first;
        } else {
            for (const auto& [node, cnt]: outcnt) {
                if (cnt > 0) { start = node; break; }
            }
        }        
        if (start == INF) {
            cerr << "ERROR: no start node found\n";
            exit(1);
        }
        VINT tour, curpath;
        curpath.push_back(start);
        
        while(!curpath.empty()) {
            INT u = curpath.back();
            if (outcnt[u] > 0) {
                bool moved = false;
                auto it_map = dmap.find(u);
                if (it_map != dmap.end() && !it_map->second.empty()) {
                    INT v = it_map->second.back();
                    it_map->second.pop_back();
                    curpath.push_back(v);
                    --outcnt[u];
                    used_dummy_set.insert({u, v});
                    moved = true;
                }
                if (moved) {
                    progress(++m, M, "Finding Eulerian tour");
                    continue;
                }

                for (int b = 0; b < 4; ++b) {
                    if (adj[u] & (1 << (b + 4))) {
                        INT mask = (1ULL << (2 * (K - 2))) - 1;
                        INT v = ((u & mask) << 2) | b;
                        curpath.push_back(v);
                        --outcnt[u]; // used one oedge
                        adj[u] &= ~(1 << (b + 4)); // erase used edge
                        moved = true;
                        progress(++m, M, "Finding Eulerian tour");
                        break;
                    }
                }
                if (!moved) {
                    cerr << "ERROR: stuck at "
                         << decode_kmer(u, K-1)
                         << " outcnt=" << outcnt[u]
                         << " remaining dummies_for_u=";
                    auto it = dmap.find(u);
                    if (it != dmap.end())   cerr << it->second.size();
                    else                    cerr << 0;
                    cerr << "\n";
                    exit(1);
                }
            } else {
                tour.push_back(u);
                curpath.pop_back();
            }
        }
        finished("Finding Eulerian tour");

        reverse(tour.begin(), tour.end());
        cout << "Eulerian cycle found: " << tour.size() << " nodes\n\n";

        return tour;
    }

    string spell_out(const VINT& tour) {
        string txt;
        size_t n = tour.size();
        if (n < 2) { cerr << "Tour too short\n"; exit(1); }

        for (size_t i = 0; i < n - 1; ++i) {
            pair<INT, INT> edge{tour[i], tour[i + 1]};
            auto it = used_dummy_set.find(edge);
            bool is_dummy = (it != used_dummy_set.end());
            if(is_dummy) {
                used_dummy_set.erase(it);
                if (i > 0) txt += ",";
                txt += decode_kmer(edge.second, K - 1);
            } else {
                txt += decode_base(edge.second & 3);
            }
        }
        if (!txt.empty() && txt.back() == ',') {
            txt.pop_back();
        }
        return txt;
    }

    void write(const string& rep) {
        // --- Write .txt file ---
        if (rep.empty()) {
            cerr << "Note: txt is empty.\n";
        }
    
        string txtfilename = remove_extension(filename, ".fa") + ".ue" + to_string(K) + ".txt";
        ofstream txtfile(txtfilename);
        if (!txtfile) {
            cerr << "Error: Could not open file " << txtfilename << " for writing.\n";
            return;
        }
        txtfile << rep;
        txtfile.close();
    
        fs::path txtfilepath = txtfilename;
        auto txtfilesize = fs::file_size(txtfilepath);
        cout << "\nFile created: " << path_to_filename(txtfilename, '/')
             << " (File Size: " << txtfilesize << " B)\n";
    
        // --- Log file ---
        ofstream logfile(logfilename, ios::app);
        if (!logfile) {
            cerr << "Error: Could not open file " << logfilename << " for writing.\n";
            return;
        }
        logfile << "\nFile created: " << path_to_filename(txtfilename, '/')
                << " (File Size: " << txtfilesize << " B)\n";
        logfile.close();
    }
};

class NDBG{
public:
    string filename, logfilename;
    INT K;
    INT option;
    size_t pool_size, reserve_size;
    
    pmr::monotonic_buffer_resource pool;
    using UMAP = pmr::unordered_map<INT, INT>;
    UMAP kmers;
    UMAP heads;

    NDBG(string _filename, INT _K, INT _option)
    : filename(_filename), K(_K), option(_option),
      pool_size(estimate_pool_size_from_file(_filename, _K)),
      reserve_size(pool_size / 32),
      pool(pool_size),
      kmers(&pool), heads(&pool)
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
        get_kmers(kmers, kmerv);
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> get_kmers_time = end_time - start_time;

        VINT match_u, match_v;

        start_time = chrono::high_resolution_clock::now();
        INT M = hopcroft_karp(kmers, kmerv, match_u, match_v);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> hopcroft_karp_time = end_time - start_time;

        VVINT cycles, paths;

        start_time = chrono::high_resolution_clock::now();
        decompose(match_u, match_v, cycles, paths);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> decompose_time = end_time - start_time;

        cout << "\n# of (k-mers, matching, cycles, paths)\n= (" 
            << kmers.size() << ", " <<  M << ", " << cycles.size() << ", " << paths.size() << ")\n";

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
            rep = unsorted(kmers, kmerv, cycles, paths);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else if (option == 2) {
            // sorted
            start_time = chrono::high_resolution_clock::now();
            rep = sorted(kmers, kmerv, cycles, paths);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else if (option == 3) {
            // tree (BP)
            start_time = chrono::high_resolution_clock::now();
            rep = forest(kmers, kmerv, cycles, paths);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else {
            cerr << "Invalid option !\n";
            exit(1);
        }

        // display benchmarks
        print_table_header(logfile);
        log_time("get_kmers", get_kmers_time, logfile);
        log_time("hopcroft_karp", hopcroft_karp_time, logfile);
        log_time("decompose", decompose_time, logfile);
        log_time((option == 0 ? "plain" : option == 1 ? "unsorted" : option == 2 ? "sorted" : "tree (BP)"),
                            align_time, logfile);
        logfile << string(30, '-');
        cout << string(30, '-');

        auto total_time = get_kmers_time.count() + 
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

    void get_kmers(UMAP& kmers, VINT& kmerv) {
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
        cout << "\n";
    }

    INT step(const UMAP& kmers, const VINT& kmerv, INT id, char c, bool is_forward) {
        INT hash = kmerv[id]; 
        if (hash == INF) return INF;

        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return INF;

        if (is_forward) {
            INT mask = (1ULL << (2 * (K - 1))) - 1;
            hash = (hash & mask) << 2;
            if (c == 'C')       hash |= 1;
            else if (c == 'G')  hash |= 2;
            else if (c == 'T')  hash |= 3;
        } else {
            hash >>= 2;
            if (c == 'C')       hash |= (1ULL << (2 * (K - 1)));
            else if (c == 'G')  hash |= (2ULL << (2 * (K - 1)));
            else if (c == 'T')  hash |= (3ULL << (2 * (K - 1)));
        }

        auto it = kmers.find(hash);
        if (it != kmers.end() && it->second != id)
            return it->second;
        
        return INF; // no branch to c
    }

    bool bfs (const UMAP& kmers, const VINT& kmerv, const VINT& match_u, const VINT& match_v, VINT& dist) {
        INT N = (INT)kmers.size();
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
            for (auto c: base) {
                auto v = step(kmers, kmerv, u, c, true);
                if (v == INF) continue;

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

    bool dfs (const UMAP& kmers, const VINT& kmerv, VINT& match_u, VINT& match_v, 
            VINT& dist, INT u) {
        for (auto c: base) {
            auto v = step(kmers, kmerv, u, c, true);
            if (v == INF) continue;

            if (match_v[v] == INF || (dist[match_v[v]] == dist[u] + 1 && dfs(kmers, kmerv, match_u, match_v, dist, match_v[v]))) {
                match_u[u] = v;
                match_v[v] = u;
                return true;
            }
        }
        dist[u] = INF; // not found any augpath
        return false;
    }

    INT hopcroft_karp(const UMAP& kmers, const VINT& kmerv, VINT& match_u, VINT& match_v) {
        INT N = (INT)kmers.size();
        VINT dist(N); // distance of bfs
        INT M = 0;
        match_u.assign(N, INF);
        match_v.assign(N, INF);
        while (bfs(kmers, kmerv, match_u, match_v, dist)) {
            for (INT u = 0; u < N; ++u) {
                if (match_u[u] == INF && dfs(kmers, kmerv, match_u, match_v, dist, u))
                    M++; // found an augpath
            }
            progress(M, N, "Maximum matching");
        }
        finished("Maximum matching");
        cout << "Found maximum matching of size " << M << "\n\n";

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

    void decompose(const VINT& match_u, const VINT& match_v, VVINT& cycles, VVINT& paths) {
        INT N = (INT)kmers.size();
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
                    cerr << "NODE" << node << " IS DUPLICATED\n";
                    exit(1);
                }
                seen[node] = 1;
                ++base_count;
            }
        }
        for (const auto& path: paths) {
            for (const auto& node: path) {
                if (seen[node]) {
                    cerr << "NODE" << node << " IS DUPLICATED\n";
                    exit(1);
                }
                seen[node] = 1;
                ++base_count;
            }
        }
        if (base_count != N) {
            cerr << "WRONG NUMBER OF BASES IN CYCLES AND PATHS\n";
            exit(1);
        }
        cout << "No duplicates found. Total number of bases is correct.\n";
        cout << "Found " << cycles.size() << " cycles and " 
                         << paths.size() << " paths.\n";
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

    REP unsorted(const UMAP& kmers, const VINT& kmerv, const VVINT& cycles, const VVINT& paths) {
        INT N = (INT)kmers.size();
        INT C = (INT)cycles.size();
        INT P = (INT)paths.size(); 
        INT S = N + C + P;

        USET root_paths;
        for (INT i = 0; i < P; ++i){
            bool is_root_path = true;
            INT head = paths[i][0];
            for (auto c: base) {
                if (!is_root_path) continue;
                if (step(kmers, kmerv, head, c, false) != INF)
                    is_root_path = false;
            }
            if (is_root_path) root_paths.insert(i);
            else heads[paths[i][0]] = i;
        }

        INT pos = 0;
        VINT pnt(P, -1);
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                for (const auto& c: base) {
                    auto next = step(kmers, kmerv, node, c, true);
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
                    auto next = step(kmers, kmerv, node, c, true);
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
        for (const auto& pid: root_paths) {
            pnt[pid] = pos;
            auto path = paths[pid];
            for (const auto& node: path) {
                for (const auto& c: base) {
                    auto next = step(kmers, kmerv, node, c, true);
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
        finished("Pointing");
        end:
        if (heads.empty()) 
            cout << "All paths are pointed.\n";
        else cout << heads.size() << " paths are not pointed.\n";

        // generate representation
        string txt;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                auto c = decode_base(kmerv[node] % 4);
                txt += c;
            } txt += "$";
        } for (INT i = 0; i < P; ++i) {
            auto path = paths[i];
            if (root_paths.find(i) != root_paths.end())
                txt += "$" + decode_kmer(kmerv[path[0]], K).substr(0, K - 1);
            for (const auto& node: path) {
                auto c = decode_base(kmerv[node] % 4);
                txt += c;
            } txt += "$";
        } if(!txt.empty()) txt.pop_back();

        return {txt, pnt};
    }

    VINT find_new_cycle(const UMAP& kmers, const VINT& kmerv, const VVINT& paths, Vint& in_pord,  Vint& is_leaf, Vint& visited, INT start) {
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
                is_leaf[top.current] = 1; // memo as a leaf
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
                INT next_node = step(kmers, kmerv, cur_node, c, true);
                if (next_node == INF) continue;
                auto it = heads.find(next_node);
                if (it == heads.end()) continue;
                INT next_pid = it->second;
                if (in_pord[next_pid]) continue;
                if (next_pid == start) return cycle;
                // skip visited path
                if (visited[next_pid]) continue;
                if (is_leaf[next_pid]) continue; // skip leaves

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

    REP sorted(const UMAP& kmers, const VINT& kmerv, VVINT& cycles, VVINT& paths) {
        INT N = (INT)kmers.size();
        INT C = (INT)cycles.size();
        INT P = (INT)paths.size();
        INT S = N + C + P, ncycs = 0;

        USET roots;
        for (INT i = 0; i < P; ++i) {
            bool is_root_path = true;
            INT head = paths[i][0];
            for (auto c: base) {
                if (!is_root_path) continue;
                if (step(kmers, kmerv, head, c, false) != INF)
                    is_root_path = false;
            }            
            if (is_root_path) roots.insert(i);
            else heads[paths[i][0]] = i;
        }
        cout << "\n# of (non-root paths, root paths) = (" << heads.size() << ", "
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
                        auto next = step(kmers, kmerv, node, c, true);
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
                    auto next = step(kmers, kmerv, node, c, true);
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
                        auto next = step(kmers, kmerv, node, c, true);
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

        Vint is_leaf(P, 0);

        // find a new cycle and bfs from it
        while (INT res = count(in_ord.begin(), in_ord.end(), 0)) {
            cout << "\n\n#residue paths = " << res << "\n"
                 << "Searching for a new cycle...";
            for (INT i = 0; i < P; ++i) {
                if (in_ord[i]) continue;
                if (is_leaf[i]) continue; // skip leaves
                Vint vis(P, 0);
                VINT pids = find_new_cycle(kmers, kmerv, paths, in_ord, is_leaf, vis, i);
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
                                auto next = step(kmers, kmerv, node, c, true);
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
                            auto next = step(kmers, kmerv, node, c, true);
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
                                auto next = step(kmers, kmerv, node, c, true);
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
        if (res == 0) cout << "All paths are pointed.\n\n";
        else {
            cout << res << " paths are not pointed.\n";
            for (INT i = 0; i < P; ++i) if (!in_ord[i])
                cout << "path" << i << " not added\n";
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

    string write_subtree(INT root_pid,
                                const VVINT& paths,
                                const UMAP& kmers,
                                const VINT& kmerv,
                                const UMAP& heads,
                                Vint& added) {
        struct Frame {
            INT pid;
            size_t node_idx    = 0;
            size_t base_idx    = 0;
            bool   printed    = false;
        };

        std::ostringstream oss;
        std::stack<Frame> st;

        added[root_pid] = 1;
        st.push({root_pid, 0, 0, false});

        while (!st.empty()) {
            Frame &fr = st.top();
            const auto &path = paths[fr.pid];

            if (fr.node_idx >= path.size()) {
                st.pop();
                if (!st.empty()) oss << ')';
                continue;
            }

            INT node = path[fr.node_idx];
            if (!fr.printed) {
                oss << decode_base(kmerv[node] % 4);
                fr.printed = true;
                fr.base_idx = 0;
            }

            bool pushed = false;
            for (; fr.base_idx < 4; ++fr.base_idx) {
                char c = base[fr.base_idx];
                INT nxt = step(kmers, kmerv, node, c, true);
                if (nxt == INF) continue;
                auto it = heads.find(nxt);
                if (it == heads.end()) continue;
                INT child_pid = it->second;
                if (added[child_pid]) continue;
                added[child_pid] = 1;
                ++fr.base_idx;
                oss << '(';
                st.push({child_pid, 0, 0, false});
                pushed = true;
                break;
            }

            if (!pushed) {
                fr.node_idx++;
                fr.printed = false;
            }
        }
        return oss.str();
    }

    REP forest(const UMAP& kmers, const VINT& kmerv, VVINT& cycles, VVINT& paths) {
        string txt;
        INT C = (INT)cycles.size();
        INT P = (INT)paths.size();
        INT S = C + P, n = 0, ncycs = 0;

        USET roots;        
        for (INT i = 0; i < P; ++i) {
            bool is_root_path = true;
            INT head = paths[i][0];
            for (auto c: base) {
                if (!is_root_path) continue;
                if (step(kmers, kmerv, head, c, false) != INF)
                    is_root_path = false;
            }
            if (is_root_path) roots.insert(i);
            else heads[paths[i][0]] = i;
        }
        cout << "\n(#non-root paths, #root paths) = (" << heads.size() << ", "
             << roots.size() << ")\n";

        Vint added(P, 0); // 1 iff path already embedded

        for (const auto& pid: roots) {
            txt += decode_kmer(kmerv[paths[pid][0]], K).substr(0, K - 1)
                   + write_subtree(pid, paths, kmers, kmerv, heads, added)
                   + ",";
            ++n; progress(n, S, "Finding pseudoforest");
        }

        for (const auto& cycle: cycles) {
            txt += "*";
            for (const auto& node: cycle) {
                txt += decode_base(kmerv[node] % 4);
                for (const auto& c: base) {
                    auto next = step(kmers, kmerv, node, c, true);
                    if (next == INF) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    auto npid = it->second;
                    if (added[npid]) continue;
                    txt += "(" + write_subtree(npid, paths, kmers, kmerv, heads, added) + ")";
                }
            }
            txt += ","; ++n; progress(n, S, "Finding pseudoforest");
        }
        cout << "\n";

        Vint is_leaf(P, 0);

        while (INT res = count(added.begin(), added.end(), 0)) {
            cout << "\r#residue paths = " << res << ", searching for a new cycle..." << flush;
            for (INT i = 0; i < P; ++i) {
                if (added[i]) continue;
                if (is_leaf[i]) continue;
                Vint vis(P, 0);
                VINT pids = find_new_cycle(kmers, kmerv, paths, added, is_leaf, vis, i);
                if (!pids.empty()) {
                    cout << "\r#residue paths = " << count(added.begin(), added.end(), 0) << ", found a cycle of " << pids.size() << " paths          " << flush;
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
                                auto next = step(kmers, kmerv, node, c, true);
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
                            auto next = step(kmers, kmerv, node, c, true);
                            if (next == INF) continue;
                            auto it = heads.find(next);
                            if (it == heads.end()) continue;
                            auto npid = it->second;
                            if (added[npid]) continue;
                            txt += "(" + write_subtree(npid, paths, kmers, kmerv, heads, added) + ")";
                        }
                    }
                    txt += ","; ++n; progress(n, S, "Finding pseudoforest");
                    break;
                }
            }
            if ((INT)count(added.begin(), added.end(), 0) == res) {
                cerr << "ERROR: no new cycle found in remaining paths\n";
                exit(1);
            }
        }
        finished("Finding pseudoforest");

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

void print_usage(const string& cmdname) {
    cerr << "Usage: " << cmdname << " [0(node-centric) / 1(edge-centric)] [****.fa] [k] ([0 / 1 / 2 / 3])\n"
        << "Output options(node-centric) are...\n"
        << "0: explicit representation without pointers (SPSS)\n"
        << "1: representation with unsorted pointers\n"
        << "2: representation with sorted pointers\n"
        << "3: tree representation (BP)\n"
        << "fna or fasta files are also ok as an input\n"
        << "if no extension, \".fa\" is automatically added\n";
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
        string rep = edbg.process();
        edbg.write(rep);
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