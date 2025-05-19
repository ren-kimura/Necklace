#include <sdsl/bit_vectors.hpp>
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

using namespace sdsl;
using namespace std;
namespace fs = std::filesystem;
namespace chrono = std::chrono;

using VSTR = vector<string>;
using V8 = vector<uint8_t>;
using VT = vector<size_t>;
using VVT = vector<VT>;
using REP = pair<string, VT>;
using USET = unordered_set<size_t>;

vector<char> base = {'A', 'C', 'G', 'T'};
const size_t INF = SIZE_MAX;

size_t get_peak_memory_kb() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // in kilobytes
}

void progress(size_t now, size_t total, const string& task) {
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

size_t encode_kmer(const string& kmer) {
    size_t val = 0;
    for (char c : kmer) {
        if (c == 'A')      val = (val << 2) | 0; // A -> 00
        else if (c == 'C') val = (val << 2) | 1; // C -> 01
        else if (c == 'G') val = (val << 2) | 2; // G -> 10
        else if (c == 'T') val = (val << 2) | 3; // T -> 11
        else return INF;
    }
    return val;
}

string decode_kmer(size_t hash, size_t K) {
    string kmer;
    for (size_t i = 0; i < K; ++i) {
        int base = hash & 3;  // take rightmost 2 bits
        if (base == 0)        kmer = 'A' + kmer;  // 00 -> A
        else if (base == 1)   kmer = 'C' + kmer;  // 01 -> C
        else if (base == 2)   kmer = 'G' + kmer;  // 10 -> G
        else if (base == 3)   kmer = 'T' + kmer;  // 11 -> T
        hash >>= 2;
    }
    return kmer;
}

char decode_base(size_t val) {
    return "ACGT"[val];
}

size_t rolling_valid_kmer_start(const string& read, size_t start, size_t K) {
    size_t len = read.size();
    size_t count = 0;
    for (size_t j = start; j < len; ++j) {
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
    size_t K;
    bit_vector B0;
    bit_vector B1;
    rank_support_v<> b1rank1; // hash -> id
    select_support_mcl<> b1select1; // id -> hash
    size_t N; // total # of kmers
    unordered_map<size_t, int> outcnt;
    unordered_map<size_t, VT> dummy_map;

    struct pair_hash {
        size_t operator()(pair<size_t, size_t> const &p) const {
            return hash<size_t>()(p.first) ^ (hash<size_t>()(p.second) << 1);
        }
    };
    pmr::unordered_multiset<pair<size_t, size_t>, pair_hash> used_dummy_set;

    EDBG(string _filename, size_t _K)
        : filename(_filename), K(_K)
    {
        logfilename = remove_extension(filename, ".fa") + ".ue" + to_string(K) + ".bv.log";
        size_t U = size_t(1) << (2 * K);
        B0 = bit_vector(U, 0);
        B1 = bit_vector(U >> 2, 0);
    };

    inline size_t i2h(size_t id) const {
        return b1select1(id + 1);
    }

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
        util::init_support(b1rank1, &B1);
        util::init_support(b1select1, &B1);
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> get_kmers_time = end_time - start_time;
        cout << "Collected " << N << " unique k-mers\n\n";

        start_time = chrono::high_resolution_clock::now();
        balance_graph();
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> balance_graph_time = end_time - start_time;

        start_time = chrono::high_resolution_clock::now();
        VT tour = find_eulerian_cycle();
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> find_eulerian_cycle_time = end_time - start_time;

        start_time = chrono::high_resolution_clock::now();
        string txt = spell_out(tour);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> spell_out_time = end_time - start_time;

        print_table_header(logfile);
        log_time("get_kmers", get_kmers_time, logfile);
        log_time("balance_graph", balance_graph_time, logfile);
        log_time("find_eulerian_cycle", find_eulerian_cycle_time, logfile);
        log_time("spell_out", spell_out_time, logfile);
        logfile << string(30, '-');
        cout << string(30, '-');
        auto total_time = get_kmers_time.count()
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

        size_t tlen = 0;
        to_uppercase(reads);
        for (string& read: reads)
            tlen += read.size();
        
        size_t footing = 0;
        size_t mask = (1ULL << (2 * (K - 1))) - 1; // 2K-bits mask

        for (size_t i = 0; i < reads.size(); ++i) {
            const string& read = reads[i];
            size_t len = read.size();
            if (len < K) continue;
            size_t j = rolling_valid_kmer_start(read, 0, K);
            if (j == INF) continue; // if not any valid k-mer, go to next read
            size_t hash = encode_kmer(read.substr(j, K));
            if (B0[hash] == 0) {
                B0[hash] = 1; ++N; // the k-mer exists
                B1[hash >> 2] = 1; // pre k-1-mer
                B1[hash & mask] = 1; // suf k-1-mer
                ++outcnt[hash >> 2];
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
                    B1[hash >> 2] = 1; // pre k-1-mer
                }
                if (B0[hash] == 0) {
                    B0[hash] = 1; ++N; // the k-mer exists
                    B1[hash & mask] = 1; // suf k-1-mer
                    ++outcnt[hash >> 2];
                }
                progress(footing + j, tlen, "Detecting k-mers");
            }
            footing += len;
        }
        finished("Detecting k-mers");
    }

    size_t step(size_t hash, char c, bool is_forward) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return INF;
        if (is_forward) {
            hash = hash << 2;
            if (c == 'C')       hash |= 1;
            else if (c == 'G')  hash |= 2;
            else if (c == 'T')  hash |= 3;
            if (B0[hash]) {
                size_t mask = (1ULL << (2 * (K - 1))) - 1;
                return hash & mask;
            }
        } else {
            if (c == 'C')       hash |= (1ULL << (2 * (K - 1)));
            else if (c == 'G')  hash |= (2ULL << (2 * (K - 1)));
            else if (c == 'T')  hash |= (3ULL << (2 * (K - 1)));
            if (B0[hash]) return hash >> 2;
        }
        return INF;
    }

    void balance_graph() {
        VT sources, sinks;
        size_t M = b1rank1(B1.size());

        for (size_t i = 0; i < M; ++i) {
            auto hash = i2h(i);
            int odeg = 0, ideg = 0;
            for (const auto c: base) {
                if (step(hash, c, true) != INF) ++odeg;
                if (step(hash, c, false) != INF) ++ideg;
            }
            if (odeg > ideg) {
                int diff = odeg - ideg;
                for (int j = 0; j < diff; ++j) sources.emplace_back(hash);
            } else if (ideg > odeg) {
                int diff = ideg - odeg;
                for (int j = 0; j < diff; ++j) sinks.emplace_back(hash);
            }
            progress(i, M, "Counting imbalances");
        }
        finished("Counting imbalances"); cout << "\n";

        if (sources.size() != sinks.size()) {
            cerr << "ERROR: mismatch between sources and sinks!\n";
            exit(1);
        }
        auto I = sources.size();
        for (size_t i = 0; i < I; ++i) {
            size_t u = sinks[i];
            size_t v = sources[i];
            ++outcnt[u];
            dummy_map[u].push_back(v);
            progress(i, I, "Adding dummy edges");
        }
        finished("Adding dummy edges");
        cout << "Added " << sources.size() << " dummy edges to balance the graph\n\n";
    }

    VT find_eulerian_cycle() {
        size_t m = 0, M = 0; // M: total moves
        for (auto& [node, cnt]: outcnt) M += cnt;
        auto dmap = dummy_map;

        size_t start = INF;
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
        VT tour, curpath;
        curpath.push_back(start);
        
        while(!curpath.empty()) {
            size_t u = curpath.back();
            if (outcnt[u] > 0) {
                bool moved = false;
                auto it_map = dmap.find(u);
                if (it_map != dmap.end() && !it_map->second.empty()) {
                    size_t v = it_map->second.back();
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
                for (const auto c: base) {
                    auto v = step(u, c, true);
                    if (v == INF) continue;
                    curpath.push_back(v);
                    --outcnt[u]; // used one edge
                    auto mhash = c == 'A' ? u << 2 :
                                 c == 'C' ? u << 2 | 1 :
                                 c == 'G' ? u << 2 | 2 :
                                            u << 2 | 3 ;
                    B0[mhash] = 0; // erase used edge
                    moved = true;
                    progress(++m, M, "Finding Eulerian tour");
                    break;
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

    string spell_out(const VT& tour) {
        string txt;
        size_t n = tour.size();
        if (n < 2) { cerr << "Tour too short\n"; exit(1); }

        for (size_t i = 0; i < n - 1; ++i) {
            pair<size_t, size_t> edge{tour[i], tour[i + 1]};
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
    
        string txtfilename = remove_extension(filename, ".fa") + ".ue" + to_string(K) + ".bv.txt";
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
    size_t K;
    size_t option;
    bit_vector B;
    rank_support_v<> rank1; // hash -> id
    select_support_mcl<> select1; // id -> hash 
    size_t N; // total # of kmers
    unordered_map<size_t, size_t> heads;

    NDBG(string _filename, size_t _K, size_t _option)
    : filename(_filename), K(_K), option(_option)
    {
        if (option != 0 && option != 1 && option != 2 && option != 3) {
            cerr << "\nInvalid option value\n";
            exit(1);
        }
        logfilename = remove_extension(filename, ".fa") + ".o" + to_string(K) + "." + to_string(option) + ".bv.log";
        size_t U = size_t(1) << (2 * K);
        B = bit_vector(U, 0);
    };

    inline size_t h2i(size_t hash) const {
        return rank1(hash); // rank1(i)=\sum_{k=0}^{k<i}B[i] (NOT INCLUDING B[i])
        // see https://simongog.github.io/assets/data/sdsl-cheatsheet.pdf
    }

    inline size_t i2h(size_t id) const {
        return select1(id + 1); // select1(i) = min{j|rank1(j+1)=i}
    }

    REP process() {
        ofstream logfile(logfilename, ios::trunc);
        if (!logfile) {
            cerr << "Error: Could not open " << logfilename << " for writing.\n";
            exit(1);
        }
        string timestamp = get_current_timestamp();
        logfile << "Timestamp: " << timestamp << "\n";
        cout << "Timestamp: " << timestamp << "\n";
        
        auto start_time = chrono::high_resolution_clock::now();        
        get_kmers();
        util::init_support(rank1, &B);
        util::init_support(select1, &B);
        N = rank1(B.size());
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> get_kmers_time = end_time - start_time;

        // // debug
        // cout << "-----------------\n";
        // for (size_t i = 0; i < N; ++i)
        //     cout << i << " : " << decode_kmer(i2h(i), K) << "\n";
        // cout << "-----------------\n";
        // for (size_t i = 0; i < N; ++i) {
        //     cout << i << " : ";
        //     for (const auto& c: base) {
        //         auto j = step(i, c, true);
        //         if (j != INF) cout << j << " ";
        //     }
        //     cout << "\n";
        // }
        // cout << "-----------------\n";

        VT match_u, match_v;

        start_time = chrono::high_resolution_clock::now();
        size_t M = hopcroft_karp(match_u, match_v);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> hopcroft_karp_time = end_time - start_time;

        VVT cycles, paths;

        start_time = chrono::high_resolution_clock::now();
        decompose(match_u, match_v, cycles, paths);
        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> decompose_time = end_time - start_time;

        cout << "\n# of (k-mers, matching, cycles, paths)\n= (" 
            << N << ", " <<  M << ", " << cycles.size() << ", " << paths.size() << ")\n";

        REP rep;
        chrono::duration<double> align_time;
        if (option == 0) {
            // plain
            start_time = chrono::high_resolution_clock::now();
            rep = plain(cycles, paths);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else if (option == 1) {
            // unsorted
            start_time = chrono::high_resolution_clock::now();
            rep = unsorted(cycles, paths);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else if (option == 2) {
            // sorted
            start_time = chrono::high_resolution_clock::now();
            rep = sorted(cycles, paths);
            end_time = chrono::high_resolution_clock::now();
            align_time = end_time - start_time;
        } else if (option == 3) {
            // tree (BP)
            start_time = chrono::high_resolution_clock::now();
            rep = forest(cycles, paths);
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

        size_t tlen = 0;
        to_uppercase(reads);
        for (string& read: reads)
            tlen += read.size();
        
        size_t footing = 0;
        size_t mask = (1ULL << (2 * (K - 1))) - 1; // 2K-bits mask

        for (size_t i = 0; i < reads.size(); ++i) {
            const string& read = reads[i];
            size_t len = read.size();
            if (len < K) continue;
            size_t j = rolling_valid_kmer_start(read, 0, K);
            if (j == INF) continue; // if not any valid k-mer, go to next read
            size_t hash = encode_kmer(read.substr(j, K));
            B[hash] = 1; // the k-mer exists

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
                B[hash] = 1; // the k-mer exists
                progress(footing + j, tlen, "Detecting k-mers");
            }
            footing += len;
        }
        finished("Detecting k-mers");
        cout << "\n";
    }

    size_t step(size_t id, char c, bool is_forward) {
        size_t hash = i2h(id); 
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return INF;
        if (is_forward) {
            size_t mask = (1ULL << (2 * (K - 1))) - 1;
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
        if (hash == i2h(id)) return INF; // rule out self loops
        auto b = B[hash];
        if (b) return h2i(hash);
        else return INF;
    }

    bool bfs (const VT& match_u, const VT& match_v, VT& dist) {
        queue<size_t> q;
        bool aug = false;
        for (size_t u = 0; u < N; ++u) {
            if (match_u[u] == INF) {
                dist[u] = 0;
                q.push(u);
            } else dist[u] = INF;
        }
        while (!q.empty()) {
            size_t u = q.front();
            q.pop();
            for (auto c: base) {
                auto v = step(u, c, true);
                if (v == INF) continue;
                if (match_v[v] == INF) {
                    aug = true;
                } else if (dist[match_v[v]] == INF) {
                    dist[match_v[v]] = dist[u] + 1;
                    q.push(match_v[v]);
                }
            }
        }
        return aug;
    }

    bool dfs (VT& match_u, VT& match_v, VT& dist, size_t u) {
        for (auto c: base) {
            auto v = step(u, c, true);
            if (v == INF) continue;
            if (match_v[v] == INF || (dist[match_v[v]] == dist[u] + 1 && dfs(match_u, match_v, dist, match_v[v]))) {
                match_u[u] = v;
                match_v[v] = u;
                return true;
            }
        }
        dist[u] = INF; // augpath not found
        return false;
    }

    size_t hopcroft_karp(VT& match_u, VT& match_v) {
        VT dist(N); // distance of bfs
        size_t M = 0;
        match_u.assign(N, INF);
        match_v.assign(N, INF);
        while (bfs(match_u, match_v, dist)) {
            for (size_t u = 0; u < N; ++u) {
                if (match_u[u] == INF && dfs(match_u, match_v, dist, u))
                    M++; // found an augpath
            }
            progress(M, N, "Maximum matching");
        }
        finished("Maximum matching");
        cout << "Found maximum matching of size " << M << "\n\n";
        return M;
    }

    VT dfs_forward(const VT& match_u, const size_t& start, V8& seen_u, V8& seen_v) {
        VT trl;
        size_t crnt = start;
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

    VT dfs_backward(const VT& match_v, const size_t& start, V8& seen_u, V8& seen_v) {
        VT lrt;
        size_t crnt = start;
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

    void decompose(const VT& match_u, const VT& match_v, VVT& cycles, VVT& paths) {
        V8 seen_u(N, 0);
        V8 seen_v(N, 0);

        for (size_t u = 0; u < N; ++u) {
            if (seen_u[u]) continue;
            VT trl = dfs_forward(match_u, u, seen_u, seen_v);
            if (trl.size() > 1 && match_u[trl.back()] == trl.front())
                cycles.emplace_back(trl);
            else {
                VT lrt = dfs_backward(match_v, trl.front(), seen_u, seen_v);
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
        V8 seen(N, 0); size_t base_count = 0;
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

        // // debug
        // cout << "-------------------------------\n";
        // for (size_t i = 0; i < cycles.size(); ++i) {
        //     cout << "C" << i << ": ";
        //     auto c = cycles[i];
        //     for (const auto& node: c) cout << decode_kmer(i2h(node), K) << "(" << node << ") ";
        //     cout << "\n";
        // }
        // cout << "-------------------------------\n";
        // for (size_t i = 0; i < paths.size(); ++i) {
        //     cout << "P" << i << ": ";
        //     auto p = paths[i];
        //     for (const auto& node: p) cout << decode_kmer(i2h(node), K) << "(" << node << ") ";
        //     cout << "\n";
        // }
        // cout << "-------------------------------\n";
    }

    REP plain(const VVT& cycles, const VVT& paths) {
        string txt;
        VT pnt;
        cout << "Generating SPSS...";
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                auto c = decode_base(i2h(node) % 4);
                txt += c;
            } txt += "$";
        } txt += "$";
        for (const auto& path: paths) {
            for (const auto& node: path) {
                auto hash = i2h(node);
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

    REP unsorted(const VVT& cycles, const VVT& paths) {
        size_t C = cycles.size();
        size_t P = paths.size(); 
        size_t S = N + C + P;

        USET roots;
        for (size_t i = 0; i < P; ++i){
            bool is_root = true;
            size_t head = paths[i][0];
            for (auto c: base) {
                if (!is_root) continue;
                if (step(head, c, false) != INF)
                    is_root = false;
            }
            if (is_root) roots.insert(i);
            else heads[paths[i][0]] = i;
        }

        size_t pos = 0;
        VT pnt(P, -1);
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                for (const auto& c: base) {
                    auto next = step(node, c, true);
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
                    auto next = step(node, c, true);
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
        for (const auto& pid: roots) {
            pnt[pid] = pos;
            auto path = paths[pid];
            for (const auto& node: path) {
                for (const auto& c: base) {
                    auto next = step(node, c, true);
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
                auto c = decode_base(i2h(node) % 4);
                txt += c;
            } txt += "$";
        } for (size_t i = 0; i < P; ++i) {
            auto path = paths[i];
            if (roots.find(i) != roots.end())
                txt += "$" + decode_kmer(i2h(path[0]), K).substr(0, K - 1);
            for (const auto& node: path) {
                auto c = decode_base(i2h(node) % 4);
                txt += c;
            } txt += "$";
        } if(!txt.empty()) txt.pop_back();

        return {txt, pnt};
    }

    VT find_new_cycle(const VVT& paths, V8& in_pord,  V8& is_leaf, V8& visited, size_t start) {
        struct Frame {
            size_t current;
            size_t idx;
            size_t b_idx;
        };
        vector<Frame> stack;
        VT cycle;

        visited[start] = 1;
        stack.push_back({start, 0, 0});

        while (!stack.empty()) {
            Frame &top = stack.back();
            const VT &cur_path = paths[top.current];

            // rm from stack if explored all the nodes in cur_path
            if (top.idx >= cur_path.size()) {
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
                
                size_t cur_node = cur_path[top.idx];
                size_t next_node = step(cur_node, c, true);
                if (next_node == INF) continue;
                auto it = heads.find(next_node);
                if (it == heads.end()) continue;
                size_t next_pid = it->second;
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

    REP sorted(VVT& cycles, VVT& paths) {
        size_t C = cycles.size();
        size_t P = paths.size();
        size_t S = N + C + P, ncycs = 0;

        USET roots;
        for (size_t i = 0; i < P; ++i) {
            bool is_root = true;
            size_t head = paths[i][0];
            for (auto c: base) {
                if (!is_root) continue;
                if (step(head, c, false) != INF)
                    is_root = false;
            }            
            if (is_root) roots.insert(i);
            else heads[paths[i][0]] = i;
        }
        cout << "\n# of (non-root paths, root paths) = (" << heads.size() << ", "
             << roots.size() << ")\n";

        queue<size_t> q;
        vector<pair<uint8_t, size_t>> ord; // sorted pairs (0 or 1, id), where 0 is cycle, 1 is path
        V8 in_ord(P, 0); // 0: path is in ord, 1: not in ord yet
        VT pnt;

        size_t pos = 0;
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
                        auto next = step(node, c, true);
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
        for (size_t i = 0; i < C; ++i) {
            auto cycle = cycles[i];
            ord.emplace_back(0, i);
            ++pos; // need flag "*" for a cycle
            for (const auto& node: cycle) {
                for (const auto& c: base) {
                    auto next = step(node, c, true);
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
                        auto next = step(node, c, true);
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

        V8 is_leaf(P, 0);

        // find a new cycle and bfs from it
        while (size_t res = count(in_ord.begin(), in_ord.end(), 0)) {
            cout << "\n\n#residue paths = " << res << "\n"
                 << "Searching for a new cycle...";
            for (size_t i = 0; i < P; ++i) {
                if (in_ord[i]) continue;
                if (is_leaf[i]) continue; // skip leaves
                V8 vis(P, 0);
                VT pids = find_new_cycle(paths, in_ord, is_leaf, vis, i);
                if (!pids.empty()) {
                    cout << "\rFound a cycle of " << pids.size() << " paths          \n";
                    ++ncycs;
                    VT cycle;
                    size_t p = pids.size();
                    for (size_t j = 0; j < p; ++j) {
                        auto& path = paths[pids[j]];
                        auto h = path.size();
                        for (size_t k = 0; k < h; ++k) {
                            auto node = path[k];
                            cycle.emplace_back(node);
                            for (auto& c: base) {
                                auto next = step(node, c, true);
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
                    ord.emplace_back(0, cycles.size() - 1);
                    ++pos; // need flag "*" for a cycle
                    for (const auto& node: cycle) {
                        for (const auto& c: base) {
                            auto next = step(node, c, true);
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
                                auto next = step(node, c, true);
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
        size_t res = count(in_ord.begin(), in_ord.end(), 0);
        if (res == 0) cout << "All paths are pointed.\n\n";
        else {
            cout << res << " paths are not pointed.\n";
            for (size_t i = 0; i < P; ++i) if (!in_ord[i])
                cout << "path" << i << " not added\n";
        }
        cout << "#new cycles = " << ncycs << "\n";

        // take difference of pnt
        for (size_t i = pnt.size() - 1; i > 0; --i) {
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
                    txt += "$" + decode_kmer(i2h(path[0]), K).substr(0, K - 1);
                for (const auto& node: path)
                    txt += decode_base(i2h(node) % 4);
            } else {
                auto cycle = cycles[id];
                txt += "*";
                for (const auto& node: cycle)
                    txt += decode_base(i2h(node) % 4);
            }
            txt += "$";
        }
        if (!txt.empty()) txt.pop_back();
        return {txt, pnt};
    }

    string write_subtree(size_t root_pid,
                                const VVT& paths,
                                V8& added) {
        struct Frame {
            size_t pid;
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

            size_t node = path[fr.node_idx];
            if (!fr.printed) {
                oss << decode_base(i2h(node) % 4);
                fr.printed = true;
                fr.base_idx = 0;
            }

            bool pushed = false;
            for (; fr.base_idx < 4; ++fr.base_idx) {
                char c = base[fr.base_idx];
                size_t nxt = step(node, c, true);
                if (nxt == INF) continue;
                auto it = heads.find(nxt);
                if (it == heads.end()) continue;
                size_t child_pid = it->second;
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

    REP forest(VVT& cycles, VVT& paths) {
        string txt;
        size_t C = cycles.size();
        size_t P = paths.size();
        size_t S = C + P, n = 0, ncycs = 0;

        USET roots;        
        for (size_t i = 0; i < P; ++i) {
            bool is_root = true;
            size_t head = paths[i][0];
            for (auto c: base) {
                if (!is_root) continue;
                if (step(head, c, false) != INF)
                    is_root = false;
            }
            if (is_root) roots.insert(i);
            else heads[paths[i][0]] = i;
        }
        cout << "\n(#non-root paths, #root paths) = (" << heads.size() << ", "
             << roots.size() << ")\n";

        V8 added(P, 0); // 1 iff path already embedded

        for (const auto& pid: roots) {
            txt += decode_kmer(i2h(paths[pid][0]), K).substr(0, K - 1)
                   + write_subtree(pid, paths, added)
                   + ",";
            ++n; progress(n, S, "Finding pseudoforest");
        }

        for (const auto& cycle: cycles) {
            txt += "*";
            for (const auto& node: cycle) {
                txt += decode_base(i2h(node) % 4);
                for (const auto& c: base) {
                    auto next = step(node, c, true);
                    if (next == INF) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    auto npid = it->second;
                    if (added[npid]) continue;
                    txt += "(" + write_subtree(npid, paths, added) + ")";
                }
            }
            txt += ","; ++n; progress(n, S, "Finding pseudoforest");
        }
        cout << "\n";

        V8 is_leaf(P, 0);

        while (size_t res = count(added.begin(), added.end(), 0)) {
            cout << "\r#residue paths = " << res << ", searching for a new cycle..." << flush;
            for (size_t i = 0; i < P; ++i) {
                if (added[i]) continue;
                if (is_leaf[i]) continue;
                V8 vis(P, 0);
                VT pids = find_new_cycle(paths, added, is_leaf, vis, i);
                if (!pids.empty()) {
                    cout << "\r#residue paths = " << count(added.begin(), added.end(), 0) << ", found a cycle of " << pids.size() << " paths          " << flush;
                    ++ncycs;
                    VT cycle;
                    size_t p = pids.size();
                    for (size_t j = 0; j < p; ++j) {
                        auto& path = paths[pids[j]];
                        auto h = path.size();
                        for (size_t k = 0; k < h; ++k) {
                            auto node = path[k];
                            cycle.emplace_back(node);
                            for (auto& c: base) {
                                auto next = step(node, c, true);
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
                            auto next = step(node, c, true);
                            if (next == INF) continue;
                            auto it = heads.find(next);
                            if (it == heads.end()) continue;
                            auto npid = it->second;
                            if (added[npid]) continue;
                            txt += "(" + write_subtree(npid, paths, added) + ")";
                        }
                    }
                    txt += ","; ++n; progress(n, S, "Finding pseudoforest");
                    break;
                }
            }
            if ((size_t)count(added.begin(), added.end(), 0) == res) {
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
    
        string txtfilename = remove_extension(filename, ".fa") + ".o" + to_string(K) + "." + to_string(option) + ".bv.txt";
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
            pntfilename = remove_extension(filename, ".fa") + ".o" + to_string(K) + ".bv." + to_string(option) + "p.txt";
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

    size_t _K = atoi(argv[3]);

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
        size_t _option = atoi(argv[4]);

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