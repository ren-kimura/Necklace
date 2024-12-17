#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <tuple> // sada
#include <array> // sada
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <cstring>
#include <chrono>
#include <cmath>

using namespace std;
using INT = int64_t;
using Kmers = unordered_map<string,INT>;
using VINT = vector<INT>;
using VVINT = vector<VINT>;
using DINT = deque<INT>;

struct pair_hash {
    template <class T1, class T2>
    INT operator ( )(const pair<T1, T2> &p ) const {
        auto h1 = hash<T1>{}(p.first);
        auto h2 = hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

void print_progress_bar(INT current_kmer, INT total_kmers, const string& task) {
    double progress = static_cast<double>(current_kmer) / static_cast<double>(total_kmers);
    int bar_width = 50;
    int pos = static_cast<int>(bar_width * progress);
    cout << "\r" << task << " [";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos)        cout << "=";
        else if (i == pos)  cout << ">";
        else                cout << " ";
    }
    cout << "] " << static_cast<int>(progress * 100.0) << "%";
    cout.flush();
}

// get cumulative length up to each head of cycle in rep considering delimiters
VINT get_chead(const VVINT& cycles) {
    VINT cumls = {0};
    INT cuml = 0;
    for (const auto& cycle : cycles) {
        cuml += cycle.size() + 1;
        cumls.push_back(cuml);
    }
    return cumls;
}

// get cumulative length of current sorted arr2
INT get_phead(INT& path_id, VINT& sorted_path, const VVINT& paths) {
    INT cuml = 0;
    for (auto it = sorted_path.begin(); *it != path_id; ++it) {
        cuml += paths[*it].size() + 1;
    }
    return cuml;
}

class DeBruijnGraph{
public:
    string in_filename;
    INT K;
    bool isNodeCentric; // node-centric or edge-centric dBG?
    array<uint8_t, 4> Alphabet = {'A', 'C', 'G', 'T'};
    const INT INF = INT64_MAX;

    string sq; // stores the concatenated sequence
    Kmers kmers; // pair <kmer string of length k, its unique ID in 0..N-1>, where N = #distinct kmers in sequence
    VINT idPos; // for each ID a position i such that ID corresponds to kmer sequence.substr(i, K)
    VVINT adj, inv_adj;

    VVINT cycles, paths;
    string rep; // representation string for a cover
    VINT pointers; // pointer from and to a certain pos in rep

    DeBruijnGraph(string _in_filename, INT _K, bool _isNodeCentric): in_filename(_in_filename), K(_K), isNodeCentric(_isNodeCentric)
    {};

    void to_uppercase(string& str) {
        transform(str.begin(), str.end(), str.begin(),
                    [](unsigned char c) {return toupper(c);});
    }

    void getKmers() {
        // read input DNA sequence in FASTA format .fa
        ifstream inputFile(in_filename);
        if (!inputFile) {
            cerr << "Error opening input file." << endl;
            exit(0);
        }
        string line;
        while (getline(inputFile, line)) {
            if (line[0] != '>') { // skip header line
                sq += line;
            }
        }
        inputFile.close();

        // lower to upper
        to_uppercase(sq);

        // erase letters other than 'A, C, G, T'
        sq.erase(remove_if(sq.begin(), sq.end(), [](char x)
            {return (x != 'A' && x != 'C' && x != 'G' && x != 'T');}
            ), sq.end());

        // check every k-mer in sq
        INT total_kmers = sq.length() - K + 1;
        for (INT i = 0; i <= total_kmers - 1; ++i) {
            auto kmer = sq.substr(i, K);
            if (kmers.find(kmer) == kmers.end()) { // newly found kmer
                idPos.push_back(i);
                kmers[kmer] = idPos.size() - 1;
            }
            if (total_kmers <= 100 || i % (total_kmers / 100) == 0)
                print_progress_bar(i, total_kmers, "Detecting k-mers");
        }
        cout << endl;

        adj.resize(kmers.size()); // resize adjacency list
        inv_adj.resize(kmers.size());

    }

    INT forward(INT id, uint8_t c) {
        if (id >= static_cast<INT>(idPos.size())) return -1;
        auto pos = idPos[id];

        if (pos + 1 >= static_cast<INT>(sq.length())) return -1;
 
        string suffix = sq.substr(pos + 1, K - 1);
        auto next = suffix + static_cast<char>(c);
        // next exists AND not current
        if (kmers.find(next) != kmers.end() && kmers[next] != id)
            return kmers[next];
        else
            return -1; // no branch to c
    }

    void addEdges() {
        INT total_kmers = kmers.size();
        INT current_kmer = 0;

        // add edges in adjacency list format
        for (const auto& entry : kmers) {
            string kmer = entry.first;
            auto id = entry.second;
            for (auto const c : Alphabet) {
                auto next_id = forward(id, c);
                if (next_id != -1) {
                    // add next_id to adj of current node
                    adj[id].push_back(next_id);
                    inv_adj[next_id].push_back(id);
                }
            }
            if (total_kmers <= 100 || ++current_kmer % (total_kmers / 100) == 0)
                print_progress_bar(current_kmer, total_kmers, "Constructing de Bruijn graph");
        }
        cout << endl;
    }

    // bfs for an augmenting path
    bool bfs (const VVINT& adj, VINT& dist, VINT& match_u, VINT& match_v, INT n) {
        queue<INT> q;
        bool found_augpath = false;

        // add unmatched u in U to queue as a starting point of bfs
        for (INT u = 0; u < n; ++u) {
            if (match_u[u] == -1) {
                dist[u] = 0;
                q.push(u);
            } else {
                dist[u] = INF;
            }
        }

        while (!q.empty()) {
            INT u = q.front();
            q.pop();

            // bfs from U to V
            for (auto v : adj[u]) {
                if (match_v[v] == -1) {
                    found_augpath = true; // unmatched v gives an augpath
                } else if (dist[match_v[v]] == INF) {
                    dist[match_v[v]] = dist[u] + 1;
                    q.push(match_v[v]); // proceed to next search step
                }
            }
        }
        return found_augpath;
    }

    // update M by searching for augpath using DFS
    bool dfs (const VVINT& adj, VINT& dist, VINT& match_u, VINT& match_v, INT u) {
        for (auto v : adj[u]) {
            if (match_v[v] == -1 || (dist[match_v[v]] == dist[u] + 1 && dfs(adj, dist, match_u, match_v, match_v[v]))) {
                // update matching
                match_u[u] = v;
                match_v[v] = u;
                return true;
            }
        }
        dist[u] = INF; // not found any augpath
        return false;
    }

    // solve maximum bipartite matching
    INT hopcroft_karp() {
        INT n = kmers.size();
        VINT match_u(n, -1); // matching of U
        VINT match_v(n, -1); // matching of V
        VINT dist(n); // distance of bfs being executed

        // update M as long as augpath exists
        INT matching_size = 0;
        INT itr = 0;
        while (bfs(adj, dist, match_u, match_v, n)) {
            for (INT u = 0; u < n; ++u) {
                if (match_u[u] == -1 && dfs(adj, dist, match_u, match_v, u)) {
                    matching_size++; // update M when found an augpath
                }
            }
            if (n < 1000 || itr % 1000 == 0)
                print_progress_bar(++itr, static_cast<INT>(sqrt(n)), "Finding maximum matching");
        }
        cout << "\rFound maximum matching of size " << matching_size << endl; itr = 0;

        // decompose the matching into cycles and paths
        DINT active;
        vector<uint8_t> vis(n, 0);
        for (INT u = 0; u < n; ++u) {
            if (vis[u]) continue;
            active.push_back(u);
            vis[u] = 1;
            auto utmp = match_u[u]; //sada, next to u
            while (utmp != -1 && !vis[utmp]) {
                active.push_back(utmp);
                vis[utmp] = 1;
                utmp = match_u[utmp];
            }
            if (utmp != -1 && vis[utmp]) {
                // cycle detected
                VINT ctmp(active.begin(), active.end());
                cycles.push_back(ctmp);
            } else {
                utmp = match_v[u];
                while (utmp != -1 && !vis[utmp]) {
                    active.push_front(utmp);
                    vis[utmp] = 1;
                    utmp = match_v[utmp];
                }
                VINT ptmp(active.begin(), active.end());
                paths.push_back(ptmp);
            }
            active.clear();

            if (n <= 100 || ++itr % (n / 100) == 0)
                print_progress_bar(u, n, "Matching -> cycles and paths");
        }
        cout << endl;

        // add unmatched nodes to paths
        // for (INT u = 0; u < n; u++) {
        //     if (visited[u]) continue; // skip already matched nodes
        //     paths.push_back({u});
        // }

        // debug
        INT id = 0;
        for (const auto& cycle: cycles) {
            cout << "cycle" << id++ << ": ";
            for (const auto& node : cycle) {
                cout << sq[idPos[node] + K - 1];
            } cout << endl;
        } cout << endl;

        // point from paths to paths and cycles
        vector<tuple<INT, INT, INT>> pc;
        unordered_map<pair<INT, INT>, INT, pair_hash> pp;
        for (INT i = 0; i < static_cast<INT>(paths.size()); ++i) {
            auto head = paths[i][0];
            if (head == 0 && match_v[head] == -1 && inv_adj[head].empty()) {
                pp[make_pair(i, 0)] = i;
                goto label;
            }
            for (INT j = 0; j < static_cast<INT>(cycles.size()); ++j) {
                for (INT k = 0; k < static_cast<INT>(cycles[j].size()); ++k) {
                    for (const auto& inneigh : inv_adj[head]) {
                        if (cycles[j][k] == inneigh) {
                            pc.push_back(make_tuple(i, j, k));
                            goto label;
                        }
            }}}
            for (INT j = 0; j < static_cast<INT>(paths.size()); ++j) {
                for (INT k = 0; k < static_cast<INT>(paths[j].size()); ++k) {
                    for (const auto& inneigh : inv_adj[head]) {
                        if (paths[j][k] == inneigh) {
                            pp[make_pair(j, k)] = i;
                            goto label;
                        }
            }}}
            label:
            if (paths.size() < 100 || i % (paths.size() / 100) == 0)
                print_progress_bar(i, paths.size(), "Pointing paths to its preds");
        }
        
        // sort pc by cycles' index-wise ASC order
        sort(pc.begin(), pc.end(), [](const tuple<INT, INT, INT>& a, const tuple<INT, INT, INT>& b) {
            if (get<1>(a) == get<1>(b))
                return get<2>(a) < get<2>(b);
            return get<1>(a) < get<1>(b);
        });

        // sort paths
        const VINT cumls = get_chead(cycles);  // compute cumulative length up to every cycle in rep
        VINT sorted_paths;
        vector<uint8_t> added(paths.size(), 0);
        // step 1: add paths pointing to cycles
        for (INT i = 0; i < static_cast<INT>(pc.size()); ++i) {
            INT path_id = get<0>(pc[i]);
            sorted_paths.push_back(path_id);
            pointers.push_back(cumls[get<1>(pc[i])] + get<2>(pc[i]));
            added[path_id] = 1;
        }
        // step 2: process paths pointing to other paths
        INT offset = cumls.back();
        for (INT i = 0; i < static_cast<INT>(pc.size()); ++i) {
            INT path_id = get<0>(pc[i]);
            for (INT j = 0; j < static_cast<INT>(paths[path_id].size()); ++j) {
                auto it = pp.find(make_pair(path_id, j));                
                if (it != pp.end()) {
                    sorted_paths.push_back(it->second);
                    pointers.push_back(offset + get_phead(path_id, sorted_paths, paths) + j);
                    added[it->second] = 1;
                }
        }}
        // step 3: add any unadded paths
        for (const auto& _pp : pp) {
            auto path_id = _pp.second;
            for (INT j = 0; j < static_cast<INT>(paths[path_id].size()); ++j) {
                auto it = pp.find(make_pair(path_id, j));
                if (it != pp.end() && added[it->second] == 0) {
                    sorted_paths.push_back(it->second);
                    pointers.push_back(offset + get_phead(path_id, sorted_paths, paths) + j);
                    added[it->second] = 1;
                }
        }}
        INT false_count = count(added.begin(), added.end(), 0);
        if (false_count >= 1) {
            if (false_count > 1) {
                cerr << "exist" << false_count << "unadded path" << endl;
            } else {
                INT idx = distance(added.begin(), find(added.begin(), added.end(), 0));
                sorted_paths.push_back(idx);
                pointers.push_back(offset + get_phead(idx, sorted_paths, paths));
            }
        }

        // debug
        cout << "Added complete ?" <<endl;
        for (INT i = 0; i < static_cast<INT>(added.size()); ++i) {
            cout << "path" << i << ": " << (added[i] ? "added" : "not added") << endl;
        } cout << endl;

        // debug
        cout << "sorted_paths = ";
        for (const auto& sorted_path : sorted_paths) {
            cout << sorted_path << " ";
        } cout << endl;

       for (const auto& cycle : cycles) {
        for (const auto& node : cycle) {
            rep += sq[idPos[node] + K - 1];
        } rep += "$";
       }
       for (const auto& sorted_path : sorted_paths) {
        auto it = pp.find(make_pair(sorted_path, 0));
        if (it != pp.end() && it->second == sorted_path) rep += "$";
        for (const auto& node : paths[sorted_path]) {
            rep += sq[idPos[node] + K - 1];
        } rep += "$";
       } rep.pop_back();

        return matching_size;
    }

};

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " [K] [input.fa] [output((no extension)]" << endl;
        return 1;
    }
    INT _K = atoi(argv[1]);
    string _in_filename = argv[2];
    string _out_filename = argv[3];

    // build node-centric dBG
    DeBruijnGraph ncdbg = DeBruijnGraph(_in_filename, _K, true);
    ncdbg.getKmers();
    ncdbg.addEdges();

    auto start_time = chrono::high_resolution_clock::now();
    ncdbg.hopcroft_karp();
    auto end_time = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time);
    cout << "Hopcroft-Karp alg. completed in " << duration.count() << " sec" << endl;

    // write out txt file
    string txt_filename = _out_filename + ".txt";
    ofstream txtFile(txt_filename);
    if (!txtFile) {
        cerr << "Error: Could not open file " << txt_filename << " for writing.\n";
        return 1;
    }
    if (ncdbg.rep.empty()) cerr << "Warning: rep is empty.\n";
    if (ncdbg.pointers.empty()) cerr << "Warning: pointers are empty.\n";

    txtFile << ncdbg.rep << "\n";
    for (const auto& pointer : ncdbg.pointers) {
        txtFile << pointer << " ";
    }
    txtFile << "\n";
    txtFile.close();
    cout << "Text file created: " << txt_filename << endl;

    return 0;
}