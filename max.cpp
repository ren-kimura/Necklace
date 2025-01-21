#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <vector>
#include <array> // sada
#include <tuple> // sada
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <chrono>
#include <math.h>

using namespace std;
using INT = int64_t;
using Kmers = unordered_map<string,INT>;
using VINT = vector<INT>;
using VVINT = vector<VINT>;

struct pair_hash {
    template <class T1, class T2>
    size_t operator ()(const pair<T1, T2>& p) const {
        auto h1 = hash<T1>{}(p.first);
        auto h2 = hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

void progress(INT now, INT total, const string& task) {
    double progress = static_cast<double>(now) / static_cast<double>(total);
    int width = 30;
    int pos = static_cast<int>(width * progress);
    cout << "\r" << task << " [";
    for (int i = 0; i < width; ++i) {
        if (i < pos)        cout << "=";
        else if (i == pos)  cout << ">";
        else                cout << " ";
    }
    cout << "] " << static_cast<int>(progress * 100.0) << "%";
    cout.flush();
}

class DeBruijnGraph{
public:
    string in_filename;
    INT K, n;
    bool isNodeCentric;
    array<char, 4> Alphabet = {'A', 'C', 'G', 'T'};
    const INT INF = INT64_MAX;

    string sq;
    Kmers kmers; // pair <kmer string of length k, its unique ID in 0..N-1>, where N = #distinct kmers in sequence
    VINT idPos; // for each ID a position i such that ID corresponds to kmer sequence.substr(i, K)
    VVINT adj, inv_adj;
    VINT match_u, match_v;

    VVINT cycles, paths;
    string rep;
    VINT pointers;

    DeBruijnGraph(string _in_filename, INT _K, bool _isNodeCentric): in_filename(_in_filename), K(_K), isNodeCentric(_isNodeCentric)
    {};

    void to_uppercase(string& str) {
        transform(str.begin(), str.end(), str.begin(),
                    [](unsigned char c) {return toupper(c);});
    }

    void getKmers() {
        ifstream inputFile(in_filename);
        if (!inputFile) {
            cerr << "Error opening input file." << "\n";
            exit(0);
        }
        string line;
        while (getline(inputFile, line))
            if (line[0] != '>') sq += line; // skip header line
        inputFile.close();

        to_uppercase(sq);

        sq.erase(remove_if(sq.begin(), sq.end(), [this](char x)
            {return find(Alphabet.begin(), Alphabet.end(), x) == Alphabet.end();}
            ), sq.end());

        INT total_kmers = sq.length() - K + 1;
        for (INT i = 0; i <= total_kmers - 1; ++i) {
            auto kmer = sq.substr(i, K);
            if (kmers.find(kmer) == kmers.end()) { // newly found kmer
                idPos.push_back(i);
                kmers[kmer] = idPos.size() - 1;
            }
            if (total_kmers <= 100 || i % (total_kmers / 100) == 0)
                progress(i, total_kmers, "Detecting k-mers");
        }
        cout << "\n";
        n = kmers.size();
        adj.resize(n);
        inv_adj.resize(n);
    }

    INT forward(INT id, char c) {
        if (id >= static_cast<INT>(idPos.size())) return -1;
        auto pos = idPos[id];
        if (pos + 1 >= static_cast<INT>(sq.length())) return -1;
 
        string suffix = sq.substr(pos + 1, K - 1);
        auto next = suffix + c;
        // next exists AND not current
        if (kmers.find(next) != kmers.end() && kmers[next] != id)
            return kmers[next];
        return -1; // no branch to c
    }

    void addEdges() {
        INT cnt = 0;
        for (const auto& entry : kmers) {
            string kmer = entry.first;
            auto id = entry.second;
            for (auto const c : Alphabet) {
                auto next_id = forward(id, c);
                if (next_id != -1 && next_id != id) {
                    adj[id].push_back(next_id);
                    inv_adj[next_id].push_back(id);
                }
            }
            if (n <= 100 || ++cnt % (n / 100) == 0)
                progress(cnt, n, "Constructing de Bruijn graph");
        }
        cout << "\n";
    }

    bool bfs (VINT& dist) {
        queue<INT> q;
        bool found_augpath = false;

        for (INT u = 0; u < n; ++u) {
            if (match_u[u] == -1) {
                dist[u] = 0;
                q.push(u);
            } else dist[u] = INF;
        }

        while (!q.empty()) {
            INT u = q.front();
            q.pop();
            for (auto v : adj[u]) {
                if (match_v[v] == -1) {
                    found_augpath = true;
                } else if (dist[match_v[v]] == INF) {
                    dist[match_v[v]] = dist[u] + 1;
                    q.push(match_v[v]);
                }
            }
        }
        return found_augpath;
    }

    // update M by searching for augpath using DFS
    bool dfs (VINT& dist, INT u) {
        for (auto v : adj[u]) {
            if (match_v[v] == -1 || (dist[match_v[v]] == dist[u] + 1 && dfs(dist, match_v[v]))) {
                match_u[u] = v;
                match_v[v] = u;
                return true;
            }
        }
        dist[u] = INF; // not found any augpath
        return false;
    }

    // maximum bipartite matching
    INT hopcroft_karp() {
        VINT dist(n); // distance of bfs being executed
        match_u.assign(n, -1);
        match_v.assign(n, -1);

        INT m = 0;
        while (bfs(dist)) {
            for (INT u = 0; u < n; ++u) {
                if (match_u[u] == -1 && dfs(dist, u))
                    m++; // update M when found an augpath
            }
            if (n < 1000 || m % 1000 == 0)
                progress(m, static_cast<INT>(sqrt(n)), "Maximum matching");
        }
        cout << "\rFound maximum matching of size " << m << "\n";

        return m;
    }

    void dfs_forward(INT start, vector<int8_t>& seen, VINT& trail) {
        INT current = start;
        while(1) {
            if (seen[current]) break;
            seen[current] = 1;
            trail.push_back(current);
            if (match_u[current] == -1) break;
            current = match_u[current];
        }
    }

    void dfs_backward(INT start, vector<int8_t>& seen, VINT& back_trail) {
        INT current = start;
        while(1) {
            if (seen[current]) break;
            seen[current] = 1;
            back_trail.push_back(current);
            if (match_v[current] == -1) break;
            current = match_v[current];
        }
    }

    void decompose() {
        vector<int8_t> seen(n, 0);
        for (INT u = 0; u < n; ++u) {
            if (seen[u]) continue;
            VINT trail;
            dfs_forward(u, seen, trail);
            if (trail.size() > 1 && match_u[trail.back()] == trail.front()) cycles.push_back(trail);
            else {
                VINT back_trail;
                dfs_backward(trail.front(), seen, back_trail);
                if (!back_trail.empty()) {
                    reverse(back_trail.begin(), back_trail.end());
                    back_trail.pop_back();
                    trail.insert(trail.begin(), back_trail.begin(), back_trail.end());
                }
                paths.push_back(trail);
            }
            if (n < 1000 || u % 1000 == 0)
                progress(u, n, "Decomposing");
        }
        if (has_dup(cycles, paths))
            cout << "Duplicates found.\n";
        else
            cout << "No duplicates found.\n";
    }

    bool has_dup(const VVINT& cycles, const VVINT& paths) {
        unordered_set<INT> seen;
        auto check = [&](const VVINT& v) {
            for (const auto& vec: v) {
                for (auto x: vec) {
                    if (!seen.insert(x).second) return true;
                }
            } return false;
        };
        return check(cycles) || check(paths);
    }

    void cp_bare() {
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                rep += sq[idPos[node] + K - 1];
            } rep += "$";
        } rep += "$";
        for (const auto& path: paths) {
            for (const auto& node: path) {
                if (node == path.front())
                    rep += sq.substr(idPos[node], K - 1);
                rep += sq[idPos[node] + K - 1];
            } rep += "$";
        } if (!rep.empty()) rep.pop_back();
    }

    void cp_unsorted() {
        INT p = paths.size(), z = cycles.size() + p + n;
        pointers.resize(paths.size(), -1);
        unordered_map<INT, INT> heads;
        INT offset;
        // record paths' heads
        for (INT i = 0; i < static_cast<INT>(paths.size()); ++i)
            heads[paths[i][0]] = i;
        // scan nodes in cycles one by one to find paths pointing to them
        INT pos = 0;
        bool self = false;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                for (const auto& c: Alphabet) {
                    auto next = forward(node, c);
                    if (next == -1) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pointers[heads[next]] = pos;
                    heads.erase(it);
                    if (heads.empty()) goto end;
                }
                pos++;
                if (z < 1000 || pos % 1000 == 0)
                    progress(pos, z, "Pointing");
            }
            pos++;
        } offset = pos;
        for (const auto& path: paths) {
            for (const auto& node: path) {
                for (const auto& c: Alphabet) {
                    auto next = forward(node, c);
                    if (next == -1) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pointers[heads[next]] = pos;
                    heads.erase(it);
                    if (heads.empty()) goto end;
                }
                pos++;
                if (z < 1000 || pos % 1000 == 0)
                    progress(pos, z, "Pointing");
            }
            pos++;
        }
        if (heads.size() > 0) {
            if (heads.size() > 1) {
                cerr << "Error: More than one self-pointing path.\n";
                return;
            }
            auto it = heads.begin();
            if (it->first != 0) {
                cerr << "Error: Invalid self-pointing path.\n";
                return;
            }
            self = true;
            pointers[heads[0]] = offset;
            for (auto& to: pointers) 
                if (to > offset) to++;
            heads.erase(it);
        }
        end:;
        if (heads.empty()) 
            cout << "All paths are pointed.\n";
        // generate representation
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                rep += sq[idPos[node] + K - 1];
            } rep += "$";
        } for (const auto& path: paths) {
            for (const auto& node: path) {
                if (node == 0 && self) rep += "$";
                rep += sq[idPos[node] + K - 1];
            } rep += "$";
        } if(!rep.empty()) rep.pop_back();
        // take difference of pointers
        to_diff(pointers);
    }

    void cp_sorted() {
        VINT pord;
        unordered_map<INT, INT> heads;
        // record paths' heads
        for (INT i = 0; i < static_cast<INT>(paths.size()); ++i)
            heads[paths[i][0]] = i;
        // scan nodes in cycles one by one to find paths pointing to them
        INT pos = 0;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                for (const auto& c: Alphabet) {
                    auto next = forward(node, c);
                    if (next == -1) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pord.push_back(it->second);
                    pointers.push_back(pos);
                    heads.erase(it);
                    if (heads.empty()) goto end;
                }
                pos++;
            }
            pos++;
        }
        // handle a self-pointing path if exists
        if (pord.size() == 0) pord.push_back(pos);
        // add remaining paths
        for (const auto& id: pord) {
            for(const auto& node: paths[id]) {
                for (const auto& c: Alphabet) {
                    auto next = forward(node, c);
                    if (next == -1) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pord.push_back(it->second);
                    pointers.push_back(pos);
                    heads.erase(it);
                    if (heads.empty()) goto end;
                }
                pos++;
            }
            pos++;
        }
        // handle a self-pointing path if exists
        /*パスだけならノード0はパス0に含まれているはず*/
        if (heads.size() > 0) {
            if (heads.size() > 1) {
                cerr << "Error: More than one self-pointing path.\n";
                return;
            }
            auto it = heads.begin();
            if (it->first != 0) {
                cerr << "Error: Invalid self-pointing path.\n";
                return;
            }
            pord.push_back(it->second);
            pointers.push_back(pos);
        }
        end:;
        // generate representation
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                rep += sq[idPos[node] + K - 1];
            } rep += "$";
        } for (const auto& id: pord) {
            for (const auto& node: paths[id]) {
                if (node == 0) rep += "$";
                rep += sq[idPos[node] + K - 1];
            } rep += "$";
        } if(!rep.empty()) rep.pop_back();
    }

    void to_diff(VINT& v) {
        if (v.size() < 2) return;
        INT prev = v[0];
        for (INT i = 1; i < static_cast<INT>(v.size()); ++i) {
            INT tmp = v[i];
            v[i] -= prev;
            prev = tmp;
        }
    }
};

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " [K] [input.fa] [output((no extension)]" << "\n";
        return 1;
    }
    INT _K = atoi(argv[1]);
    string _in_filename = argv[2];
    string _out_filename = argv[3];

    // build node-centric dBG
    DeBruijnGraph ncdbg = DeBruijnGraph(_in_filename, _K, true);
    auto start0 = chrono::high_resolution_clock::now();
    ncdbg.getKmers();
    ncdbg.addEdges();

    auto start1 = chrono::high_resolution_clock::now();
    INT m = ncdbg.hopcroft_karp();
    auto end1 = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration_cast<chrono::milliseconds>(end1 - start1);
    cout << "Hopcroft-Karp alg. completed in " << duration1.count() << " ms" << "\n";

    auto start2 = chrono::high_resolution_clock::now();
    ncdbg.decompose();
    auto end2 = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration_cast<chrono::milliseconds>(end2 - start2);
    cout << "Decomposition completed in " << duration2.count() << " ms" << "\n";

    // ncdbg.cp_bare();
    ncdbg.cp_unsorted();
    // ncdbg.cp_sorted();
    auto end0 = chrono::high_resolution_clock::now();
    auto duration0 = chrono::duration_cast<chrono::milliseconds>(end0 - start0);
    cout << "Total time elapsed: " << duration0.count() << " ms" << "\n";

    // write out string representation
    string str_filename = _out_filename + "_str.txt";
    ofstream strFile(str_filename);
    if (!strFile) {
        cerr << "Error: Could not open file " << str_filename << " for writing.\n";
        return 1;
    }
    if (ncdbg.rep.empty()) cerr << "Warning: rep is empty.\n";
    strFile << ncdbg.rep << "\n";
    strFile.close();
    cout << "Text file created: " << str_filename << "\n";

    string pnt_filename = _out_filename + "_pnt.txt";
    ofstream pntFile(pnt_filename);
    if (!pntFile) {
        cerr << "Error: Could not open file " << pnt_filename << " for writing.\n";
        return 1;
    }
    if (ncdbg.pointers.empty()) cerr << "Warning: pointers are empty.\n";
    for (const auto& pointer : ncdbg.pointers) {
        pntFile << pointer << " ";
    }
    pntFile.close();
    cout << "Text file created: " << pnt_filename << "\n";
    cout << "(#kmers, #cycles, #paths, #matching) = (" << ncdbg.n << ", " << ncdbg.cycles.size() << ", " << ncdbg.paths.size() << ", " << m << ")\n";

    return 0;
}