#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple> // sada
#include <array> // sada
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <cstring>

using namespace std;
using INT = int64_t;
using Kmers = unordered_map<string,INT>;
using VecInt = vector<INT>;
using VecVecInt = vector<VecInt>;

class DeBruijnGraph{
public:
    string filename;
    INT K;
    bool isNodeCentric; // node-centric or edge-centric dBG?
    array<uint8_t, 4> Alphabet = {'A', 'C', 'G', 'T'};
    const INT INF = INT64_MAX;

    string sq; // stores the concatenated sequence
    Kmers kmers; // pair <kmer string of length k, its unique ID in 0..N-1>, where N = #distinct kmers in sequence
    VecInt idPos; // for each ID a position i such that ID corresponds to kmer sequence.substr(i, K)
    VecVecInt adj, inv_adj;

    string rep; // representation string for a cover
    VecInt pointers; // pointer from and to a certain pos in rep

    DeBruijnGraph(string _filename, INT _K, bool _isNodeCentric): filename(_filename), K(_K), isNodeCentric(_isNodeCentric)
    {};

    void to_uppercase(string& str) {
        transform(str.begin(), str.end(), str.begin(),
                    [](unsigned char c) {return toupper(c);});
    }

    void getKmers() {
        // read input DNA sequence in FASTA format .fa
        ifstream inputFile(filename);
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
        for (INT i = 0; i <= sq.length() - K; i++) {
            auto kmer = sq.substr(i, K);
            if (kmers.find(kmer) == kmers.end()) { // newly found kmer
                idPos.push_back(i);
                kmers[kmer] = idPos.size() - 1;
            }
        }
        adj.resize(kmers.size()); // resize adjacency list
        inv_adj.resize(kmers.size());

    }

    INT forward(INT id, uint8_t c) {
        if (id >= idPos.size()) return -1;
        auto pos = idPos[id];

        if (pos + 1 >= sq.length()) return -1;
 
        string suffix = sq.substr(pos + 1, K - 1);
        auto next = suffix + static_cast<char>(c);
        // next exists AND not current
        if (kmers.find(next) != kmers.end() && kmers[next] != id)
            return kmers[next];
        else
            return -1; // no branch to c
    }

    void addEdges() {
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
        }
    }

    void printGraph() {
        // map current node to next nodes
        for (auto const &entry : kmers) {
            auto kmer = entry.first;
            auto id = entry.second;
            cout << kmer << " >>> ";

            for (auto const &c : Alphabet) {
                // rule out edges of k-mers not occuring in sq (ecdBG)
                if (!isNodeCentric && (sq.find(kmer + static_cast<char>(c)) == string::npos)) continue;

                auto nextid = forward(id, c);
                if (nextid != -1) {
                    auto nextpos = idPos[nextid];
                    cout << sq.substr(nextpos, K) << " ";
                }
            }
            cout << endl;
        }
    }

    // bfs for an augmenting path
    bool bfs (const VecVecInt& adj, VecInt& dist, VecInt& match_u, VecInt& match_v, INT n) {
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
    bool dfs (const VecVecInt& adj, VecInt& dist, VecInt& match_u, VecInt& match_v, INT u) {
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
        VecInt match_u(n, -1); // matching of U
        VecInt match_v(n, -1); // matching of V
        VecInt dist(n); // distance of bfs being executed
        
        INT matching_size = 0;

        // update M as long as augpath exists
        while (bfs(adj, dist, match_u, match_v, n)) {
            for (INT u = 0; u < n; ++u) {
                if (match_u[u] == -1 && dfs(adj, dist, match_u, match_v, u)) {
                    matching_size++; // update M when found an augpath
                }
            }
        }

        // store cycles and paths
        VecVecInt cycles;
        VecVecInt paths;
        VecInt current_segment;
        vector<bool> visited(n, false);

        // decompose the matching into cycles and paths
        for (INT u = 0; u < n; u++) {
            if (visited[u]) continue;
            auto utmp = u; //sada
            while (!visited[utmp] && utmp != -1) {
                current_segment.push_back(utmp);
                visited[utmp] = true;
                utmp = match_u[utmp];
            }
            if (utmp == -1) {
                // path ends here (not part of a cycle)
                paths.push_back(current_segment);
            } else {
                // cycle detected
                cycles.push_back(current_segment);
            }
            current_segment.clear();
        }

        // add unmatched nodes to paths
        for (INT u = 0; u < n; u++) {
            if (visited[u]) continue; // skip already matched nodes
            paths.push_back({u});
        }

        INT current_startpos = 0;
        VecVecInt sorted_paths;
        // then, for each cycle, check if it has pointers from paths
        for (const auto& cycle: cycles) {
            for (INT i = 0; i < cycle.size(); i++) {
                const auto& adjs = adj[cycle[i]];
                // check if cycle[i] has pointer from paths
                auto it = paths.begin();
                while (it != paths.end()) {
                    auto head = (*it)[0];
                    if (find(adjs.begin(), adjs.end(), head) != adjs.end()) {
                        sorted_paths.push_back(*it);
                        pointers.push_back(current_startpos + i);
                        it = paths.erase(it);
                    } else {
                        ++it;
                    }
                }
            }
            current_startpos += cycle.size() + 1;
        }

        // for each path, check if it has pointers from paths
        INT idx = 0;
        while (idx < sorted_paths.size()) {
            auto& sorted_path = sorted_paths[idx];
            for (INT i = 0; i < sorted_path.size(); i++) {
                const auto& adjs = adj[sorted_path[i]];
                // check if path[i] has pointer from paths
                auto it = paths.begin();
                while (it != paths.end()) {
                    auto head = (*it)[0];
                    if (find(adjs.begin(), adjs.end(), head) != adjs.end()) {
                        sorted_paths.push_back(*it);
                        pointers.push_back(current_startpos + i);
                        it = paths.erase(it);
                    } else {
                        ++it;
                    }
                }
            }
            current_startpos += sorted_path.size() + 1;
            idx++;
        }

        // if the very first k-mer doesn't repeat, add the path to the end of sorted_paths
        bool is_selfpointer = false;
        if (!paths.empty() && paths[0][0] == 0 && match_v[paths[0][0]] == -1) {
            sorted_paths.push_back(paths[0]);
            pointers.push_back(current_startpos + 1);
            is_selfpointer = true;
        }

        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                rep += sq[idPos[node] + K - 1];
            }
            rep += "$";
        }

        for (auto it = sorted_paths.begin(); it != sorted_paths.end(); ++it) {
            const auto& sorted_path = *it;
            if (is_selfpointer && it == prev(sorted_paths.end()))
                rep += "$";
            for (const auto& node : sorted_path) {
                rep += sq[idPos[node] + K - 1];
            }
            rep += "$";
        }
        
        rep.pop_back();        
        return matching_size;
    }
};

int main(int argc, char *argv[]) {
    string _filename;
    INT _K = 3;
    if (argc >= 2) {
      _filename = argv[1];
      if (argc >= 3) {
        _K = atoi(argv[2]);
      }
    } else {
      cout << "Input file: ";
      cin >> _filename;
      cout << "K = ";
      cin >> _K;
      cout << endl;
    }

    /*
    // build edge-centric dBG
    DeBruijnGraph ecdbg = DeBruijnGraph(_filename, _K - 1, false);
    ecdbg.getKmers();
    cout << "Edge-centric De Bruijn Graph:" << endl;
    ecdbg.printGraph();
    */

    // build node-centric dBG
    DeBruijnGraph ncdbg = DeBruijnGraph(_filename, _K, true);
    ncdbg.getKmers();
    ncdbg.addEdges();
    // cout << "Node-centric De Bruijn Graph:" << endl;
    // ncdbg.printGraph();

    ncdbg.hopcroft_karp();
    cout << "Processed sequence: " << endl;
    cout << ncdbg.rep << endl;
    cout << "Pointer array: " << endl;
    for (INT i = 0; i < ncdbg.pointers.size(); ++i) {
        cout << ncdbg.pointers[i] << " ";
    }
    cout << endl;

    return 0;
}