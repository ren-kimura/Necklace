#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <vector>
#include <queue>
#include <array> // sada
#include <tuple> // sada
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <chrono>

using namespace std;
using INT = int64_t;
using Kmers = unordered_map<string,INT>;
using VSTR = vector<string>;
using Vint = vector<int8_t>;
using VINT = vector<INT>;
using PINT = pair<INT, INT>;
using VPINT = vector<PINT>;
using VTINT = vector<tuple<INT, INT, INT, INT>>;
using VVINT = vector<VINT>;
using REP = pair<string, VINT>;

unordered_set<char> base = {'A', 'C', 'G', 'T'};
const INT INF = INT64_MAX;

void progress(INT now, INT total, const string& task) {
    if (total >= 100 && now % (total / 100) != 0) return;
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

void remove_suffix(string& str, const string& suffix) {
    if (str.size() >= suffix.size() && str.rfind(suffix) == str.size() - suffix.size()) {
        str.erase(str.size() - suffix.size());
    }
}

class DeBruijnGraph{
public:
    string filename;
    INT K;
    INT option;
    bool is_node_centric;

    DeBruijnGraph(string _filename, INT _K, INT _option, bool _is_node_centric): filename(_filename), K(_K), option(_option), is_node_centric(_is_node_centric)
    {};

    REP process() {
        VSTR reads;
        Kmers kmers;
        VPINT idpos;        
        INT N = get_kmers(reads, kmers, idpos);

        VVINT adj(N), inv_adj(N);
        INT E = add_edges(reads, kmers, idpos, N, adj, inv_adj);

        VINT match_u, match_v;
        INT M = hopcroft_karp(N, adj, match_u, match_v);

        VVINT cycles, paths;
        PINT CnP = decompose(N, match_u, match_v, cycles, paths);
        // debug
        cout << "(#k-mers, #edges, #matching, #cycles, #paths) = (" 
            << N << ", " << E << ", " << M << ", "
            << CnP.first << ", " << CnP.second << ")\n";
        if (option == 0) return plain(reads, idpos, cycles, paths);
        else if (option == 1) return unsorted(reads, kmers, N, idpos, cycles, paths, CnP);
        else if (option == 2) return sorted(reads, kmers, N, idpos, cycles, paths, CnP);
        else cerr << "Error: Invalid option value.\n";
        return {"",{}};
    }

    void to_uppercase(VSTR& strs) {
        for (auto& str: strs)
            transform(str.begin(), str.end(), str.begin(),
                    [](unsigned char c) {return toupper(c);});
    }

    INT get_kmers(VSTR& reads, Kmers& kmers, VPINT& idpos) {
        ifstream inputFile(filename, ios::in | ios::binary);
        if (!inputFile) {
            cerr << "Error opening input file." << "\n";
            exit(1);
        }
        const size_t BUF_SIZE = 64 * 1024 * 1024; // 64MB
        char *buffer = new char[BUF_SIZE];
        inputFile.rdbuf()->pubsetbuf(buffer, BUF_SIZE);
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
        for (string& read: reads) {
            read.erase(remove_if(read.begin(), read.end(), 
                [this](char x) {return base.find(x) == base.end();}), read.end());
            tlen += read.size();
        }
        INT id = 0, footing = 0;
        for (auto i = 0; i < static_cast<INT>(reads.size()); ++i) {
            auto read = reads[i];
            INT len = read.size();
            if (len < K) continue;
            INT t = len - K + 1;
            string crnt = read.substr(0, K);
            if (kmers.find(crnt) == kmers.end()) {
                idpos.emplace_back(i, 0);
                kmers[crnt] = static_cast<INT>(id++);
            }
            for (INT j = 1; j < t; ++j) {
                crnt.erase(0, 1);
                crnt.push_back(read[j + K - 1]);
                if (kmers.find(crnt) == kmers.end()) {
                    idpos.emplace_back(i, j);
                    kmers[crnt] = static_cast<INT>(id++);
                }
                progress(footing + j, tlen, "Detecting k-mers");
            }
            footing += len;
        }
        cout << "\n";
        INT N = kmers.size();
        cout << "Total number of k-mers: " << N << "\n";
        return N;
    }

    INT forward(const VSTR& reads, Kmers& kmers, 
                const VPINT& idpos, INT id, char c) {
        auto pos = idpos[id]; 
        string suffix = reads[pos.first].substr(pos.second + 1, K - 1);
        auto next = suffix + c;
        if (kmers.find(next) != kmers.end() && kmers[next] != id)
            return kmers[next];
        return -1; // no branch to c
    }

    INT add_edges(const VSTR& reads, Kmers& kmers, 
                const VPINT& idpos, const INT& N, VVINT& adj, VVINT& inv_adj) {
        INT cnt = 0, E = 0;
        adj.resize(N);
        inv_adj.resize(N);
        for (const auto& entry : kmers) {
            string kmer = entry.first;
            auto id = entry.second;
            for (auto const c : base) {
                auto next_id = forward(reads, kmers, idpos, id, c);
                if (next_id != -1 && next_id != id) {
                    adj[id].emplace_back(next_id);
                    inv_adj[next_id].emplace_back(id);
                    ++E;
                }
            }
            progress(++cnt, N, "Constructing de Bruijn graph");
        }
        cout << "\nTotal number of edges: " << E << "\n";
        return E;
    }

    bool bfs (const INT& N, const VVINT& adj, 
            const VINT& match_u, const VINT& match_v, VINT& dist) {
        queue<INT> q;
        bool found_augpath = false;
        for (INT u = 0; u < N; ++u) {
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

    bool dfs (const VVINT& adj, VINT& match_u, VINT& match_v, 
            VINT& dist, INT u) {
        for (auto v : adj[u]) {
            if (match_v[v] == -1 || (dist[match_v[v]] == dist[u] + 1 && dfs(adj, match_u, match_v, dist, match_v[v]))) {
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
        match_u.assign(N, -1);
        match_v.assign(N, -1);
        while (bfs(N, adj, match_u, match_v, dist)) {
            for (INT u = 0; u < N; ++u) {
                if (match_u[u] == -1 && dfs(adj, match_u, match_v, dist, u))
                    M++; // update M when found an augpath
            }
            progress(M, N, "Maximum matching");
        }
        cout << "\nFound maximum matching of size " << M << "\n";
        return M;
    }

    VINT dfs_forward(const VINT& match_u, const INT& start, Vint& seen) {
        VINT trl;
        INT crnt = start;
        while(1) {
            if (seen[crnt]) break;
            seen[crnt] = 1;
            trl.emplace_back(crnt);
            if (match_u[crnt] == -1) break;
            crnt = match_u[crnt];
        }
        return trl;
    }

    VINT dfs_backward(const VINT& match_v, const INT& start, Vint& seen) {
        VINT lrt;
        INT crnt = start;
        while(1) {
            if (seen[crnt]) break;
            seen[crnt] = 1;
            lrt.emplace_back(crnt);
            if (match_v[crnt] == -1) break;
            crnt = match_v[crnt];
        }
        return lrt;
    }

    PINT decompose(const INT& N, const VINT& match_u, const VINT& match_v,
                VVINT& cycles, VVINT& paths) {
        Vint seen(N, 0);
        for (INT u = 0; u < N; ++u) {
            if (seen[u]) continue;
            VINT trl = dfs_forward(match_u, u, seen);
            if (trl.size() > 1 && match_u[trl.back()] == trl.front())
                cycles.emplace_back(trl);
            else {
                VINT lrt = dfs_backward(match_v, trl.front(), seen);
                if (!lrt.empty()) {
                    reverse(lrt.begin(), lrt.end());
                    lrt.pop_back();
                    trl.insert(trl.begin(), lrt.begin(), lrt.end());
                }
                paths.emplace_back(trl);
            }
            progress(u, N, "Decomposing");
        }
        if (has_dup(cycles, paths))
            cout << "\nDuplicates found.\n";
        else
            cout << "\nNo duplicates found.\n";
        INT C = cycles.size(), P = paths.size();
        cout << "Found " << C << " cycles and " << P << " paths.\n";
        return {C, P};
    }

    bool has_dup(const VVINT& cycles, const VVINT& paths) {
        unordered_set<INT> seen;
        auto check = [&](const VVINT& v) {
            for (const auto& vec: v)
                for (const auto& x: vec)
                    if (!seen.insert(x).second)
                        return true;
            return false;
        };
        return check(cycles) || check(paths);
    }

    REP plain(const VSTR& reads, const VPINT& idpos,
                const VVINT& cycles, const VVINT& paths) {
        string txt;
        VINT pnt;
        cout << "Generating SPSS...";
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                auto x = idpos[node];
                txt += reads[x.first][x.second + K - 1];
            } txt += "$";
        } txt += "$";
        for (const auto& path: paths) {
            for (const auto& node: path) {
                auto x = idpos[node];
                if (node == path.front())
                        txt += reads[x.first].substr(x.second, K - 1);
                txt += reads[x.first][x.second + K - 1];
            } txt += "$";
        } if (!txt.empty()) txt.pop_back();
        cout << "\rCompleted generating SPSS.\n";
        return {txt, {}};
    }

    REP unsorted(const VSTR& reads, Kmers& kmers, const INT& N, const VPINT& idpos,
                const VVINT& cycles, const VVINT& paths, const PINT& CnP) {
        INT P = CnP.second;  
        INT S = CnP.first + P + N;      
        unordered_map<INT, INT> heads;
        for (INT i = 0; i < P; ++i)
            heads[paths[i][0]] = i; // record paths' heads
        INT offset = 0, pos = 0;
        bool self = false;
        VINT pnt(P, -1);
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                for (const auto& c: base) {
                    auto next = forward(reads, kmers, idpos, node, c);
                    if (next == -1) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pnt[heads[next]] = pos;
                    heads.erase(it);
                    if (heads.empty()) goto end;
                } pos++;
                progress(pos, S, "Pointing");
            } pos++;
        } offset = pos;
        for (const auto& path: paths) {
            for (const auto& node: path) {
                for (const auto& c: base) {
                    auto next = forward(reads, kmers, idpos, node, c);
                    if (next == -1) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pnt[heads[next]] = pos;
                    heads.erase(it);
                    if (heads.empty()) goto end;
                } pos++;
                progress(pos, S, "Pointing");
            } pos++;
        }
        if (heads.size() > 0) {
            auto it = heads.begin();
            self = true;
            pnt[heads[0]] = offset;
            for (auto& to: pnt) 
                if (to > offset) to++;
            heads.erase(it);
        }
        end:
        if (heads.empty()) 
            cout << "All paths are pointed.\n";
        // generate representation
        string txt;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                auto x = idpos[node];
                txt += reads[x.first][x.second + K - 1];
            } txt += "$";
        } for (const auto& path: paths) {
            for (const auto& node: path) {
                if (node == 0 && self) txt += "$";
                auto x = idpos[node];
                txt += reads[x.first][x.second + K - 1];
            } txt += "$";
        } if(!txt.empty()) txt.pop_back();        
        to_diff(pnt); // take difference of pointers
        return {txt, pnt};
    }

    REP sorted(const VSTR& reads, Kmers& kmers, const INT& N, const VPINT& idpos,
                VVINT& cycles, VVINT& paths, const PINT& CnP) {
        INT P = CnP.second;
        INT S = CnP.first + P + N, from = 0;
        bool self = false; // exists pointer from node 0 to node 0 ?
        unordered_set<INT> pointees;
        unordered_map<INT, INT> heads;
        for (INT i = 0; i < P; ++i)
            heads[paths[i][0]] = i; // record paths' heads
        VINT pord;

        // check each node in cycles if it can be pointed from paths
        INT pos = 0;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                for (const auto& c: base) {
                    auto next = forward(reads, kmers, idpos, node, c);
                    if (next == -1) continue;
                    auto it = heads.find(next);
                    if (it == heads.end()) continue;
                    pord.emplace_back(it->second);
                    pointees.insert(node);
                    heads.erase(it);
                    if (heads.empty()) goto end;
                } ++pos;
                progress(pos, S, "Pointing");
            } ++pos;
        } 

        while (from < P && !heads.empty()) {
            // check each node in paths in "pord" if it can be pointed from other paths
            auto bound = pord.size();
            for (auto i = from; i < min(static_cast<INT>(bound), P); ++i) {
                for (const auto& node: paths[pord[i]]) {
                    for (const auto& c: base) {
                        auto next = forward(reads, kmers, idpos, node, c);
                        if (next == -1) continue;
                        auto it = heads.find(next);
                        if (it == heads.end()) continue;
                        pord.emplace_back(it->second);
                        ++bound;
                        pointees.insert(node);
                        heads.erase(it);
                        if (heads.empty()) goto end;
                    } ++pos;
                    progress(pos, S, "Pointing");
                } ++pos;
            } 
            from = pord.size();

            // resolve unpointed paths
            for (auto it1 = heads.begin(); it1 != heads.end(); ++it1) {
                VINT new_cycle;
                VTINT cands; // candidates of (pointee, path to modify, pos to delete, head to delete)
                auto start = it1->second;
                auto crnt = start;
                do {
                    auto& path = paths[crnt];
                    for (auto i = 0; i < static_cast<INT>(path.size()); ++i) {
                        auto node = path[i];
                        new_cycle.emplace_back(node);
                        for (const auto& c: base) {
                            auto next = forward(reads, kmers, idpos, node, c);
                            if (next == -1) continue;
                            auto it2 = heads.find(next);
                            if (it2 == heads.end()) continue;
                            cands.emplace_back(node, crnt, i, path[0]);
                            crnt = it2->second;
                            goto next;
                        }
                    }
                    new_cycle.clear();
                    break;
                    next:;
                } while (crnt != start);

                if (new_cycle.empty()) continue;
                cycles.emplace_back(new_cycle);
                for (const auto& cand: cands) {
                    auto pointee = get<0>(cand);
                    auto path_id = get<1>(cand);
                    auto path_to_modify = paths[path_id];
                    auto pos_to_delete = get<2>(cand);
                    auto head_to_delete = get<3>(cand);

                    pointees.insert(pointee);
                    pord.emplace_back(path_id);
                    path_to_modify.erase(path_to_modify.begin(), path_to_modify.begin() + pos_to_delete);
                    auto it2 = heads.find(head_to_delete);
                    if (it2 != heads.end()) heads.erase(it2);
                    if (heads.empty()) goto end;
                }
            }
            if (heads.size() == 1 && heads.begin()->first == 0 && heads.begin()->second == 0) {
                auto it = heads.begin();
                self = true;
                pord.emplace_back(it->second);
                pointees.insert(0);
                heads.erase(it);
                goto end;
            }
        }
        end:
        if (heads.empty()) 
            cout << "All paths are pointed.\n";
        else cout << heads.size() << " paths are not pointed.\n";

        // generate representation
        string txt;
        VINT pnt;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                auto x = idpos[node];
                txt += reads[x.first][x.second + K - 1];
            } txt += "$";
        } for (const auto& path: paths) {
            for (const auto& node: path) {
                if (node == 0 && self) txt += "$";
                auto x = idpos[node];
                txt += reads[x.first][x.second + K - 1];
            } txt += "$";
        } if(!txt.empty()) txt.pop_back();
        pos = 0;
        for (const auto& cycle: cycles) {
            for (const auto& node: cycle) {
                if (pointees.find(node) != pointees.end())
                    pnt.emplace_back(pos);
                ++pos;
            } ++pos;
        } for (const auto& pid: pord) {
            for (const auto& node: paths[pid]) {
                if (pointees.find(node) != pointees.end())
                    pnt.emplace_back(pos);
                ++pos;
            } ++pos;
        }      
        to_diff(pnt); // take difference of pointers
        return {txt, pnt};
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

    void write(const REP& rep) {
        if (rep.first.empty()) {
            cerr << "Error: txt is empty.\n";
            return;
        }
        remove_suffix(filename, ".fa");
        string txtfilename = filename + ".ours.k" + to_string(K) + ".txt";
        ofstream txtfile(txtfilename);
        if (!txtfile) {
            cerr << "Error: Could not open file " << txtfilename << " for writing.\n";
            return;
        }
        txtfile << rep.first;
        txtfile.close();
        cout << "\nFile created: " << txtfilename << "\n";

        if (rep.second.empty()) {
            cerr << "Note: pnt is empty.\n";
            return;
        }
        string pntfilename = filename + ".ours.k" + to_string(K) + ".pnt.txt";
        ofstream pntfile(pntfilename);
        if (!pntfile) {
            cerr << "Error: Could not open file " << pntfilename << " for writing.\n";
            return;
        }
        for (const auto& p: rep.second)
            pntfile << p << " ";
        pntfile.close();
        cout << "\nFile created: " << pntfilename << "\n";
    }
};

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " [****.fa] [k] [0 / 1 / 2]\n"
             << "Output options are...\n"
             << "0: explicit representation without pointers (SPSS)\n"
             << "1: representation with unsorted pointers\n"
             << "2: representation with sorted pointers\n";
        return 1;
    }
    string _filename = argv[1];
    INT _K = atoi(argv[2]);
    INT _option = atoi(argv[3]);

    DeBruijnGraph ncdbg = DeBruijnGraph(_filename, _K, _option, true);
    REP rep = ncdbg.process();
    ncdbg.write(rep);
    return 0;
}