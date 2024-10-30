#include <iostream>
#include <vector>
#include <unordered_set>
#include <queue>
#include <cstdint> // uint64_t
#include <utility> // std::pair
#include <optional> // std::optional

using namespace std;

using Vertex = uint64_t;
using Vertices = unordered_set<Vertex>;
using Parent = vector<optional<Vertex>>;
using Queue = queue<Vertex>;

// Define a hash function for std::pair<Vertex, Vertex>
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ (std::hash<T2>()(pair.second) << 1);
    }
};

using Edge = unordered_set<pair<Vertex, Vertex>, pair_hash>;

class Bipartite {
public:
    Vertex size_u, size_v;
    Edge edges;
    Parent parent;

    Bipartite(Vertex _u, Vertex _v, Edge _edges) 
        : size_u(_u), size_v(_v), edges(_edges), parent(size_u + size_v, nullopt) {};

    Edge asm1() {
        Edge m;
        Queue q;

        for (Vertex i = 0; i < size_u; i++) {
            // Clear the queue and push the current vertex
            while (!q.empty()) q.pop();
            q.push(i);
            Vertices s = {i};

            optional<Vertex> bestV = nullopt;
            uint64_t minDeg = UINT64_MAX;

            // Reset parent before BFS
            fill(parent.begin(), parent.end(), nullopt);

            while (!q.empty()) {
                Vertex w = q.front();
                q.pop();
                Vertices neighs;

                // Find neighbors based on matching status
                if (w < size_u) {
                    neighs = neighbors(m, w, false);
                } else {
                    neighs = neighbors(m, w, true);
                    if (degM(m, w) < minDeg) {
                        bestV = w;
                        minDeg = degM(m, w);
                    }
                }

                // Update parents and add neighbors to the queue
                for (const auto& neigh : neighs) {
                    if (!parent[neigh].has_value()) {
                        parent[neigh] = w;
                        q.push(neigh);
                    }
                }
            }

            // Update the semi-matching M if a better vertex was found
            if (bestV) {
                Vertex v = *bestV;
                Vertex u = parent[v].value();

                // Debug output for matching process
                cout << "Matching vertex " << i << " to " << v << endl;

                m.emplace(u, v);
                while (u != i && u != -1) {
                    v = parent[u].value();
                    m.erase(make_pair(u, v));
                    u = parent[v].value();
                    m.emplace(u, v);
                }
            } else {
                cout << "No match found for vertex " << i << endl;  // Debug output for unmatched case
            }
        }
        return m;
    }

    Vertices neighbors(const Edge& m, Vertex w, bool isWantMatched) {
        Vertices neighs;
        // wが左側の頂点か右側の頂点かによって処理を分岐
        if (w < size_u) {  // 左側の頂点の場合
            for (Vertex v = 0; v < size_v; v++) {
                if (edges.find(make_pair(w, v)) != edges.end()) {
                    // マッチングされていない頂点を探す
                    bool isMatched = m.find(make_pair(w, v)) != m.end();
                    if ((isMatched && isWantMatched) || (!isMatched && !isWantMatched)) {
                        neighs.insert(v);
                    }
                }
            }
        } else {  // 右側の頂点の場合
            for (Vertex u = 0; u < size_u; u++) {
                if (edges.find(make_pair(u, w)) != edges.end()) {
                    bool isMatched = false;
                    for (const auto& edge : m) {
                        if (edge.first == u && edge.second == w) {
                            isMatched = true;
                            break;
                        }
                    }
                    if ((isMatched && isWantMatched) || (!isMatched && !isWantMatched)) {
                        neighs.insert(u);
                    }
                }
            }
        }
        return neighs;
    }


    uint64_t degM(const Edge& m, Vertex w) {
        uint64_t count = 0;
        for (const auto& edge : m) {
            if (edge.second == w) {  // wが右側の頂点の場合
                count++;
            }
        }
        return count;
    }

    void print(const Edge& m) {
        if (m.empty()) {
            cout << "No matching found" << endl; // Debug output when no matching exists
            return;
        }

        for (const auto& edge : m) {
            cout << "(" << edge.first << ", " << edge.second << ")" << endl;
        }
    }
};

int main() {
    Vertex _u = 6; // Number of vertices in U
    Vertex _v = 4; // Number of vertices in V
    Edge _edges = {
        {0, 0},
        {0, 2},
        {1, 1},
        {2, 0},
        {2, 2},
        {2, 3},
        {3, 0},
        {4, 1},
        {5, 3}
    };

    // Create an instance of the Bipartite class
    Bipartite bg = Bipartite(_u, _v, _edges);
    Edge semiM = bg.asm1(); // Compute semi-matching
    bg.print(semiM); // Print the matching result

    return 0;
}
