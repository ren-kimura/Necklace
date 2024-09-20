#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <deque>
#include <variant>

#define Alphabet {'A', 'C', 'G', 'T'}

using namespace std;
using AdjList = vector<vector<int64_t>>;
using Kmers = unordered_map<string,int64_t>;
using NodeId = vector<int64_t>;

using Path = pair<deque<pair<int64_t, int64_t>>, uint8_t>;
using Paths = pair<vector<pair<deque<pair<int64_t, int64_t>>, uint8_t>>, int64_t>;


class DeBruijnGraph {
public:
    NodeId idPosition;
    Kmers kmers;
    int K;
    bool isNodeCentric;
    string sequence;  // To store the concatenated sequence


    DeBruijnGraph(int _K, bool _isNodeCentric) : K(_K), isNodeCentric(_isNodeCentric)
    {};

    void getKmers(string filename) {

        // cout << "Designate an input file(****.fa): ";
        // cin >> filename;

        // Convert FASTA into a string
        ifstream inputFile(filename);
        if (!inputFile) {
            cerr << "Error opening input file." << endl;
            exit(0);
        }

        string line;
        sequence = "";  // To store the concatenated sequence

        while (getline(inputFile, line)) {
            if (line[0] != '>') {  // Skip header lines
                sequence += line;
            }
        }

        inputFile.close();

        // Erase letters other than "a,t,g,c"
        sequence.erase(remove_if(sequence.begin(), sequence.end(), [](char x) {return (x != 'a' && x != 't' && x != 'g' && x != 'c' && x != 'A' && x != 'T' && x != 'G' && x != 'C');}), sequence.end());

        // Output the result as a string
        // cout << "Processed Sequence: " << sequence << "\n" << endl;
        
        // Set k
        // cout << "Set k to: ";
        // cin >> k;

        // Check every k-mer in "sequence"
        for (int i = 0; i < sequence.length() - K + 1; i++){
            auto kmer = sequence.substr(i, K);
            if (kmers.find(kmer) == kmers.end()){
                kmers[kmer] = idPosition.size();
                idPosition.push_back(i);
            }
        }

        
        // Output all k-mers
        // cout << "Distinct k-mers: " << endl;
        // for (const string& kmer : kmers) {
        //     cout << kmer << " ";
        // }
        // cout << endl;

        // Number of distinct k-mers
        cout << "# k-mers: " << kmers.size() << "\n" << endl;
    }


    // Edges of the graph (current node, successive nodes)
    AdjList outNeighbors;
    // Edges of the graph (current node, previous nodes)
    AdjList inNeighbors;
    
    // Construct node-centric De Bruijn Graph from the given UPPERCASE k-mer sets
    void buildNodeCentricGraph( ) { 
        outNeighbors.resize(kmers.size());
        inNeighbors.resize(kmers.size());
        // Map current node to next nodes
        for (auto const &x : kmers {
            auto kmer = x.first;
            auto id = x.second;
            auto position = idPostion[id];

            // Define node as a current looking k-mer, then determine the edges using the last letter
            string kmerSuffix = sequence.substr(position + 1, K-1); // k-mer's last (k-1) letters
            for (auto c : Alphabet){
                auto next = kmerSuffix.append(c);
                if (kmers.find(next) != kmers.end())
                    auto nextid = kmers[next];
                    outNeighbors[id].push_back(nextid);
                    inNeighbors[nextid].push_back(id);           
            }

    }

    uint64_t forward(uint64_t id, uint8_t c){
        auto position = idPostion[id];
        string kmerSuffix = sequence.substr(position + 1, K-1); 
        auto next = kmerSuffix.append(c);
        if (kmers.find(next) != kmers.end())
            return kmers[next];
        else 
            return -1; // no branching exists
    }

    void buildEdgeCentricGraph(const string S, const Kmers& kmers, const NodeId& idPos) {  // Rough invariant: kmers.first == S.substr(kmers.second,k)
    }

    // Show the graph
    void printGraph(const vector<string>& kmers) {
        for (int i = 0; i < kmers.size(); i++) {
            cout << "Node: " << kmers[i] << " has edges to: ";
            for (int j = 0; j < nc_outAdjacencyList[i].size(); j++) {
                cout << kmers[nc_outAdjacencyList[i][j]] << " ";
            }
            cout << endl;
        }
    }

    Path greedyPath(const int& start, const NcAdjList& adjList, const NcAdjList& revAdjList, vector<uint8_t>& visited, vector<uint8_t>& running){
        Path path; // pair<deque<pair<int, int>>, uint8_t>
        path.first = {};
        path.second = 0;

        int current = start;
        path.first.push_back(make_pair(current, -1));
        visited[current] = 1;
        running[current] = 1;
        bool tmp = 0;

        // Forward tracking
        while (true) {
            bool tmp1 = false;
            for (int i = 0; i < adjList[current].size(); i++) {
                if (!visited[adjList[current][i]]) {
                    visited[adjList[current][i]] = 1;
                    running[adjList[current][i]] = 1;
                    path.first.push_back(make_pair(adjList[current][i], -1));
                    current = adjList[current][i];
                    tmp1 = true;
                    break;
                }
                if (running[adjList[current][i]]) {
                    // Update running[all kmers in path] = 0.
                    for (int j = 0; j < path.first.size(); j++) {
                        running[path.first[j].first] = 0;
                    }
                    // Update visited[outside cycle] = 0, and erase them.
                    auto it = find(path.first.begin(), path.first.end(), make_pair(adjList[current][i], -1));
                    auto delete_before = distance(path.first.begin(), it);
                    for (int j = 0; j < delete_before; j++) {
                        visited[path.first[0].first] = 0;
                        path.first.pop_front();
                    }
                    // Mark the path as a "cycle"
                    path.second = 1;
                    // Skip the reverse search because we already found a cycle.
                    tmp = true;
                    break;
                }
            }
            if(!tmp1) break;
        }
        // Jump <current> back to <start>
        current = start;

        // Backward tracking
        while (true) {
            if(tmp) break;
            bool tmp2 = false;
            for (int i = 0; i < revAdjList[current].size(); i++) {
                if (!visited[revAdjList[current][i]]) {
                    visited[revAdjList[current][i]] = 1;
                    running[revAdjList[current][i]] = 1;
                    path.first.push_front(make_pair(revAdjList[current][i], -1));
                    current = revAdjList[current][i];
                    tmp2 = true;
                    break;
                }
                if (running[revAdjList[current][i]]) {
                    // Update running[all kmers in <path>] = 0.
                    for (int j = 0; j < path.first.size(); j++) {
                        running[path.first[j].first] = 0;
                    }
                    // Update visited[outside cycle] = 0, and erase them.
                    auto it = find(path.first.begin(), path.first.end(), make_pair(revAdjList[current][i], -1));
                    auto delete_from = distance(path.first.begin(), it); 
                    auto path_size = path.first.size();
                    for (int j = path_size - 1; j > delete_from; j--) {
                        visited[path.first[path.first.size() - 1].first] = 0;
                        path.first.pop_back();
                    }
                    // Mark the path as a "cycle"
                    path.second = 1;
                    break;
                }
            }
            if(!tmp2) break;
        }
        // Update running[all kmers in <path>] = 0.
        for (int i = 0; i < path.first.size(); i++) {
            running[path.first[i].first] = 0;
        }
        return path;
    }

    Paths findPaths(const NcAdjList& adjList, const NcAdjList& revAdjList) {
        vector<uint8_t> visited(adjList.size(), 0);
        vector<uint8_t> running(adjList.size(), 0);
        Paths paths; // pair<vector<pair<deque<pair<int, int>>, uint8_t>>, int>
        paths.first = {};
        paths.second = 0;

        // Start greedy search from an unvisited node
        for (int i = 0; i < adjList.size(); i++) {
            if (!visited[i]) {
                Path path = greedyPath(i, adjList, revAdjList, visited, running);
                paths.first.push_back(path);
                if (path.second) paths.second++;
            }
        }

        // Add isolated kmers to paths as "Length-one path"s.
        for (int i = 0; i < adjList.size(); i++) {
            if (!visited[i]) {
                paths.first.push_back(make_pair(deque<pair<int, int>>{make_pair(i, -1)}, 0));
            }
        }

        // Delete empty paths.
        paths.first.erase(
            remove_if(paths.first.begin(), paths.first.end(), 
                [](const Path& path) {
                    return path.first.empty();
                }
            ), 
            paths.first.end()
        );

        return paths;
    }

    Paths attachPendants(Paths& paths, const vector<string>& kmers) {
    // List of already attached length-1 paths
    vector<int> attached; 
    
    // Attach length-1 paths to an adjacent path
    for (int i = 0; i < paths.first.size(); i++) {
        if (paths.first[i].first.size() != 1) continue; // Skip non-length-1 paths
        bool attachedToMainPath = false;
        
        for (int j = 0; j < paths.first.size(); j++) {
            if (paths.first[j].first.size() <= 1) continue; // Skip short paths
            
            for (int k = 0; k < paths.first[j].first.size(); k++) {
                // Get k-mer strings
                string pendant = kmers[paths.first[i].first[0].first];
                string attachee = kmers[paths.first[j].first[k].first];
                
                // Check if pendant can be attached to attachee
                if (pendant.substr(0, attachee.size() - 1) == attachee.substr(0, attachee.size() - 1)) {
                    // Attach the pendant
                    paths.first[j].first[k].second = paths.first[i].first[0].first;
                    attached.push_back(paths.first[i].first[0].first);
                    attachedToMainPath = true;
                    break;
                }
            }
            if (attachedToMainPath) break;
        }
    }

    // Remove length-1 paths that were attached
    paths.first.erase(
        remove_if(paths.first.begin(), paths.first.end(), [&attached](const Path& path) {
            // Check if the path is a length-1 path and if it was attached
            return (path.first.size() == 1 && find(attached.begin(), attached.end(), path.first[0].first) != attached.end());
        }),
        paths.first.end()
    );
    
    return paths;
}

};


class EdgeCentric_DeBruijnGraph {
public:
    // Edges of the graph (current node -> successive nodes)
    EcAdjList ec_outAdjacencyList;
    // Edges of the graph (Current node -> previous nodes)
    EcAdjList ec_inAdjacencyList;

    // Construct edge-centric De Bruijn Graph from the given k-mer sets
    void buildGraph(const vector<string>& kmers) {
        for (const string& kmer : kmers) {
            //retrieve k-1 prefix and suffix
            string prefix = kmer.substr(0, kmer.size() - 1);
            string suffix = kmer.substr(1, kmer.size() - 1);

            // Append the edges (current node -> successive nodes)
            ec_outAdjacencyList[prefix].push_back(suffix);
            // Append the edges (current node -> previous nodes)
            ec_inAdjacencyList[suffix].push_back(prefix);
        }
    }

    // Count the number of edges to make ecdbg eulerian.
    int countEdgesToAdd() {
        int edgesToAdd = 0;
        for (const auto& pair : ec_outAdjacencyList) {
            int outMinusIn = ec_outAdjacencyList[pair.first].size() - ec_inAdjacencyList[pair.first].size();
            if (outMinusIn <= 0) continue;
            edgesToAdd += outMinusIn;
        }

        return edgesToAdd;
    }

    // Show the graph
    void printGraph() const {
        for (const auto& pair : ec_outAdjacencyList) {
            const string& node = pair.first;
            const vector<string>& neighbors = pair.second;
            cout << "Node: " << node << " has edges to: ";
            for (const string& neighbor : neighbors) {
                cout << neighbor << " ";
            }
            cout << endl;
        }
    }
};

class NodeCentric_DeBruijnGraph {
public:
    // Edges of the graph (current node, successive nodes)
    NcAdjList nc_outAdjacencyList;
    // Edges of the graph (current node, previous nodes)
    NcAdjList nc_inAdjacencyList;

    // Construct node-centric De Bruijn Graph from the given k-mer sets
    void buildGraph(const vector<string>& kmers) {
        nc_outAdjacencyList.resize(kmers.size());
        // Map current node to next nodes
        for (int i = 0; i < kmers.size(); i++) {
            // Define node as a current looking k-mer, then determine the edges using the last letter
            string kmerSuffix = kmers[i].substr(1); // k-mer's last (k-1) letters
            
            // Search for the next nodes
            for (int j = 0; j < kmers.size(); j++) {
                if (kmers[i] == kmers[j]) continue; // Rule out the edges such as "aaa -> aaa"
                string nextKmerPrefix = kmers[j].substr(0, kmers[j].size() - 1); // next node's first (k-1) letters

                if (kmerSuffix == nextKmerPrefix) {
                    nc_outAdjacencyList[i].push_back(j);
                }
            }
        }
    }

    // Construct the reverse adjacency list of the above
    void buildReverseGraph(const vector<string>& kmers) {
        nc_inAdjacencyList.resize(kmers.size());
        // Map current node to previous nodes
        for (int i = 0; i < kmers.size(); i++) {
            // Define node as a current looking k-mer, then determine the edges using the first letter
            string kmerPrefix = kmers[i].substr(0, kmers[i].size() - 1); // k-mer's first (k-1) letters
            
            // Search for the next nodes
            for (int j = 0; j < kmers.size(); j++) {
                if (kmers[i] == kmers[j]) continue; // Rule out the edges such as "aaa <- aaa"
                string previousKmerSuffix = kmers[j].substr(1); // previous node's last (k-1) letters

                if (kmerPrefix == previousKmerSuffix) {
                    nc_inAdjacencyList[i].push_back(j);
                }
            }
        }
    }

    // Show the graph
    void printGraph(const vector<string>& kmers) {
        for (int i = 0; i < kmers.size(); i++) {
            cout << "Node: " << kmers[i] << " has edges to: ";
            for (int j = 0; j < nc_outAdjacencyList[i].size(); j++) {
                cout << kmers[nc_outAdjacencyList[i][j]] << " ";
            }
            cout << endl;
        }
    }

    Path greedyPath(const int& start, const NcAdjList& adjList, const NcAdjList& revAdjList, vector<uint8_t>& visited, vector<uint8_t>& running){
        Path path; // pair<deque<pair<int, int>>, uint8_t>
        path.first = {};
        path.second = 0;

        int current = start;
        path.first.push_back(make_pair(current, -1));
        visited[current] = 1;
        running[current] = 1;
        bool tmp = 0;

        // Forward tracking
        while (true) {
            bool tmp1 = false;
            for (int i = 0; i < adjList[current].size(); i++) {
                if (!visited[adjList[current][i]]) {
                    visited[adjList[current][i]] = 1;
                    running[adjList[current][i]] = 1;
                    path.first.push_back(make_pair(adjList[current][i], -1));
                    current = adjList[current][i];
                    tmp1 = true;
                    break;
                }
                if (running[adjList[current][i]]) {
                    // Update running[all kmers in path] = 0.
                    for (int j = 0; j < path.first.size(); j++) {
                        running[path.first[j].first] = 0;
                    }
                    // Update visited[outside cycle] = 0, and erase them.
                    auto it = find(path.first.begin(), path.first.end(), make_pair(adjList[current][i], -1));
                    auto delete_before = distance(path.first.begin(), it);
                    for (int j = 0; j < delete_before; j++) {
                        visited[path.first[0].first] = 0;
                        path.first.pop_front();
                    }
                    // Mark the path as a "cycle"
                    path.second = 1;
                    // Skip the reverse search because we already found a cycle.
                    tmp = true;
                    break;
                }
            }
            if(!tmp1) break;
        }
        // Jump <current> back to <start>
        current = start;

        // Backward tracking
        while (true) {
            if(tmp) break;
            bool tmp2 = false;
            for (int i = 0; i < revAdjList[current].size(); i++) {
                if (!visited[revAdjList[current][i]]) {
                    visited[revAdjList[current][i]] = 1;
                    running[revAdjList[current][i]] = 1;
                    path.first.push_front(make_pair(revAdjList[current][i], -1));
                    current = revAdjList[current][i];
                    tmp2 = true;
                    break;
                }
                if (running[revAdjList[current][i]]) {
                    // Update running[all kmers in <path>] = 0.
                    for (int j = 0; j < path.first.size(); j++) {
                        running[path.first[j].first] = 0;
                    }
                    // Update visited[outside cycle] = 0, and erase them.
                    auto it = find(path.first.begin(), path.first.end(), make_pair(revAdjList[current][i], -1));
                    auto delete_from = distance(path.first.begin(), it); 
                    auto path_size = path.first.size();
                    for (int j = path_size - 1; j > delete_from; j--) {
                        visited[path.first[path.first.size() - 1].first] = 0;
                        path.first.pop_back();
                    }
                    // Mark the path as a "cycle"
                    path.second = 1;
                    break;
                }
            }
            if(!tmp2) break;
        }
        // Update running[all kmers in <path>] = 0.
        for (int i = 0; i < path.first.size(); i++) {
            running[path.first[i].first] = 0;
        }
        return path;
    }

    Paths findPaths(const NcAdjList& adjList, const NcAdjList& revAdjList) {
        vector<uint8_t> visited(adjList.size(), 0);
        vector<uint8_t> running(adjList.size(), 0);
        Paths paths; // pair<vector<pair<deque<pair<int, int>>, uint8_t>>, int>
        paths.first = {};
        paths.second = 0;

        // Start greedy search from an unvisited node
        for (int i = 0; i < adjList.size(); i++) {
            if (!visited[i]) {
                Path path = greedyPath(i, adjList, revAdjList, visited, running);
                paths.first.push_back(path);
                if (path.second) paths.second++;
            }
        }

        // Add isolated kmers to paths as "Length-one path"s.
        for (int i = 0; i < adjList.size(); i++) {
            if (!visited[i]) {
                paths.first.push_back(make_pair(deque<pair<int, int>>{make_pair(i, -1)}, 0));
            }
        }

        // Delete empty paths.
        paths.first.erase(
            remove_if(paths.first.begin(), paths.first.end(), 
                [](const Path& path) {
                    return path.first.empty();
                }
            ), 
            paths.first.end()
        );

        return paths;
    }

    Paths attachPendants(Paths& paths, const vector<string>& kmers) {
    // List of already attached length-1 paths
    vector<int> attached; 
    
    // Attach length-1 paths to an adjacent path
    for (int i = 0; i < paths.first.size(); i++) {
        if (paths.first[i].first.size() != 1) continue; // Skip non-length-1 paths
        bool attachedToMainPath = false;
        
        for (int j = 0; j < paths.first.size(); j++) {
            if (paths.first[j].first.size() <= 1) continue; // Skip short paths
            
            for (int k = 0; k < paths.first[j].first.size(); k++) {
                // Get k-mer strings
                string pendant = kmers[paths.first[i].first[0].first];
                string attachee = kmers[paths.first[j].first[k].first];
                
                // Check if pendant can be attached to attachee
                if (pendant.substr(0, attachee.size() - 1) == attachee.substr(0, attachee.size() - 1)) {
                    // Attach the pendant
                    paths.first[j].first[k].second = paths.first[i].first[0].first;
                    attached.push_back(paths.first[i].first[0].first);
                    attachedToMainPath = true;
                    break;
                }
            }
            if (attachedToMainPath) break;
        }
    }

    // Remove length-1 paths that were attached
    paths.first.erase(
        remove_if(paths.first.begin(), paths.first.end(), [&attached](const Path& path) {
            // Check if the path is a length-1 path and if it was attached
            return (path.first.size() == 1 && find(attached.begin(), attached.end(), path.first[0].first) != attached.end());
        }),
        paths.first.end()
    );
    
    return paths;
}

};



void printResult(const Paths necklaces, const vector<string>& kmers, const int countEdgesToAdd) {
    // Output the resulting paths
    int kmers_covered = 0;
    for (int i = 0; i < necklaces.first.size(); i++) {
        cout << "Path" << i + 1 << " : ";
        for (int j = 0; j < necklaces.first[i].first.size(); j++) {
            if (necklaces.first[i].first[j].second == -1) {
                cout << kmers[necklaces.first[i].first[j].first] << " ";
                kmers_covered += 1;
            }
            else {
                cout << "(" << kmers[necklaces.first[i].first[j].first] << ", "
                << kmers[necklaces.first[i].first[j].second] << ") ";
                kmers_covered += 2;
            }
        }
        if (necklaces.first[i].second) cout << "*";
        cout << endl;
    }
    cout << "Found " << necklaces.first.size() << " necklaces in total! (Letter \"*\" at the end means it is a cycle)" << endl;

    if (kmers_covered == kmers.size()) {
        cout << "Successfully covered all kmers!" << endl;
    }
    else {
        cout << "There are " << kmers.size() - kmers_covered << " nodes uncovered..." << endl;
    }

    cout << endl;
    cout << string(60, '-') << endl;
    cout << "# Edges to add to make Eulerian: " << countEdgesToAdd << endl;
    cout << "# Open necklaces in the obtained Necklace cover: " << necklaces.first.size() - necklaces.second << endl;
    cout << string(60, '-') << endl;
}

int main() {
    // Construct a vector of k-mers from an input FASTA file(.fa).
    vector<string> kmers = getKmers();
    int kmers_size = kmers.size();

    // Construct and show the edge-centric De Bruijn Graph
    EdgeCentric_DeBruijnGraph ecdbg;
    ecdbg.buildGraph(kmers);
    int countEdgesToAdd = ecdbg.countEdgesToAdd();
    cout << "Edge-centric De Bruijn Graph:" << endl;
    ecdbg.printGraph();
    cout << endl;
    
    // Construct and show the node-centric De Bruijn Graph
    NodeCentric_DeBruijnGraph ncdbg;
    ncdbg.buildGraph(kmers);
    ncdbg.buildReverseGraph(kmers);
    cout << "Node-centric De Bruijn Graph:" << endl;
    ncdbg.printGraph(kmers);
    cout << endl;

    // Find paths from a given adjacency list
    Paths paths = ncdbg.findPaths(ncdbg.nc_outAdjacencyList, ncdbg.nc_inAdjacencyList);
    
    // Attach length-1 path to an adjacent path as a pendant
    Paths necklaces = ncdbg.attachPendants(paths, kmers);

    // Output a necklace cover and a comparison with Eulertigs
    printResult(necklaces, kmers, countEdgesToAdd);

    return 0;
}
