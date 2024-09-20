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

using Path = pair<deque<pair<int64_t, int64_t>>, bool>;
// Path = deque<int64_t> {(-N, ..., -1, 1, ..., N)^(#nodes in the path including pendants)[,0]}. Pendants are represented as a negative of the corresponding positive <id>. 0 in the end means the path is a cycle.
using Paths = pair<vector<pair<deque<pair<int64_t, int64_t>>, bool>>, int64_t>;
// Paths = 


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
    void buildNodeCentricGraph( ){ 
        outNeighbors.resize(kmers.size());
        inNeighbors.resize(kmers.size());
        // Map current node to next nodes
        for (auto const &x : kmers) {
            auto kmer = x.first;
            auto id = x.second;
            auto position = idPosition[id];

            // Define node as a current looking k-mer, then determine the edges using the last letter
            auto kmerSuffix = sequence.substr(position + 1, K-1); // k-mer's last (k-1) letters
            for (auto& c : Alphabet){
                auto next = kmerSuffix.append(to_string(c));
                if (kmers.find(next) != kmers.end()) {
                    auto nextid = kmers[next];
                    outNeighbors[id].push_back(nextid);
                    inNeighbors[nextid].push_back(id);
                }       
            }
        }
    }

    int64_t forward(uint64_t id, uint8_t c){
        auto position = idPosition[id];
        string kmerSuffix = sequence.substr(position + 1, K-1); 
        auto next = kmerSuffix.append(to_string(c));
        if (kmers.find(next) != kmers.end())
            return kmers[next];
        else 
            return -1; // no outgoing branching exists for c
    }

    int64_t backward(uint64_t id, uint8_t c){
        auto position = idPosition[id];
        string kmerPrefix = sequence.substr(position, K-1);
        auto previous = to_string(c) + kmerPrefix;
        if (kmers.find(previous) != kmers.end())
            return kmers[previous];
        else
            return -1; // no incoming branching exists for c
    }

    void buildEdgeCentricGraph(const string S, const Kmers& kmers, const NodeId& idPos) {  // Rough invariant: kmers.first == S.substr(kmers.second,k)
    }

    // Show the graph
    void printGraph() {
        for (auto const &x : kmers) { 
            cout << "Node: " << x.first << " -> ";  // kmer string
            for (auto const & idNeigh : outNeighbors[x.second]) {
                auto position = idPosition[idNeigh];
                cout << sequence.substr(position, K) << " ";
            }
            cout << endl;
        }
    }

    Path greedyPath(const int64_t start, bool visited[], bool running[]){
        Path path; // Path = pair<deque<pair<int64_t, int64_t>>, bool>;
        path.first = {};
        path.second = false;

        int current = start;
        path.first.push_back(make_pair(current, -1));
        visited[current] = true;
        running[current] = true;
        bool tmp = 0;

        // Forward tracking
        while (true) {
            bool tmp1 = false;
            for (int i = 0; i < outNeighbors[current].size(); i++) {
                if (!visited[outNeighbors[current][i]]) {
                    visited[outNeighbors[current][i]] = true;
                    running[outNeighbors[current][i]] = true;
                    path.first.push_back(make_pair(outNeighbors[current][i], -1));
                    current = outNeighbors[current][i];
                    tmp1 = true;
                    break;
                }
                if (running[outNeighbors[current][i]]) {
                    // Update running[all kmers in path] = 0.
                    for (int j = 0; j < path.first.size(); j++) {
                        running[path.first[j].first] = false;
                    }
                    // Update visited[outside cycle] = 0, and erase them.
                    auto it = find(path.first.begin(), path.first.end(), make_pair(outNeighbors[current][i], -1));
                    auto delete_before = distance(path.first.begin(), it);
                    for (int j = 0; j < delete_before; j++) {
                        visited[path.first[0].first] = false;
                        path.first.pop_front();
                    }
                    // Mark the path as a "cycle"
                    path.second = true;
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
            for (int i = 0; i < inNeighbors[current].size(); i++) {
                if (!visited[inNeighbors[current][i]]) {
                    visited[inNeighbors[current][i]] = true;
                    running[inNeighbors[current][i]] = true;
                    path.first.push_front(make_pair(inNeighbors[current][i], -1));
                    current = inNeighbors[current][i];
                    tmp2 = true;
                    break;
                }
                if (running[inNeighbors[current][i]]) {
                    // Update running[all kmers in <path>] = 0.
                    for (int j = 0; j < path.first.size(); j++) {
                        running[path.first[j].first] = false;
                    }
                    // Update visited[outside cycle] = 0, and erase them.
                    auto it = find(path.first.begin(), path.first.end(), make_pair(inNeighbors[current][i], -1));
                    auto delete_from = distance(path.first.begin(), it); 
                    auto path_size = path.first.size();
                    for (int j = path_size - 1; j > delete_from; j--) {
                        visited[path.first[path.first.size() - 1].first] = false;
                        path.first.pop_back();
                    }
                    // Mark the path as a "cycle"
                    path.second = true;
                    break;
                }
            }
            if(!tmp2) break;
        }
        // Update running[all kmers in <path>] = 0.
        for (int i = 0; i < path.first.size(); i++) {
            running[path.first[i].first] = false;
        }
        return path;
    }

    deque<int64_t> greedyPath_rev(const int64_t start, bool visited[], bool running[]){
        deque<int64_t> path;
        auto current = start;
        path.push_back(current);
        visited[current] = true;
        running[current] = true;
        
        // Forward search
        for (auto c : Alphabet) { 
            auto next = forward(current, c);
            if(next == -1) continue; // skip if there isn't an edge from current to next
            if(visited[next]){
                if(!running[next]) break; // break if it is not a cycle (------CASE 1------)
                // when found a cycle (------CASE 2------)
                for (auto kmer : path) { // update running[every kmer in the path] = false
                    running[kmer] = false;
                }
                auto itr = path.finds(next);
                auto delete_before = distance(path.begin(), itr);
                while (delete_before) { // update visited[kmers outside the cycle] = false, and erase them
                    visited[0] = false;
                    path.pop_front();
                    delete_before--;
                }
                path.push_back(0); // add zero meaning cycle in the end of the path 
                return path;
            }
            // when next is not visited (------CASE 3------)
            visited[next] = true;
            running[next] = true;
            path.push_back(next);
            current = next;
        }
        // Backward search
        for (auto c : Alphabet) { 
            auto previous = backward(current, c);
            if(previous == -1) continue; // skip if there isn't an edge from current to previous
            if(visited[previous]){
                if(!running[previous]) break; // break if it is not a cycle (------CASE 4------)
                // when found a cycle (------CASE 5------)
                for (auto kmer : path) { // update running[every kmer in the path] = false
                    running[kmer] = false;
                }
                auto itr = path.finds(previous);
                auto delete_from = distance(path.begin(), itr);
                auto dlt = path.size() - delete_from // how many k-mers do we delete
                while (dlt) { // update visited[kmers outside the cycle] = false, and erase them
                    visited[path.size()] = false;
                    path.pop_back();
                    dlt--;
                }
                path.push_front(0); // add zero meaning cycle in the end of the path 
                return path;
            }
            // when previous is not visited (------CASE 6------)
            visited[previous] = true;
            running[previous] = true;
            path.push_front(previous);
            current = previous;
        }
    }

    Paths findPaths() {
        Paths paths; // Paths = pair<vector<pair<deque<pair<int64_t, int64_t>>, bool>>, int64_t>;
        paths.first = {};
        paths.second = 0;
        bool visited[] = {false};
        bool running[] = {false};

        // Start greedy search from an unvisited node
        for (int i = 0; i < kmers.size(); i++) {
            if (!visited[i]) {
                Path path = greedyPath(i, visited, running);
                paths.first.push_back(path);
                if (path.second) paths.second++;
            }
        }

        // Add isolated kmers to paths as "Length-one path"s.
        for (auto i = 0; i < kmers.size(); i++) {
            if (!visited[i])
                paths.first.push_back(make_pair(deque<pair<int64_t, int64_t>>{make_pair(i, -1)}, false));
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

    Paths attachPendants(Paths& paths) {
        // List of already attached length-1 paths
        vector<int64_t> attached; 
        
        // Attach length-1 paths to an adjacent path
        for (int i = 0; i < paths.first.size(); i++) {
            if (paths.first[i].first.size() != 1) continue; // Skip non-length-1 paths
            bool attachedToMainPath = false;
            
            for (int j = 0; j < paths.first.size(); j++) {
                if (paths.first[j].first.size() <= 1) continue; // Skip short paths
                
                for (int k = 0; k < paths.first[j].first.size(); k++) {
                    // Get k-mer strings
                    string pendant = sequence.substr(idPosition[paths.first[i].first[0].first], K);
                    string attachee = sequence.substr(idPosition[paths.first[j].first[k].first], K);
                    
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

    void printResult(const Paths& necklaces, int64_t countEdgesToAdd) {
        // Output the resulting paths
        int kmers_covered = 0;
        for (int i = 0; i < necklaces.first.size(); i++) {
            cout << "Path" << i + 1 << " : ";
            for (int j = 0; j < necklaces.first[i].first.size(); j++) {
                if (necklaces.first[i].first[j].second == -1) {
                    cout << sequence.substr(idPosition[necklaces.first[i].first[j].first], K) << " ";
                    kmers_covered += 1;
                }
                else {
                    cout << "(" << sequence.substr(idPosition[necklaces.first[i].first[j].first], K) << ", "
                    << sequence.substr(idPosition[necklaces.first[i].first[j].second], K) << ") ";
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
};

int main() {
    // Construct a vector of k-mers from an input FASTA file(.fa).
    int64_t _K;
    string filename;
    cout << "Set K to : ";
    cin >> _K;
    cout << "Input file: ";
    cin >> filename;
    
    // Construct and show the node-centric De Bruijn Graph
    DeBruijnGraph ncdbg(_K, true);
    ncdbg.getKmers(filename);
    ncdbg.buildNodeCentricGraph();
    cout << "Node-centric De Bruijn Graph:" << endl;
    ncdbg.printGraph();
    cout << endl;

    // Find paths from a given adjacency list
    Paths paths = ncdbg.findPaths();
    
    // Attach length-1 path to an adjacent path as a pendant
    Paths necklaces = ncdbg.attachPendants(paths);

    // Output a necklace cover and a comparison with Eulertigs
    ncdbg.printResult(necklaces, -1);

    return 0;
}