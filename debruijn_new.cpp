// Ren Kimura, 2024
// compile as g++ -std=c++17 debruijn_new.cpp 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <deque>
#include <variant>
#include <cctype>

using namespace std;
using AdjList = vector<vector<int64_t>>;
using Kmers = unordered_map<string,int64_t>;
using NodeId = vector<int64_t>;

using Path = deque<int64_t>;
// {(-N, ..., -1, 1, ..., N)^(#nodes in the path including pendants)[,0]}. 
// Pendants are represented as a negative of the corresponding positive ID, where 0 in the end means the path is a cycle.
using Paths =  vector<deque<int64_t>>;

class DeBruijnGraph {
public:
    int K;
    bool isNodeCentric; // node-centric or edge-centing dBG?
    vector<uint8_t> Alphabet = {'A', 'C', 'G', 'T'};

    string sequence;  // To store the concatenated sequence
    Kmers kmers;  // pair <kmer string of length k, its unique ID in 1..N>, where N = #distinct kmers in sequence
    NodeId idPosition; // for each ID a position i such that ID corresponds to kmer sequence.substr(i, K)
    unordered_set<int64_t> pdCands;

    int64_t makeEulerian = 0; // Number of edges to make Eulerian
    int64_t countOpenNecklaces = 0; // Number of open necklaces;

    DeBruijnGraph(int _K, bool _isNodeCentric) : K(_K), isNodeCentric(_isNodeCentric)
    {};

    __inline void to_uppercase(string& str) {
        transform(str.begin(), str.end(), str.begin(),
                    [](unsigned char c) { return toupper(c); });
    }

    void getKmers(string filename) {
        // Read input DNA sequence in FASTA format .fa
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

        // Change lower into capital case
        to_uppercase(sequence);

        // Erase letters other than "A, T, G, C"
        sequence.erase(remove_if(sequence.begin(), sequence.end(), [](char x) 
            {return (x != 'A' && x != 'T' && x != 'G' && x != 'C');}
            ), sequence.end());

        // Check every k-mer in "sequence"
        for (int i = 0; i < sequence.length() - K; i++){
            auto kmer = sequence.substr(i, K);
            if (kmers.find(kmer) == kmers.end()){ // newly found kmer
                kmers[kmer] = idPosition.size() + 1; // "+ 1" because we want kmers ID starting from 1 to N
                idPosition.push_back(i);
            }
        }

        for (auto const &x : kmers) {
            cout << x.second << ": " << x.first << " at position " << idPosition[x.second - 1] << endl;
        }
    }

    int64_t forward(int64_t id, uint8_t c){
        auto position = idPosition[id - 1]; // "id - 1" because we have kmers ID starting from 1
        string suffix = sequence.substr(position + 1, K - 1); 
        auto next = suffix + char(c);
        if (kmers.find(next) != kmers.end()) // kmer exists
            return kmers[next];
        else 
            return 0; // no outgoing branching exists for c
    }

    int64_t backward(int64_t id, uint8_t c){
        auto position = idPosition[id - 1]; // "id - 1" because we have kmers ID starting from 1
        string prefix = sequence.substr(position, K - 1);
        auto previous = char(c) + prefix;
        if (kmers.find(previous) != kmers.end()) // kmer exists
            return kmers[previous];
        else
            return 0; // no incoming branching exists for c
    }

    // Edges of the graph (current node, outgoing nodes)
    AdjList outNeighbors;
    // Edges of the graph (current node, incoming nodes)
    AdjList inNeighbors;
    
    // Construct node-centric De Bruijn Graph from the given UPPERCASE k-mer sets
    void buildNodeCentricGraph( ){ 
        outNeighbors.resize(kmers.size());
        inNeighbors.resize(kmers.size());
        // Map current node to next nodes
        for (auto const &x : kmers) {
            auto kmer = x.first;
            auto id = x.second;
            auto position = idPosition[id - 1]; // "id - 1" because we have kmers ID starting from 1

            auto kmerSuffix = sequence.substr(position + 1, K - 1); // k-mer's last (k-1) letters
            for (auto& c : Alphabet){
                auto next = kmerSuffix + char(c);
                if (kmers.find(next) != kmers.end()) { // kmer exists, so its node
                    auto nextid = kmers[next];
                    outNeighbors[id - 1].push_back(nextid - 1);
                    inNeighbors[nextid - 1].push_back(id - 1);
                }       
            }
        }
    }

    void buildEdgeCentricGraph(const string S, const Kmers& kmers, const NodeId& idPos) {  // Rough invariant: kmers.first == S.substr(kmers.second - 1, k)
    }

    // Display the graph
    void printGraph() {
        for (auto const &x : kmers) { 
            cout << x.first << " -> ";  // kmer string
            for (auto const& idNeigh_minus1 : outNeighbors[x.second - 1]) {
                auto position = idPosition[idNeigh_minus1]; // "idNeigh - 1" because we have kmers ID starting from 1
                cout << sequence.substr(position, K) << " ";
            }
            cout << endl;
        }
    }

    Path greedyPath(const int64_t start, bool visited[], bool running[]){  // to be initialized as visited = running = { false }
        Path path;
        auto current = start;
        path.push_back(current);
        visited[current] = true;
        running[current] = true;
        
        // [suggestion] from current, how to choose the best neighbor: 
        // - first choice is to get a cycle
        // - second choice is to extend the path
        // - third choice is CASE 1 when neither of the first two choices above occur

        // forward search
        int64_t i = 0;
        while (i < Alphabet.size())
        {
            auto c = Alphabet[i];
            auto next = forward(current, c);
            if(next == 0) continue; // skip if there isn't an edge from current to next
            if(visited[next]){
                if(!running[next]) break; // break if it is not a cycle (------CASE 1------)
                // when found a cycle (------CASE 2------)
                for (auto kmerId : path) { // update running[every kmer in the path] = false
                    running[kmerId] = false;
                }
                auto itr = find(path.begin(), path.end(), next);
                auto delete_before = distance(path.begin(), itr);
                while (delete_before--) { // update visited[kmers outside the cycle] = false, and erase them
                    visited[path.front()] = false;
                    path.pop_front();
                }
                path.push_back(0); // add zero meaning cycle in the end of the path 
                return path; // return a cycle
            }
            // when next is not visited (------CASE 3------)
            visited[next] = true;
            running[next] = true;
            path.push_back(next);
            current = next;
            i = 0; // need this to restart the scan of the neighbors of the new current
        }
        // jump current back to start
        current = start;
        // backward search
        i = 0;
        while (i < Alphabet.size())
        {
            auto c = Alphabet[i];
            auto previous = backward(current, c);
            if(previous == 0) continue; // skip if there isn't an edge from current to previous
            if(visited[previous]){
                if(!running[previous]) break; // break if it is not a cycle (------CASE 4------)
                // when found a cycle (------CASE 5------)
                for (auto kmer : path) { // update running[every kmer in the path] = false
                    running[kmer] = false;
                }
                auto itr = find(path.begin(), path.end(), previous);
                auto delete_before = distance(path.begin(), itr);
                auto delete_from = path.size() - delete_before; // how many k-mers do we delete
                while (delete_from--) { // update visited[kmers outside the cycle] = false, and erase them
                    visited[path.back()] = false;
                    path.pop_back();
                }
                path.push_back(0); // add zero meaning cycle in the end of the path 
                return path; // return a cycle
            }
            // when previous is not visited (------CASE 6------)
            visited[previous] = true;
            running[previous] = true;
            path.push_front(previous);
            current = previous;
            i = 0; // need this to restart the scan of the neighbors of the new current
        }
        countOpenNecklaces += 1;
        // update running[every kmer in the path] = false
        for (auto kmer : path) {
            running[kmer] = false;
        } 
        return path; // return a non-cyclic path
    }

    Path greedyPathRev(const int64_t start, vector<bool>& visited, vector<bool>& running) {
        Path path;
        auto current = start;
        path.push_back(current);
        visited[current] = true;
        running[current] = true;

        // forward search
        for (auto c : Alphabet) {
            auto next = forward(current, c);
            if (next == 0) continue;
            if (!visited[next]) {
                visited[next] = true;
                running[next] = true;
                path.push_back(next);
                current = next;
                continue;
            }
            if (!running[next]) continue;
            for (auto v : path) {
                running[v] = false;
            }
            auto itr = find(path.begin(), path.end(), next);
            auto delete_before = distance(path.begin(), itr);
            while (delete_before) {
                visited[path[0]] = false;
                path.pop_front();
            }
            path.push_back(0);
            return path;
        }
        current = start;
        // backward search
        for (auto c : Alphabet) {
            auto prev = backward(current, c);
            if (prev == 0) continue;
            if (!visited[prev]) {
                visited[prev] = true;
                running[prev] = true;
                path.push_front(prev);
                current = prev;
                continue;
            }
            if (!running[prev]) continue;
            for (auto v : path) {
                running[v] = false;
            }
            auto itr = find(path.begin(), path.end(), prev);
            auto delete_from = distance(path.begin(), itr);
            auto path_size = path.size();
            while (delete_from < path_size) {
                visited[path[path.size() - 1]] = false;
                path.pop_back();
                delete_from++;
            }
            path.push_back(0);
            return path;
        }
        countOpenNecklaces += 1;
        for (auto v : path) {
            visited[v] = false;
        }
        return path;
    }

    Path greedyPath2(const int64_t start, vector<bool>& visited, vector<bool>& running) {
        Path path;
        auto current = start;
        path.push_back(current);
        visited[current] = true;
        running[current] = true;
        uint8_t c_extend = 0;

        for (int8_t i = 0; i < 4; i++) {
            uint8_t c = Alphabet[i];
            int64_t next = forward(current, c);
            if (next && visited[next] && running[next]) {
                for (auto id : path) {
                    running[id] = false;
                }
                auto itr = find(path.begin(), path.end(), next);
                auto delete_before = distance(path.begin(), itr);
                while (delete_before--) {
                    visited[path[0]] = false;
                    path.pop_front();
                }
                path.push_back(0);
                return path;   
            }             
            else if (next && !visited[next] && !running[next]) { // when extendable
                if(i < 3) {c_extend = c; continue;} // keep c & pend an extension
                visited[next] = running[next] = true; // for the last alphabet, extend
                path.push_back(next);
                current = next;
                c_extend = 0;
            }     
            else if (!next || (next && visited[next] && !running[next])) { // "no edge to c" or "occupied by another path"
                if (i < 3) continue; // go on the next char c
                if (!c_extend) break; // proceed to backward search
                next = forward(current, c_extend);
                visited[next] = running[next] = true;
                path.push_back(next);
                current = next;
                c_extend = 0;
            }       
        }
        current = start;
        for (int8_t i = 0; i < 4; i++) {
            uint8_t c = Alphabet[i];
            int64_t prev = backward(current, c);
            if (prev && visited[prev] && running[prev]){
                for (auto id : path) {
                    running[id] = false;
                }
                auto itr = find(path.begin(), path.end(), prev);
                auto delete_from = distance(path.begin(), itr);
                auto dlt = path.size() - delete_from - 1;
                while (dlt--) {
                    visited[path[path.size() - 1]] = false;
                    path.pop_back();
                }
                path.push_back(0);
                return path;   
            }
            else if (prev && !visited[prev] && !running[prev]) { // when extendable
                if(i < 3) {c_extend = c; continue;} // keep c & pend an extension
                visited[prev] = running[prev] = true; // for the last alphabet, extend
                path.push_front(prev);
                current = prev;
                c_extend = 0;
            }
            else if (!prev || (prev && visited[prev] && !running[prev])) { // "no edge to c" or "occupied by another path"
                if (i < 3) continue; // go on the next char c
                if (!c_extend) break; // proceed to the last part
                prev = backward(current, c_extend);
                visited[prev] = running[prev] = true;
                path.push_front(prev);
                current = prev;
                c_extend = 0;
            }
        }
        if (path.size() > 2) countOpenNecklaces += 1;
        else pdCands.insert(current); // here it is a len-1 path which means to be a pdCand
        // update running[every kmer in the path] = false
        for (auto id : path) {
            running[id] = false;
        }
        return path;
    }

    Paths findPaths() {
        Paths paths;
        vector<bool> visited(kmers.size() + 1);
        vector<bool> running(kmers.size() + 1);

        // Start greedy search from an unvisited node
        for (auto const x : kmers) {
            auto id = x.second;
            if (!visited[id]) {
                Path path = greedyPath2(id, visited, running);
                if (path.size() > 1) // not adding pdCand to paths
                    paths.push_back(path);
            }
        }

        // [suggestion] maybe you can decide in the loop below where to attach isolated kmers (instead of having attachPendants() )

        // attach pendants
        for (auto& path : paths) {
            for (auto& id : path) {
                if (id <= 0) continue;
                for (auto const &c : Alphabet) {
                    auto next = forward(id, c);
                    if (next == 0) continue;
                    if (pdCands.find(next) != pdCands.end()) {
                        auto it = find(path.begin(), path.end(), id);
                        if (it != path.end()) {
                            path.insert(it + 1, -next);
                            pdCands.erase(next);
                        }
                    }
                }
            }
        }
        return paths;
    }

    void printResult(const Paths& paths) {
        for (auto const& path : paths) {
            for (auto const& id : path) {
                cout << id << " ";
            }
            cout << endl;
        }

        int64_t kmers_covered = 0;
        int64_t count = 1;
        for (auto const& path : paths) {
            cout << "Path" << count << ": ";
            for (auto const& id : path) {
                if (id > 0) {
                    cout << sequence.substr(idPosition[id - 1], K) << " ";
                    kmers_covered += 1;
                }
                else if (id < 0) {
                    cout << "+" << sequence.substr(idPosition[-id - 1], K) << " ";
                    kmers_covered += 1;
                }
                else { // id == 0
                    cout << "*";
                }
            }
            cout << endl;
            count++;
        }
        cout << "Found " << paths.size() 
             << " necklaces in total!\n" 
             << "\"*\" at the end means it is a cycle.\n"
             << "\"+(k-mer)\" are pendants." << endl;

        if (kmers_covered == kmers.size()) {
            cout << "Successfully covered all kmers!" << endl;
        }
        else {
            cout << "There are " << kmers.size() - kmers_covered << " nodes uncovered..." << endl;
        }

        cout << endl;
        cout << string(60, '-') << endl;
        cout << "# Edges to add to make Eulerian: " << makeEulerian << endl;
        cout << "# Open necklaces: " << countOpenNecklaces << endl;
        cout << string(60, '-') << endl;
    }
};

int main() {
    // Construct a vector of k-mers from an input FASTA file(.fa).
    int64_t _K;
    string filename;
    cout << "Input file: ";
    cin >> filename;
    cout << "Set K to : ";
    cin >> _K;
    
    // Construct and show the node-centric De Bruijn Graph
    DeBruijnGraph ncdbg = DeBruijnGraph(_K, true);
    ncdbg.getKmers(filename);
    ncdbg.buildNodeCentricGraph();
    cout << "Node-centric De Bruijn Graph:" << endl;
    ncdbg.printGraph();
    cout << endl;

    // Find paths from a given adjacency list
    Paths paths = ncdbg.findPaths();

    // Output a necklace cover and a comparison with Eulertigs
    ncdbg.printResult(paths);

    return 0;
}