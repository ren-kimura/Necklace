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
    string filename;
    int64_t K;
    bool isNodeCentric; // node-centric or edge-centing dBG?
    vector<uint8_t> Alphabet = {'A', 'C', 'G', 'T'};

    string sequence;  // To store the concatenated sequence
    Kmers kmers;  // pair <kmer string of length k, its unique ID in 1..N>, where N = #distinct kmers in sequence
    NodeId idPosition; // for each ID a position i such that ID corresponds to kmer sequence.substr(i, K)
    unordered_set<int64_t> pdCands;

    int64_t makeEulerian = 0; // # edges to make ecdBG eulerian
    int64_t countOpenNecklaces = 0; // # open necklaces

    DeBruijnGraph(string _filename, int64_t _K, bool _isNodeCentric) : filename(_filename), K(_K), isNodeCentric(_isNodeCentric)
    {};

    __inline void to_uppercase(string& str) {
        transform(str.begin(), str.end(), str.begin(),
                    [](unsigned char c) { return toupper(c); });
    }

    void getKmers() {
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
        for (int64_t i = 0; i < (int64_t)(sequence.length() - K + 1); i++){
            auto kmer = sequence.substr(i, K);
            if (kmers.find(kmer) == kmers.end()){ // newly found kmer
                idPosition.push_back(i);
                kmers[kmer] = idPosition.size();
            }
        }
    }

    int64_t forward(int64_t id, uint8_t c){
        auto position = idPosition[id - 1]; // "id - 1" because we have kmers ID starting from 1
        string suffix = sequence.substr(position + 1, K - 1); 
        auto next = suffix + static_cast<char>(c);
        if (kmers.find(next) != kmers.end()) // kmer exists
            return kmers[next];
        else 
            return 0; // no outgoing branching exists for c
    }

    int64_t backward(int64_t id, uint8_t c){
        auto position = idPosition[id - 1]; // "id - 1" because we have kmers ID starting from 1
        string prefix = sequence.substr(position, K - 1);
        auto previous = static_cast<char>(c) + prefix;
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
    void buildGraph(){ 
        outNeighbors.resize(kmers.size());
        inNeighbors.resize(kmers.size());
        // Map current node to next nodes
        for (auto const &x : kmers) {
            auto kmer = x.first;
            auto id = x.second;
            auto position = idPosition[id - 1]; // "id - 1" because we have kmers ID starting from 1

            auto kmerSuffix = sequence.substr(position + 1, K - 1); // k-mer's last (k-1) letters
            for (auto& c : Alphabet){
                auto next = kmerSuffix + static_cast<char>(c);
                // rule out edges of k-mers not occuring in sequence (Edge centric)
                if (!isNodeCentric && sequence.find(kmer + static_cast<char>(c)) == string::npos) continue;

                if (kmers.find(next) != kmers.end()) { // kmer exists, so its node
                    auto nextid = kmers[next];
                    outNeighbors[id - 1].push_back(nextid - 1);
                    inNeighbors[nextid - 1].push_back(id - 1);
                }       
            }
        }
        
        // count #unbalanced edges(indegree != outdegree) in edge centric dBG
        if (!isNodeCentric) {
            for (auto const &x: kmers) {
                auto id = x.second;
                int64_t outMinusIn = outNeighbors[id - 1].size() - inNeighbors[id - 1].size();
                if (outMinusIn > 0) makeEulerian += outMinusIn;
            }
        }
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

    Path greedyPath(const int64_t start, vector<bool>& visited, vector<bool>& running) {
        Path path;
        auto current = start;
        path.push_back(current);
        visited[current] = true;
        running[current] = true;
        uint8_t c_extend = 0;

        for (int i = 0; i < 4; i++) {
            uint8_t c = Alphabet[i];
            int64_t next = forward(current, c);
            if (next == current) continue;
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
        for (int i = 0; i < 4; i++) {
            uint8_t c = Alphabet[i];
            int64_t prev = backward(current, c);
            if (prev == current) continue;
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
        for (auto const &x : kmers) {
            auto id = x.second;
            if (!visited[id]) {
                Path path = greedyPath(id, visited, running);
                if (path.size() > 1) // not adding pdCand to paths
                    paths.push_back(path);
            }
        }

        // attach pendants
        for (auto &path : paths) {
            for (auto &id : path) {
                if (id <= 0) continue;
                for (auto const &c : Alphabet) {
                    auto next = forward(id, c);
                    if (next == 0 || next == id) continue;
                    if (pdCands.find(next) == pdCands.end()) continue;
                    auto it = find(path.begin(), path.end(), id);
                    if (it != path.end()) {
                        path.insert(it + 1, -next);
                        pdCands.erase(next);
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
        
        int64_t uncovered = kmers.size() - kmers_covered; 
        if (!uncovered) {
            cout << "Successfully covered all kmers!" << endl;
        }
        else {
            cout << "There are " << uncovered << " nodes uncovered..." << endl;
        }
        cout << endl;
    }
};

int main() {
    // Construct a vector of k-mers from an input FASTA file(.fa).
    string _filename;
    int64_t _K;
    cout << "Input file: ";
    cin >> _filename;
    cout << "Set K to : ";
    cin >> _K;
    cout << endl;
    
    // Construct and show the edge-centric De Bruijn Graph
    DeBruijnGraph ecdbg = DeBruijnGraph(_filename, _K - 1, false);
    ecdbg.getKmers();
    ecdbg.buildGraph();
    cout << "Edge-centric De Bruijn Graph:" << endl;
    ecdbg.printGraph();
    cout << endl;

    // Construct and show the node-centric De Bruijn Graph
    DeBruijnGraph ncdbg = DeBruijnGraph(_filename , _K, true);
    ncdbg.getKmers();
    ncdbg.buildGraph();
    cout << "Node-centric De Bruijn Graph:" << endl;
    ncdbg.printGraph();
    cout << endl;

    // Find paths from a given adjacency list
    Paths paths = ncdbg.findPaths();

    // Output a necklace cover and a comparison with Eulertigs
    ncdbg.printResult(paths);

    cout << string(60, '-') << endl;
    cout << "# Edges to add to make Eulerian: " << ecdbg.makeEulerian << endl;
    cout << "# Open necklaces: " << ncdbg.countOpenNecklaces << endl;
    cout << string(60, '-') << endl;

    return 0;
}