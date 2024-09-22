#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
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
    NodeId idPosition = {0}; // for each ID a position i such that ID corresponds to kmer sequence.substr(i, K)
    NodeId pdCands; // candidates of pendants <=> length-1 found paths

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
        for (int i = 0; i < sequence.length() - K + 1; i++){
            auto kmer = sequence.substr(i, K);
            if (kmers.find(kmer) == kmers.end()){ // newly found kmer
                kmers[kmer] = idPosition.size() + 1; // "+ 1" because we want kmers ID starting from 1 to N
                idPosition.push_back(i);
            }
        }
    }

    int64_t forward(int64_t id, int8_t c){
        auto position = idPosition[id - 1]; // "id - 1" because we have kmers ID starting from 1
        string kmerSuffix = sequence.substr(position + 1, K - 1); 
        auto next = kmerSuffix + to_string(c);
        if (kmers.find(next) != kmers.end()) // kmer exists
            return kmers[next];
        else 
            return 0; // no outgoing branching exists for c
    }

    int64_t backward(int64_t id, int8_t c){
        auto position = idPosition[id - 1]; // "id - 1" because we have kmers ID starting from 1
        string kmerPrefix = sequence.substr(position, K - 1);
        auto previous = to_string(c) + kmerPrefix;
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
                auto next = kmerSuffix + to_string(c);
                if (kmers.find(next) != kmers.end()) { // kmer exists, so its node
                    auto nextid = kmers[next];
                    outNeighbors[id].push_back(nextid);
                    inNeighbors[nextid].push_back(id);
                }       
            }
        }
    }

    void buildEdgeCentricGraph(const string S, const Kmers& kmers, const NodeId& idPos) {  // Rough invariant: kmers.first == S.substr(kmers.second - 1, k)
    }

    // Display the graph
    void printGraph() {
        for (auto const &x : kmers) { 
            cout << "Node: " << x.first << " -> ";  // kmer string
            for (auto const& idNeigh : outNeighbors[x.second]) {
                auto position = idPosition[idNeigh - 1]; // "idNeigh - 1" because we have kmers ID starting from 1
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
        
        // suggestion: from current choose the best neighbor: first choice is cycle, second choice is extend path, and third choice CASE 1 if neither of the first two choices occur.

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
        int64_t i = 0;
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

    Paths findPaths() {
        Paths paths;
        bool visited[] = {false};
        bool running[] = {false};

        // Start greedy search from an unvisited node
        for (auto const &x : kmers) {
            auto id = x.second;
            if (!visited[id]) {
                Path path = greedyPath(id, visited, running);
                paths.push_back(path);
            }
        }

        // Add isolated kmers to paths as "Length-one path"s.
        for (auto const &x : kmers) {
            auto id = x.second;
            if (!visited[id]) {
                paths.push_back(deque<int64_t>{id});
                pdCands.push_back(id);
            }
        }

        // // Delete empty paths.
        // paths.erase(
        //     remove_if(paths.begin(), paths.end(), 
        //         [](const Path& path) {
        //             return path.empty();
        //         }
        //     ), 
        //     paths.end()
        // );

        return paths;
    }

    Paths attachPendants(Paths paths) {
        // Attach length-1 paths to an adjacent path
        for (auto id_cand : pdCands) {
            for (auto path : paths) {
                if (path.size() == 1) continue; // Skip pendant candidates 

                for (auto const& id : path) {
                    if (!id) continue; // skip the already attached pendants
                    string pendantPrefix = sequence.substr(idPosition[id_cand - 1], K - 1);
                    string attacheePrefix = sequence.substr(idPosition[id - 1], K - 1);
                    
                    // Check if pendant can be attached to attachee
                    if (pendantPrefix != attacheePrefix) continue; // skip kmer which cannot be attached the pendant to 
                    // Attach the pendant
                    auto after_this = find(path.begin(), path.end(), id);
                    path.insert(after_this + 1, -id_cand);
                    goto next_id_cand;
                }
            }
            next_id_cand:
            id_cand *= -1;
        }

        // Remove length-1 paths that were attached
        paths.erase(
            remove_if(
                paths.begin(), paths.end(), [this](const Path& path){return is_attached(path, pdCands);}
            ), paths.end()
        );
        
        return paths;
    }

    bool is_attached(const Path& path, const NodeId& pdCands) {
        return (path.size() == 1 && find(pdCands.begin(), pdCands.end(), -path[0]) != pdCands.end());
    }

    void printResult(const Paths& paths) {
        // Output the resulting paths
        int64_t kmers_covered = 0;
        int64_t count = 1;
        for (auto const& path : paths) {
            cout << "Path" << count << " : ";
            for (auto const& id : path) {
                if (id > 0) {
                    cout << sequence.substr(idPosition[abs(id) - 1], K) << " ";
                    kmers_covered += 1;
                }
                else if (id < 0) {
                    cout << "+" << sequence.substr(idPosition[abs(id) - 1], K) << " ";
                    kmers_covered += 2;
                }
                else { // id == 0
                    cout << "*" << endl;
                }
            }
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
    cout << "Set K to : ";
    cin >> _K;
    cout << "Input file: ";
    cin >> filename;
    
    // Construct and show the node-centric De Bruijn Graph
    DeBruijnGraph ncdbg = DeBruijnGraph(_K, true);
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
    ncdbg.printResult(necklaces);

    return 0;
}