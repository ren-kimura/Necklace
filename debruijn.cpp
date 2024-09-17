#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <deque>

using namespace std;
using AdjList = unordered_map<string, vector<string>>;

class EdgeCentric_DeBruijnGraph {
public:
    // Edges of the graph (previous node -> successive node)
    AdjList ec_adjacenceyList;

    // Construct edge-centric De Bruijn Graph from the given k-mer sets
    void buildGraph(const vector<string>& kmers) {
        for (const string& kmer : kmers) {
            //retrieve k-1 prefix and suffix
            string prefix = kmer.substr(0, kmer.size() - 1);
            string suffix = kmer.substr(1, kmer.size() - 1);

            // Append the edge (previous node -> successive node)
            ec_adjacenceyList[prefix].push_back(suffix);
        }
    }

    // Show the graph
    void printGraph() const {
        for (const auto& pair : ec_adjacenceyList) {
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
    AdjList nc_outAdjacencyList;
    // Edges of the graph (current node, previous nodes)
    AdjList nc_inAdjacencyList;

    // Construct node-centric De Bruijn Graph from the given k-mer sets
    void buildGraph(const vector<string>& kmers) {
        // Map current node to next nodes
        for (const string& kmer : kmers) {
            // Define node as a current looking k-mer, then determine the edges using the last letter
            string kmerSuffix = kmer.substr(1); // k-mer's last (k-1) letters
            
            // Search for the next nodes
            for (const string& nextKmer : kmers) {
                string nextKmerPrefix = nextKmer.substr(0, nextKmer.size() - 1); // next node's first (k-1) letters

                if (kmerSuffix == nextKmerPrefix) {
                    nc_outAdjacencyList[kmer].push_back(nextKmer);
                }
            }
        }
    }

    // Construct the reverse adjacency list of the above
    void buildReverseGraph(const vector<string>& kmers) {
        // Map current node to previous nodes
        for (const string& kmer : kmers) {
            // Define node as a current looking k-mer, then determine the edges using the first letter
            string kmerPrefix = kmer.substr(0, kmer.size() - 1); // k-mer's first (k-1) letters
            
            // Search for the next nodes
            for (const string& previousKmer : kmers) {
                string previousKmerSuffix = previousKmer.substr(1); // previous node's last (k-1) letters

                if (kmerPrefix == previousKmerSuffix) {
                    nc_inAdjacencyList[kmer].push_back(previousKmer);
                }
            }
        }
    }

    // Show the graph
    void printGraph() const {
        for (const auto& pair : nc_outAdjacencyList) {
            const string& node = pair.first;
            const vector<string>& neighbors = pair.second;
            cout << "Node: " << node << " has edges to: ";
            for (const string& neighbor : neighbors) {
                cout << neighbor << " ";
            }
            cout << endl;
        }
    }

    deque<string> greedyPath(const string& start, AdjList& adjList, AdjList& revAdjList, unordered_map<string, bool>& visited, unordered_map<string, bool>& running){
        // Define a deque of string "path" that contains k-mers in a path
        deque<string> path;
        string current = start;
        path.push_back(current);
        visited[current] = true;
        running[current] = true;
        bool tmp = false;

        while (true) {
            bool tmp1 = false;
            for (const string& neighbor : adjList[current]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    running[neighbor] = true;
                    path.push_back(neighbor);
                    current = neighbor;
                    tmp1 = true;
                    break;
                }
                if (running[neighbor]) {
                    // Update running[all kmers in <path>] = false.
                    for (const string& kmer : path) {
                        running[kmer] = false;
                    }
                    // Update visited[outside cycle] = false, and erase them.
                    auto it = find(path.begin(), path.end(), neighbor);
                    auto delete_before = distance(path.begin(), it);
                    for (auto i = 0; i < delete_before; i++) {
                        visited[path[0]] = false;
                        path.pop_front();
                    }
                    // Skip the reverse search because we already found a cycle.
                    tmp = true;
                    break;
                }
            }
            if(!tmp1) break;
        }
        while (true) {
            if(tmp) break;
            bool tmp2 = false;
            for (const string& neighbor : revAdjList[current]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    running[neighbor] = true;
                    path.push_front(neighbor);
                    current = neighbor;
                    tmp2 = true;
                    break;
                }
                if (running[neighbor]) {
                    // Update running[all kmers in <path>] = false.
                    for (const string& kmer : path) {
                        running[kmer] = false;
                    }
                    // Update visited[outside cycle] = false, and erase them.
                    auto it = find(path.begin(), path.end(), neighbor);
                    auto delete_from = distance(path.begin(), it); 
                    auto path_size = path.size();
                    for (auto i = path_size - 1; i > delete_from; i--) {
                        visited[path[path.size() - 1]] = false;
                        path.pop_back();
                    }
                    break;
                }
            }
            if(!tmp2) break;
        }
        return path;
    }

    vector<deque<string>> findPaths(AdjList& adjList, AdjList& revAdjList) {
        unordered_map<string, bool> visited;
        unordered_map<string, bool> running;
        vector<deque<string>> paths;

        // Initialize the visiting states of all nodes
        for (const auto& pair : adjList) {
            visited[pair.first] = false;
        }

        // Start greedy search from an unvisited node
        for (const auto& pair : adjList) {
            if (!visited[pair.first]) {
                deque<string> path = greedyPath(pair.first, adjList, revAdjList, visited, running);
                paths.push_back(path);
            }
        }

        // Add isolated kmers to paths as "Length-one path"s.
        for (const auto& pair : adjList) {
            if (!visited[pair.first]) {
                deque<string> tmp_dq = {pair.first};
                paths.push_back(tmp_dq);
            }
        }

        // Delete empty paths.
        paths.erase(
            remove_if(paths.begin(), paths.end(), 
                [](const deque<string>& dq) {
                    return dq.empty();
                }
            ), 
            paths.end()
        );

        return paths;
    }

    // vector<vector<string>> attachPendants(vector<vector<string>>& paths) {
    //     // for (size_t i = 0; i < paths.size(); ) {
    //     //     if (paths[i].size() != 1) {
    //     //         i++;
    //     //         continue;
    //     //     }
    //     //     for (size_t j = 0; j < paths.size(); j++) {
    //     //         if (paths[j].size() == 1) continue;
    //     //         for (size_t k = 0; k < paths[j].size(); k++) {
    //     //             if (paths[i][0].substr(0, paths[j][k].size() - 1) != paths[j][k].substr(0, paths[j][k].size() - 1)) continue;
    //     //             paths[j][k] = "(" + paths[j][k] + ", " + paths[i][0] + ")";
    //     //             // attached.push_back(paths[i][0]);
    //     //             break;
    //     //         }
    //     //         break;
    //     //     }
    //     //     vector<int>::iterator it = paths.begin() + i;
    //     //     paths.erase(it);
    //     //     // paths.erase(remove_if(paths.begin(), paths.end(), attachedLen1(paths[i][0], attached)), paths.end());
    //     // }

    //     // list of already attached length-1 path
    //     vector<string> attached;
    //     auto attachedLen1 = [](string kmer, vector<string> attached){return (find(attached.begin(), attached.end(), kmer) != attached.end());};
        
    //     // attach length-1 path to an adjacent path
    //     for (vector<string>& path : paths) {
    //         if (path.size() == 1) {
    //             for (vector<string>& mainPath : paths) {
    //                 if (mainPath.size() > 2) {
    //                     for (string& kmer : mainPath) {
    //                         if (path[0].substr(0, kmer.size() - 1) == kmer.substr(0, kmer.size() - 1)) {
    //                             kmer = "(" + kmer + ", " + path[0] + ")";
    //                             attached.push_back(path[0]);
    //                         }
    //                     }
    //                 }
    //             }
    //             // auto itr = paths.begin();
    //             // while (itr != paths.end()) {
    //             //     if((*itr) == path){
    //             //         itr = paths.erase(itr);
    //             //     }
    //             //     else itr++;
    //             // }
    //         }
    //     }
    //     paths.erase(remove_if(paths.begin(), paths.end(), attachedLen1(path, attached)), paths.end());
        
    //     return paths;
    // }
};



int main() {
    // Define k and a filename
    int k;
    string filename;

    cout << "Designate an input file(****.fa): ";
    cin >> filename;

    // Convert FASTA into a string
    ifstream inputFile(filename);
    if (!inputFile) {
        cerr << "Error opening input file." << endl;
        return 1;
    }

    string line;
    string sequence = "";  // To store the concatenated sequence

    while (getline(inputFile, line)) {
        if (line[0] != '>') {  // Skip header lines
            sequence += line;
        }
    }

    inputFile.close();

    // Erase letters other than "a,t,g,c"
    sequence.erase(remove_if(sequence.begin(), sequence.end(), [](char x) {return (x != 'a' && x != 't' && x != 'g' && x != 'c' && x != 'A' && x != 'T' && x != 'G' && x != 'C');}), sequence.end());

    // Output the result as a string
    cout << "Processed Sequence: " << sequence << "\n" << endl;
    
    // Set k
    cout << "Set k to: ";
    cin >> k;

    // Define dynamic array that has strings as an element
    vector<string> rep_kmers;

    // Append every k-mer in "sequence"
    for (int i = 0; i < sequence.length() - k; i++){
        rep_kmers.push_back(sequence.substr(i, k));
    }

    // Delete the repetitions
    set<string> uniqueSet(rep_kmers.begin(), rep_kmers.end());
    vector<string> kmers(uniqueSet.begin(), uniqueSet.end());
    
    // Output all k-mers
    cout << "Distinct k-mers: " << endl;
    for (const string& kmer : kmers) {
        cout << kmer << " ";
    }
    cout << endl;

    // Number of distinct k-mers
    cout << "# k-mers: " << kmers.size() << "\n" << endl;

    // Construct and show the edge-centric De Bruijn Graph
    EdgeCentric_DeBruijnGraph ecdbg;
    ecdbg.buildGraph(kmers);
    cout << "Edge-centric De Bruijn Graph:" << endl;
    ecdbg.printGraph();
    cout << endl;
    
    // Construct and show the node-centric De Bruijn Graph
    NodeCentric_DeBruijnGraph ncdbg;
    ncdbg.buildGraph(kmers);
    ncdbg.buildReverseGraph(kmers);
    cout << "Node-centric De Bruijn Graph:" << endl;
    ncdbg.printGraph();
    cout << endl;

    // Find paths from a given adjacency list
    vector<deque<string>> paths = ncdbg.findPaths(ncdbg.nc_outAdjacencyList, ncdbg.nc_inAdjacencyList);
    
    // Attach length-1 path to an adjacent path as a pendant
    // paths = ncdbg.attachPendants(paths);
    
    // Output the resulting paths
    int i = 1, kmers_covered = 0;
    for (const auto& path : paths) {
        cout << "Path" << i << " : ";
        for (const string& node : path) {
            cout << node << " ";
            if (node[0] == '(') kmers_covered += 1;
        }
        cout << endl;
        i++;
        kmers_covered += path.size();
    }
    cout << "Found " << i - 1 << " paths in total!" << endl;

    if (kmers_covered == kmers.size()) {
        cout << "Successfully covered all kmers!" << endl;
    }
    else {
        cout << "There are " << kmers.size() - kmers_covered << " nodes uncovered..." << endl;
    }

    return 0;
}