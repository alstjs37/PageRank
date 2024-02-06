#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <chrono>

using namespace std;
using namespace chrono;

struct Edge {
    int src;
    int dest;
};

class Graph {
public:
    vector<Edge> edges;

    void addEdge(int src, int dest) {
        edges.push_back({src, dest});
    }

    int getNumNodes() const {
        unordered_set<int> nodes;
        for (const Edge& edge : edges) {
            nodes.insert(edge.src);
            nodes.insert(edge.dest);
        }
        return nodes.size();
    }

    int getNumEdges() const {
        return edges.size();
    }

    unordered_map<int, double> pageRank(double tolerance, double dampingFactor) {
        unordered_map<int, double> nodeScores;

        int numNodes = getNumNodes();

        // initialize pagerank value
        for (int node = 0; node < numNodes; ++node) {
            nodeScores[node] = 1.0 / numNodes;
        }

        double diff;
        int iteration = 0;

        do {
            unordered_map<int, double> newScores;

            // Calculate pagerank
            for (const Edge& edge : edges) {
                newScores[edge.dest] += (1.0 - dampingFactor) * (nodeScores[edge.src] / getOutDegree(edge.src));
            }

            // Calculate diff with before value
            diff = 0.0;
            for (int node = 0; node < numNodes; ++node) {
                diff += abs(newScores[node] - nodeScores[node]);
            }

            nodeScores = newScores;
            ++iteration;

            // print calculate step
            cout << "Calculate: " << iteration << " step" << endl;
        } while (diff > tolerance);

        cout << "Converged after " << iteration << " iterations." << endl;

        return nodeScores;
    }

private:
    int getOutDegree(int node) const {
        int outDegree = 0;
        for (const Edge& edge : edges) {
            if (edge.src == node) {
                ++outDegree;
            }
        }
        return outDegree;
    }
};

Graph readGraphFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    Graph graph;
    int src, dest;
    while (file >> src >> dest) {
        graph.addEdge(src, dest);
    }

    file.close();
    return graph;
}

int main() {

    string filename = "./facebook_combined.txt";

    Graph graph = readGraphFromFile(filename);

    int numNodes = graph.getNumNodes();
    int numEdges = graph.getNumEdges();

    cout << "Number of nodes: " << numNodes << endl;
    cout << "Number of edges: " << numEdges << endl;

    // Pagerank Algorithm parameter
    double tolerance = 0.00001;
    double dampingFactor = 0.85;

    // start time
    auto start = high_resolution_clock::now();

    // pagerank function
    unordered_map<int, double> nodeScores = graph.pageRank(tolerance, dampingFactor);

    // end time
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // sort result
    vector<pair<int, double>> sortedScores(nodeScores.begin(), nodeScores.end());
    sort(sortedScores.begin(), sortedScores.end(), [](const auto& a, const auto& b) -> bool {
        return a.second > b.second;
    });

    cout << endl;
    // print top 5 pagerank value
    cout << "Top 5 PageRank Scores:" << endl;
    int count = 0;
    for (const auto& pair : sortedScores) {
        cout << "Node " << pair.first << ": " << pair.second << endl;
        count++;
        if (count == 5) {
            break;
        }
    }
    // print Calculate time -> microseconds
    cout << "PageRank calculation time: " << duration.count() << " microseconds" << endl;

    return 0;
}
