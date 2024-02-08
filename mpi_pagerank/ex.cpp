#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <mpi.h>

using namespace std;
using namespace chrono;

int ROOT = 0;

// Pagerank Algorithm parameter
double tolerance = 0.001;
double dampingFactor = 0.85;

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
        } while (diff > tolerance);

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

int main(int argc, char *argv[]) {
    // MPI initialize
    int rank, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (comm_size < 2) {
        cerr << "Error: This program requires at least 2 MPI processes." << endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Master process reads graph data
    if (rank == ROOT) {
        string filename = "../graph_data/facebook_combined.txt";
        Graph graph = readGraphFromFile(filename);

        int numNodes = graph.getNumNodes();
        int numEdges = graph.getNumEdges();

        // Distribute graph data to MPI processes
        vector<int> numEdgesPerProcess(comm_size, numEdges / comm_size);
        for (int i = 0; i < numEdges % comm_size; ++i) {
            numEdgesPerProcess[i]++;
        }

        vector<int> displacements(comm_size);
        displacements[0] = 0;
        for (int i = 1; i < comm_size; ++i) {
            displacements[i] = displacements[i - 1] + numEdgesPerProcess[i - 1];
        }

        vector<Edge> edgeData(numEdges);
        for (int i = 0; i < numEdges; ++i) {
            edgeData[i] = graph.edges[i];

            
        }

        MPI_Scatterv(edgeData.data(), numEdgesPerProcess.data(), displacements.data(), MPI_INT, edgeData.data(), numEdgesPerProcess[rank], MPI_INT, ROOT, MPI_COMM_WORLD);

        // PageRank calculation
        auto start = high_resolution_clock::now();
        unordered_map<int, double> nodeScores = graph.pageRank(tolerance, dampingFactor);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);

        // Gather results from MPI processes
        vector<pair<int, double>> results;
        for (int i = 0; i < comm_size; ++i) {
            results.emplace_back(i, nodeScores[i]);
        }

        // Output results
        cout << "PageRank scores:" << endl;
        for (const auto& result : results) {
            cout << "Node " << result.first << ": " << result.second << endl;
        }

        cout << "PageRank calculation time: " << duration.count() << " seconds" << endl;
    } else {
        // Worker processes receive graph data
        int numEdges;
        MPI_Scatter(NULL, 0, MPI_INT, &numEdges, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

        vector<Edge> edgeData(numEdges);
        MPI_Scatter(NULL, numEdges * 2, MPI_INT, edgeData.data(), numEdges * 2, MPI_INT, ROOT, MPI_COMM_WORLD);

        // PageRank calculation (worker processes don't calculate)
    }

    MPI_Finalize();
    return 0;
}
