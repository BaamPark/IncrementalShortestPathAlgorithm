#include <iostream>
#include <vector>
#include <limits>
#include <cassert>
#include <cstdlib>
#include <chrono>
#include <fstream>

using namespace std;

const double INF = numeric_limits<double>::infinity();

vector<vector<double>> generate_random_graph(int V, int E) {
    vector<vector<double>> graph(V, vector<double>(V, INF));

    for (int i = 0; i < E; ++i) {
        int u = rand() % V;
        int v = rand() % V;
        double w = rand() % 100 + 1;

        // Ensure we don't set a node to itself
        if (u != v) {
            graph[u][v] = w;
            graph[v][u] = w;  // Ensure symmetry
        }
    }

    // Diagonal should be 0
    for (int i = 0; i < V; ++i) {
        graph[i][i] = 0;
    }

    return graph;
}

/* vector<vector<double>> floyd_warshall(const vector<vector<double>>& graph) {
    int V = graph.size();
    vector<vector<double>> distance = graph;

    for (int k = 0; k < V; ++k) {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (distance[i][k] < INF && distance[k][j] < INF) {
                    distance[i][j] = min(distance[i][j], distance[i][k] + distance[k][j]);
                }
            }
        }
    }

    return distance;
} */

// Optimised version of the floyd warshall function 
vector<vector<double>> floyd_warshall_optimized(const vector<vector<double>>& graph) {
    int V = graph.size();
    vector<vector<double>> distance = graph;

    for (int k = 0; k < V; ++k) {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (distance[i][k] != INF && distance[k][j] != INF) {
                    distance[i][j] = min(distance[i][j], distance[i][k] + distance[k][j]);
                }
            }
        }
    }

    return distance;
}

int main() {
    int V = 20;  // Number of vertices
    int E = 25;  // Number of edges

    // Generate random graph
    vector<vector<double>> graph = generate_random_graph(V, E);

    // Validate graph
    bool valid = true;
    for (int i = 0; i < V; ++i) {
        if (graph[i][i] != 0) {
            valid = false;
            break;
        }
        for (int j = i + 1; j < V; ++j) {
            if (graph[i][j] != graph[j][i]) {
                valid = false;
                break;
            }
        }
    }

    if (!valid) {
        cout << "Error: Graph is not symmetric." << endl;
        cout << "Dataset validation failed!" << endl;
        return 1;
    }

    // Measure execution time
    auto start = chrono::high_resolution_clock::now();
    
    // Run Floyd-Warshall algorithm
    vector<vector<double>> distance = floyd_warshall_optimized(graph);

     // Save graph to CSV
    ofstream graph_file("graph.csv");
    for (const auto& row : graph) {
        for (double val : row) {
            graph_file << val << ",";
        }
        graph_file << endl;
    }
    graph_file.close();

    // Save distances to CSV
    ofstream distance_file("distance.csv");
    for (const auto& row : distance) {
        for (double val : row) {
            distance_file << val << ",";
        }
        distance_file << endl;
    }
    distance_file.close();

    auto stop = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

    cout << "Dataset validation passed!" << endl;

    // Optionally, print the distances matrix
    for (const auto& row : distance) {
        for (double val : row) {
            if (val == INF) {
                cout << "INF\t";
            } else {
                cout << val << "\t";
            }
        }
        cout << endl;
    }

    return 0;
}
