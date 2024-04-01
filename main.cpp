#include <iostream>
#include <vector>
#include <limits>

using namespace std;

// Define infinity as the maximum value for a double
const double INF = numeric_limits<double>::infinity();

// Function to run the Floyd-Warshall algorithm and return the distance matrix
vector<vector<double>> floyd_warshall(const vector<vector<double>>& graph) {
    int V = graph.size();
    vector<vector<double>> distance = graph; // Start with the adjacency matrix as a base for the distance matrix

    // Initialize the distance matrix
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (i != j && graph[i][j] == 0)
                distance[i][j] = INF;
        }
    }

    // Floyd-Warshall algorithm
    for (int k = 0; k < V; ++k) {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (distance[i][k] < INF && distance[k][j] < INF)
                    distance[i][j] = min(distance[i][j], distance[i][k] + distance[k][j]);
            }
        }
    }

    return distance;
}

// Function for dynamic shortest path algorithms by Ramalingam and Rep
// time complexity O(V+E)
vector<vector<double>> incremental_apsp(const vector<vector<double>>& graph, vector<vector<double>>& distance,
                                        int source, int u, int v, double w) {
    int V = graph.size();
    vector<vector<double>> new_graph = graph;
    new_graph[u][v] = w; // Update the graph with the new weight w for edge u->v

    vector<bool> vis(V, false);
    vector<int> Q;

    // Update the distance from source to v if the new path is shorter
    if (distance[source][v] > distance[source][u] + w) {
        distance[source][v] = distance[source][u] + w;
        Q.push_back(v);
        vis[v] = true;
    }

    // Truncated BFS (Breadth-First Search)
    while (!Q.empty()) {
        int y = Q.front();
        Q.erase(Q.begin());
        distance[source][y] = distance[source][u] + w + distance[v][y];
        for (int w_index = 0; w_index < V; ++w_index) {
            if (new_graph[y][w_index] > 0 && !vis[w_index] && distance[source][w_index] > distance[source][u] + w + distance[v][w_index]) {
                vis[w_index] = true;
                Q.push_back(w_index);
            }
        }
    }

    return distance;
}

int main() {
    int V = 4;
    vector<vector<double>> graph(V, vector<double>(V, 0));

    // Generate graph data
    graph[0][1] = 5;
    graph[0][3] = 10;
    graph[1][2] = 3;
    graph[2][3] = 1;

    // Run Floyd-Warshall algorithm and return the distance matrix
    vector<vector<double>> distance = floyd_warshall(graph);

    // Print the distance matrix
    cout << "Initial Distance Matrix:\n";
    for (const auto& row : distance) {
        for (double val : row) {
            if (val == INF)
                cout << "INF ";
            else
                cout << val << " ";
        }
        cout << "\n";
    }

    // Set source, u, v, and w
    int source = 0, u = 1, v = 2;
    double w = 2;

    // Run the incremental APSP algorithm
    distance = incremental_apsp(graph, distance, source, u, v, w);

    // Print the updated distance row from the source to all nodes
    cout << "\nUpdated distances from source " << source << " to all nodes after applying incremental APSP:\n";
    for (double dist : distance[source]) {
        if (dist == INF)
            cout << "INF ";
        else
            cout << dist << " ";
    }
    cout << "\n";

    return 0;
}
