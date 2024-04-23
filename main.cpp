//Written by Beomseok Park
#include <iostream>
#include <vector>
#include <limits>
#include <queue>

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


vector<int> find_affected_source(vector<vector<double>>& distance, int u, int v, double w) {
    int n = distance.size();

    vector<bool> vis(n, false);
    vector<int> affected_sources;

    if (distance[u][v] > w) {
        queue<int> Q;
        Q.push(u);
        vis[u] = true;

        while (!Q.empty()) {
            int x = Q.front();
            Q.pop();

            for (int z = 0; z < n; ++z) {
                if (!vis[z] && distance[z][v] > distance[z][u] + w) {
                    Q.push(z);
                    vis[z] = true;
                    affected_sources.push_back(z);
                }
            }
        }

        // Reset vis vector to false for all vertices.
        fill(vis.begin(), vis.end(), false);
    }

    return affected_sources;
}

// PR algorithm proposed by Ramalingam and Rep
// time complexity O(V+E)
vector<vector<double>> incremental_apsp_pr(vector<vector<double>>& distance, int u, int v, double w) {

    vector<int> sources = find_affected_source(distance, u, v, w);
    int V = distance.size();
    for (int s : sources) {
        vector<bool> vis(V, false);
        queue<int> Q;

        distance[u][v] = w;
        // Update the distance from source to v if the new path is shorter
        if (distance[s][v] > distance[s][u] + w) {
            distance[s][v] = distance[s][u] + w;
            Q.push(v);
            vis[v] = true;


            // Truncated BFS (Breadth-First Search)
            while (!Q.empty()) {
                int y = Q.front();
                Q.pop();
                distance[s][y] = distance[s][u] + w + distance[v][y];
                for (int w_index = 0; w_index < V; ++w_index) {
                    if (!vis[w_index] && distance[s][w_index] > distance[s][u] + w + distance[v][w_index]) {
                        vis[w_index] = true;
                        Q.push(w_index);
                    }
                }
            }
        }
    }
    return distance;
}

//QUINCA algorithm proposed by Slobbe, Bergamini, and Meyerhenke
vector<vector<double>> incremental_apsp_quinca(vector<vector<double>>& distances, int u, int v, double w_prime) {
    int n = distances.size();
    vector<int> affected_sources = find_affected_source(distances, u, v, w_prime);

    if (distances[u][v] > w_prime) {
        queue<int> Q;
        vector<bool> visited(n, false);
        distances[u][v] = w_prime; // Update the distance for the edge (u,v)

        Q.push(v);
        visited[v] = true;

        while (!Q.empty()) {
            int y = Q.front();
            Q.pop();

            for (int x : affected_sources) {
                if (distances[x][y] > distances[x][u] + w_prime + distances[v][y]) {
                    distances[x][y] = distances[x][u] + w_prime + distances[v][y];
                }
            }

            for (int w = 0; w < n; ++w) {
                if (!visited[w] && distances[u][w] > w_prime + distances[v][w]) {
                    Q.push(w);
                    visited[w] = true;
                }
            }
        }
    }

    return distances; // Return the updated distance matrix
}

void print_distance(vector<vector<double>> distance) {
    for (const auto& row : distance) {
        for (double val : row) {
            if (val == INF)
                cout << "INF ";
            else
                cout << val << " ";
        }
        cout << "\n";
    }
}

int main() {
    int V = 6;
    vector<vector<double>> graph(V, vector<double>(V, 0));

    graph[0][2] = 1;
    graph[1][2] = 2;
    graph[3][4] = 2;
    graph[3][5] = 1;
    graph[4][0] = 3;
    graph[5][1] = 3;

    // Run Floyd-Warshall algorithm and return the distance matrix
    vector<vector<double>> distance = floyd_warshall(graph);
    vector<vector<double>> distance_for_PR = distance;
    vector<vector<double>> distance_for_QUINCA = distance;
    // Print the distance matrix
    cout << "Initial Distance Matrix:\n";
    print_distance(distance);

    int u = 2, v = 3;
    double w = 1;
    vector<int> sources = find_affected_source(distance_for_PR, u, v, w);

    distance_for_PR = incremental_apsp_pr(distance_for_PR, u, v, w);
    cout << "Updated Distance Matrix by PR algorithm:\n";
    print_distance(distance_for_PR);

    distance_for_QUINCA = incremental_apsp_quinca(distance_for_QUINCA, u, v, w);
    cout << "Updated Distance Matrix by PR algorithm:\n";
    print_distance(distance_for_QUINCA);

    return 0;
}
