//Written by Beomseok Park
#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <chrono>
#include <cstdlib>

using namespace std::chrono;
using namespace std;

// Define infinity as the maximum value for a double
const double INF = numeric_limits<double>::infinity();

// Generate a random directed weighted graph
vector<vector<double>> generate_random_graph(int V, int E) {
    vector<vector<double>> graph(V, vector<double>(V, INF));

    for (int i = 0; i < E; ++i) {
        int u = rand() % V;
        int v = rand() % V;
        double w = rand() % 100 + 1;  // Edge weight between 1 and 100

        // Ensure we don't set a node to itself and only set the directed edge from u to v
        if (u != v) {
            graph[u][v] = w;  // Only set the edge from u to v
        }
    }

    // Diagonal should be 0, indicating no self-loops
    for (int i = 0; i < V; ++i) {
        graph[i][i] = 0;
    }

    return graph;
}


//this function doesn't modify graph but return u,v, and w for edge update
tuple<int, int, double> pick_random_edge_for_update(const vector<vector<double>>& graph) {
    int V = graph.size();
    vector<tuple<int, int, double>> edges;

    // Collect all edges with a set weight
    for (int u = 0; u < V; ++u) {
        for (int v = 0; v < V; ++v) {
            if (graph[u][v] != INF && graph[u][v] != 0 && u != v) {
                edges.push_back(make_tuple(u, v, graph[u][v]));
            }
        }
    }

    // Randomly select one edge and suggest a new weight less than the current weight
    if (!edges.empty()) {
        int index = rand() % edges.size();
        int u = get<0>(edges[index]);
        int v = get<1>(edges[index]);
        double current_weight = get<2>(edges[index]);
        double new_weight = (rand() % static_cast<int>(current_weight)) + 1;  // Ensure new weight is less and > 0

        return make_tuple(u, v, new_weight);
    } else {
        throw std::runtime_error("No valid edges found in the graph for updating.");
    }
}


// Function to run the Floyd-Warshall algorithm and return the distance matrix
//O(V3)
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

//O(V3)
vector<vector<double>> update_and_floyd_warshall(vector<vector<double>>& graph, int u, int v, double w) {
    // Update the graph with the new weight for the edge between u and v
    graph[u][v] = w;

    // Run the Floyd-Warshall algorithm on the updated graph
    return floyd_warshall(graph);
}

//O(m)
vector<int> find_affected_source(vector<vector<double>>& distance, int u, int v, double w) {
    int n = distance.size();

    vector<bool> vis(n, false);
    vector<int> affected_sources;

    if (distance[u][v] > w) {
        queue<int> Q;
        Q.push(u);
        affected_sources.push_back(u);
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
// time complexity O(V^2), Θ(∣∣S(v)∣∣+∑∣∣T(x)∣∣).
vector<vector<double>> incremental_apsp_pr(vector<vector<double>>& distance, int u, int v, double w) {

    vector<int> sources = find_affected_source(distance, u, v, w);
    int V = distance.size();
    for (int s : sources) {
        vector<bool> vis(V, false);
        queue<int> Q;

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
//O(V^2) Θ(∣∣S(v)∣∣+∣∣T(u)∣∣+∑∣S(P(y))∣)
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


void printVector(const std::vector<int>& sources) {
    for (int num : sources) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
}

int main() {
    int V = 500;  // Number of vertices
    int E = 50;  // Number of edges

    vector<vector<double>> graph = generate_random_graph(V, E);
    auto [u, v, w] = pick_random_edge_for_update(graph);
    std::cout << "Updating edge (" << u << ", " << v << ") from weight " << graph[u][v] << " to " << w << std::endl;

    // Run Floyd-Warshall algorithm and return the distance matrix
    vector<vector<double>> distance = floyd_warshall(graph);
    vector<vector<double>> distance_for_PR = distance;
    vector<vector<double>> distance_for_QUINCA = distance;
    vector<vector<double>> distance_naive;

    //Measure naive floyd warshall algorithm
    auto start_naive = high_resolution_clock::now();
    distance_naive = update_and_floyd_warshall(graph, u, v, w);
    auto end_naive = high_resolution_clock::now();
    auto duration_naive = duration_cast<microseconds>(end_naive - start_naive);
    cout << "Time taken by Naive Floyd Warshall algorithm: " << duration_naive.count() << " microseconds\n";

    // Measure PR algorithm
    auto start_pr = high_resolution_clock::now();
    distance_for_PR = incremental_apsp_pr(distance_for_PR, u, v, w);
    auto end_pr = high_resolution_clock::now();
    auto duration_pr = duration_cast<microseconds>(end_pr - start_pr);
    cout << "Time taken by PR algorithm: " << duration_pr.count() << " microseconds\n";

    // Measure QUINCA algorithm
    auto start_quinca = high_resolution_clock::now();
    distance_for_QUINCA = incremental_apsp_quinca(distance_for_QUINCA, u, v, w);
    auto end_quinca = high_resolution_clock::now();
    auto duration_quinca = duration_cast<microseconds>(end_quinca - start_quinca);
    cout << "Time taken by QUINCA algorithm: " << duration_quinca.count() << " microseconds\n";

    return 0;
}
