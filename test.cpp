#include <iostream>
#include <assert.h>
#include <vector>
#include <limits>
#include <queue>
#include <cstdlib>
using namespace std;

const double INF = numeric_limits<double>::infinity();

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


vector<vector<double>> update_and_floyd_warshall(vector<vector<double>>& graph, int u, int v, double w) {
    // Update the graph with the new weight for the edge between u and v
    graph[u][v] = w;

    // Run the Floyd-Warshall algorithm on the updated graph
    return floyd_warshall(graph);
}


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
// time complexity O(V+E)
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

void test_floyd() {
    // Test case 1: Graph with 3 vertices in a straight line
    vector<vector<double>> graph1 = {
            {0, 5, 0, 10},
            {0, 0, 3, 0},
            {0, 0, 0, 1},
            {0, 0, 0, 0}
    };

    //generate file
    vector<vector<double>> expected1 = {
            {0, 5, 8, 9},
            {INF, 0, 3, 4},
            {INF, INF, 0, 1},
            {INF, INF, INF, 0}
    };

    vector<vector<double>> graph2 = {
            {0, 5.1, 0, 2.5},
            {0, 0, 3.2, 0},
            {0, 0, 0, 1.8},
            {0, 0, 0, 0}
    };
    vector<vector<double>> expected2 = {
            {0, 5.1, 8.3, 2.5},
            {INF, 0, 3.2, 5},
            {INF, INF, 0, 1.8},
            {INF, INF, INF, 0}
    };

    vector<vector<double>> graph3 = {
            {0, 0},
            {0, 0}
    };

    vector<vector<double>> expected3 = {
            {0, INF},
            {INF, 0}
    };

    vector<vector<double>> graph4 = {
            {0, 2500, 0, 0, 3100},
            {2500, 0, 1700, 0, 0},
            {0, 1700, 0, 2200, 0},
            {0, 0, 2200, 0, 1100},
            {3100, 0, 0, 1100, 0}
    };

    vector<vector<double>> expected4 = {
            {0, 2500, 4200, 4200, 3100},
            {2500, 0, 1700, 3900, 5000},
            {4200, 1700, 0, 2200, 3300},
            {4200, 3900, 2200, 0, 1100},
            {3100, 5000, 3300, 1100, 0}
    };

    vector<vector<double>> graph5 = {
            {0}
    };

    vector<vector<double>> expected5 = {
            {0}
    };

    
    assert(floyd_warshall(graph1) == expected1); //pass int graph
    assert(floyd_warshall(graph2) == expected2); //pass double graph
    assert(floyd_warshall(graph3) == expected3); // pass zero entry 2x2 graph
    assert(floyd_warshall(graph4) == expected4); // pass large number graph
    assert(floyd_warshall(graph5) == expected5); // pass one vertex graph
    cout << "All tests passed for FW algorithm" << endl;
}

//The return value of update_and_floyd_warshall is considered to expected value
void test_pr() {
    // Set up for different tests
    int vertices[5] = {5, 10, 20, 30, 40};  // Number of vertices for each test
    int edges[5] = {1, 2, 3, 4, 5};    // Number of edges for each test

    for (int i = 0; i < 5; ++i) {
        int V = vertices[i];
        int E = edges[i];
        auto graph = generate_random_graph(V, E);


        auto [u, v, w] = pick_random_edge_for_update(graph);
        auto original_distance = floyd_warshall(graph);

        assert(incremental_apsp_pr(original_distance, u, v, w) == update_and_floyd_warshall(graph, u, v, w));

    }

    cout << "All tests passed for PR algorithm" << endl;
}


//The return value of update_and_floyd_warshall is considered to expected value
void test_quinca() {
    // Set up for different tests
    int vertices[5] = {5, 10, 20, 30, 40};  // Number of vertices for each test
    int edges[5] = {1, 2, 3, 4, 5};    // Number of edges for each test

    for (int i = 0; i < 5; ++i) {
        int V = vertices[i];
        int E = edges[i];
        auto graph = generate_random_graph(V, E);

        auto [u, v, w] = pick_random_edge_for_update(graph);
        auto original_distance = floyd_warshall(graph);

        assert(incremental_apsp_quinca(original_distance, u, v, w) == update_and_floyd_warshall(graph, u, v, w));

    }

    cout << "All tests passed for QUINCA algorithm" << endl;
}


int main() {
    test_floyd();
    test_pr();
    test_quinca();
}