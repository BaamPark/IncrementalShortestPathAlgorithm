#include <iostream>
#include <vector>
#include <limits>

using namespace std;

// Define infinity as the maximum value for a double
const double INF = numeric_limits<double>::infinity();

//global variable V and graph
int V = 4;
vector<vector<double> > graph(V, vector<double>(V, 0));

// Function to run the Floyd-Warshall algorithm and return the distance matrix
vector<vector<double> > floyd_warshall(const vector<vector<double> >& graph) {
    int V = graph.size();
    vector<vector<double> > distance = graph; // Start with the adjacency matrix as a base for the distance matrix

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
vector<vector<double> > incremental_cg(const vector<vector<double> >& graph, vector<vector<double> >& distance,
                                        int source, int u, int v, double w) {
    int V = graph.size();
    vector<vector<double> > new_graph = graph;
    new_graph[u][v] = w; // Update the graph with the new weight w for edge u->v
    //todo: implement the CG algorithm here.
    return distance;
}
void print_graph(const vector<vector<double> >& graph){
    for (const auto& row : graph) {
        for (double val : row) {
            if (val == INF)
                cout << "INF   ";
            else
                cout << val << "   ";
        }
        cout << "\n";
    }
    cout << "\n";
}
int main() {
    // Generate graph data
    graph[0][1] = 5.1;
    graph[0][3] = 2.5;
    graph[1][2] = 3.2;
    graph[2][3] = 1.8;

    printf("Matrix of graph\n");
    print_graph(graph);
    // Run Floyd-Warshall algorithm and return the distance matrix
    vector<vector<double> > distance = floyd_warshall(graph);

    printf("Initial distance:\n");
    print_graph(distance);

    // Set source, u, v, and w
    int source = 0, u = 1, v = 2;
    double w = 2;

    printf("add an edge with value %f from %d to %d",w, u,v);
    // Run the incremental APSP algorithm
    distance = incremental_cg(graph, distance, source, u, v, w);

    // Print the updated distance row from the source to all nodes
    cout << "\nUpdated distances from source " << source << " to all nodes after applying incremental APSP:\n";
    for (double dist : distance[source]) {
        if (dist == INF)
            cout << "INF   ";
        else
            cout << dist << "   ";
    }
    cout << "\n";

    return 0;
}
