#include <iostream>
#include <assert.h>
#include <vector>
#include <limits>

using namespace std;

const double INF = numeric_limits<double>::infinity();

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
}

int main() {
    test_floyd();
}