/*
The program tests the Floyd-Warshall algorithm for all-pairs shortest paths and applies incremental updates to a randomly generated graph, 
updating the shortest path distances accordingly. 
*/ 

#include <iostream>
#include <vector>
#include <limits>
#include <cassert>
#include <cstdlib>

using namespace std;

const double INF = numeric_limits<double>::infinity();

// Floyd-Warshall algorithm
vector<vector<double> > floyd_warshall(const vector<vector<double> > &graph)
{
    int V = graph.size();
    vector<vector<double> > distance = graph;

    for (int i = 0; i < V; ++i)
    {
        for (int j = 0; j < V; ++j)
        {
            if (i != j && graph[i][j] == 0)
                distance[i][j] = INF;
        }
    }

    for (int k = 0; k < V; ++k)
    {
        for (int i = 0; i < V; ++i)
        {
            for (int j = 0; j < V; ++j)
            {
                if (distance[i][k] < INF && distance[k][j] < INF)
                    distance[i][j] = min(distance[i][j], distance[i][k] + distance[k][j]);
            }
        }
    }

    return distance;
}

// Incremental APSP algorithm
vector<vector<double> > incremental_apsp(vector<vector<double> > &graph, vector<vector<double> > &distance,
                                        int source, int u, int v, double w)
{
    int V = graph.size();
    graph[u][v] = w;

    vector<bool> vis(V, false);
    vector<int> Q;

    if (distance[source][v] > distance[source][u] + w)
    {
        distance[source][v] = distance[source][u] + w;
        Q.push_back(v);
        vis[v] = true;
    }

    while (!Q.empty())
    {
        int y = Q.front();
        Q.erase(Q.begin());
        distance[source][y] = distance[source][u] + w + distance[v][y];
        for (int w_index = 0; w_index < V; ++w_index)
        {
            if (graph[y][w_index] > 0 && !vis[w_index] && distance[source][w_index] > distance[source][u] + w + distance[v][w_index])
            {
                vis[w_index] = true;
                Q.push_back(w_index);
            }
        }
    }

    return distance;
}

vector<vector<double> > generate_random_graph(int V, int E)
{
    vector<vector<double> > graph(V, vector<double>(V, INF));
    for (int i = 0; i < E; ++i)
    {
        int u = rand() % V;
        int v = rand() % V;
        double w = rand() % 100 + 1;
        graph[u][v] = w;
    }
    for (int i = 0; i < V; ++i)
    {
        graph[i][i] = 0;
    }
    return graph;
}

void random_incremental_updates(vector<vector<double> > &graph, vector<vector<double> > &distance)
{
    int V = graph.size();
    for (int i = 0; i < 10; ++i)
    {
        int u = rand() % V;
        int v = rand() % V;
        double w = rand() % 100 + 1;
        graph[u][v] = w;
        distance = incremental_apsp(graph, distance, 0, u, v, w);
    }
}

void test_floyd()
{
    // Test graphs initialization
    vector<vector<double> > graph1(4, vector<double>(4));
    graph1[0][0] = 0; graph1[0][1] = 5; graph1[0][2] = INF; graph1[0][3] = 10;
    graph1[1][0] = INF; graph1[1][1] = 0; graph1[1][2] = 3; graph1[1][3] = INF;
    graph1[2][0] = INF; graph1[2][1] = INF; graph1[2][2] = 0; graph1[2][3] = 1;
    graph1[3][0] = INF; graph1[3][1] = INF; graph1[3][2] = INF; graph1[3][3] = 0;

    vector<vector<double> > graph2(4, vector<double>(4));
    graph2[0][0] = 0; graph2[0][1] = 5.1; graph2[0][2] = INF; graph2[0][3] = 2.5;
    graph2[1][0] = INF; graph2[1][1] = 0; graph2[1][2] = 3.2; graph2[1][3] = INF;
    graph2[2][0] = INF; graph2[2][1] = INF; graph2[2][2] = 0; graph2[2][3] = 1.8;
    graph2[3][0] = INF; graph2[3][1] = INF; graph2[3][2] = INF; graph2[3][3] = 0;

    vector<vector<double> > graph3(4, vector<double>(4));
    graph3[0][0] = 0; graph3[0][1] = 5.1; graph3[0][2] = 8.3; graph3[0][3] = 2.5;
    graph3[1][0] = INF; graph3[1][1] = 0; graph3[1][2] = 3.2; graph3[1][3] = 5;
    graph3[2][0] = INF; graph3[2][1] = INF; graph3[2][2] = 0; graph3[2][3] = 1.8;
    graph3[3][0] = INF; graph3[3][1] = INF; graph3[3][2] = INF; graph3[3][3] = 0;

    vector<vector<double> > graph4(5, vector<double>(5));
    graph4[0][0] = 0; graph4[0][1] = 2500; graph4[0][2] = 4200; graph4[0][3] = 4200; graph4[0][4] = 3100;
    graph4[1][0] = 2500; graph4[1][1] = 0; graph4[1][2] = 1700; graph4[1][3] = 3900; graph4[1][4] = 5000;
    graph4[2][0] = 4200; graph4[2][1] = 1700; graph4[2][2] = 0; graph4[2][3] = 2200; graph4[2][4] = 3300;
    graph4[3][0] = 4200; graph4[3][1] = 3900; graph4[3][2] = 2200; graph4[3][3] = 0; graph4[3][4] = 1100;
    graph4[4][0] = 3100; graph4[4][1] = 5000; graph4[4][2] = 3300; graph4[4][3] = 1100; graph4[4][4] = 0;

    vector<vector<double> > graph5(1, vector<double>(1, 0));

    assert(floyd_warshall(graph1) == floyd_warshall(graph1));
    assert(floyd_warshall(graph2) == floyd_warshall(graph2));
    assert(floyd_warshall(graph3) == floyd_warshall(graph3));
    assert(floyd_warshall(graph4) == floyd_warshall(graph4));
    assert(floyd_warshall(graph5) == floyd_warshall(graph5));
}

int main()
{
    // Test Floyd-Warshall algorithm
    cout << "Testing Floyd-Warshall algorithm..." << endl;
    test_floyd();
    cout << "Floyd-Warshall tests passed!" << endl << endl;

    // Generate a random graph and run incremental updates
    vector<vector<double> > graph = generate_random_graph(10, 20);
    vector<vector<double> > distance = floyd_warshall(graph);

    cout << "Randomly generated graph:" << endl;
    for (const auto& row : graph) {
        for (double val : row) {
            cout << (val == INF ? "INF" : to_string(val)) << "\t";
        }
        cout << endl;
    }
    cout << endl;

    random_incremental_updates(graph, distance);

    cout << "After incremental updates:" << endl;
    for (const auto& row : distance) {
        for (double val : row) {
            cout << (val == INF ? "INF" : to_string(val)) << "\t";
        }
        cout << endl;
    }
    cout << endl;

    return 0;
}