#include <iostream>
#include <vector>
#include <limits>
#include <map>

using namespace std;

// Define infinity as the maximum value for a double
const double INF = numeric_limits<double>::infinity();

//global variable V and graph
#define V 4
vector<vector<double> > graph(V, vector<double>(V, 0));

vector<int> right(vector<int> path)
{
    path.erase(path.begin());
    return path;
}
vector<int> left(vector<int> path)
{
    path.erase(path.end()-1);
    return path;
}
map< vector<int>, vector< vector<int> > > L;
map< vector<int>, vector< vector<int> > > R;
map< vector<int>, vector< vector<int> > > Lstar;
map< vector<int>, vector< vector<int> > > Rstar;
vector<vector<int> > P[V][V];
vector<vector<int> > Pstar[V][V];

bool set_init()
{
    //todo:init needed sets
}

// Function to run the Floyd-Warshall algorithm and return the distance matrix
vector<vector<double> > floyd_warshall(const vector<vector<double> >& graph) {
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

// Function for ISP proposed by DEMETRESCU and ITALIANO in 2004
vector<vector<double> > incremental_2004(vector<vector<double> >& distance,
                                        int source, int u, int v, double w) {
    
    graph[u][v] = w; 
    //Step 1: Clean up
    vector< vector<int> > Q;
    vector<int> path_u;
    path_u.push_back(u);
    Q.push_back(path_u);
    while(!Q.empty())
    {
        vector<int> pi = Q.back();
        Q.pop_back();
        vector< vector<int> > Leftlist = L[pi];
        vector< vector<int> > Rightlist = R[pi];
        vector< vector<int> > LR = Leftlist;
        LR.insert(LR.end(), Rightlist.begin(), Rightlist.end());
        for (int i = 0;i < LR.size();i++)
        {
            vector<int> pi_xy = LR[i];
            Q.push_back(pi_xy);
            int x = pi_xy.front();
            int y = pi_xy.back();
            for(int i = 0; i < P[x][y].size(); i++)
            {
                if(P[x][y][i] == pi_xy)
                {
                    P[x][y].erase(P[x][y].begin()+i);
                    break;
                }
            }
            for(int i = 0; i < Pstar[x][y].size(); i++)
            {
                if(Pstar[x][y][i] == pi_xy)
                {
                    Pstar[x][y].erase(Pstar[x][y].begin()+i);
                    break;
                }
            }
            vector< vector<int> > Lrpi = L[right(pi_xy)];
            for(int i = 0;i < Lrpi.size();i++)
            {
                if (Lrpi[i] == pi_xy)
                {
                    Lrpi.erase(Lrpi.begin()+i);
                    L[right(pi_xy)] = Lrpi;
                    break;
                }
            }
            vector< vector<int> > Rlpi = R[left(pi_xy)];
            for(int i = 0;i < Rlpi.size();i++)
            {
                if (Rlpi[i] == pi_xy)
                {
                    Rlpi.erase(Rlpi.begin()+i);
                    R[left(pi_xy)] = Rlpi;
                    break;
                }
            }
            vector< vector<int> > Lstar_rpi = Lstar[right(pi_xy)];
            for(int i = 0;i < Lstar_rpi.size();i++)
            {
                if (Lstar_rpi[i] == pi_xy)
                {
                    Lstar_rpi.erase(Lstar_rpi.begin()+i);
                    Lstar[right(pi_xy)] = Lstar_rpi;
                    break;
                }
            }
            vector< vector<int> > Rstar_rpi = Rstar[left(pi_xy)];
            for(int i = 0;i < Rstar_rpi.size();i++)
            {
                if (Rstar_rpi[i] == pi_xy)
                {
                    Rstar_rpi.erase(Rstar_rpi.begin()+i);
                    Rstar[left(pi_xy)] = Lstar_rpi;
                    break;
                }
            }
        }
    }
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
    distance = incremental_2004(distance, source, u, v, w);

    // Print the updated distance row from the source to all nodes
    cout << "\nUpdated distances after applying incremental APSP:\n";
    print_graph(graph);
    vector<double> path;

    return 0;
}
