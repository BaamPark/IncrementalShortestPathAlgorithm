//Written by Beomseok Park
#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <chrono>
#include <cstdlib>
#include <map>
using namespace std::chrono;
using namespace std;

//The Data structure required by the LDSP algorithm
class DS{
    public:
        int num_V;
        map< vector<int>, vector< vector<int> > > L;
        map< vector<int>, vector< vector<int> > > R;
        map< vector<int>, vector< vector<int> > > Lstar;
        map< vector<int>, vector< vector<int> > > Rstar;
        map< vector<int>, double >  pathweight;
        //P and Pstar is priority queue here. Need to do make_heap after each access.
        vector< vector<int> > **P;
        vector< vector<int> > **Pstar;
        bool **extracted;
        DS(int v){
            num_V = v;
            P = new vector<vector<int>>*[v];
            for(int i = 0; i < v;i++)
            {
                P[i] = new vector<vector<int>> [v];
            }
            Pstar = new vector<vector<int>>*[v];
            for(int i = 0; i < v;i++)
            {
                Pstar[i] = new vector<vector<int>> [v];
            }
            extracted = new bool*[v];
            for(int i = 0; i < v;i++)
            {
                extracted[i] = new bool[v];
                memset(extracted[i],false,v*sizeof(bool));
            }
            for (int i = 0; i < v; i++)
            {
                vector<int> u;
                u.push_back(i); // path to itself
                pathweight[u] = 0;
                P[i][i].push_back(u);
                Pstar[i][i].push_back(u);
            }
        }
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
};

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

// Function for ISP proposed by DEMETRESCU and ITALIANO in 2004
vector<vector<double>> incremental_LDSP(vector<vector<double>> &graph, DS *state, int u, int v, double w) {
    //printf("add an edge with value %f from %d to %d\n",w,u,v);
    auto d = graph;
    struct pathwithweight{
    vector<int> path;
    double weight;
    };
    //This algorithm supports update all the edges incident to a vertex.
    //Here we only want to test it together with other algorithms, so only consider single update.
    graph[u][v] = w; 
    //Step 1: Clean up
    vector< vector<int> > Q;
    vector<int> path_u;
    path_u.push_back(u);
    Q.push_back(path_u);
    while(!Q.empty())
    {
        //cout << Q.size() << endl;
        vector<int> pi = Q.back();
        Q.pop_back();
        vector< vector<int> > Leftlist = state->L[pi];
        vector< vector<int> > Rightlist = state->R[pi];
        vector< vector<int> > LR = Leftlist;
        LR.insert(LR.end(), Rightlist.begin(), Rightlist.end());
        for (int i = 0;i < LR.size();i++)
        {
            vector<int> pi_xy = LR[i];
            Q.push_back(pi_xy);
            int x = pi_xy.front();
            int y = pi_xy.back();
            for(int i = 0; i < state->P[x][y].size(); i++)
            {
                if(state->P[x][y][i] == pi_xy)
                {
                    state->P[x][y].erase(state->P[x][y].begin()+i);
                    break;
                }
            }
            for(int i = 0; i < state->Pstar[x][y].size(); i++)
            {
                if(state->Pstar[x][y][i] == pi_xy)
                {
                    state->Pstar[x][y].erase(state->Pstar[x][y].begin()+i);
                    break;
                }
            }
            vector< vector<int> > Lrpi = state->L[state->right(pi_xy)];
            for(int i = 0;i < Lrpi.size();i++)
            {
                if (Lrpi[i] == pi_xy)
                {
                    Lrpi.erase(Lrpi.begin()+i);
                    state->L[state->right(pi_xy)] = Lrpi;
                    break;
                }
            }
            vector< vector<int> > Rlpi = state->R[state->left(pi_xy)];
            for(int i = 0;i < Rlpi.size();i++)
            {
                if (Rlpi[i] == pi_xy)
                {
                    Rlpi.erase(Rlpi.begin()+i);
                    state->R[state->left(pi_xy)] = Rlpi;
                    break;
                }
            }
            vector< vector<int> > Lstar_rpi = state->Lstar[state->right(pi_xy)];
            for(int i = 0;i < Lstar_rpi.size();i++)
            {
                if (Lstar_rpi[i] == pi_xy)
                {
                    Lstar_rpi.erase(Lstar_rpi.begin()+i);
                    state->Lstar[state->right(pi_xy)] = Lstar_rpi;
                    break;
                }
            }
            vector< vector<int> > Rstar_rpi = state->Rstar[state->left(pi_xy)];
            for(int i = 0;i < Rstar_rpi.size();i++)
            {
                if (Rstar_rpi[i] == pi_xy)
                {
                    Rstar_rpi.erase(Rstar_rpi.begin()+i);
                    state->Rstar[state->left(pi_xy)] = Lstar_rpi;
                    break;
                }
            }
        }
    }
    //cout << "cleanupend" << endl;
    //end of clean up
    //Step 2: fix up
    //This algorithm supports multiple updates. Here we only test single modification.
    //u->v
    graph[u][v] = w;
    //cout << "considered edges:" << endl;
    for(int i = 0; i < state->num_V; i++)
    {
        if (graph[u][i] < INF && i != u)
        {
            //cout << u << "," << i << endl;
            vector<int> path_ui;
            path_ui.push_back(u);
            path_ui.push_back(i);
            
            state->pathweight[path_ui] = graph[u][i];
            state->L[state->right(path_ui)].push_back(path_ui);
            state->R[state->left(path_ui)].push_back(path_ui);
            state->P[u][i].push_back(path_ui);
        }
        if(graph[i][u] < INF && i != u)
        {
            //cout << i << "," << u << endl;
            vector<int> path_iu;
            path_iu.push_back(i);
            path_iu.push_back(u);
                        
            state->pathweight[path_iu] = graph[i][u];
            state->L[state->right(path_iu)].push_back(path_iu);
            state->R[state->left(path_iu)].push_back(path_iu);
            state->P[i][u].push_back(path_iu);
        }
    } 
    struct Compareweight {
    bool operator()(pathwithweight const& p1, pathwithweight const& p2)
    {
        return p1.weight > p2.weight;
    }
    };
    //cout << "phase 1 end" << endl;
    //phase 2   
    priority_queue<pathwithweight, vector<pathwithweight> , Compareweight> H;
    
    //cout << "from P to H:" << endl;
    for(int i = 0;i<state->num_V;i++)
    {
        for(int j = 0;j< state->num_V; j++)
        {
            pathwithweight path;
            double minimum_weight = INF;
            for(int k = 0;k < state->P[i][j].size();k++)
            {
                if(state->pathweight[state->P[i][j][k]] < minimum_weight)
                {
                    path.weight = state->pathweight[state->P[i][j][k]];
                    path.path = state->P[i][j][k];
                    minimum_weight = path.weight;
                }
            }
            
            if(state->P[i][j].size() != 0 && path.path.size() > 1)
            {   
                /*
                for (int i = 0;i < path.path.size();i++)
                    {
                    cout << path.path[i] << " ";
                    }
                cout << endl;
                */
                H.push(path);
            }
        }
    }   
    //cout << "---------------"<<endl;
    //cout << "phase 2 end" << endl;
    //phase 3
    while(!H.empty())
    {
        //cout << H.size() << endl;
        pathwithweight p = H.top();
        H.pop();
        int x = p.path.front();
        int y = p.path.back();
        if(state->extracted[x][y]==1)continue;
        state->extracted[x][y] = 1;

        /*
        cout << "poping:" <<endl;
        for (int i = 0;i < p.path.size();i++)
        {
            cout << p.path[i] << " ";
        }
        */

        for(int i = 0; i < state->Pstar[x][y].size();i++)
        {
            if(state->Pstar[x][y][i] == p.path)
            {
                continue;
            }
        }
        state->Pstar[x][y].push_back(p.path);
        vector<vector<int> > modified_Lstarvalue = state->Lstar[state->right(p.path)];
        modified_Lstarvalue.push_back(p.path);
        state->Lstar[state->right(p.path)] = modified_Lstarvalue;

        vector<vector<int> > modified_Rstarvalue = state->Rstar[state->left(p.path)];
        modified_Rstarvalue.push_back(p.path);
        state->Rstar[state->left(p.path)] = modified_Rstarvalue;

        
        vector<vector<int> > Lstarl = state->Lstar[state->left(p.path)];
        //cout << "in Ll:" << endl;
        for(int i = 0;i < Lstarl.size();i++)
        {
            vector<int> pi_xp_b = Lstarl[i];
            int xp = pi_xp_b.front();
            bool iwantcontinue = false;
            for(int k = 0; k < p.path.size();k++)
            {
                if(xp == p.path[k])
                {
                    iwantcontinue = true;
                }
            }
            if(iwantcontinue)
            {
                continue;
            }
            /*
            for(int k = 0; k < Lstarl[i].size();k++)
            {
                cout << Lstarl[i][k] << " ";
            }
            cout << ",";
            */
            
            vector<int> pi_xp_y = p.path;
            pi_xp_y.insert(pi_xp_y.begin(),xp);
            state->pathweight[pi_xp_y] = state->pathweight[p.path] + graph[xp][x];
            state->P[xp][y].push_back(pi_xp_y);
            vector<vector<int> > temp = state->L[p.path];
            temp.push_back(pi_xp_y);
            state->L[p.path] = temp;

            temp = state->R[pi_xp_b];
            temp.push_back(pi_xp_y);
            state->R[pi_xp_b] = temp;
            pathwithweight p_w_xp_y;
            p_w_xp_y.path = pi_xp_y;
            p_w_xp_y.weight = state->pathweight[pi_xp_y];
            H.push(p_w_xp_y);

        }

        vector<int> rpath = state->right(p.path);
        vector<vector<int> > Rstarr = state->Rstar[state->right(p.path)];
        for(int i = 0;i < Rstarr.size();i++)
        {
            vector<int> pi_a_yp = Rstarr[i];
            int yp = pi_a_yp.back();
            bool iwantcontinue = false;
            for(int k = 0; k < p.path.size();k++)
            {
                if(yp == p.path[k])
                {
                    iwantcontinue = true;
                }
            }
            if(iwantcontinue)
            {
                continue;
            }
            vector<int> pi_x_yp = p.path;
            pi_x_yp.push_back(yp);
            state->pathweight[pi_x_yp] = state->pathweight[p.path] + graph[y][yp];
            state->P[x][yp].push_back(pi_x_yp);
            vector<vector<int> > temp = state->L[pi_a_yp];
            temp.push_back(pi_x_yp);
            state->L[pi_a_yp] = temp;

            temp = state->R[p.path];
            temp.push_back(pi_x_yp);
            state->R[p.path] = temp;
            pathwithweight p_w_x_yp;
            p_w_x_yp.path = pi_x_yp;
            p_w_x_yp.weight = state->pathweight[pi_x_yp];
            H.push(p_w_x_yp);
        }
    }
    for(int i = 0; i < state->num_V;i++)
    {
        memset(state->extracted[i],false,state->num_V);
    }
    //update distance 
    for(int i = 0;i < state->num_V; i++)
    {
        for(int j = 0;j < state->num_V; j++)
        {
            
            vector<vector<int> > shortestpath = state->Pstar[i][j];
            if(shortestpath.empty()) {
                d[i][j] = INF;
            }
            else {
                int imin=0;
                double min_weight = INF;
                for(int k = 0; k < shortestpath.size();k++)
                {
                    if (state->pathweight[shortestpath[k]] < min_weight)
                    {
                        imin = k;
                        min_weight = state->pathweight[shortestpath[k]];
                    }
                }
                d[i][j] = state->pathweight[shortestpath[imin]];
            }
        }           
    }
    //print the distance 
    /*
    printf("distance:\n");
    for (int i = 0; i < state->num_V; i++) 
    {
        for (int j = 0; j < state->num_V; j++) 
        {
            double val = d[i][j];
            
            if (val  >= 900000)
                cout << "INF   ";
            else
                cout << val << "   ";
            
        }
        cout << "\n";
    }
    */
    return d;
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
    int V = 100;  // Number of vertices
    int E = 80;  // Number of edges

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

    // 
    // Measure (my) LDSP algorithm
    DS state(V);
    //catch up the current graph
    vector< vector<double> > d;
    for(int i = 0;i<V;i++)
    {
        for(int j=0;j<V;j++)
        {
            if (i==j) continue;
            if(graph[i][j] < INF)
            { 
                    d = incremental_LDSP(graph, &state, i, j, graph[i][j]);
                    break;
            }
        }
    }
    auto start_LDSP = high_resolution_clock::now();
    d = incremental_LDSP(graph, &state, u, v, w);
    auto end_LDSP= high_resolution_clock::now();
    auto duration_LDSP = duration_cast<microseconds>(end_LDSP - start_LDSP);
    cout << "Time taken by LDSP algorithm: " << duration_LDSP.count() << " microseconds\n";
    return 0;
}
