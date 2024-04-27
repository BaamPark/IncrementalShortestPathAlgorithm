import random
import sys
from collections import deque
import time
import numpy as np
import copy

INF = float('inf')
class path:

    def __init__(self,*vertices):
        self.path=[]
        self.weight = 0
        for v in vertices:
            self.path.append(v)
    def right(self):
        return tuple(self.path[1:])
    def left(self):
        return tuple(self.path[:-1])
    def tupleself(self):
        return tuple(self.path)
class pathqueue:
    def __init__(self,*paths):
        self.set = []
        for p in paths:
            self.set.append(p)
    def sort(self):
        self.set = sorted(self.set, key=lambda x: x.weight)
    def top(self):
        self.set = sorted(self.set, key=lambda x: x.weight)
        return self.set[0]
    def pop(self):
        self.set = sorted(self.set, key=lambda x: x.weight)
        smallest = self.set[0]
        self.set.pop(0)
        return smallest
    def add(self,path):
        self.set.append(path)
        return self
    def delete(self,path):
        for pp in self.set:
            if path.path == pp.path and path.weight == pp.weight:
                self.set.remove(pp)
                return self
        return self
    def find(self,path):
        for pp in self.set:
            if path.path == pp.path and path.weight == pp.weight:
                return True
        return False
class State:
    def __init__(self,V):
        self.V = V
        self.P = [[pathqueue() for _ in range(V)] for _ in range(V)]
        self.Pstar = [[pathqueue() for _ in range(V)] for _ in range(V)]
        self.weight = dict([])
        self.L = dict([])
        self.R = dict([])
        self.Lstar = dict([])
        self.Rstar = dict([])
def LDSP_init(graph,V):
    state = State(V)
    for u in range(V):
        state,distance = incremental_apsp_LDSP(state, graph, u)
    return state,distance

def main():
    V = 10  # Number of vertices
    E = 5 # Number of edges


    graph = generate_random_graph(V, E)

    u, v, w = pick_random_edge_for_update(graph)
    #graph = graphdebug()
    #u,v,w = edge_degug()
    print(f"Updating edge ({u}, {v}) from weight {graph[u][v]} to {w}")
    # Run Floyd-Warshall algorithm and return the distance matrix
    original_distance = floyd_warshall(graph)
    state,LDSP_distance = LDSP_init(graph.copy(),V)
    # Update graph for subsequent tests
    graph[u][v] = w

    # Get updated distance matrices from different algorithms
    naive_distance = floyd_warshall(graph)
    pr_distance = incremental_apsp_pr(original_distance.copy(), u, v, w)
    quinca_distance = incremental_apsp_quinca(original_distance.copy(), u, v, w)
    state,LDSP_distance = incremental_apsp_LDSP(state,graph.copy(),u)
    # Compare the results
    are_same_naive_pr = np.array_equal(naive_distance, pr_distance)
    are_same_naive_quinca = np.array_equal(naive_distance, quinca_distance)
    are_same_naive_LDSP = np.array_equal(naive_distance, LDSP_distance)
    are_same_pr_quinca = np.array_equal(pr_distance, quinca_distance)

    print(f"Naive vs PR results are the same: {are_same_naive_pr}")
    print(f"Naive vs QUINCA results are the same: {are_same_naive_quinca}")
    print(f"Naive vs LDSP results are the same: {are_same_naive_LDSP}")
    print(f"PR vs QUINCA results are the same: {are_same_pr_quinca}")


def generate_random_graph(V, E):
    # Initialize the graph with INF indicating no direct path between nodes
    graph = [[INF] * V for _ in range(V)]

    for _ in range(E):
        u = random.randint(0, V - 1)
        v = random.randint(0, V - 1)
        w = random.randint(1, 100)  # Edge weight between 1 and 100

        # Ensure we don't set a node to itself and only set the directed edge from u to v
        if u != v:
            graph[u][v] = w

    # Diagonal should be 0, indicating no self-loops
    for i in range(V):
        graph[i][i] = 0

    return graph
def pick_random_edge_for_update(graph):
    V = len(graph)
    edges = []

    # Collect all edges with a set weight
    for u in range(V):
        for v in range(V):
            if graph[u][v] != float('inf') and graph[u][v] != 0 and u != v:
                edges.append((u, v, graph[u][v]))

    # Randomly select one edge and suggest a new weight less than the current weight
    if edges:
        index = random.randint(0, len(edges) - 1)
        u, v, current_weight = edges[index]
        new_weight = random.randint(1, int(current_weight))  # Ensure new weight is less and > 0
        return u, v, new_weight
    else:
        raise ValueError("No valid edges found in the graph for updating.")


def floyd_warshall(graph):
    V = len(graph)
    # Initialize the distance matrix with the values from the graph
    distance = [row[:] for row in graph]  # Make a deep copy of the graph

    # Adjust the distance matrix for initial direct paths
    for i in range(V):
        for j in range(V):
            if i != j and graph[i][j] == 0:
                distance[i][j] = float('inf')

    # Floyd-Warshall algorithm
    for k in range(V):
        for i in range(V):
            for j in range(V):
                if distance[i][k] < float('inf') and distance[k][j] < float('inf'):
                    distance[i][j] = min(distance[i][j], distance[i][k] + distance[k][j])

    return distance


def update_and_floyd_warshall(graph, u, v, w):
    # Update the graph with the new weight for the edge between u and v
    graph[u][v] = w

    # Run the Floyd-Warshall algorithm on the updated graph
    return floyd_warshall(graph)


def find_affected_source(distance, u, v, w):
    n = len(distance)
    vis = [False] * n
    affected_sources = []

    if distance[u][v] > w:
        Q = deque([u])
        affected_sources.append(u)
        vis[u] = True

        while Q:
            x = Q.popleft()

            for z in range(n):
                if not vis[z] and distance[z][v] > distance[z][u] + w:
                    Q.append(z)
                    vis[z] = True
                    affected_sources.append(z)

        # No need to reset vis to False in Python, as we don't reuse this list
    return affected_sources


def incremental_apsp_pr(distance, u, v, w):
    sources = find_affected_source(distance, u, v, w)
    V = len(distance)
    for s in sources:
        vis = [False] * V
        Q = deque()

        # Update the distance from source to v if the new path is shorter
        if distance[s][v] > distance[s][u] + w:
            distance[s][v] = distance[s][u] + w
            Q.append(v)
            vis[v] = True

            # Truncated BFS (Breadth-First Search)
            while Q:
                y = Q.popleft()
                for w_index in range(V):
                    if not vis[w_index] and distance[s][w_index] > distance[s][u] + w + distance[v][w_index]:
                        distance[s][w_index] = distance[s][u] + w + distance[v][w_index]
                        vis[w_index] = True
                        Q.append(w_index)

    return distance


def incremental_apsp_quinca(distances, u, v, w_prime):
    n = len(distances)
    affected_sources = find_affected_source(distances, u, v, w_prime)

    if distances[u][v] > w_prime:
        Q = deque()
        visited = [False] * n
        distances[u][v] = w_prime  # Update the distance for the edge (u, v)

        Q.append(v)
        visited[v] = True

        while Q:
            y = Q.popleft()

            # Update distances for affected sources using the new path through (u, v)
            for x in affected_sources:
                if distances[x][y] > distances[x][u] + w_prime + distances[v][y]:
                    distances[x][y] = distances[x][u] + w_prime + distances[v][y]

            # Check for further vertices that could use the updated path for shorter distances
            for w in range(n):
                if not visited[w] and distances[u][w] > w_prime + distances[v][w]:
                    Q.append(w)
                    visited[w] = True

    return distances  # Return the updated distance matrix

def incremental_apsp_LDSP(state,graph,u):
    LDSP_distance = copy.deepcopy(graph)
    #step1: cleanup
    pi = path(u)
    Q = pathqueue()
    Q.add(pi)
    while(len(Q.set)):
        pixy = Q.pop()

        if pixy.tupleself() in state.L:
            for p in state.L[pixy.tupleself()].set:
                Q.add(p)
                state.P[p.path[0]][p.path[-1]].delete(p)
                if (p.right() in state.L):
                    state.L[p.right()] = state.L[p.right()].delete(p)
                if (p.left() in state.R):
                    state.R[p.left()] = state.R[p.left()].delete(p)
                if state.Pstar[p.path[0]][p.path[-1]].find(p):
                    state.Pstar[p.path[0]][p.path[-1]].delete(p)
                    if (p.right() in state.Lstar):
                        state.Lstar[p.right()] = state.Lstar[p.right()].delete(p)
                    if (p.left() in state.Rstar):
                        state.Rstar[p.left()] = state.Rstar[p.left()].delete(p)
        if pixy.tupleself() in state.R:
            for p in state.R[pixy.tupleself()].set:
                Q.add(p)
                state.P[p.path[0]][p.path[-1]].delete(p)
                if (p.right() in state.L):
                    state.L[p.right()] = state.L[p.right()].delete(p)
                if (p.left() in state.R):
                    state.R[p.left()] = state.R[p.left()].delete(p)
                if state.Pstar[p.path[0]][p.path[-1]].find(p):
                    state.Pstar[p.path[0]][p.path[-1]].delete(p)
                    if (p.right() in state.Lstar):
                        state.Lstar[p.right()] = state.Lstar[p.right()].delete(p)
                    if (p.left() in state.Rstar):
                        state.Rstar[p.left()] = state.Rstar[p.left()].delete(p)
    #step2: fixup
    #phase1: check all the edges indicate to u
    for v in range(state.V):
        if u != v:
            if graph[u][v] < INF:
                p_uv = path(u,v)
                state.weight[p_uv] = graph[u][v]
                p_uv.weight = graph[u][v]
                state.P[u][v].add(p_uv)
                if p_uv.right() in state.L:
                    state.L[p_uv.right()] = state.L[p_uv.right()].add(p_uv)
                else:
                    state.L[p_uv.right()] = pathqueue(p_uv)
                if (p_uv.left() in state.R):
                    state.R[p_uv.left()] = state.R[p_uv.left()].add(p_uv)
                else:
                    state.R[p_uv.left()] = pathqueue(p_uv)
            if graph[v][u] < INF:
                p_vu = path(v,u)
                state.weight[p_vu] = graph[v][u]
                p_vu.weight = graph[v][u]
                state.P[v][u].add(p_vu)

                if (p_vu.right() in state.L):
                    state.L[p_vu.right()] = state.L[p_vu.right()].add(p_vu)
                else:
                    state.L[p_vu.right()] = pathqueue(p_vu)
                if (p_vu.left() in state.R):
                    state.R[p_vu.left()] = state.R[p_vu.left()].add(p_vu)
                else:
                    state.R[p_vu.left()] = pathqueue(p_vu)
    #phase2: add path to H.
    H = pathqueue()
    for x in range(state.V):
        for y in range(state.V):
            if len(state.P[x][y].set)!= 0:
                H.add(state.P[x][y].top())
    tpath = path(2,0)
    #phase3: process H

    extracted =  [[0 for _ in range(state.V)] for _ in range(state.V)]
    while(len(H.set)!=0):
        pixy = H.pop()
        x = pixy.path[0]
        y = pixy.path[-1]
        if(extracted[x][y]==0):
            extracted[x][y] = 1
            if not (state.Pstar[x][y].find(pixy)):
                state.Pstar[x][y].add(pixy)
                if (pixy.right() in state.Lstar):
                    state.Lstar[pixy.right()] = state.Lstar[pixy.right()].add(pixy)
                else:
                    state.Lstar[pixy.right()] = pathqueue(pixy)
                if (pixy.left() in state.Rstar):
                    state.Rstar[pixy.left()] = state.Rstar[pixy.left()].add(pixy)
                else:
                    state.Rstar[pixy.left()] = pathqueue(pixy)


                if(pixy.left() in state.Lstar):
                    for p_xp_b in state.Lstar[pixy.left()].set:

                        xp = p_xp_b.path[0]
                        if xp in pixy.path:
                            continue
                        p_xp_y = copy.deepcopy(p_xp_b)
                        p_xp_y.path.append(y)
                        state.weight[p_xp_y] = state.weight[pixy] + graph[xp][x]
                        p_xp_y.weight = state.weight[pixy] + graph[xp][x]
                        state.P[xp][y].add(p_xp_y)
                        H.add(p_xp_y)
                        if (pixy.tupleself() in state.L):
                            state.L[pixy.tupleself()] = state.L[pixy.tupleself()].add(p_xp_y)
                        else:
                            state.L[pixy.tupleself()] = pathqueue(p_xp_y)
                        if (p_xp_b.tupleself() in state.R):
                            state.R[p_xp_b.tupleself()] = state.R[p_xp_b.tupleself()].add(p_xp_y)
                        else:
                            state.R[p_xp_b.tupleself()] = pathqueue(p_xp_y)
                if (pixy.right() in state.Rstar):
                    for p_a_yp in state.Rstar[pixy.right()].set:
                        yp = p_a_yp.path[-1]
                        if yp in pixy.path:
                            continue
                        p_x_yp = copy.deepcopy(pixy)
                        p_x_yp.path.append(yp)

                        state.weight[p_x_yp] = state.weight[pixy] + graph[y][yp]
                        p_x_yp.weight = state.weight[pixy] + graph[y][yp]
                        state.P[x][yp] = state.P[x][yp].add(p_x_yp)

                        H.add(p_x_yp)
                        if (p_a_yp.tupleself() in state.L):
                            state.L[p_a_yp.tupleself()] = state.L[p_a_yp.tupleself()].add(p_x_yp)
                        else:
                            state.L[p_a_yp.tupleself()] = pathqueue(p_x_yp)
                        if (pixy.tupleself() in state.R):
                            state.R[pixy.tupleself()] = state.R[pixy.tupleself()].add(p_x_yp)
                        else:
                            state.R[pixy.tupleself()] = pathqueue(p_x_yp)
    extracted = [[0 for _ in range(state.V)] for _ in range(state.V)]
    for u in range(state.V):
        for v in range(state.V):
            if len(state.Pstar[u][v].set)!=0:
                LDSP_distance[u][v] = state.Pstar[u][v].top().weight
    return state,LDSP_distance
if __name__ == "__main__":
    main()
