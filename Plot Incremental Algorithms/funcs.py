import random
import sys
from collections import deque
import time
import numpy as np

INF = float('inf')


def main():
    V = 10  # Number of vertices
    E = 5  # Number of edges

    graph = generate_random_graph(V, E)
    u, v, w = pick_random_edge_for_update(graph)
    print(f"Updating edge ({u}, {v}) from weight {graph[u][v]} to {w}")

    # Run Floyd-Warshall algorithm and return the distance matrix
    original_distance = floyd_warshall(graph)
    
    # Update graph for subsequent tests
    graph[u][v] = w

    # Get updated distance matrices from different algorithms
    naive_distance = floyd_warshall(graph)
    pr_distance = incremental_apsp_pr(original_distance.copy(), u, v, w)
    quinca_distance = incremental_apsp_quinca(original_distance.copy(), u, v, w)

    # Compare the results
    are_same_naive_pr = np.array_equal(naive_distance, pr_distance)
    are_same_naive_quinca = np.array_equal(naive_distance, quinca_distance)
    are_same_pr_quinca = np.array_equal(pr_distance, quinca_distance)

    print(f"Naive vs PR results are the same: {are_same_naive_pr}")
    print(f"Naive vs QUINCA results are the same: {are_same_naive_quinca}")
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

if __name__ == "__main__":
    main()