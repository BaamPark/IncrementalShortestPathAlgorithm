from funcs import generate_random_graph, pick_random_edge_for_update, floyd_warshall, update_and_floyd_warshall, incremental_apsp_pr, incremental_apsp_quinca
import time
import matplotlib.pyplot as plt

def measure_performance(V, E):
    graph = generate_random_graph(V, E)
    u, v, w = pick_random_edge_for_update(graph)

    # Measure naive Floyd-Warshall algorithm
    start_naive = time.perf_counter()
    update_and_floyd_warshall(graph, u, v, w)
    end_naive = time.perf_counter()
    duration_naive = (end_naive - start_naive) * 1e6  # Convert to microseconds

    # Measure PR algorithm
    distance_for_PR = floyd_warshall(graph)
    start_pr = time.perf_counter()
    incremental_apsp_pr(distance_for_PR, u, v, w)
    end_pr = time.perf_counter()
    duration_pr = (end_pr - start_pr) * 1e6  # Convert to microseconds

    # Measure QUINCA algorithm
    distance_for_QUINCA = floyd_warshall(graph)
    start_quinca = time.perf_counter()
    incremental_apsp_quinca(distance_for_QUINCA, u, v, w)
    end_quinca = time.perf_counter()
    duration_quinca = (end_quinca - start_quinca) * 1e6  # Convert to microseconds

    return duration_naive, duration_pr, duration_quinca

def main():
    V_range = range(10, 101, 5)  # Adjust range and step as needed
    E = 10  # Example number of edges, adjust based on the graph density you want
    naive_times, pr_times, quinca_times = [], [], []

    for V in V_range:
        t_naive, t_pr, t_quinca = measure_performance(V, E)
        naive_times.append(t_naive)
        pr_times.append(t_pr)
        quinca_times.append(t_quinca)

    plt.figure(figsize=(10, 5))
    plt.plot(V_range, naive_times, label='Naive Floyd-Warshall')
    plt.plot(V_range, pr_times, label='PR Algorithm')
    plt.plot(V_range, quinca_times, label='QUINCA Algorithm')
    plt.xlabel('Number of Vertices V')
    plt.ylabel('Time (microseconds)')
    plt.title('Performance Comparison of APSP Algorithms')
    plt.legend()
    plt.grid(True)
    plt.yscale('log')  # Set y-axis to logarithmic scale
    plt.savefig('performance_comparison.png')
    plt.show()
if __name__ == "__main__":
    main()
