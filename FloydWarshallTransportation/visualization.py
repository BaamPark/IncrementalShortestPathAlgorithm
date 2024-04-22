import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data from CSV files
graph = pd.read_csv('graph.csv', header=None).values
distance = pd.read_csv('distance.csv', header=None).values

# Plotting the original graph
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.imshow(graph, cmap='viridis', interpolation='nearest')
plt.title('Original Graph')
plt.colorbar(label='Edge Weight')

# Plotting the computed distances
plt.subplot(1, 2, 2)
plt.imshow(distance, cmap='viridis', interpolation='nearest')
plt.title('Computed Shortest Distances')
plt.colorbar(label='Distance')

plt.tight_layout()

plt.savefig('visualizatio_plot.png')
plt.show()