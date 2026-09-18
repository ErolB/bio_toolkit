import numpy as np
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
import copy

from alignment import *
from utils import *

def merge_clusters(clusters, node_index1, node_index2):
    print('merge clusters')
    new_clusters = []
    for i, c in enumerate(clusters):
        if i in [node_index1, node_index2]:
            continue
        new_clusters.append(c)
    new_clusters.append(clusters[node_index1]+clusters[node_index2])
    return new_clusters

def cluster_distances(distance_matrix, cluster1, cluster2):
    distances = 0
    for i, node1 in enumerate(cluster1):
        for j, node2 in enumerate(cluster2):
            distances += distance_matrix[node1][node2]
    return distances / (len(cluster1)*len(cluster2))

def neighbor_joining(distance_matrix):
    distances = copy.copy(distance_matrix)
    current_clusters = [[i] for i in range(len(distance_matrix[0]))]
    all_clusters = copy.copy(current_clusters)
    tree_matrix = []  # used to build dendrogram
    while len(current_clusters) > 1:
        print('start of loop')
        print(current_clusters)
        smallest_distance = np.min(np.array(distances))
        closest = np.where(np.array(distances)==smallest_distance)
        neighbor1 = closest[0][0]
        neighbor2 = closest[1][0]
        # merge neighbors
        print(distances)
        d = cluster_distances(distance_matrix, current_clusters[neighbor1], current_clusters[neighbor2])
        true_index_n1 = all_clusters.index(current_clusters[neighbor1])
        true_index_n2 = all_clusters.index(current_clusters[neighbor2])
        tree_matrix.append([true_index_n1, true_index_n2, d, len(current_clusters[-1])])
        # recalculate distance matrix
        current_clusters = merge_clusters(current_clusters, neighbor1, neighbor2)
        all_clusters.append(current_clusters[-1])
        new_distances = []
        for i, c1 in enumerate(current_clusters):
            new_row = []
            for j, c2 in enumerate(current_clusters):
                if i == j:
                    new_row.append(10000)
                else:
                    new_row.append(cluster_distances(distance_matrix, c1, c2))
            new_distances.append(new_row)
        distances = copy.copy(new_distances)
    return(tree_matrix)
            

if __name__ == '__main__':
    sequences = list(read_fasta('covid.fasta').values())
    print(needleman_wunsch_score(sequences[0], sequences[1]))
    matrix = generate_matrix(sequences)
    tree = neighbor_joining(matrix)
    dendrogram(tree)
    plt.show()

