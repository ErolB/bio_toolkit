import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from matplotlib import pyplot as plt
import copy
import pickle as pkl
import os
import warnings

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

def neighbor_joining(distance_matrix, algorithm='upgma'):
    distances = np.array(copy.copy(distance_matrix), dtype=float)
    np.fill_diagonal(distances, 0.0)  # row sums below must stay finite
    current_clusters = [[i] for i in range(len(distance_matrix[0]))]
    all_clusters = copy.copy(current_clusters)
    tree_matrix = []  # used to build dendrogram
    while len(current_clusters) > 1:
        print('start of loop')
        n = len(current_clusters)
        r = distances.sum(axis=1)
        q = (n - 2) * distances - r[:, None] - r[None, :]
        np.fill_diagonal(q, np.inf)
        np.fill_diagonal(distances, np.inf)
        # both matrices are symmetric with an inf diagonal, so a plain argmin
        # lands on an off-diagonal pair (np.triu would zero the lower triangle,
        # and those zeros win the argmin for non-negative distances)
        if algorithm == 'nj':
            n1, n2 = np.unravel_index(np.argmin(q), q.shape)
        elif algorithm == 'upgma':
            n1, n2 = np.unravel_index(np.argmin(distances), distances.shape)
        neighbor1 = current_clusters[n1]
        neighbor2 = current_clusters[n2]
        print(current_clusters)
        # merge neighbors
        current_clusters = merge_clusters(current_clusters, n1, n2)
        all_clusters.append(current_clusters[-1])
        true_index_n1 = all_clusters.index(neighbor1)
        true_index_n2 = all_clusters.index(neighbor2)
        d = cluster_distances(distance_matrix, neighbor1, neighbor2)
        tree_matrix.append([true_index_n1, true_index_n2, d, len(current_clusters[-1])])
        # recalculate distance matrix
        new_distances = []
        for i, c1 in enumerate(current_clusters):
            new_row = []
            for j, c2 in enumerate(current_clusters):
                if i == j:
                    new_row.append(0.0)
                else:
                    new_row.append(cluster_distances(distance_matrix, c1, c2))
            new_distances.append(new_row)
        distances = np.array(new_distances)
    return(tree_matrix)


def to_phylo(tree_matrix, labels=None):
    # convert the linkage matrix into a Bio.Phylo tree; a child's branch length
    # is how far it sits below its parent, so leaves end up at height zero
    leaf_count = len(tree_matrix) + 1
    if labels is None:
        labels = [str(i) for i in range(leaf_count)]
    clades = {i: Clade(name=labels[i]) for i in range(leaf_count)}
    heights = {i: 0.0 for i in range(leaf_count)}
    clamped = 0
    for step, row in enumerate(tree_matrix):
        height = float(row[2])
        children = []
        for child_index in (int(row[0]), int(row[1])):
            child = clades.pop(child_index)
            length = height - heights[child_index]
            if length < 0:  # non-monotonic merge heights, see README
                clamped += 1
                length = 0.0
            child.branch_length = length
            children.append(child)
        node_index = leaf_count + step
        clades[node_index] = Clade(clades=children)
        heights[node_index] = height
    if clamped:
        warnings.warn(f'{clamped} branch length(s) clamped to zero: the merge '
                      'heights are not monotonic, so the drawn tree understates '
                      'those branches')
    return Tree(root=clades[leaf_count + len(tree_matrix) - 1], rooted=True)


if __name__ == '__main__':
    records = read_fasta('covid.fasta')
    sequences = list(records.values())
    #sequences = [item[:100] for item in sequences]
    score, space = needleman_wunsch_score(sequences[0], sequences[5])
    print(space.shape)
    print(len(sequences[0]), len(sequences[5]))
    align1, align2 = backtrack(sequences[0], sequences[5], space)
    print(align1[:100])
    print(''.join(['|' if align1[i]!=align2[i] and align1[i]!='-' and align2[i]!='-' else ' ' for i in range(100)]))
    print(align2[:100])
    '''if os.path.exists('C:\\Users\\Owner\\Documents\\bio_toolkit\\nw_matrix.pkl'):
        print('loading matrix')
        with open('nw_matrix.pkl', 'rb') as matrix_file:
            matrix = pkl.load(matrix_file)
    else:
        print('regenerating matrix')
        matrix = generate_matrix(sequences)
        with open('nw_matrix.pkl', 'wb') as matrix_file:
            pkl.dump(matrix, matrix_file)
    print(matrix)
    tree = to_phylo(neighbor_joining(matrix, algorithm='nj'), labels=list(records))
    Phylo.draw_ascii(tree)
    Phylo.draw(tree, do_show=False)
    plt.show()'''

