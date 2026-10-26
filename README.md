# bio_toolkit

Custom implementations of common bioinformatics algorithms.

## Introduction

`bio_toolkit` is a small, dependency-light collection of sequence analysis
algorithms written from scratch in Python. The goal is to make the mechanics
visible: each algorithm is implemented directly rather than delegated to
Biopython or scikit-bio, so the dynamic programming tables, distance updates
and tree-building steps can be read and modified.


### Modules

| File | Contents |
| --- | --- |
| `utils.py` | `read_fasta` — FASTA parser handling wrapped lines and legacy `;` comments |
| `alignment.py` | Needleman–Wunsch scoring and all-pairs distance matrix construction |
| `clustering.py` | Agglomerative tree building and dendrogram plotting |

### Requirements

Python 3 with `numpy`, `scipy` and `matplotlib`.

### Example

`covid.fasta` contains seven SARS-CoV-2 spike protein sequences sampled
between January 2020 and February 2024. Running the clustering module end to
end builds a distance matrix from those sequences and plots the tree:

```bash
python clustering.py
```

## Alignment algorithms

### `needleman_wunsch_score(seq1, seq2, mismatch_penalty=3, indel_penalty=5)`

Global alignment by dynamic programming. The scoring table has an extra row
and column representing the empty prefix of each sequence, initialised to
multiples of the indel penalty. Every interior cell takes the cheapest of
three moves — an insertion, a deletion, or a substitution — and the bottom
right cell holds the optimal alignment cost.

This is a cost-minimising formulation: identical characters cost nothing,
mismatches cost `mismatch_penalty` and gaps cost `indel_penalty` per position,
so **lower scores mean more similar sequences**. That makes the score usable
directly as a distance. Runtime and memory are both O(mn).

Only the score is returned; the traceback that would recover the aligned
sequences themselves is not implemented, and neither is affine gap scoring
(a gap of length *k* costs `k * indel_penalty`). The `Alignment` class in
`alignment.py` is a placeholder for that work.

```python
from alignment import needleman_wunsch_score

needleman_wunsch_score('AATGCTCAAAA', 'AATCCTCGAAA')  # -> 6.0
```

### `generate_matrix(seqs, distance_function=needleman_wunsch_score)`

Applies a pairwise distance function across every pair in `seqs` and returns
an n x n `numpy` matrix. The diagonal is set to `np.inf` rather than `0` so
that self-comparisons are never selected as a nearest pair downstream. Any
function with the signature `(seq1, seq2) -> float` can be substituted.

Note that all n^2 pairs are computed, including both orders of each pair, so
the matrix costs twice what symmetry requires.

## Clustering algorithms

### `neighbor_joining(distance_matrix)`

Builds a tree from a distance matrix by repeated agglomeration. Each round
finds the smallest entry in the current matrix, merges that pair of clusters,
and recomputes the distances from the *original* matrix — recomputing from
source rather than from the collapsed matrix avoids accumulating rounding
drift across rounds.

Returns a linkage matrix in SciPy's format: one row per merge, holding
`[cluster_id_1, cluster_id_2, distance, size]`, where ids `0..n-1` are the
input sequences and each merge creates a new id. That matrix can be passed
straight to `scipy.cluster.hierarchy.dendrogram`.

```python
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt

tree = neighbor_joining(matrix)
dendrogram(tree)
plt.show()
```

Supporting functions:

- `cluster_distances(distance_matrix, cluster1, cluster2)` — average linkage:
  the mean of all `|cluster1| * |cluster2|` pairwise distances between the two
  clusters.
- `merge_clusters(clusters, node_index1, node_index2)` — returns a new cluster
  list with the two named clusters replaced by their union, appended last.

