# bio_toolkit

Custom implementations of common bioinformatics algorithms.

## Introduction

`bio_toolkit` is a small, dependency-light collection of sequence analysis
algorithms written from scratch in Python. The goal is to make the mechanics
visible: each algorithm is implemented directly rather than delegated to
Biopython or scikit-bio, so the dynamic programming tables, index structures,
distance updates and graph traversals can be read and modified. Biopython is
used only to represent and draw the finished trees.


### Modules

| File | Contents |
| --- | --- |
| `utils.py` | `read_fasta` — FASTA parser handling wrapped lines and legacy `;` comments |
| `alignment.py` | Needleman–Wunsch scoring and traceback, all-pairs distance matrices, and an FM index for exact search |
| `clustering.py` | UPGMA / neighbor-joining tree building and conversion to a `Bio.Phylo` tree |
| `assembly.py` | k-mer composition, de Bruijn graph construction and Eulerian path reconstruction |
| `testing.py` | End-to-end check that a sequence can be rebuilt from its own k-mers |

### Requirements

Python 3 with `numpy`, `pandas`, `matplotlib` and `biopython`.

### Data

- `covid.fasta` — seven SARS-CoV-2 spike **protein** sequences sampled between
  January 2020 and February 2024.
- `sars_cov_2.fna` — the SARS-CoV-2 reference genome (NC_045512.2,
  Wuhan-Hu-1) as nucleotides.
- `nw_matrix.pkl` — a cached Needleman–Wunsch distance matrix for
  `covid.fasta`, so the tree can be rebuilt without realigning every pair.

### Example

Each module has a `__main__` block that exercises it on the bundled data:

```bash
python alignment.py   # FM index search for a spike peptide
python clustering.py  # align two spike sequences and print the alignment
python assembly.py    # rebuild a short test sequence from its 5-mers
python testing.py     # rebuild a spike sequence from its 7-mers and compare
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
directly as a distance. Runtime and memory are both O(mn). Affine gap scoring
is not implemented (a gap of length *k* costs `k * indel_penalty`).

Returns a `(score, table)` tuple; the table is what `backtrack` needs to
recover the alignment itself.

```python
from alignment import needleman_wunsch_score

score, table = needleman_wunsch_score('AATGCTCAAAA', 'AATCCTCGAAA')  # score -> 6.0
```

### `backtrack(seq1, seq2, alignment_space)`

Walks the table from the bottom right corner back to the top left and returns
the two aligned strings, with `-` marking gaps.

```python
backtrack('AATGCTCAAAA', 'AATCCTCGAAA', table)
# -> ('AATGCTCAAAA', 'AATCCTCGAAA')
```

Known limitations:

- At each step it moves to whichever neighbouring cell is smallest, rather than
  the cell that actually produced the current value (neighbour plus the move's
  penalty). On sequences of unequal length this can produce a non-optimal
  alignment: `backtrack('GAAC', 'AC', ...)` gives `A-C-` instead of `--AC`.
- The walk stops as soon as either index reaches zero, so any leading
  characters still left in the other sequence are not emitted as gaps.

### `generate_matrix(seqs, distance_function=needleman_wunsch_score)`

Applies a pairwise distance function across every pair in `seqs` and returns
an n x n `numpy` matrix with a zero diagonal. The diagonal must be finite
because neighbor joining sums whole rows. Any function with the signature
`(seq1, seq2) -> float` can be substituted. All n^2 pairs are computed,
including both orders of each pair, so the matrix costs twice what symmetry
requires.

Because `needleman_wunsch_score` now returns a tuple, the default
`distance_function` raises `ValueError: setting an array element with a
sequence`. Until the default is updated, pass a wrapper that keeps only the
score:

```python
from alignment import generate_matrix, needleman_wunsch_score

matrix = generate_matrix(seqs, lambda a, b: needleman_wunsch_score(a, b)[0])
```

### `FMIndex(reference_sequence)`

Builds a Burrows–Wheeler-based index of the reference: it appends a `$`
terminator, sorts every suffix, and keeps the first (`F`) and last (`L`)
columns of the BWT along with the suffix array and its inverse.

`search(query)` returns the 0-based start position of every exact match of
`query`. It starts at each row whose suffix begins with the query's last
character and extends the match one character to the left at a time, checking
the `L` column and stepping to the preceding suffix through the suffix array.

```python
from alignment import FMIndex

FMIndex('BANANABANANA').search('ANA')  # -> [1, 3, 7, 9] in some order
```

Limitations:

- Construction materialises every suffix as a string in a `pandas` frame, so
  memory is O(n^2). That is fine for a protein but not for a whole genome.
- Search uses the suffix array rather than occurrence counts (LF mapping), so
  it tries each candidate row one by one. That is correct, but it doesn't have
  the O(|query|) speed of a real FM index.
- `mismatch_limit` and `indel_limit` are accepted but ignored: only exact
  matches are found.
- The constructor prints the whole table and its build time.

## Clustering algorithms

### `neighbor_joining(distance_matrix, algorithm='upgma')`

Builds a tree from a distance matrix by repeated agglomeration. Each round
picks a pair of clusters to merge, records the merge, and recomputes the
cluster distances from the *original* matrix. Recomputing from source rather
than from the collapsed matrix avoids accumulating rounding drift across
rounds.

The `algorithm` argument controls only how the pair is chosen:

- `'upgma'` (default) merges the pair with the smallest average-linkage
  distance.
- `'nj'` merges the pair that minimises the neighbor-joining Q matrix,
  `Q(i, j) = (n - 2) d(i, j) - r(i) - r(j)`, where `r` is the row sum.

In both modes the recorded merge height is the average-linkage distance
between the two clusters. These are not true neighbor-joining branch
lengths.

Returns a linkage matrix in SciPy's format: one row per merge, holding
`[cluster_id_1, cluster_id_2, distance, size]`, where ids `0..n-1` are the
input sequences and each merge creates a new id.

Supporting functions:

- `cluster_distances(distance_matrix, cluster1, cluster2)` — average linkage:
  the mean of all `|cluster1| * |cluster2|` pairwise distances between the two
  clusters.
- `merge_clusters(clusters, node_index1, node_index2)` — returns a new cluster
  list with the two named clusters replaced by their union, appended last.

### `to_phylo(tree_matrix, labels=None)`

Converts the linkage matrix into a rooted `Bio.Phylo` tree. Each branch's
length is the parent's merge height minus the child's, so the leaves sit at
height zero. `labels` names the leaves in input order; it defaults to their
indices.

```python
from Bio import Phylo
from matplotlib import pyplot as plt

tree = to_phylo(neighbor_joining(matrix, algorithm='nj'), labels=list(records))
Phylo.draw_ascii(tree)
Phylo.draw(tree, do_show=False)
plt.show()
```

**Non-monotonic merge heights.** UPGMA merge heights never decrease, but with
`algorithm='nj'` a later merge can sit *lower* than an earlier one. The Q
matrix can choose a pair whose average-linkage distance is smaller than a
merge already made. That would give a negative branch length, so `to_phylo`
clamps it to zero and issues a warning saying how many branches were
affected. The drawn tree then understates those branches.

## Assembly algorithms

### `kmer_composition(original_sequence, k)`

Returns every length-*k* substring of the sequence in sorted order. Duplicates
are kept, because a k-mer that occurs several times needs that many edges in
the graph.

### `DeBruijnGraph(kmers)`

Builds a de Bruijn graph where each k-mer is an **edge** from its prefix
(first k-1 characters) to its suffix (last k-1 characters), and each distinct
(k-1)-mer is a **node**. `Node` keeps its incoming and outgoing edges, and
`Edge.set_start_node` / `set_end_node` attach an edge to both of its nodes.

`is_balanced()` returns `True` when every node has equal in- and out-degree.
For the k-mers of a single linear sequence that is normally **not** the case,
and that is expected rather than a bug. The first (k-1)-mer has one extra
outgoing edge and the last has one extra incoming edge. The graph is balanced
only when those two (k-1)-mers are the same, i.e. the sequence wraps around,
as in a circular genome.

### `get_eulerian_path()`

Reconstructs a sequence by walking every edge exactly once, using Hierholzer's
algorithm. If the graph is balanced it starts at a random node. Otherwise it
starts at the node with more outgoing than incoming edges. It walks forward
until stuck and adds each node to the path once all its edges are used, so
any side loops are spliced in automatically. The path's first node
contributes its full (k-1)-mer and each later node contributes its last
character.

```python
from assembly import kmer_composition, DeBruijnGraph

seq = 'AAGCTGCCGAATGCATTCGAGCTTAACTGAGAG'
DeBruijnGraph(kmer_composition(seq, k=5)).get_eulerian_path()  # -> seq
```

Things to keep in mind:

- **The answer is not always unique.** If a (k-1)-mer is repeated, several
  Eulerian paths can share the same k-mers, and the random edge choice picks
  one of them. With k=3, `AAGCTGCCGAATGCATTCG` has 12 different valid
  reconstructions. A larger *k* resolves shorter repeats.
- **Balanced graphs start at a random node**, so a circular sequence comes
  back as a rotation of the original.
- The start/end check only looks for *a* node with unequal degrees. It doesn't
  confirm there is exactly one start and one end, or that the graph is
  connected. On a graph with no Eulerian path it can return a partial
  sequence instead of failing. In that case it prints `ERROR` and returns
  `None` only if no start or end node exists at all.
- Used edges are tracked in a list, so each step costs O(E) and the whole walk
  is O(E^2). A set or a per-node index of unused edges would make it O(E).

Real sequencing reads add problems this graph doesn't handle yet: missing
k-mers, sequencing errors, unknown repeat counts, reverse-complement strands
and multiple chromosomes. Practical assemblers clean the graph and then
output contigs (maximal non-branching paths) instead of one Eulerian path.
