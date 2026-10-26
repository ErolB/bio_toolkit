from alignment import *
from utils import *
from assembly import *

if __name__ == '__main__':
    print('reading file...')
    genome = list(read_fasta('covid.fasta').values())[0]
    print('generating k-mers...')
    genome_kmers = kmer_composition(genome, k=7)
    print('building graph...')
    graph = DeBruijnGraph(genome_kmers)
    print('rebuilding sequence...')
    reconstruction = graph.get_eulerian_path()
    print('scoring...')
    print(len(genome), len(reconstruction))
    print(genome == reconstruction)