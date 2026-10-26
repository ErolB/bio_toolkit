import random

# utility functions

def kmer_composition(original_sequence, k):
    kmers = []
    for i in range(len(original_sequence)-k+1):
        kmers.append(original_sequence[i:i+k])
    return sorted(kmers)

def get_prefix(sequence):
    return sequence[:-1]

def get_suffix(sequence):
    return sequence[1:]


# classes for graph structure

class Node:
    def __init__(self, sequence):
        self.sequence = sequence
        self.in_edges = []
        self.out_edges = []

    def connect_in(self, edge):
        self.in_edges.append(edge)

    def connect_out(self, edge):
        self.out_edges.append(edge)

    def in_degree(self):
        return len(self.in_edges)

    def out_degree(self):
        return len(self.out_edges)


class Edge:
    def __init__(self, sequence):
        self.sequence = sequence
        self.start_node = None
        self.end_node = None

    def set_start_node(self, node):
        self.start_node = node
        node.connect_out(self)

    def set_end_node(self, node):
        self.end_node = node
        node.connect_in(self)

    def get_start(self):
        return self.start_node

    def get_destination(self):
        return self.end_node

class DeBruijnGraph:
    def __init__(self, seqs):
        self.edges = [Edge(s) for s in seqs]
        self.nodes = []
        for e in self.edges:
            prefix = get_prefix(e.sequence)
            suffix = get_suffix(e.sequence)
            node_seqs = {n.sequence: n for n in self.nodes}
            if prefix in node_seqs.keys():
                e.set_start_node(node_seqs[prefix])
            else:
                new_node = Node(prefix)
                self.nodes.append(new_node)
                e.set_start_node(new_node)
                node_seqs[prefix] = new_node
            
            if suffix in node_seqs.keys():
                e.set_end_node(node_seqs[suffix])
            else:
                new_node = Node(suffix)
                self.nodes.append(new_node)
                e.set_end_node(new_node)
                node_seqs[suffix] = new_node

    def is_balanced(self):
        for n in self.nodes:
            if n.in_degree() != n.out_degree():
                return False
        return True

    def get_eulerian_path(self):
        visited_edges = []
        if self.is_balanced():
            start_node = random.choice(self.nodes)
            end_node = start_node
        else:
            start_node = None
            for n in self.nodes:
                if n.in_degree() < n.out_degree():
                    start_node = n
                    break
            if start_node is None:
                print('ERROR')
                return
            end_node = None
            for n in self.nodes:
                if n.in_degree() > n.out_degree():
                    end_node = n
                    break
            if end_node is None:
                print('ERROR')
                return
        # traverse graph
        path = []
        stack = [start_node]
        while stack:
            n = stack[-1]
            outgoing_edges = n.out_edges
            unused_edges = [item for item in outgoing_edges if item not in visited_edges]
            if unused_edges:
                selected_edge = random.choice(unused_edges)
                stack.append(selected_edge.get_destination())
                visited_edges.append(selected_edge)
            else:
                stack.pop()
                path.append(n)
        path_string = ''
        for i, n in enumerate(reversed(path)):
            if i == 0:
                path_string += n.sequence
            else:
                path_string += n.sequence[-1]
        return path_string
            



if __name__ == '__main__':
    kmers = kmer_composition('AAGCTGCCGAATGCATTCGAGCTTAACTGAGAG', k=5)
    graph = DeBruijnGraph(kmers)
    print(graph.get_eulerian_path())
