import numpy as np
import scipy as sc


class Alignment:
    def __init__(self, seq1, seq2):
        pass


def needleman_wunsch_score(seq1, seq2, mismatch_penalty=3, indel_penalty=5):
    # the extra row and column represent the empty prefix of each sequence
    alignment_space = np.zeros((len(seq1)+1, len(seq2)+1))
    alignment_space[0,:] = np.arange(len(seq2)+1) * indel_penalty
    alignment_space[:,0] = np.arange(len(seq1)+1) * indel_penalty
    for i, char1 in enumerate(seq1, start=1):
        for j, char2 in enumerate(seq2, start=1):
            # try both indels and substitution and use the lowest value
            up_score = alignment_space[i-1,j] + indel_penalty
            left_score = alignment_space[i,j-1] + indel_penalty
            if char1 == char2:
                diag_score = alignment_space[i-1,j-1]
            else:
                diag_score = alignment_space[i-1,j-1] + mismatch_penalty
            alignment_space[i,j] = min([up_score, left_score, diag_score])
    # identify path with lowest cost
    score = alignment_space[-1,-1]
    return score


def generate_matrix(seqs, distance_function=needleman_wunsch_score):
    distance_matrix = np.zeros((len(seqs), len(seqs)))
    for i, s1 in enumerate(seqs):
        for j, s2 in enumerate(seqs):
            if i == j:
                distance_matrix[i,j] = np.inf
            else:
                distance_matrix[i,j] = distance_function(s1, s2)
    return distance_matrix


if __name__ == '__main__':
    sequences = [
        'AATGCTCAAAA',
        'AATCCTCGAAA',
        'ATTGCTCAAAA',
        'AATTGCAAAA',
        'AATGCCAAAA'
    ]
    print(needleman_wunsch_score(sequences[0], sequences[1]))
    matrix = generate_matrix(sequences)

            
            
            
