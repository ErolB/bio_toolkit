import numpy as np
import pandas as pd
import time

from utils import *


class Alignment:
    def __init__(self, seq1, seq2):
        pass

class FMIndex:
    def __init__(self, reference_sequence):
        # generate vectors based on BWT
        start = time.time()
        seq = reference_sequence + '$'
        self.size = len(seq)
        suffixes = [seq[i:] for i in range(len(seq)-1, -1, -1)]
        suffix_array = list(range(len(seq)-1, -1, -1))
        f = [seq[i] for i in range(len(seq)-1, -1, -1)]
        l = [seq[i-1] for i in range(len(seq)-1, 0, -1)] + ['$']
        frame = pd.DataFrame([suffixes, suffix_array, f, l]).transpose()
        frame.columns = ['suffix', 'order', 'F', 'L']
        frame = frame.sort_values('suffix').reset_index(drop=True)
        self.F = frame['F']
        self.L = frame['L']
        self.suffix_array = frame['order']
        self.inverse_suffix_array = {pos: row for row, pos in enumerate(self.suffix_array)}
        end = time.time()
        print(frame)
        print('FM Index generated in %s seconds' % str(end-start))

    def exact_search(self, query):
        # helper function
        def fm_search(whole_query, match_len, current_index):
            if match_len == len(whole_query):
                return current_index
            preceeding_char = whole_query[len(whole_query)-match_len-1]
            if self.L[current_index] == preceeding_char:
                new_index = self.inverse_suffix_array[self.suffix_array[current_index] - 1]
                return fm_search(whole_query, match_len+1, new_index)
            else:
                return None
        # initiate
        matches = []
        for i, char in enumerate(self.F):
            if char == query[-1]:
                result = fm_search(query, 1, i)
                if result is not None:
                    matches.append(self.suffix_array[result])
        return matches
        

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
    return score, alignment_space


def backtrack(seq1, seq2, alignment_space):
    final_seq1 = []
    final_seq2 = []
    if len(seq1)+1 != alignment_space.shape[0] or len(seq2)+1 != alignment_space.shape[1]:
        print('mismatch')
        return
    i = len(seq1)
    j = len(seq2)
    while i > 0 and j > 0:
        up_score = alignment_space[i-1,j]
        left_score = alignment_space[i, j-1]
        diag_score = alignment_space[i-1, j-1]
        min_score = min([up_score, left_score, diag_score])
        if diag_score == min_score:
            final_seq1.insert(0, seq1[i-1])
            final_seq2.insert(0, seq2[j-1])
            i -= 1
            j -= 1
        elif up_score == min_score:
            final_seq1.insert(0, seq1[i-1])
            final_seq2.insert(0, '-')
            i -= 1
        elif left_score == min_score:
            final_seq1.insert(0, '-')
            final_seq2.insert(0, seq2[j-1])
            j -= 1
    return ''.join(final_seq1), ''.join(final_seq2)
        

def generate_matrix(seqs, distance_function=needleman_wunsch_score):
    distance_matrix = np.zeros((len(seqs), len(seqs)))
    for i, s1 in enumerate(seqs):
        for j, s2 in enumerate(seqs):
            if i == j:
                # zero, not inf: neighbor_joining sums whole rows to build the
                # Q matrix, so a non-finite diagonal poisons every row sum
                distance_matrix[i,j] = 0.0
            else:
                distance_matrix[i,j] = distance_function(s1, s2)
    return distance_matrix


if __name__ == '__main__':
    records = read_fasta('covid.fasta')
    sequence = list(records.values())[0]
    query = 'YYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDL'
    print(sequence)
    print(query)
    fm = FMIndex(sequence)
    print(fm.exact_search(query))
    '''sequences = [
        'AATGCTCAAAA',
        'AATCCTCGAAA',
        'ATTGCTCAAAA',
        'AATTGCAAAA',
        'AATGCCAAAA'
    ]
    print(needleman_wunsch_score(sequences[0], sequences[1]))
    matrix = generate_matrix(sequences)'''

            
            
            
