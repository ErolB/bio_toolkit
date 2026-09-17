import numpy as np
import scipy as sc


class Alignment:
    def __init__(self, seq1, seq2):
        pass


def needleman_wunsch(seq1, seq2, mismatch_penalty=3, indel_penalty=5):
    alignment_space = np.zeros((len(seq1), len(seq2)))
    for i, char1 in enumerate(seq1):
        for j, char2 in enumerate(seq2):
            if i == 0:
                if j == 0:
                    continue  # leave top left corner as zero
                alignment_space[0,j] = alignment_space[0,j-1] + indel_penalty
                continue
            if j == 0:
                alignment_space[i,0] = alignment_space[i-1,0] + indel_penalty
                continue
            # try both indels and substitution and use the lowest value
            up_score = alignment_space[i-1,j] + indel_penalty
            left_score = alignment_space[i,j-1] + indel_penalty
            if char1 == char2:
                diag_score = alignment_space[i-1,j-1]
            else:
                diag_score = alignment_space[i-1,j-1] + mismatch_penalty
            alignment_space[i,j] = min([up_score, left_score, diag_score])
    # identify path with lowest cost
    score = min(alignment_space[-1,:])


if __name__ == '__main__':
    sequence1 = 'AAGATCCC'
    sequence2 = 'AACTTTTC'
    print(needleman_wunsch(sequence1, sequence2))
            
            
            
