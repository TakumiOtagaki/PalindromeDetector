# Local Alignment: Smith-Waterman Algorithm with multiple matches


# import pandas as pd
import sys
from .util import *

INT_MAX = sys.maxsize

# Sequence is 1-indexed. matrix F is 0 padded.
# the restriction is as follows:
# 1. the band width is D.
# 2. the band is centered at the diagonal.
# In other words, the loop is n-D < i + j <= n


def palindrome_detector_with_band_overlap_restriction_in_linearspace(seq1, match, mismatch, gap, Threshold, band_width):
    seq2 = inverse_order_comlementary(seq1)
    gap = abs(gap)
    mismatch = - abs(mismatch)
    D = band_width
    n, m = len(seq1), len(seq2)

    def s(a, b):
        return (a == b) * match + (a != b) * mismatch

    def top_j(i):
        return max(0, n - D - i)

    def bottom_j(i):
        # return min(m, n - i)
        return n - i

    def traceback_local_palindrome(traceback):
        (i, j) = n-1, traceback[n][0]
        # alignments = [(x_start, x_end, y_start, y_end, x_string, y_string) for alignment in all alignments]
        # x_start, x_end, y_start, y_end are 0-indexed.
        alignments = []
        X, Y, score = "", "", 0
        x_start, y_start = INT_MAX, INT_MAX  # infty
        x_end, y_end = None, None
        while True:
            if j == top_j(i):
                # match region completed.
                if X != "" or Y != "":
                    y_end, y_start = (m - 1 - y_start) + 1, (m - 1 - y_end) + 1
                    alignments.append(
                        (score, x_start, x_end, y_start, y_end, X, Y))
                    X, Y, score = "", "", 0
                trace = traceback[i][j - top_j(i)]

                if i == 0:
                    break
                elif trace == 0:  # gap. from left.
                    i, j = i-1, top_j(i-1)
                else:
                    i, j = i-1, trace
                    x_end, y_end = i, j
                    x_start, y_start = INT_MAX, INT_MAX  # infty

            else:  # j != top_j(i)
                # 0: jump to zero, 1: match/mismatch, 2: gap, 3: gap
                trace = traceback[i][j - top_j(i)]
                if trace == 0:
                    i, j = i, top_j(i)
                elif trace == 1:
                    X, Y = seq1[i-1] + X, seq2[j-1] + Y
                    if seq1[i-1] == seq2[j-1]:
                        # matching region starts from this match.
                        x_start, y_start = min(x_start, i-1), min(y_start, j-1)
                    i, j = i-1, j-1
                    score += s(seq1[i], seq2[j])

                elif trace == 2:
                    X, Y = seq1[i-1] + X, "-" + Y
                    score -= gap
                    i, j = i-1, j

                elif trace == 3:
                    X, Y = "-" + X, seq2[j-1] + Y
                    i, j = i, j-1
                    score -= gap

                else:
                    print("invalid traceback matrix in palindrome LinearSpace")
                    sys.exit(1)

        return alignments

    # Initialize the scoring matrix
    F = zeros((n+1, D+1))
    F[0][0] = 0
    # F[0, ..., n][0, ..., m];
    for k in range(1, D+1):
        F[0][k] = max(0, F[0][k - 1])
    # Initialize the traceback matrix
    traceback = zeros((n+1, D+1))

    # Recurrence:
    for i in range(1, n+1):
        max_candidate = [F[i-1][0]] + [F[i-1][j - top_j(i-1)] -
                                       Threshold for j in range(1 + top_j(i-1), 1 + bottom_j(i-1))]
        F[i][0] = max(max_candidate)
        traceback[i][0] = argmax(max_candidate) + top_j(i-1)
        # for j in range(1, m+1):
        for j in range(top_j(i) + 1, bottom_j(i) + 1):
            max_candidate = [
                F[i][0],  ((j-1) - top_j(i-1) >= 0) * (F[i-1][(j-1) - top_j(i-1)] +
                                                       s(seq1[i-1], seq2[j-1])),
                F[i-1][j - top_j(i-1)] - gap, F[i][(j-1) - top_j(i)] - gap]

            F[i][j - top_j(i)] = max(max_candidate)
            # Update the traceback matrix
            traceback[i][j - top_j(i)] = argmax(max_candidate)


    alignments = traceback_local_palindrome(traceback)
    return alignments

def main():
    return


if __name__ == "__main__":
    main()
