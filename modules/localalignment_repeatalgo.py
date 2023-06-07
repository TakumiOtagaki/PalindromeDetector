# Local Alignment: Smith-Waterman Algorithm with multiple matches

import pandas as pd
import sys
from .util import *

INT_MAX = sys.maxsize


def local_alignment(seq1, seq2, match, mismatch, gap, Threshold):
    gap = abs(gap)
    mismatch = -abs(mismatch)

    def s(a, b):
        if a == b:
            return match
        else:
            return mismatch
    n, m = len(seq1), len(seq2)

    def traceback_local(traceback):
        (i, j) = n-1, traceback[n][0]
        # alignments = [(x_start, x_end, y_start, y_end, x_string, y_string) for alignment in all alignments]
        # x_start, x_end, y_start, y_end are 0-indexed.
        alignments = []
        X, Y, score = "", "", 0
        x_end, y_end = i-1, j-1
        x_start, y_start = INT_MAX, INT_MAX  # infty
        while True:
            # print("i, j = ", i, j)
            if j == 0:
                # match region completed.
                if X != "" or Y != "":
                    # print((score, x_start, x_end, y_start, y_end, X, Y))
                    alignments.append(
                        (score, x_start, x_end, y_start, y_end, X, Y))
                    X, Y, score = "", "", 0
                trace = traceback[i][j]
                if i == 0:
                    break
                elif trace == 0:  # gap. from left.
                    i, j = i-1, 0
                else:
                    i, j = i-1, trace
                    x_end, y_end = i-1, j-1
                    x_start, y_start = INT_MAX, INT_MAX  # infty
                    # print("x_end, y_end = ", x_end, y_end)

            else:  # j != 0
                # 0: jump to zero, 1: match/mismatch, 2: gap, 3: gap
                trace = traceback[i][j]
                # print("trace = ", trace)
                if trace == 0:
                    i, j = i, 0
                elif trace == 1:
                    i, j = i-1, j-1
                    X, Y = seq1[i] + X, seq2[j] + Y
                    if seq1[i] == seq2[j]:
                        x_start, y_start = min(x_start, i), min(y_start, j)
                    score += s(seq1[i], seq2[j])

                elif trace == 2:
                    i, j = i-1, j
                    X, Y = seq1[i] + X, "-" + Y
                    score -= gap
                elif trace == 3:
                    i, j = i, j-1
                    X, Y = "-" + X, seq2[j] + Y
                    score -= gap
                else:
                    print("invalid traceback matrix in Local alignment")
                    # print(pd.DataFrame(F, index=list(" "+seq1)).T)
                    # print(pd.DataFrame(traceback, index=list(" "+seq1)).T)
                    sys.exit(1)
        return alignments

    # Initialize the scoring matrix
    F = zeros((n+1, m+1))
    # F[0, ..., n][0, ..., m];
    for j in range(m+1):
        F[0][j] = max(0, F[0][j-1])
    # Initialize the traceback matrix
    traceback = zeros((n+1, m+1))
    # traceback = [["" for j in range(m+1)] for i in range(n+1)]
    # arrows = ["↑", "↖", "←", "↑"]

    # Recurrence:
    for i in range(1, n+1):
        max_candidate = [F[i-1][0]] + [F[i-1][j] -
                                       Threshold for j in range(1, m+1)]
        F[i][0] = max(max_candidate)
        traceback[i][0] = argmax(max_candidate)
        for j in range(1, m+1):
            # Calculate the score for each possible move
            max_candidate = [F[i][0], F[i-1][j-1] +
                             s(seq1[i-1], seq2[j-1]), F[i-1][j] - gap, F[i][j-1] - gap]
            F[i][j] = max(max_candidate)
            # Update the traceback matrix
            # traceback[i][j] = arrows[np.argmax(max_candidate)]
            traceback[i][j] = argmax(max_candidate)
    alignments = traceback_local(traceback)

    return F, traceback, alignments


def main():
    seq1, seq2 = "ATTTTTAAAAAATTCCCCCGG", "ATATATAGGGGGGGTTTTTTTTTTAAAAAAAAAA"
    match, mismatch, gap = 10, -5, -4
    Threshold = 20
    F, traceback, alignments = local_alignment(
        seq1, seq2, match, mismatch, gap, Threshold)
    print(pd.DataFrame(F, columns=list(" "+seq2), index=list(" "+seq1)).T)
    print(pd.DataFrame(traceback, columns=list(" "+seq2), index=list(" "+seq1)).T)
    print(alignments)
    return


if __name__ == "__main__":
    main()
