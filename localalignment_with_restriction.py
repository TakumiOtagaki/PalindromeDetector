# Local Alignment: Smith-Waterman Algorithm with multiple matches

import numpy as np
import pandas as pd
import sys
import time
from modules.util import *

INT_MAX = sys.maxsize

# Sequence is 1-indexed. matrix F is 0 padded.
# the restriction is as follows:
# 1. the band width is D.
# 2. the band is centered at the diagonal.
# In other words, the loop is n-D < i + j <= n


def OSD_detector_with_band_overlap_restriction(seq1, match, mismatch, gap, Threshold, band_width):
    seq2 = inverse_order_comlementary(seq1)
    gap = np.abs(gap)
    D = band_width
    n, m = len(seq1), len(seq2)
    D = min(D, n)

    def traceback_local_OSD(traceback, F):
        (i, j) = n-1, traceback[n][0]
        # alignments = [(x_start, x_end, y_start, y_end, x_string, y_string) for alignment in all alignments]
        # x_start, x_end, y_start, y_end are 0-indexed.
        alignments = []
        X, Y, score = "", "", 0
        x_start, y_start = INT_MAX, INT_MAX  # infty
        while True:
            if j == top_j(i):
                # match region completed.
                if X != "" or Y != "":
                    y_end, y_start = (m - 1 - y_start) + 1, (m - 1 - y_end) + 1
                    alignments.append(
                        (score, x_start, x_end, y_start, y_end, X, Y))
                    X, Y, score = "", "", 0
                trace = traceback[i][j]
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
                trace = traceback[i][j]
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
                    print("invalid traceback matrix.")
                    pass
        return alignments

    def s(a, b):
        return (a == b) * match + (a != b) * mismatch

    def top_j(i):
        return max(0, n - D - i)

    def bottom_j(i):
        return min(m, n - i)

    # Initialize the scoring matrix
    F = np.zeros((n+1, m+1), dtype=int)
    F[n-D+1][0] = 0
    # F[0, ..., n][0, ..., m];
    for j in range(n-D+2, n+1):
        F[0][j] = np.max([0, F[0][j-1]])
    # Initialize the traceback matrix
    traceback = np.zeros((n+1, m+1), dtype=int)
    # traceback = [["" for j in range(m+1)] for i in range(n+1)]
    # arrows = ["⇧", "↖", "←", "↑"]

    # Recurrence:
    for i in range(1, n+1):
        max_candidate = [F[i-1][top_j(i-1)]] + [F[i-1][j] -
                                                Threshold for j in range(top_j(i-1) + 1, bottom_j(i-1) + 1)]
        F[i][top_j(i)] = np.max(max_candidate)
        traceback[i][top_j(i)] = np.argmax(max_candidate) + top_j(i-1)
        # for j in range(1, m+1):
        for j in range(top_j(i) + 1, bottom_j(i) + 1):
            max_candidate = [
                F[i][top_j(i)], (F[i-1][j-1] + s(seq1[i-1], seq2[j-1])), F[i-1][j] - gap, F[i][j-1] - gap]
            F[i][j] = np.max(max_candidate)
            # Update the traceback matrix
            traceback[i][j] = np.argmax(max_candidate)
            # print("i, j = ", i, j)
            # print(f"F[{i-1}][{(j-1)}] = {F[i-1][j-1]}")
            # traceback[i][j] = arrows[np.argmax(max_candidate)]
    print(pd.DataFrame(F, columns=list(" "+seq2), index=list(" "+seq1)).T)
    print(pd.DataFrame(traceback, columns=list(" "+seq2), index=list(" "+seq1)).T)
    alignments = traceback_local_OSD(traceback, F)
    return F, traceback, alignments


def main():
    # seq1 = "TAATTTCCCCCCCCCCCAAAAAAAATTTTTTTT"
    # seq1 = "AAAGAAACCTTTTCA"
    seq1 = "ATTCCATAGGGGGAATCCTAGGTGACTGAACTC"
    match, mismatch, gap = 10, -5, -4
    Threshold = 20
    D = 10
    seq2 = inverse_order_comlementary(seq1)
    F, traceback, alignments = OSD_detector_with_band_overlap_restriction(
        seq1, match, mismatch, gap, Threshold, D)
    print("Threshold = ", Threshold)
    print("______________________________________")
    # print("alignments = ", alignments)
    alignments = sorted(alignments, key=lambda x: x[0], reverse=True)
    for alignment in alignments:
        score, x_start, x_end, y_start, y_end, X, Y = alignment
        print(alignment)
        print("score = ", score)
        print(
            f"x_start, x_end = {x_start}, {x_end}, y_start, y_end = {y_start}, {y_end}")
        print(
            f"seq1[{x_start}:{x_end}) = {seq1[x_start:x_end]} and seq1[{y_start}:{y_end}) = {seq1[y_start:y_end]}")
        print("X, Y = ", X, Y)
        print("______________________________________")
    return


if __name__ == "__main__":
    main()
