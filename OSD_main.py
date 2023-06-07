# to execute: python3 main.py inputmode
# inputmode = 1: input from command line
# inputmode = 0: No input. Use the default values.
# python3 main.py: use the default values.

import sys
from util import *
from OSDDetector_with_restriction_in_linearspace import *


def main():
    args = sys.argv
    if len(args) == 1:
        inputmode = 0
    elif len(args) == 2:
        inputmode = int(args[1])
    else:
        print("Error: invalid number of arguments.")
        return
    # seq1 = "TAATTTCCCCCCCCCCCAAAAAAAATTTTTTTT"
    # seq1 = "AAAGAAACCTTTTCA"
    if inputmode == 1:
        seq1 = input("Enter the sequence: ")
        Threshold = int(input("Enter the Threshold: "))
    else:
        print("use the default sequence.")
        seq1 = "ATTCCATAGGGGGAATCCTAGGTGACTGAACTC"
        Threshold = 20
    match, mismatch, gap = 10, -5, -4
    D = 10
    F, traceback, alignments = OSD_detector_with_band_overlap_restriction_in_linearspace(
        seq1, match, mismatch, gap, Threshold, D)
    print("Threshold = ", Threshold)
    print("______________________________________")
    seq2 = inverse_order_comlementary(seq1)
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
