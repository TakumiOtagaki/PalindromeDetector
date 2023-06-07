# 1. cd Frith研究室/opposite_dup_detector
# 2. python3 main.py input_file_path
# Or you can: python3 main.py input_file_path >> output_file_path

import sys
from modules.util import *
from modules.PalindromeDetector_with_restriction_in_linearspace import *


def main():
    args = sys.argv
    if len(args) == 2:
        input_file = args[1]
    else:
        print("Error: invalid number of arguments.")
        return
    # seq = "TAATTTCCCCCCCCCCCAAAAAAAATTTTTTTT"
    # seq = "AAAGAAACCTTTTCA"

    # else:
    #     print("use the default sequence.")
    #     seq = "ATTCCATAGGGGGAATCCTAGGTGACTGAACTC"
    #     Threshold = 20

    seq = read_single_fasta_file(input_file)
    seq = seq[:10000]
    # Threshold, D, match, mismatch, gap = 20, 10, 10, -5, -4
    print("Parameters:")
    print("Threshold, D, match, mismatch, gap")
    print(get_parameters("./parameters.txt"))
    Threshold, D, match, mismatch, gap = get_parameters("./parameters.txt")


    alignments = palindrome_detector_with_band_overlap_restriction_in_linearspace(seq, match, mismatch, gap, Threshold, D)
    print(f"\n\nResults: alignments which score >= {Threshold}")
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
            f"seq[{x_start}:{x_end}) = {seq[x_start:x_end]} and seq[{y_start}:{y_end}) = {seq[y_start:y_end]}")
        print("X, Y = ", X, Y)
        print("______________________________________")
    return


if __name__ == "__main__":
    main()
