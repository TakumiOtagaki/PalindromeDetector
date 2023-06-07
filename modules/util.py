def get_parameters(parameters_file):
    f = open(parameters_file, mode='r')
    data = f.readlines()[1]
    T, D, match, mismatch, gap = tuple(map(int, data.split(",")))

    return T, D, match, mismatch, gap



def inverse_order_comlementary(seq):
    """
    Returns the reverse order of the complementary sequence.
    """
    return seq[::-1].translate(str.maketrans('ATCG-', 'TAGC-'))


def abs(x):
    return x if x >= 0 else -x


def read_single_fasta_file(file_path):
    seq = ""
    multi = 0
    with open(file_path) as f:
        for line in f:
            if line[0] == ">":
                multi += 1
                continue
            if multi > 1:
                print("Error: multi fasta file.")
                print("Please use single fasta file, which contains only one DNA sequence.")
                return
            seq += line.strip().upper()
    return seq


def zeros(size):
    i, j = size
    return [[0 for _ in range(j)] for _ in range(i)]


def argmax(list_):
    return max(range(len(list_)), key=lambda x: list_[x])


def print_seq(seq):
    line = ""
    for i in range(1, 1 + len(seq)):
        line += seq[i - 1]
        if i % 10 == 0:
            line += " "
        if i % 50 == 0:
            line += "\n"
    print(line)
    return


def print_TSD_and_TIR(seq, TIR_and_TSDs, output_file_path):
    rev_seq = inverse_order_comlementary(seq)

    with open(output_file_path, mode='w') as f:
        print("----------------------------------------")
        f.write("----------------------------------------\n")
        if len(TIR_and_TSDs) == 0:
            print("No TIR and TSD found.")
            f.write("No TIR and TSD found.\n")
        for TIR_and_TSD in TIR_and_TSDs:
            TIR, TSD = TIR_and_TSD["TIR"], TIR_and_TSD["TSD"]
            print(
                f"TIR: (score = {TIR[0]}) \n \t Forward = {TIR[-2]} (seq[{TIR[1]}:{TIR[2]})) \n \t RevComp = {TIR[-1]} (seq[{TIR[3]}:{TIR[4]}))")
            f.write(
                f"TIR: (score = {TIR[0]}) \n \t Forward = {TIR[-2]} (seq[{TIR[1]}:{TIR[2]})) \n \t RevComp = {TIR[-1]} (seq[{TIR[3]}:{TIR[4]}))\n")

            print(
                f"TSD: (score = {TSD[0]}) \n \t Left  = {TSD[-2]} (seq[{TSD[1]}:{TSD[2]})) \n \t Right = {TSD[-1]} (seq[{TSD[3]}:{TSD[4]}))")
            f.write(
                f"TSD: (score = {TSD[0]}) \n \t Left  = {TSD[-2]} (seq[{TSD[1]}:{TSD[2]})) \n \t Right = {TSD[-1]} (seq[{TSD[3]}:{TSD[4]}))\n")

            print("----------------------------------------")
            f.write("----------------------------------------\n")
