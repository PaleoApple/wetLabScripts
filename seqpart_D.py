#!/usr/bin/env python3

import argparse
import sys


# --------------------------------------------------
# Properties
# --------------------------------------------------

CHARGE = {"E": -1, "D": -1, "R": 1, "K": 1}
AROMATICS = {"W": 1, "Y": 1, "F": 1}
PROLINES = {"P": 1}

MIN_SEGMENT = 1


# --------------------------------------------------
# FASTA Reader
# --------------------------------------------------

def read_fasta(filename):
    headers = []
    sequences = []

    with open(filename) as fh:
        current_seq = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                headers.append(line)
                sequences.append("")
                current_seq = []
            else:
                sequences[-1] += line.replace(" ", "")

    return headers, sequences


# --------------------------------------------------
# Scoring
# --------------------------------------------------

def score(sequence, start, stop, property_dict):
    total = 0
    for aa in sequence[start:stop + 1]:
        total += property_dict.get(aa, 0)
    return total


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Finds sequence partitioning point with maximal difference "
                    "in amino acid feature composition."
    )

    parser.add_argument("-fasta", required=True, help="Protein FASTA file")
    parser.add_argument("-charge", action="store_true",
                        help="Partition by charge (KR vs ED)")
    parser.add_argument("-aroma", action="store_true",
                        help="Partition by aromatics (FYW)")
    parser.add_argument("-pro", action="store_true",
                        help="Partition by proline (P)")
    parser.add_argument("-v", action="store_true", help="Verbose mode")
    parser.add_argument("-d", action="store_true", help="Debug mode")

    args = parser.parse_args()

    # Validate mode selection
    modes = [args.charge, args.aroma, args.pro]
    if sum(modes) != 1:
        sys.exit("Provide exactly one analysis mode: -charge OR -aroma OR -pro")

    if args.charge:
        property_dict = CHARGE
    elif args.aroma:
        property_dict = AROMATICS
    else:
        property_dict = PROLINES

    headers, sequences = read_fasta(args.fasta)

    for header, sequence in zip(headers, sequences):

        sequence_size = len(sequence)
        if sequence_size < 2 * MIN_SEGMENT:
            continue

        max_score = -100
        percent_score = percentN = percentC = 0
        score_map = {}

        for i in range(MIN_SEGMENT, sequence_size - MIN_SEGMENT):

            scoreN = score(sequence, 0, i - 1, property_dict)
            scoreC = score(sequence, i + 1, sequence_size - 1, property_dict)

            diff = abs(scoreN - scoreC)
            score_map[i] = diff

            if diff > max_score:
                max_score = diff

                percent_score = int((max_score / sequence_size) * 100)

                percentN = int(abs(scoreN / (i + 1) * 100))
                percentC = int(abs(scoreC / (sequence_size - i - 1) * 100))

        print(f"\n{header}")
        print(
            f"Max score difference: {max_score};\t"
            f"%score= {percent_score}\t"
            f"%scoreN= {percentN}\t"
            f"%scoreC= {percentC}\t\tResidues:\t",
            end=""
        )

        for i, diff in score_map.items():
            if diff == max_score:
                print(f"{sequence[i]}{i+1} ", end="")

        print("\n")


if __name__ == "__main__":
    main()
