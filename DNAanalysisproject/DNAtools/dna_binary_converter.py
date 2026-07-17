#!/usr/bin/env python3
"""
dna_binary_converter.py
========================
Converts between binary and DNA using the 2-bit encoding:

    A = 00   C = 01   G = 10   T = 11

Usage
-----
    python dna_binary_converter.py --to-dna 0001101100
    python dna_binary_converter.py --to-binary ACGT
"""

import argparse

BASE_TO_BITS = {"A": "00", "C": "01", "G": "10", "T": "11"}
BITS_TO_BASE = {v: k for k, v in BASE_TO_BITS.items()}


def binary_to_dna(bits: str) -> str:
    bits = bits.strip().replace(" ", "")
    if len(bits) % 2 != 0:
        raise ValueError("Binary string length must be a multiple of 2.")
    dna = []
    for i in range(0, len(bits), 2):
        pair = bits[i:i + 2]
        if pair not in BITS_TO_BASE:
            raise ValueError(f"Invalid bit pair: '{pair}' (must be 00, 01, 10, or 11)")
        dna.append(BITS_TO_BASE[pair])
    return "".join(dna)


def dna_to_binary(seq: str) -> str:
    seq = seq.strip().upper()
    bits = []
    for base in seq:
        if base not in BASE_TO_BITS:
            raise ValueError(f"Invalid base: '{base}' (must be A, C, G, or T)")
        bits.append(BASE_TO_BITS[base])
    return "".join(bits)


def write_fasta(seq: str, path: str, header: str = "converted_sequence", wrap: int = 70) -> None:
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


def main():
    parser = argparse.ArgumentParser(description="Convert between binary and DNA (A=00 C=01 G=10 T=11).")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--to-dna", metavar="BITS", help="Binary string to convert to DNA.")
    group.add_argument("--to-binary", metavar="SEQ", help="DNA sequence to convert to binary.")
    parser.add_argument("--out", metavar="PATH", help="Write result to a file instead of printing it. "
                         "If converting --to-dna, writes FASTA format; otherwise writes plain text.")
    parser.add_argument("--header", default="converted_sequence",
                         help="FASTA header/ID to use when writing DNA output (default: converted_sequence).")
    args = parser.parse_args()

    if args.to_dna:
        result = binary_to_dna(args.to_dna)
        if args.out:
            write_fasta(result, args.out, header=args.header)
            print(f"Wrote {args.out} ({len(result)} bp)")
        else:
            print(result)
    else:
        result = dna_to_binary(args.to_binary)
        if args.out:
            with open(args.out, "w") as fh:
                fh.write(result + "\n")
            print(f"Wrote {args.out}")
        else:
            print(result)


if __name__ == "__main__":
    main()
