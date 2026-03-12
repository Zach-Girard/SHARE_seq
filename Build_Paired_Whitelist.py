#!/usr/bin/env python3

import argparse


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build STARsolo paired-end whitelist by generating all 3-barcode "
            "combinations from an 8bp barcode list."
        )
    )
    parser.add_argument(
        "--barcodes",
        required=True,
        help="Input file with one 8bp barcode per line.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output whitelist file (24bp barcodes, one per line).",
    )
    parser.add_argument(
        "--reverse-complement",
        action="store_true",
        help="Reverse complement each 8bp barcode before generating combinations.",
    )

    args = parser.parse_args()

    # Load 8bp barcodes
    with open(args.barcodes) as f:
        barcodes = [line.strip() for line in f if line.strip()]

    if args.reverse_complement:
        barcodes = [revcomp(bc) for bc in barcodes]

    if not barcodes:
        raise SystemExit("No barcodes loaded from input file.")

    with open(args.output, "w") as out:
        for b1 in barcodes:
            for b2 in barcodes:
                for b3 in barcodes:
                    out.write(f"{b1}{b2}{b3}\n")


if __name__ == "__main__":
    main()

