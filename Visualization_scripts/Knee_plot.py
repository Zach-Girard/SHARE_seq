#!/usr/bin/env python3

import argparse
import pathlib

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a GeneFull knee plot using STARsolo "
            "GeneFull/UMIperCellSorted.txt and an estimated cell count "
            "provided by the pipeline."
        )
    )
    parser.add_argument(
        "--genefull-dir",
        required=True,
        help="Path to STARsolo GeneFull directory for a sample (e.g. Solo.out/GeneFull)",
    )
    parser.add_argument(
        "--estimated-cells",
        type=float,
        required=True,
        help="Estimated number of cells (from STARsolo Summary.csv).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output PNG path for the knee plot.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    gene_dir = pathlib.Path(args.genefull_dir)

    umi_path = gene_dir / "UMIperCellSorted.txt"

    # Load UMI-per-cell sorted data
    data = pd.read_csv(umi_path, header=None, names=["umis"])
    data["rank"] = range(1, len(data) + 1)

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.loglog(data["rank"], data["umis"], label="Cells", color="blue", linewidth=2)

    # Add a vertical line for STAR's estimate (passed in by the pipeline)
    est = args.estimated_cells
    plt.axvline(
        x=est,
        color="red",
        linestyle="--",
        label=f"STAR Estimated Cells ({int(est)})",
    )

    plt.title("Knee Plot (GeneFull)")
    plt.xlabel("Cell Rank (Log Scale)")
    plt.ylabel("UMIs per Cell (Log Scale)")
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved as {out_path}")


if __name__ == "__main__":
    main()
