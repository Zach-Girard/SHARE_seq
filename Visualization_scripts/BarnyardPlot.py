#!/usr/bin/env python3

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate a human-mouse barnyard plot from STARsolo "
            "GeneFull/filtered/ matrix.mtx and features.tsv."
        )
    )
    p.add_argument(
        "--filtered-dir",
        required=True,
        help="Path to STARsolo GeneFull/filtered directory for a sample.",
    )
    p.add_argument(
        "--human-prefix",
        default="ENSG",
        help="Prefix for human gene IDs (default: ENSG).",
    )
    p.add_argument(
        "--mouse-prefix",
        default="ENSMUSG",
        help="Prefix for mouse gene IDs (default: ENSMUSG).",
    )
    p.add_argument(
        "--human-threshold",
        type=float,
        required=True,
        help="Minimum human fraction to classify a cell as human.",
    )
    p.add_argument(
        "--mouse-threshold",
        type=float,
        required=True,
        help="Maximum human fraction to classify a cell as mouse.",
    )
    p.add_argument(
        "--collision-low",
        type=float,
        required=True,
        help="Lower bound (exclusive) of human fraction to classify as collision.",
    )
    p.add_argument(
        "--collision-high",
        type=float,
        required=True,
        help="Upper bound (exclusive) of human fraction to classify as collision.",
    )
    p.add_argument(
        "--title",
        required=True,
        help="Plot title.",
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output PNG file for the barnyard plot.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    filtered_dir = pathlib.Path(args.filtered_dir)

    matrix_path = filtered_dir / "matrix.mtx"
    features_path = filtered_dir / "features.tsv"

    # 1. Load data
    matrix = scipy.io.mmread(matrix_path).tocsr()
    features = pd.read_csv(
        features_path, sep="\t", header=None, names=["id", "name", "type"]
    )

    # 2. Species masks
    human_mask = features["id"].astype(str).str.startswith(args.human_prefix)
    mouse_mask = features["id"].astype(str).str.startswith(args.mouse_prefix)

    # 3. Sum UMIs
    human_counts = matrix[human_mask, :].sum(axis=0).A1
    mouse_counts = matrix[mouse_mask, :].sum(axis=0).A1
    total_counts = human_counts + mouse_counts

    # 4. Classify cells by human fraction
    with np.errstate(divide="ignore", invalid="ignore"):
        h_ratio = human_counts / total_counts

    is_human = h_ratio >= args.human_threshold
    is_mouse = h_ratio <= args.mouse_threshold
    is_collision = (h_ratio > args.collision_low) & (h_ratio < args.collision_high)

    # 5. Plotting
    plt.figure(figsize=(8, 8))

    plt.scatter(
        human_counts[is_human],
        mouse_counts[is_human],
        color="red",
        s=10,
        alpha=0.5,
        label="Human",
    )
    plt.scatter(
        human_counts[is_mouse],
        mouse_counts[is_mouse],
        color="blue",
        s=10,
        alpha=0.5,
        label="Mouse",
    )
    plt.scatter(
        human_counts[is_collision],
        mouse_counts[is_collision],
        color="grey",
        s=10,
        alpha=0.5,
        label="Collision",
    )

    plt.xlabel("Human UMIs (ENSG)")
    plt.ylabel("Mouse UMIs (ENSMUSG)")
    plt.title(args.title)
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.3)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Barnyard plot saved as {out_path}")


if __name__ == "__main__":
    main()

