#!/usr/bin/env python


import pandas as pd
import scipy.io
import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Split STARsolo PDX output into Human, Mouse, and Doublet sets.')
    
    # Required Arguments
    parser.add_id = parser.add_argument('--input', type=str, required=True, 
                        help='Path to STARsolo filtered output folder (containing matrix.mtx, features.tsv, barcodes.tsv)')
    parser.add_argument('--output', type=str, required=True, 
                        help='Directory to save the split species files')
    
    # Optional Arguments
    parser.add_argument('--purity', type=float, default=0.9, 
                        help='Purity threshold for species assignment (default: 0.9)')
    parser.add_argument('--human_prefix', type=str, default='ENSG', 
                        help='Prefix for human genes (default: ENSG)')
    parser.add_argument('--mouse_prefix', type=str, default='ENSMUSG', 
                        help='Prefix for mouse genes (default: ENSMUSG)')

    args = parser.parse_args()

    # Verify input paths
    if not os.path.exists(args.input):
        print(f"Error: Input directory {args.input} does not exist.")
        sys.exit(1)

    # 1. Load the data
    print(f"--- Loading Matrix from {args.input} ---")
    try:
        # STARsolo can output .mtx or .mtx.gz; scipy handle both if named correctly
        mat = scipy.io.mmread(os.path.join(args.input, "matrix.mtx")).tocsr()
        features = pd.read_csv(os.path.join(args.input, "features.tsv"), sep='\t', header=None)
        barcodes = pd.read_csv(os.path.join(args.input, "barcodes.tsv"), sep='\t', header=None)
    except FileNotFoundError as e:
        print(f"Error: Missing essential files in input directory. {e}")
        sys.exit(1)

    # 2. Identify indices based on prefixes
    human_idx = features[features[0].str.startswith(args.human_prefix)].index
    mouse_idx = features[features[0].str.startswith(args.mouse_prefix)].index

    print(f"Detected {len(human_idx)} human genes and {len(mouse_idx)} mouse genes.")

    # 3. Calculate per-barcode sums
    # STARsolo matrices are typically (genes x barcodes)
    if mat.shape[0] == len(features):
        h_sums = mat[human_idx, :].sum(axis=0).A1
        m_sums = mat[mouse_idx, :].sum(axis=0).A1
    else:
        # Handle transposed matrix if necessary
        h_sums = mat[:, human_idx].sum(axis=1).A1
        m_sums = mat[:, mouse_idx].sum(axis=1).A1

    total_sums = h_sums + m_sums
    
    # 4. Species Assignment Logic
    # Use Series to handle division safely
    h_prop = pd.Series(h_sums) / pd.Series(total_sums)
    m_prop = pd.Series(m_sums) / pd.Series(total_sums)

    labels = []
    for h, m in zip(h_prop, m_prop):
        if h >= args.purity:
            labels.append('human')
        elif m >= args.purity:
            labels.append('mouse')
        else:
            labels.append('collision')

    # 5. Create Metadata Table
    master_df = pd.DataFrame({
        'barcode': barcodes[0],
        'label': labels,
        'human_prop': h_prop.fillna(0),
        'mouse_prop': m_prop.fillna(0)
    })
    
    if not os.path.exists(args.output): 
        os.makedirs(args.output)
    
    master_df.to_csv(os.path.join(args.output, "species_metadata.csv"), index=False)

    # 6. Export Filtered Matrices
    for species in ['human', 'mouse', 'collision']:
        subset_indices = master_df[master_df['label'] == species].index
        if len(subset_indices) == 0:
            print(f"No cells found for category: {species}")
            continue

        print(f"Writing {len(subset_indices)} cells to {species} folder...")
        spec_path = os.path.join(args.output, species)
        os.makedirs(spec_path, exist_ok=True)
        
        # Subset matrix and barcodes
        if mat.shape[0] == len(features):
            sub_mat = mat[:, subset_indices]
        else:
            sub_mat = mat[subset_indices, :]
            
        sub_barcodes = barcodes.iloc[subset_indices]

        # Save files to output directory
        scipy.io.mmwrite(os.path.join(spec_path, "matrix.mtx"), sub_mat)
        sub_barcodes.to_csv(os.path.join(spec_path, "barcodes.tsv"), sep='\t', index=False, header=False)
        features.to_csv(os.path.join(spec_path, "features.tsv"), sep='\t', index=False, header=False)

    print("--- Processing Complete ---")

if __name__ == "__main__":
    main()
