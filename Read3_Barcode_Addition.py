#!/usr/bin/env python3


import argparse
import os
import gzip
try:
    from tqdm import tqdm
except ImportError:
    # Fallback when tqdm is unavailable in the active environment.
    class tqdm:
        def __init__(self, *args, **kwargs):
            pass
        def __enter__(self):
            return self
        def __exit__(self, exc_type, exc, tb):
            return False
        def update(self, n=1):
            pass


# For parsing the input barcode positions
def parse_coords(coord_str):
    """Converts 'start-end' string into (start, end) integer tuple."""
    try:
        start, end = coord_str.split('-')
        return int(start), int(end)
    except:
        raise argparse.ArgumentTypeError("Coords must be 'start-end' (e.g., 0-9)")


# ARGPARASE 
def main():
    parser = argparse.ArgumentParser(description="Extract barcodes from R2 and prepend to R3.")
    parser.add_argument("-r3", "--read3", required=True, help="Path to Read 3 FASTQ") # Currently written to accept zipped fastq
    parser.add_argument("-r2", "--read2", required=True, help="Path to Read 2 FASTQ")
    parser.add_argument("-bc", "--barcodes", required=True, 
                        help="Three barcode coords in R2 (e.g., 15-23,53-61,91-99)")
    parser.add_argument("-umi", "--umi_coords", required=True, 
                        help="UMI coords in R3 (e.g., 0-10)")
    
    args = parser.parse_args()


# Example script call
#  ./Read3_Barcode_Addition.py -r3 SHAREseq_R3.fastq -r2 SHAREseq_R2.fastq -bc 15-23,53-61,91-99 -umi 0-10



    # Parse coordinates
    bc_coords = [parse_coords(c.strip()) for c in args.barcodes.split(',')]
    if len(bc_coords) != 3:
        print("Error: Please provide exactly three barcode coordinate ranges.")
        return



    umi_start, umi_end = parse_coords(args.umi_coords)
    

    # Creation of Output file name
    # Logic to remove .gz from the output filename
    base_name = os.path.basename(args.read3)
    if base_name.endswith('.gz'):
        output_r3 = "withBarcodes_" + base_name[:-3] # Removes '.gz'
    else:
        output_r3 = "withBarcodes_" + base_name


    if os.path.exists(output_r3):
        raise FileExistsError(f"Safety Stop: The file {output_r3} already exists. "
                          "Delete it or rename your input to avoid overwriting.")
    

    # Barcode length to adjust UMI length for STARsolo
    total_bc_length = sum([(end - start) for start, end in bc_coords])
    

    # Calculate new UMI coordinates
    new_umi_start = umi_start + total_bc_length
    new_umi_end = umi_end + total_bc_length


    

    print(f"Reading: {args.read3}")
    print(f"Writing: {output_r3}")


    # Do Work 
    with gzip.open(args.read3, 'rt') as f1, gzip.open(args.read2, 'rt') as f2, open(output_r3, 'w') as out:
        with tqdm(unit=" reads", desc="Processing Sequences", dynamic_ncols=True) as pbar:
            for i, (r3_line, r2_line) in enumerate(zip(f1, f2)):
                
                # Sequence (Line 2) and Quality Score (Line 4)
                if i % 4 == 1 or i % 4 == 3:
                    extracted = "".join([r2_line[s:e] for s, e in bc_coords]) # extract barcodes from out list of coordinates
                    out.write(extracted + r3_line) # Prepend the barcodes/quality scores to the r3 line
                else:
                    out.write(r3_line)
                
                # Update the pbar every time we finish 1 full FASTQ record (4 lines)
                # First run processed approximately 13k-14k reads/sec
                if i % 4 == 3:
                    pbar.update(1)



    print("\n" + "-"*30)
    print(f"Processing FASTQ files complete.")
    print(f"New file saved as: {output_r3}")
    print(f"Concatenated Barcode Length: {total_bc_length} bp")
    print(f"Updated UMI Coordinates: {new_umi_start}-{new_umi_end}")
    print("-" * 30)

if __name__ == "__main__":
    main()
