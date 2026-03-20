#!/usr/bin/env python3

"""
Generate synthetic demultiplexed FASTQ files for testing the pipeline.

Creates gzipped FASTQs in demux/ for a single sample 'sampleA':

  - sampleA_R1.fastq.gz : 5000 reads, 30 bp random cDNA
  - sampleA_R2.fastq.gz : 5000 reads, 99 bp layout (matches bc_coords
                          15-23,53-61,91-99 as Python half-open slices on R2):
                            * 15 bp constant linker (indices 0–14)
                            * 8 bp barcode 1  ([15:23])
                            * 30 bp linker    ([23:53])
                            * 8 bp barcode 2  ([53:61])
                            * 30 bp linker    ([61:91])
                            * 8 bp barcode 3  ([91:99])
  - sampleA_R3.fastq.gz : 5000 reads, 30 bp:
                            * 10 bp UMI (pseudo-random)
                            * 15 bp poly-T
                            * 5 bp cDNA

Barcodes are taken from barcodes_RC.txt in the project root.
We cycle through many barcodes (not just the first 3), and every 10th
read introduces 1 mismatch per 8 bp block to exercise the correction logic.
"""

import gzip
import random
from pathlib import Path


def load_barcodes(path: Path, n: int = 12):
    with path.open() as f:
        barcodes = [line.strip() for line in f if line.strip()]
    if len(barcodes) < n:
        raise RuntimeError(f"Need at least {n} barcodes in {path}, found {len(barcodes)}")
    return barcodes


def mutate_base(b: str) -> str:
    bases = ["A", "C", "G", "T"]
    b = b.upper()
    if b in bases:
        bases.remove(b)
    return random.choice(bases)


def make_r1_seq(i: int) -> str:
    # 30 bp cDNA: pseudo-random but deterministic per read index
    random.seed(100000 + i)
    return "".join(random.choice("ACGT") for _ in range(30))


def make_r3_seq(i: int) -> str:
    # 10 bp UMI: pseudo-random but deterministic per read index
    random.seed(i)
    umi = "".join(random.choice("ACGT") for _ in range(10))
    poly_t = "T" * 15
    cdna = ("GCTAAGCT" * 3)[:5]
    return umi + poly_t + cdna  # length 30


def make_r2_seq(bc1: str, bc2: str, bc3: str, mismatch: bool = False) -> str:
    """
    Build 99 bp R2 with constant linkers and three 8 bp barcodes (total 99 bp):
      0–14  : constant linker (15 bp)
      15–22 : bc1 (8 bp)   — same span as slice [15:23] in Read3_Barcode_Addition.py
      23–52 : linker (30 bp)
      53–60 : bc2 (8 bp)   — slice [53:61]
      61–90 : linker (30 bp)
      91–98 : bc3 (8 bp)   — slice [91:99]
    """
    fillerA = ("ACGTACGTACGTACG")[:15]
    fillerB = ("GGGGCCCCAAAATTTT" * 2)[:30]
    fillerC = ("TTTTGGGGCCCCAAAA" * 2)[:30]

    if mismatch:
        # Introduce 1 mismatch into exactly one of the 3 barcode blocks
        def mutate_block(bc):
            pos = random.randrange(len(bc))
            return bc[:pos] + mutate_base(bc[pos]) + bc[pos + 1 :]

        which = random.choice([1, 2, 3])
        if which == 1:
            bc1 = mutate_block(bc1)
        elif which == 2:
            bc2 = mutate_block(bc2)
        else:
            bc3 = mutate_block(bc3)

    seq = fillerA + bc1 + fillerB + bc2 + fillerC + bc3
    assert len(seq) == 99, f"R2 sequence length is {len(seq)}, expected 99"
    return seq


def main():
    base_dir = Path(__file__).resolve().parent
    demux_dir = base_dir / "demux"
    demux_dir.mkdir(exist_ok=True)

    barcodes_path = base_dir / "barcodes_RC.txt"
    bc_list = load_barcodes(barcodes_path, n=12)  # cycle through barcodes from file

    r1_path = demux_dir / "sampleA_R1.fastq.gz"
    r2_path = demux_dir / "sampleA_R2.fastq.gz"
    r3_path = demux_dir / "sampleA_R3.fastq.gz"

    n_reads = 5000

    with gzip.open(r1_path, "wt") as r1_out, gzip.open(
        r2_path, "wt"
    ) as r2_out, gzip.open(r3_path, "wt") as r3_out:
        for i in range(1, n_reads + 1):
            read_id = f"sampleA_read{i}"

            # R1: 30 bp random cDNA
            seq_r1 = make_r1_seq(i)
            qual_r1 = "I" * len(seq_r1)
            r1_out.write(f"@{read_id}\n{seq_r1}\n+\n{qual_r1}\n")

            # R3: 10 bp UMI + 15 bp poly-T + 5 bp cDNA
            seq_r3 = make_r3_seq(i)
            qual_r3 = "I" * len(seq_r3)
            r3_out.write(f"@{read_id}\n{seq_r3}\n+\n{qual_r3}\n")

            # Choose 3 barcodes for this read (cycle through list)
            idx1 = (i - 1) % len(bc_list)
            idx2 = (i * 3) % len(bc_list)
            idx3 = (i * 5) % len(bc_list)
            bc1, bc2, bc3 = bc_list[idx1], bc_list[idx2], bc_list[idx3]

            # Every 10th read: introduce 1 mismatch in 1 random barcode block
            mismatch = (i % 10 == 0)
            seq_r2 = make_r2_seq(bc1, bc2, bc3, mismatch=mismatch)
            qual_r2 = "I" * len(seq_r2)
            r2_out.write(f"@{read_id}\n{seq_r2}\n+\n{qual_r2}\n")

    print("Generated test FASTQs in demux/:")
    print(f"  {r1_path} (5000 reads, 30 bp random cDNA)")
    print(f"  {r2_path} (5000 reads, 99 bp with barcodes at default positions)")
    print(f"  {r3_path} (5000 reads, 30 bp: 10 bp UMI + 15 bp poly-T + 5 bp cDNA)")


if __name__ == "__main__":
    main()

