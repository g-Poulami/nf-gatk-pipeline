#!/usr/bin/env python3
"""
generate_test_data.py
Writes synthetic paired FASTQ files and a minimal reference genome
to test/data/ for use with the test profile and -stub-run.

Usage:
    python3 test/generate_test_data.py
"""

import gzip
import os
import random

OUT_DIR  = os.path.join(os.path.dirname(__file__), "data")
SEED     = 42
N_READS  = 400
READ_LEN = 75

os.makedirs(OUT_DIR, exist_ok=True)
random.seed(SEED)

BASES = "ACGT"
QUAL  = "I" * READ_LEN


def random_seq(length):
    return "".join(random.choices(BASES, k=length))


def write_fastq_gz(path, n_reads, read_len):
    with gzip.open(path, "wt") as fh:
        for i in range(n_reads):
            fh.write(f"@read_{i}\n{random_seq(read_len)}\n+\n{QUAL[:read_len]}\n")
    print(f"  wrote {path}")


for sample in ("sample1", "sample2"):
    for mate in ("R1", "R2"):
        write_fastq_gz(
            os.path.join(OUT_DIR, f"{sample}_{mate}.fastq.gz"),
            N_READS,
            READ_LEN,
        )

ref_path = os.path.join(OUT_DIR, "ref.fa")
with open(ref_path, "w") as fh:
    for chrom in (1, 2, 3):
        seq = random_seq(10_000)
        fh.write(f">chr{chrom}\n")
        for start in range(0, len(seq), 60):
            fh.write(seq[start : start + 60] + "\n")
print(f"  wrote {ref_path}  (3 x 10 kb chromosomes)")

print()
print("Test data ready. Run the stub test with:")
print("  nextflow run main.nf -profile test -stub-run")
