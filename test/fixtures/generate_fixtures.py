#!/usr/bin/env python3
"""Generate the 10x-structured synthetic fixtures for oxo-flow-scrnaseq.

The previous kit was 200 reads with 200 distinct barcodes: cellranger's
auto chemistry detection aborted (TXRNGR10001: minimum 10000 reads) and
cell-calling would have found 200 one-read cells. This generator emits,
per sample, 10000 cell read pairs (100 barcodes x 100 reads each) plus
2000 ambient reads across 2000 distinct empty-droplet barcodes (1-2
reads each): 10x v2 structure (R1 = 16bp cell barcode + 10bp UMI),
R2 = 75bp drawn from the reference gene intervals (refs/refdata.fa +
refdata_genes GTF exons), so chemistry detection sees >=10k reads and
cell-calling finds ~100 real cells per sample. The ambient fraction
matters: cellbender's empty-drop priors need a non-cell UMI
distribution to fit — a fixture whose non-cell barcodes all have
exactly 0 reads degenerates (live: IndexError in
_peak_density_given_cutoff on empty noncell_counts).

Regenerate with:  python3 test/fixtures/generate_fixtures.py [refdata.fa]
"""
import gzip
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
REFDATA = sys.argv[1] if len(sys.argv) > 1 else "refs/refdata.fa"
READS_PER_SAMPLE = 10000
N_BARCODES = 100
N_AMBIENT_BARCODES = 2000
AMBIENT_READS = (1, 3)  # inclusive range per empty droplet
SEED = 5

# gene intervals in refdata.fa coordinates (1-based, inclusive) — keep in
# sync with refs/refdata_genes.gtf exons: (100-500), (800-1300), (1500-1900)
EXONS = [(100, 500), (800, 1300), (1500, 1900)]
R2_LEN = 75


def load_genome():
    return "".join(l.strip() for l in open(REFDATA) if not l.startswith(">"))


def make_read(srng, genome, bc, rid):
    umi = "".join(srng.choice("ACGT") for _ in range(10))
    r1 = bc + umi  # 10x v2: 16bp barcode + 10bp UMI
    es, ee = srng.choice(EXONS)
    start = srng.randrange(es - 1, ee - R2_LEN)
    r2 = genome[start : start + R2_LEN]
    return (
        f"{rid}/1\n{r1}\n+\n{'I' * len(r1)}",
        f"{rid}/2\n{r2}\n+\n{'I' * R2_LEN}",
    )


def main():
    genome = load_genome()
    rng = random.Random(SEED)
    # REAL 10x v2 whitelist barcodes — the first 2100 of
    # 737K-august-2016.txt (cellranger's SC3Pv2 whitelist, the same
    # whitelist the chemistry detector matches against). Random 16bp
    # barcodes fail detection even with 10k reads (live: TXRNGR10002).
    barcode_file = os.path.join(HERE, "v2_barcodes.txt")
    with open(barcode_file) as fh:
        barcodes = [line.strip() for line in fh if line.strip()]
    need = N_BARCODES + N_AMBIENT_BARCODES
    if len(barcodes) < need:
        raise SystemExit(f"{barcode_file} has {len(barcodes)} barcodes, need {need}")
    os.makedirs(RAW, exist_ok=True)
    for sample in ("S1", "S2"):
        srng = random.Random(SEED + (1 if sample == "S2" else 0))
        r1_lines, r2_lines = [], []
        # cell reads: 100 barcodes x 100 reads
        for i in range(READS_PER_SAMPLE):
            bc = barcodes[i % N_BARCODES]
            r1, r2 = make_read(srng, genome, bc, f"@{sample}_barcode{i % N_BARCODES}_r{i}")
            r1_lines.append(r1)
            r2_lines.append(r2)
        # ambient reads in empty droplets (cellbender priors; see docstring)
        for i in range(N_AMBIENT_BARCODES):
            bc = barcodes[N_BARCODES + i]
            n_reads = srng.randint(*AMBIENT_READS)
            for j in range(n_reads):
                r1, r2 = make_read(srng, genome, bc, f"@{sample}_ambient{i}_r{j}")
                r1_lines.append(r1)
                r2_lines.append(r2)
        with gzip.open(os.path.join(RAW, f"{sample}_R1.fastq.gz"), "wt") as f1, gzip.open(
            os.path.join(RAW, f"{sample}_R2.fastq.gz"), "wt"
        ) as f2:
            f1.write("\n".join(r1_lines) + "\n")
            f2.write("\n".join(r2_lines) + "\n")
    print(
        f"scrnaseq fixtures regenerated: {READS_PER_SAMPLE} cell reads + "
        f"{N_AMBIENT_BARCODES} ambient droplets x 2 samples ({N_BARCODES} cell barcodes)"
    )


if __name__ == "__main__":
    main()
