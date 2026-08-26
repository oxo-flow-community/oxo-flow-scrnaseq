#!/usr/bin/env python3
"""Generate the 10x-structured synthetic fixtures for oxo-flow-scrnaseq.

The previous kit was 200 reads with 200 distinct barcodes: cellranger's
auto chemistry detection aborted (TXRNGR10001: minimum 10000 reads) and
cell-calling would have found 200 one-read cells. This generator emits,
per sample, ~150k cell reads across 800 barcodes with a realistic
log-normal count spread (median ~150, 50-2000) plus ~1.6k ambient reads
across 800 distinct empty-droplet barcodes (1-3 reads each): 10x v2
structure (R1 = 16bp cell barcode + 10bp UMI), R2 = 75bp drawn from the
reference gene intervals (refs/refdata.fa + refdata_genes GTF exons), so
chemistry detection sees >=10k reads and cell-calling finds ~800 real
cells per sample (the 800-cell / 800-empty-droplet ratio also keeps
cellbender's encoder minibatches populated with cells — 100 cells
among 2000 empties tripped 'Fewer than 4 cells passed to encoder
minibatch', live).

Both statistical shapes matter for cellbender remove-background:
- the cell-count SPREAD is essential — its empty-count estimation puts
  the cutoff INSIDE the cell-count distribution; a flat fixture (every
  cell exactly 100 reads) leaves the <=-cutoff partition empty and
  _peak_density_given_cutoff dies with an IndexError (live);
- the ambient droplets give the count distribution a non-cell tail.

Regenerate with:  python3 test/fixtures/generate_fixtures.py [refdata.fa]
"""
import gzip
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
REFDATA = sys.argv[1] if len(sys.argv) > 1 else "refs/refdata.fa"
N_BARCODES = 800
N_AMBIENT_BARCODES = 800
AMBIENT_READS = (1, 3)  # inclusive range per empty droplet
CELL_READS_LOGNORM = (math.log(150), 0.9)  # (log-median, log-sigma)
MIN_CELL_READS = 50
MAX_CELL_READS = 2000
SEED = 5

# gene exon intervals in refdata.fa coordinates (1-based, inclusive) — keep
# in sync with refs/refdata_genes.gtf: two exons per gene (grangers needs
# introns; single-exon transcripts crashed simpleaf's polars layer, live)
EXONS = [(100, 300), (400, 600), (800, 1000), (1100, 1300), (1500, 1700), (1800, 2000)]
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
    # REAL 10x v2 whitelist barcodes — the first 1600 of
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
    total_cell_reads = 0
    for sample in ("S1", "S2"):
        srng = random.Random(SEED + (1 if sample == "S2" else 0))
        r1_lines, r2_lines = [], []
        # cell reads: realistic log-normal per-cell count spread
        sample_cell_reads = 0
        for i in range(N_BARCODES):
            n_reads = max(MIN_CELL_READS, min(MAX_CELL_READS, int(srng.lognormvariate(*CELL_READS_LOGNORM))))
            for j in range(n_reads):
                r1, r2 = make_read(srng, genome, barcodes[i], f"@{sample}_barcode{i}_r{j}")
                r1_lines.append(r1)
                r2_lines.append(r2)
            sample_cell_reads += n_reads
        # ambient reads in empty droplets (non-cell tail; see docstring)
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
        print(f"  {sample}: {sample_cell_reads} cell reads, ambient droplets x{N_AMBIENT_BARCODES}")
        total_cell_reads += sample_cell_reads
    print(f"scrnaseq fixtures regenerated: {total_cell_reads} cell reads total ({N_BARCODES} barcodes/sample)")


if __name__ == "__main__":
    main()
