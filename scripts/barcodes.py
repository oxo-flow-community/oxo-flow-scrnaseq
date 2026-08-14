#!/usr/bin/env python3

# Port of nf-core/scrnaseq 4.2.0 modules/nf-core/anndata/barcodes/templates/barcodes.py
# Logic is identical; the Nextflow "${...}" substitutions are replaced by CLI arguments.
# SPDX-License-Identifier: MIT (upstream) — see NOTICE.md

import argparse

import anndata as ad
import pandas as pd

parser = argparse.ArgumentParser(description="Subset an h5ad file to a list of barcodes (CellBender filtered)")
parser.add_argument("--h5ad", required=True, help="input h5ad file")
parser.add_argument("--barcodes", required=True, help="CSV file with one barcode per line")
parser.add_argument("--output", required=True, help="output h5ad file")
args = parser.parse_args()

df = pd.read_csv(args.barcodes, header=None)
adata = ad.read_h5ad(args.h5ad)

adata = adata[df[0].values]

adata.write_h5ad(args.output)
