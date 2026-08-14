#!/usr/bin/env python

# Port of nf-core/scrnaseq 4.2.0 modules/local/mtx_to_h5ad/templates/mtx_to_h5ad_cellranger.py
# Logic is identical; the Nextflow "${...}" substitutions are replaced by CLI arguments.
# SPDX-License-Identifier: MIT (upstream) — see NOTICE.md

# Set numba cache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import argparse
import platform

import anndata
import pandas as pd
import scanpy as sc


def _mtx_to_adata(
    input: str,
    sample: str,
):
    adata = sc.read_10x_h5(input)
    adata.var["gene_symbols"] = adata.var_names
    adata.var.set_index("gene_ids", inplace=True)
    adata.obs["sample"] = sample

    # reorder columns for 10x mtx files
    adata.var = adata.var[["gene_symbols", "feature_types", "genome"]]

    return adata


def input_to_adata(
    input_data: str,
    output: str,
    sample: str,
):
    print(f"Reading in {input_data}")

    # open main data
    adata = _mtx_to_adata(input_data, sample)

    # standard format
    # index are gene IDs and symbols are a column
    adata.var["gene_versions"] = adata.var.index
    adata.var.index = adata.var["gene_versions"].str.split(".").str[0].values
    adata.var_names_make_unique()

    # write results
    adata.write_h5ad(f"{output}")
    print(f"Wrote h5ad file to {output}")

    return adata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a 10x feature-barcode h5 file to h5ad (Cell Ranger)")
    parser.add_argument("--input", required=True, help="*_feature_bc_matrix.h5 file")
    parser.add_argument("--output", required=True, help="output h5ad file")
    parser.add_argument("--sample", required=True, help="sample id (meta.id)")
    args = parser.parse_args()

    # create the directory with the sample name
    os.makedirs(args.sample, exist_ok=True)

    # input_type (raw/filtered) comes from the caller via the file name
    adata = input_to_adata(
        input_data=args.input,
        output=args.output,
        sample=args.sample,
    )
