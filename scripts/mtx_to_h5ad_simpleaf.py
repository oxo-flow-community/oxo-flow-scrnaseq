#!/usr/bin/env python

# Port of nf-core/scrnaseq 4.2.0 modules/local/mtx_to_h5ad/templates/mtx_to_h5ad_simpleaf.py
# Logic is identical; the Nextflow "${...}" substitutions are replaced by CLI arguments.
# SPDX-License-Identifier: MIT (upstream) — see NOTICE.md

# Set numba cache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import argparse
import json
import platform

import anndata
import pandas as pd
import scanpy as sc


def input_to_adata(
    input_data: str,
    output: str,
    sample: str,
):
    print(f"Reading in {input_data}")

    # open main data
    simpleaf_h5ad_path = f"{input_data}/alevin/quants.h5ad"

    # the simpleaf quant module exports an h5ad file.
    adata = sc.read_h5ad(simpleaf_h5ad_path)
    adata.obs_names = adata.obs["barcodes"].values
    adata.var_names = adata.var["gene_id"].values
    adata.obs["sample"] = sample

    # sort adata column- and row- wise to avoid positional differences
    adata = adata[adata.obs_names.sort_values(), adata.var_names.sort_values()].copy()

    # standard format
    # index are gene IDs and symbols are a column
    adata.var["gene_versions"] = adata.var["gene_id"]
    adata.var.index = adata.var["gene_versions"].str.split(".").str[0].values
    adata.var_names_make_unique()

    # sort adata column- and row- wise to avoid positional differences
    adata = adata[adata.obs_names.sort_values(), adata.var_names.sort_values()].copy()

    # Remove runtime to prevent hash changes
    simpleaf_map_info = json.loads(adata.uns["simpleaf_map_info"])
    simpleaf_map_info.pop("runtime_seconds")
    adata.uns["simpleaf_map_info"] = json.dumps(simpleaf_map_info, sort_keys=True)

    # write results
    adata.write_h5ad(f"{output}")
    print(f"Wrote h5ad file to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a simpleaf alevin quants.h5ad to a standard h5ad")
    parser.add_argument("--input", required=True, help="simpleaf quant output directory (containing alevin/quants.h5ad)")
    parser.add_argument("--output", required=True, help="output h5ad file")
    parser.add_argument("--sample", required=True, help="sample id (meta.id)")
    args = parser.parse_args()

    # create the directory with the sample name
    os.makedirs(args.sample, exist_ok=True)

    # input_type (always raw for simpleaf) comes from the caller via the file name
    input_to_adata(
        input_data=args.input,
        output=args.output,
        sample=args.sample,
    )
