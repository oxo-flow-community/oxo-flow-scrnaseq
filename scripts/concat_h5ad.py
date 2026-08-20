#!/usr/bin/env python

# Port of nf-core/scrnaseq 4.2.0 modules/local/concat_h5ad/templates/concat_h5ad.py
# Logic is identical (incl. the obs.sample label = file stem with input_type suffix, and the
# samplesheet join on `sample`); the Nextflow "${...}" substitutions are replaced by CLI arguments.
# SPDX-License-Identifier: MIT (upstream) — see NOTICE.md

# Set numba cache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def read_samplesheet(samplesheet):
    df = pd.read_csv(samplesheet)
    df.set_index("sample")

    # samplesheet may contain replicates, when it has,
    # group information from replicates and collapse with commas
    # only keep unique values using set()
    df = df.groupby(["sample"]).agg(lambda column: ",".join(set(column.astype(str))))

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate per-sample h5ad files into a combined h5ad")
    parser.add_argument("--samplesheet", required=True, help="input samplesheet CSV")
    parser.add_argument("--output", required=True, help="output h5ad file")
    parser.add_argument("--inputs", nargs="+", default=None,
                        help="per-sample h5ad files of one input_type (upstream contract)")
    parser.add_argument("--input-type", choices=["filtered", "cellbender_filter", "raw"], default=None,
                        help="when --inputs is omitted, derive per-sample paths from the samplesheet "
                             "(results/cellranger/mtx_conversions/<sample>/<sample>_<input-type>_matrix.h5ad)")
    args = parser.parse_args()

    # Open samplesheet as dataframe
    df_samplesheet = read_samplesheet(args.samplesheet)

    if args.inputs:
        inputs = args.inputs
    elif args.input_type:
        inputs = [
            f"results/cellranger/mtx_conversions/{sample}/{sample}_{args.input_type}_matrix.h5ad"
            for sample in df_samplesheet.index
        ]
    else:
        parser.error("either --inputs or --input-type is required")

    # find all h5ad and append to dict; keys are the basename minus '_matrix.h5ad'
    # (upstream uses str(path).replace("_matrix.h5ad", "") over files staged flat in the workdir)
    dict_of_h5ad = {str(Path(path).name).replace("_matrix.h5ad", ""): sc.read_h5ad(path) for path in inputs}

    # concat h5ad files
    adata = ad.concat(dict_of_h5ad, label="sample", merge="unique", index_unique="_")

    # merge with data.frame, on sample information
    adata.obs = adata.obs.join(df_samplesheet, on="sample", how="left").astype(str)
    adata.write_h5ad(args.output)

    print("Wrote h5ad file to {}".format(args.output))
