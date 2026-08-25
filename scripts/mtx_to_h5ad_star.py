#!/usr/bin/env python

# Port of nf-core/scrnaseq 4.2.0 modules/local/mtx_to_h5ad/templates/mtx_to_h5ad_star.py
# Logic is identical; the Nextflow "${...}" substitutions are replaced by CLI arguments.
# SPDX-License-Identifier: MIT (upstream) — see NOTICE.md

# Set numba cache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import argparse

import anndata
import pandas as pd
import scanpy as sc
from anndata import AnnData, concat
from scipy.sparse import csr_matrix


def _mtx_to_adata(
    input: str,
    sample: str,
):
    adata = sc.read_10x_mtx(input)
    adata.obs["sample"] = sample
    adata.layers["count"] = adata.X.copy()

    velocyto_dir = f"velocyto_{input}"
    if os.path.exists(velocyto_dir):
        barcodes = os.path.join(velocyto_dir, "barcodes.tsv.gz")
        features = os.path.join(velocyto_dir, "features.tsv.gz")

        for matrix in ["ambiguous", "spliced", "unspliced"]:
            adata_state = sc.read_mtx(os.path.join(velocyto_dir, f"{matrix}.mtx.gz")).T

            adata_state.obs_names = pd.read_csv(barcodes, header=None, sep="\t")[0].values
            adata_state.var_names = pd.read_csv(features, header=None, sep="\t")[0].values

            missing_obs = adata.obs_names[~adata.obs_names.isin(adata_state.obs_names)]
            adata_missing = AnnData(
                X=csr_matrix((len(missing_obs), adata.shape[1])),
                obs=pd.DataFrame(index=missing_obs),
                var=adata_state.var,
            )
            adata_state = concat([adata_state, adata_missing], join="outer")
            adata_state = adata_state[adata.obs_names, adata.var["gene_ids"]].copy()

            adata.layers[matrix] = adata_state.X

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
    adata.var["gene_symbol"] = adata.var.index
    adata.var["gene_versions"] = adata.var["gene_ids"]
    adata.var.index = adata.var["gene_versions"].str.split(".").str[0].values
    adata.var_names_make_unique()  # in case user does not use ensembl references, names might not be unique

    # write results
    adata.write_h5ad(f"{output}")
    print(f"Wrote h5ad file to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a STARsolo counts directory (Gene/raw or Gene/filtered) to h5ad")
    parser.add_argument("--input", required=True, help="STARsolo Solo.out counts directory (raw or filtered)")
    parser.add_argument("--output", required=True, help="output h5ad file")
    parser.add_argument("--sample", required=True, help="sample id (meta.id)")
    args = parser.parse_args()

    # create the directory with the sample name
    os.makedirs(args.sample, exist_ok=True)

    # input_type (raw/filtered) comes from the caller via the file name
    input_to_adata(
        input_data=args.input,
        output=args.output,
        sample=args.sample,
    )
