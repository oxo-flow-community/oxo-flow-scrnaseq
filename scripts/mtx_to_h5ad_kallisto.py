#!/usr/bin/env python

# Port of nf-core/scrnaseq 4.2.0 modules/local/mtx_to_h5ad/templates/mtx_to_h5ad_kallisto.py
# Logic is identical; the Nextflow "${...}" substitutions are replaced by CLI arguments.
# SPDX-License-Identifier: MIT (upstream) — see NOTICE.md

# Set numba cache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import argparse
import glob
import platform

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from anndata import concat as concat_ad
from scipy.sparse import csr_matrix


def _mtx_to_adata(
    matrix: str,
    barcodes: str,
    features: str,
):
    """Load kallisto-formatted mtx files into AnnData."""
    adata = sc.read_mtx(matrix)
    adata.obs_names = pd.read_csv(barcodes, header=None, sep="\t")[0].values
    adata.var_names = pd.read_csv(features, header=None, sep="\t")[0].values
    return adata


def _add_metadata(adata: AnnData, t2g: str, sample: str):
    """Add var and obs metadata."""
    adata.obs["sample"] = sample

    txp2gene = pd.read_table(t2g, header=None, names=["gene_id", "gene_symbol"], usecols=[1, 2])
    txp2gene = txp2gene.drop_duplicates(subset="gene_id").set_index("gene_id")
    adata.var = adata.var.join(txp2gene, how="left")

    # sanitize gene IDs into standard format
    # index are gene IDs and symbols are a column
    adata.var["gene_versions"] = adata.var.index
    adata.var.index = adata.var["gene_versions"].str.split(".").str[0].values
    adata.var_names_make_unique()  # in case user does not use ensembl references, names might not be unique


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a kallisto/bustools count output directory to h5ad")
    parser.add_argument("--input", required=True, help="kallisto counts_* output directory (containing *.mtx files)")
    parser.add_argument("--output", required=True, help="output h5ad file")
    parser.add_argument("--sample", required=True, help="sample id (meta.id)")
    parser.add_argument("--t2g", required=True, help="transcript-to-gene mapping file (t2g.txt)")
    parser.add_argument("--workflow", required=True, choices=["standard", "lamanno", "nac"], help="kb workflow mode")
    args = parser.parse_args()

    # create the directory with the sample name
    os.makedirs(args.sample, exist_ok=True)

    if args.workflow == "standard":
        adata = _mtx_to_adata(
            matrix=glob.glob(f"{args.input}/*.mtx")[0],
            barcodes=glob.glob(f"{args.input}/*.barcodes.txt")[0],
            features=glob.glob(f"{args.input}/*.genes.txt")[0],
        )

    elif args.workflow == "lamanno":
        spliced = _mtx_to_adata(
            matrix=glob.glob(f"{args.input}/spliced*.mtx")[0],
            barcodes=glob.glob(f"{args.input}/spliced*.barcodes.txt")[0],
            features=glob.glob(f"{args.input}/spliced*.genes.txt")[0],
        )
        unspliced = _mtx_to_adata(
            matrix=glob.glob(f"{args.input}/unspliced*.mtx")[0],
            barcodes=glob.glob(f"{args.input}/unspliced*.barcodes.txt")[0],
            features=glob.glob(f"{args.input}/unspliced*.genes.txt")[0],
        )

        # The barcodes of spliced / non-spliced are not necessarily the same.
        # We fill the missing barcodes with zeros
        all_barcodes = list(set(unspliced.obs_names) | set(spliced.obs_names))
        missing_spliced = list(set(all_barcodes) - set(spliced.obs_names))
        missing_unspliced = list(set(all_barcodes) - set(unspliced.obs_names))
        ad_missing_spliced = AnnData(
            X=csr_matrix((len(missing_spliced), spliced.shape[1])),
            obs=pd.DataFrame(index=missing_spliced),
            var=spliced.var,
        )
        ad_missing_unspliced = AnnData(
            X=csr_matrix((len(missing_unspliced), unspliced.shape[1])),
            obs=pd.DataFrame(index=missing_unspliced),
            var=unspliced.var,
        )

        spliced = concat_ad([spliced, ad_missing_spliced], join="outer")[all_barcodes, :]
        unspliced = concat_ad([unspliced, ad_missing_unspliced], join="outer")[all_barcodes, :]

        assert np.all(spliced.var_names == unspliced.var_names)

        adata = AnnData(
            X=spliced.X + unspliced.X,
            layers={"unspliced": unspliced.X, "spliced": spliced.X},
            obs=pd.DataFrame(index=all_barcodes),
            var=pd.DataFrame(index=spliced.var_names),
        )

    elif args.workflow == "nac":
        barcodes = glob.glob(f"{args.input}/*.barcodes.txt")[0]
        features = glob.glob(f"{args.input}/*.genes.txt")[0]

        nascent = _mtx_to_adata(
            matrix=glob.glob(f"{args.input}/*nascent.mtx")[0],
            barcodes=barcodes,
            features=features,
        )
        ambiguous = _mtx_to_adata(
            matrix=glob.glob(f"{args.input}/*ambiguous.mtx")[0],
            barcodes=barcodes,
            features=features,
        )
        mature = _mtx_to_adata(
            matrix=glob.glob(f"{args.input}/*mature.mtx")[0],
            barcodes=barcodes,
            features=features,
        )

        # The barcodes of nascent / mature / ambiguous are not necessarily the same.
        # We fill the missing barcodes with zeros
        all_barcodes = list(set(nascent.obs_names) | set(mature.obs_names) | set(ambiguous.obs_names))
        missing_nascent = list(set(all_barcodes) - set(nascent.obs_names))
        missing_mature = list(set(all_barcodes) - set(mature.obs_names))
        missing_ambiguous = list(set(all_barcodes) - set(ambiguous.obs_names))

        ad_missing_nascent = AnnData(
            X=csr_matrix((len(missing_nascent), nascent.shape[1])),
            obs=pd.DataFrame(index=missing_nascent),
            var=nascent.var,
        )
        ad_missing_ambiguous = AnnData(
            X=csr_matrix((len(missing_ambiguous), ambiguous.shape[1])),
            obs=pd.DataFrame(index=missing_ambiguous),
            var=ambiguous.var,
        )
        ad_missing_mature = AnnData(
            X=csr_matrix((len(missing_mature), mature.shape[1])),
            obs=pd.DataFrame(index=missing_mature),
            var=mature.var,
        )

        nascent = concat_ad([nascent, ad_missing_nascent], join="outer")[all_barcodes, :]
        ambiguous = concat_ad([ambiguous, ad_missing_ambiguous], join="outer")[all_barcodes, :]
        mature = concat_ad([mature, ad_missing_mature], join="outer")[all_barcodes, :]

        assert np.all(nascent.var_names == ambiguous.var_names)
        assert np.all(mature.var_names == ambiguous.var_names)

        adata = AnnData(
            X=nascent.X + ambiguous.X + mature.X,
            layers={"nascent": nascent.X, "ambiguous": ambiguous.X, "mature": mature.X},
            obs=pd.DataFrame(index=all_barcodes),
            var=pd.DataFrame(index=nascent.var_names),
        )

    #
    # out of the conditional: snippet for both standard and non-standard workflows
    #

    # finalize generated adata object
    _add_metadata(adata, t2g=args.t2g, sample=args.sample)
    adata.write_h5ad(args.output)
