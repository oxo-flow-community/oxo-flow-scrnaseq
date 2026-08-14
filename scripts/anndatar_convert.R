#!/usr/bin/env Rscript

# Port of nf-core/scrnaseq 4.2.0 modules/local/anndatar_convert/templates/anndatar_convert.R
# Logic is identical (read h5ad, convert to Seurat and SingleCellExperiment, saveRDS each);
# the Nextflow "${...}" substitutions are replaced by CLI arguments. The upstream template
# also writes a versions.yml file; the port records versions in collect_versions instead.
# SPDX-License-Identifier: MIT (upstream) — see NOTICE.md

# load libraries
library(anndataR)
library(SeuratObject)
library(SingleCellExperiment)

# parse CLI arguments (--input <h5ad> --output-seurat <rds> --output-sce <rds>)
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
    idx <- grep(paste0("^", name, "="), args)
    stopifnot(length(idx) == 1)
    sub(paste0("^", name, "="), "", args[idx])
}
h5ad <- get_arg("--input")
out_seurat <- get_arg("--output-seurat")
out_sce <- get_arg("--output-sce")

# read input
adata <- read_h5ad(h5ad)

# convert to Seurat
obj <- adata$as_Seurat()

# save files
saveRDS(obj, file = out_seurat)

# convert to SingleCellExperiment
obj <- adata$as_SingleCellExperiment()

# save files
saveRDS(obj, file = out_sce)
