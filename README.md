# oxo-flow-scrnaseq — Single-cell RNA-seq: alignment, quantification and QC

> ★ Verified · ⇄ Official port of [`nf-core/scrnaseq`](https://github.com/nf-core/scrnaseq) @ `4.2.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

[![CI](https://github.com/oxo-flow-community/oxo-flow-scrnaseq/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-scrnaseq/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Run single-cell RNA-seq analysis on 10x Genomics data from raw FASTQ reads to a
final MultiQC report: FastQC read QC, Cell Ranger reference preparation and
per-sample `cellranger count` (alignment and quantification, with optional BAM
output), conversion of the raw and filtered 10x matrices to h5ad, CellBender
ambient-RNA background removal, sample-wise h5ad concatenation, and optional
export to Seurat and SingleCellExperiment objects. The workflow covers the
Cell Ranger execution path of nf-core/scrnaseq (see the fidelity table for the
exact scope); other aligners and assays such as ATAC or multiome are not
included. All results land under `results/`.

## Installation

### 1. Install oxo-flow

This workflow requires **oxo-flow >= 0.12.0**. The recommended way is the
prebuilt release binary:

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternatively, install via conda:

```bash
conda install -c bioconda oxo-flow-cli
```

Note: the conda package may lag behind releases; binaries for other platforms
are available on the [releases page](https://github.com/Traitome/oxo-flow/releases).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-scrnaseq.git
cd oxo-flow-scrnaseq
```

### 3. Requirements

- **Reference data** (you provide; see the `[config]` section of `main.oxoflow`):
  a genome FASTA (optionally gzipped) and a gene-annotation GTF (optionally
  gzipped) — defaults `refs/refdata.fa.gz` and `refs/refdata.gtf.gz`. The Cell
  Ranger reference is built from them by the workflow into
  `refs/cellranger_reference` unless you set `build_cellranger_index = false`
  and point `transcriptome` at a pre-built index.
- **Sample data**: one raw 10x FASTQ pair per sample as
  `raw/<sample>_R1.fastq.gz` / `raw/<sample>_R2.fastq.gz` (one pair per sample),
  plus a `samplesheet.csv` with columns `sample,fastq_1,fastq_2,protocol,expected_cells`
  for the combined-h5ad step.
- **Compute**: up to 12 CPUs and 72 GB RAM per rule (Cell Ranger `mkref` /
  `count`); h5ad conversion, CellBender and concatenation rules use 6 CPUs /
  36 GB, QC/prep rules 1–2 CPUs / 6–12 GB. Per-sample rules run concurrently,
  so peak usage scales with the number of samples and `-j`.

> **Cell Ranger memory sizing.** `mkref`/`count` auto-size `--localmem`
> to 2/3 of the *actually free* physical memory (`/proc/meminfo`
> MemAvailable, 1 GB floor) — never the machine's effective total, which
> counts swap: cellranger's job manager waits forever when `--localmem`
> exceeds the free RAM (live: `Need 6 GB ... (2.6 GB available)` looped
> for hours on a 3.7 GB box). Real runs still want 10x's documented 8 GB
> floor; the tiny test fixtures run on smaller machines with swap
> absorbing overflow. Set `cellranger_localmem` to force a value.
- **Tool delivery**: containers with pinned images — every rule pins the exact
  upstream nf-core container string (no `latest`), so Docker (or Singularity)
  is required at runtime. Conda alternatives for local runs ship in `envs/`,
  except the Cell Ranger rules, which are docker-only (the upstream modules
  refuse conda profiles).

## Usage

```bash
# 1. install oxo-flow (see Installation)
# 2. prepare data: raw/<sample>_R1.fastq.gz / _R2.fastq.gz, a samplesheet.csv,
#    and reference files (refs/refdata.fa.gz, refs/refdata.gtf.gz) — a synthetic 3-gene test reference is committed under refs/ (swap in a real FASTA/GTF for real data)
# 3. preview the plan
oxo-flow dry-run main.oxoflow
# 4. run
oxo-flow run main.oxoflow -j 8
# 5. run a subset
oxo-flow run main.oxoflow -t multiqc --samples first:2
```

### Configuration

- **Default config knobs** mirror upstream params: `protocol=auto` → `--chemistry auto`,
  `save_align_intermeds=true` → `--create-bam true`, `expected_cells=""` → no
  `--expect-cells` flag.
- **`.gz` detection**: upstream decides at runtime whether to gunzip (`fasta.endsWith('.gz')`);
  oxo-flow cannot inspect file names, so `fasta_gz`/`gtf_gz` are explicit flags. When a
  plain FASTA/GTF is used, set the flag to `false` and point `fasta_prepared`/`gtf_prepared`
  at the actual file paths.
- **Pre-built index**: with `build_cellranger_index=false` and `transcriptome=<index dir>`,
  `cellranger_mkgtf`/`cellranger_mkref` are skipped and count runs directly (upstream
  `--cellranger_index` behavior).
- **Reference chain** (defaults shown): `fasta=refs/refdata.fa.gz` →
  `fasta_prepared=refs/refdata.fa` → `gtf_filtered=refs/refdata_genes.gtf` →
  `gtf_mkgtf=refs/refdata_genes.filtered.gtf` → `transcriptome=refs/cellranger_reference`.
  When enabling `gtf_source_fix`, also set `gtf_mkgtf_input` to the source-fixed file.
- **Conda alternatives**: the workflow pins the exact upstream container images; conda
  env yamls ship in `envs/` for local runs (cellranger rules are docker-only — the
  upstream modules refuse conda profiles). `python-igraph`/`leidenalg` are not pinned
  upstream and were resolved at port time (2026-08-15) to conda-forge releases.
- **`samplesheet`** (for `CONCAT_H5AD`): user CSV with the upstream columns
  (`sample,fastq_1,fastq_2,protocol,expected_cells`); see `test/fixtures/samplesheet.csv`.
  The `raw/` symlink in the repo root points at `test/fixtures/raw/` so `validate`/`dry-run`
  see the fixture reads as existing inputs; replace the fixtures with real data for runs.

### Reference genome

The pipeline expects `refs/refdata.fa.gz` + `refs/refdata.gtf.gz` and derives
the rest: gunzip → `refs/refdata.fa|.gtf`, GTF gene filter →
`refs/refdata_genes.gtf`, mkgtf → `refs/refdata_genes.filtered.gtf`, mkref →
`refs/cellranger_reference/`. The committed synthetic reference is a
3-protein-coding-gene chr1 (exons 100-500, 800-1300, 1500-1900) sized so
mkref/count run on small machines; `test/fixtures/generate_fixtures.py` draws
the test reads from those exons.

## Source

Ported from **[nf-core/scrnaseq](https://github.com/nf-core/scrnaseq)**, version
`4.2.0` (commit `3fc17b4f971a89e47c88337de71d0e777ffad8cc`, MIT). Created
2026-08-15; this workflow may lag behind upstream releases. Upstream
attribution and licensing details are in [NOTICE.md](NOTICE.md); the upstream
MIT license is retained verbatim in [LICENSE.upstream](LICENSE.upstream).

## Fidelity

Rows cover every upstream process/subworkflow involved in the default path; the
`aligner` branches not ported are listed at the bottom with reasons. Container image
strings and conda pins are copied verbatim from the upstream modules (all pinned, no
`latest`).

| Upstream process/rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| `PIPELINE_INITIALISATION` (samplesheet check) | sample source + `config.samplesheet` | — | Samplesheet (`sample, fastq_1, fastq_2, protocol, expected_cells`) maps to `[[sample_groups]]`; `expected_cells` column → `config.expected_cells`; per-sample `protocol` column is informational (chemistry comes from `--protocol`). Schema checks are enforced by the port's fixtures + README contract. |
| `FASTQC` | `fastqc` | fastqc 0.12.1 | Identical command: `printf … \| while read; ln -s` staging loop, `fastqc --quiet --threads N --memory <12G/N clamped 100-10000>`. Published under `results/fastqc/` (upstream default publishDir). `--memory` computed in-shell from the process_low 12G/2 cpus. |
| `GUNZIP` (as `GUNZIP_FASTA`) | `gunzip_fasta` | gzip 1.13 | Identical command (`gzip -cd <fasta> > <out>`). Runs only when `config.fasta_gz` (upstream decides by `.endsWith('.gz')` at runtime — port uses an explicit flag, see Gotchas). |
| `GUNZIP` (as `GUNZIP_GTF`) | `gunzip_gtf` | gzip 1.13 | Same as above for the GTF. |
| `GTF_GENE_FILTER` | `gtf_gene_filter` | python 3.9 | Same bundled script `filter_gtf_for_genes_in_genome.py`, same flags (`--gtf --fasta -o`); output name `<fasta_stem>_genes.gtf` is `config.gtf_filtered`. |
| `GAWK` (as `GTF_SOURCE_FIX`) | `gtf_source_fix` | gawk 5.3.1 | Same awk program (`FS=OFS="\t"`, source-field spaces→underscores, output suffix `gtf`). Off by default, exactly like upstream (only fires for iGenomes entries flagged `gtf_source_has_spaces`). |
| `CELLRANGER_MKGTF` | `cellranger_mkgtf` | cellranger 10.0.0 | Same command incl. the three `--attribute=gene_biotype:` filters. Runs only when `build_cellranger_index=true` (mirrors upstream `if (!cellranger_index)`). |
| `CELLRANGER_MKREF` | `cellranger_mkref` | cellranger 10.0.0 | Same command (`--genome=… --fasta=… --genes=… --localcores --localmem --nthreads`). `--genome` is `config.transcriptome` (default `refs/cellranger_reference`) instead of a bare workdir name — same reference name, path relocated to the workflow tree. |
| `CELLRANGER_COUNT` | `cellranger_count` | cellranger 10.0.0 | Same command: reads staged under Cell Ranger naming (`<sample>_S1_L001_R1/R2_001.fastq.gz`), `cellranger count --id <sample> --fastqs fastq_all --transcriptome … --localcores … --localmem … --chemistry <protocol> --create-bam <bool>` + `--expect-cells` when set. The outs tree is then relocated to `results/cellranger/count/<sample>/outs/` (upstream publishDir `outdir/cellranger/count`). Multi-lane samples (several fastq pairs per sample) are not represented — one pair per sample. |
| `MTX_TO_H5AD` | `mtx_to_h5ad_raw`, `mtx_to_h5ad_filtered` | scanpy 1.10.2 / pandas / anndata | Same template script `mtx_to_h5ad_cellranger.py` (read_10x_h5, gene_symbols, gene_ids index, version-stripped gene ids, `var_names_make_unique`), one rule per input_type (upstream runs the process once per raw/filtered channel). Outputs `<sample>_{raw,filtered}_matrix.h5ad`. |
| `CELLBENDER_REMOVEBACKGROUND` | `cellbender_removebackground` | cellbender 0.3.2 | Same command `TMPDIR=. cellbender remove-background --cpu-threads … --estimator-multiple-cpu --input … --output <sample>.h5` (no `--cuda`: GPU profile is out of scope). Full output file set moved to `results/cellranger/<sample>/cellbender_removebackground/`. |
| `ANNDATA_BARCODES` | `anndata_barcodes` | anndata 0.11.4 / pandas | Same template script (barcode CSV → subset → write), same output name `<sample>_cellbender_filter_matrix.h5ad`. |
| `CONCAT_H5AD` | `concat_h5ad_filtered`, `concat_h5ad_cellbender_filter`, `concat_h5ad_raw` | scanpy 1.10.2 | Same template script (`ad.concat(label="sample", merge="unique", index_unique="_")` + samplesheet join on `sample`). Upstream runs one process per input_type; the port has one rule per input_type. `concat_h5ad_raw` only runs when `skip_cellbender=true`, mirroring upstream's channel replacement (raw is superseded by the CellBender-filtered h5ad). |
| `ANNDATAR_CONVERT` | `anndatar_convert_{filtered,cellbender_filter,raw}` + `anndatar_convert_combined_{…}` | anndataR 1.0.2, SeuratObject 5.5.0, SingleCellExperiment 1.32.0 | Same R template (read_h5ad → `as_Seurat()`/`as_SingleCellExperiment()` → saveRDS). Six rules: per sample and per combined h5ad, per input_type; type gating mirrors the concat rules. Upstream `dir.create(<sample>)` calls and versions.yml writing dropped (output dirs are pre-created by the engine; versions are recorded in `collect_versions`). |
| `softwareVersionsToYAML` + `collectFile` | `collect_versions` | — | Writes the same file `results/pipeline_info/nf_core_scrnaseq_software_mqc_versions.yml` consumed by MultiQC. Content is the port's pinned versions (upstream collates live tool versions from a channel topic, which has no oxo-flow equivalent); since containers are pinned, the recorded versions equal the executed ones. |
| `MULTIQC` | `multiqc` | multiqc 1.34 | Same command (`multiqc --force [--title] --config <assets/multiqc_config.yml> .`) with inputs staged flat like the module's `stageAs '?/*'`. Default `assets/multiqc_config.yml` copied verbatim from upstream. |

**Not ported (with reasons):**

| Upstream branch | Reason |
|---|---|
| `aligner = simpleaf` (upstream default aligner; alevin-fry) | Not on the ported cellranger default path; requires qcatch protocol handling and a simpleaf index/whitelist matrix. |
| `aligner = kallisto` (kallisto/bustools, incl. `smartseq`/`dropseq` protocols) | Not on the ported cellranger path. |
| `aligner = star` (STARsolo, incl. `smartseq`/`dropseq` protocols) | Not on the ported cellranger path. |
| `aligner = cellrangerarc` (ATAC, `fastq_barcode` samplesheet column) | Not on the ported cellranger path. |
| `aligner = cellrangermulti` (multiome: VDJ/Ab-seq/CRO, `cellranger_multi_barcodes`) | Not on the ported cellranger path. |
| `PIPELINE_COMPLETION` (email/notification) | nf-core boilerplate, out of scope. |
| `paramsSummaryMultiqc` / methods-description MultiQC inputs | nf-core reporting boilerplate; MultiQC still aggregates FastQC + Cell Ranger + versions. |
| `skip_cellranger_renaming` (multi-lane samples) | One fastq pair per sample is supported; the staging rename hard-codes lane `L001`. |

## Test

Run the acceptance test — `validate`, `lint`, and a `dry-run` plan check with
the bundled fixtures:

```bash
bash test/run.sh
```

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md); the upstream MIT license is retained verbatim in
[LICENSE.upstream](LICENSE.upstream).

## Community

https://oxo-flow-community.github.io/
