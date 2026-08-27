# oxo-flow-scrnaseq — Single-cell RNA-seq: alignment, quantification and QC

> ★ Verified · ⇄ Official port of [`nf-core/scrnaseq`](https://github.com/nf-core/scrnaseq) @ `4.2.0` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

[![CI](https://github.com/oxo-flow-community/oxo-flow-scrnaseq/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-scrnaseq/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Run single-cell RNA-seq analysis from raw FASTQ reads to a final MultiQC
report, on **all five upstream aligner branches** of nf-core/scrnaseq 4.2.0
(see the fidelity table for the exact scope):

- `aligner = cellranger` (default): FastQC read QC, Cell Ranger reference
  preparation, per-sample `cellranger count` (with optional BAM output).
- `aligner = simpleaf` (upstream's default aligner): simpleaf/piscem index +
  `simpleaf quant` (alevin-fry), optional QCatch empty-droplet QC report.
- `aligner = kallisto`: `kb ref` + `kb count` (kallisto/bustools, `standard` /
  `lamanno` / `nac` workflows).
- `aligner = star`: STAR `genomeGenerate` + STARsolo per-sample alignment
  (legacy iGenomes index upgrade optional).
- `aligner = cellrangerarc`: cellranger-arc `mkref` + `count` for multiome
  ATAC+GEX data.

All branches converge on the same downstream path: conversion of raw/filtered
matrices to h5ad, CellBender ambient-RNA background removal (not for
cellrangerarc, exactly like upstream), sample-wise h5ad concatenation, and
optional export to Seurat and SingleCellExperiment objects. All results land
under `results/`.

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
  for the combined-h5ad step. For `aligner = cellrangerarc`, five pre-named
  files per sample are required instead (the port's naming contract replaces
  the upstream samplesheet `sample_type`/`fastq_barcode` columns):
  `raw/<sample>_gex_S1_L001_R1_001.fastq.gz`, `raw/<sample>_gex_S1_L001_R2_001.fastq.gz`,
  `raw/<sample>_atac_S1_L001_R1_001.fastq.gz`, `raw/<sample>_atac_S1_L001_R2_001.fastq.gz`,
  `raw/<sample>_atac_S1_L001_R3_001.fastq.gz`.
- **Barcode whitelists**: simpleaf and STARsolo pass a whitelist per protocol
  (upstream maps it automatically in `assets/protocols.json`; the port uses one
  explicit `config.whitelist` path). The four upstream whitelists ship under
  `assets/whitelist/10x_V{1..4}_barcode_whitelist.txt.gz`; set `whitelist` per
  protocol (e.g. `assets/whitelist/10x_V1_barcode_whitelist.txt.gz` for 10XV1).
  Protocols without a whitelist upstream (dropseq/smartseq) work with
  `whitelist = ""`.
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
  `--expect-cells` flag. For `simpleaf`/`kallisto`/`star`, `protocol = "auto"`
  fails fast with upstream's exact error ("Only cellranger supports protocol =
  'auto'…") — pick `10XV1`–`10XV4`, `dropseq` or `smartseq` (kallisto/star only).
- **`.gz` detection**: upstream decides at runtime whether to gunzip (`fasta.endsWith('.gz')`);
  oxo-flow cannot inspect file names, so `fasta_gz`/`gtf_gz` are explicit flags. When a
  plain FASTA/GTF is used, set the flag to `false` and point `fasta_prepared`/`gtf_prepared`
  at the actual file paths.
- **Pre-built index**: with `build_cellranger_index=false` and `transcriptome=<index dir>`,
  `cellranger_mkgtf`/`cellranger_mkref` are skipped and count runs directly (upstream
  `--cellranger_index` behavior). The same pattern exists per aligner:
  `simpleaf_index` (pre-built simpleaf index dir), `kallisto_index` (pre-built
  `kb_ref_out.idx`; prebuilt mode requires `txp2gene`, mirroring upstream's
  assert), `star_index` (pre-built STAR index dir), `cellrangerarc_reference`
  (pre-built cellranger-arc reference). The corresponding build rule is skipped
  and the run rules fail fast with a clear message if the path is absent.
- **Reference chain** (defaults shown): `fasta=refs/refdata.fa.gz` →
  `fasta_prepared=refs/refdata.fa` → `gtf_filtered=refs/refdata_genes.gtf` →
  `gtf_mkgtf=refs/refdata_genes.filtered.gtf` → `transcriptome=refs/cellranger_reference`.
  When enabling `gtf_source_fix`, also set `gtf_mkgtf_input` to the source-fixed file.
- **Aligners** (all switches are `--arg aligner=… --arg protocol=…`):
  `simpleaf` (index from `fasta`+`gtf`, or from `transcript_fasta` which
  requires `txp2gene`; `skip_qcatch`, `simpleaf_umi_resolution`, `remove_doublets`,
  `qcatch_n_partitions`); `kallisto` (`kb_workflow` = standard/lamanno/nac,
  `kb_t1c`/`kb_t2c` overrides); `star` (`star_feature` → `--soloFeatures`,
  `star_index_legacy` = upgrade a legacy iGenomes 2.6.x index,
  `star_ignore_sjdbgtf`, `seq_center` → `--outSAMattrRGline CN`);
  `cellrangerarc` (`cellrangerarc_config` = optional mkref config json,
  auto-generated when empty).
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

Rows cover every upstream process/subworkflow of nf-core/scrnaseq 4.2.0, on all
five aligner branches. Container image strings and conda pins are copied
verbatim from the upstream modules (all pinned, no `latest`). Deviations from
upstream mechanics are called out per row; two structural exclusions and one
multi-lane data limitation remain and are listed at the bottom with evidence.

**Live verification** (2026-08-26/27, tx-ubuntu, engine 0.15.0 + apptainer):
five configurations passed end-to-end — `aligner = cellranger`, `simpleaf`,
`kallisto`, `star` (10X) and `star` with `protocol = dropseq`. The
`cellrangerarc` branch is ported and validate/lint-clean but was not live-run
in this wave.

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
| `CELLRANGER_COUNT` | `cellranger_count` | cellranger 10.0.0 | Same command: reads staged under Cell Ranger naming (`<sample>_S1_L001_R1/R2_001.fastq.gz`), `cellranger count --id <sample> --fastqs fastq_all --transcriptome … --localcores … --localmem … --chemistry <protocol> --create-bam <bool>` + `--expect-cells` when set. The outs tree is then relocated to `results/<aligner>/count/<sample>/outs/` (upstream publishDir `outdir/cellranger/count`). Multi-lane samples (several fastq pairs per sample) are not represented — one pair per sample. |
| `SIMPLEAF_INDEX` | `simpleaf_index` | simpleaf 0.19.5, piscem 0.12.2, alevin-fry 0.11.2, salmon 1.10.3 | Same command (`simpleaf set-paths` + `simpleaf index --threads … [--ref-seq <transcript_fasta> | --fasta … --gtf …] -o simpleaf_index`; `ulimit -n 2048` and `ALEVIN_FRY_HOME` exported). Transcript-fasta mode requires `txp2gene`, mirroring upstream's assert. Output under `refs/simpleaf_index/`. |
| `SIMPLEAF_QUANT` | `simpleaf_quant` | simpleaf 0.19.5, alevin-fry 0.11.2, piscem 0.12.2, salmon 1.10.3 | Same command (`simpleaf quant [--t2g-map …] --chemistry <protocol-mapped> --index … --reads1/2 … --resolution cr-like --output simpleaf_quant --threads … --anndata-out --unfiltered-pl <whitelist>`; cell filtering hardcoded to `unfiltered-pl` upstream → input_type is always raw). Protocol→chemistry mapping and per-protocol whitelist mirror `assets/protocols.json`. Output `results/<aligner>/<sample>/simpleaf_quant/af_quant/`. |
| `QCATCH` | `qcatch` | qcatch 0.2.12 | Same command (`qcatch --input <af_quant dir> --output qcatch [--chemistry 10X_3p_v2/v3/v4] --save_filtered_h5ad --export_summary_table [--n_partitions] [--remove_doublets --visualize_doublets]`), same output renames (`QCatch_report.html` → `<sample>_qcatch_report.html`, `filtered_quants.h5ad` → `<sample>_filtered_quants.h5ad`, `summary_table.csv` → `<sample>_metrics_summary.csv`). Chemistry mapping for 10XV2-4 only, exactly like upstream. |
| `KALLISTOBUSTOOLS_REF` | `kallistobustools_ref_standard`, `kallistobustools_ref_velocity` | kb-python 0.28.2 | Same commands: standard `kb ref -i … -g … -f1 cdna.fa --workflow standard`; non-standard workflows (`lamanno`/`nac`) add `-f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow <mode>`. Mutual exclusion is a `when` on `kb_workflow` (upstream picks the command by the same variable). Outputs under `refs/kallisto/`. |
| `KALLISTOBUSTOOLS_COUNT` | `kallistobustools_count` | kb-python 0.28.2 | Same command (`kb count -t … -i … -g … [-c1 …] [-c2 …] -x <technology> --workflow <kb_workflow> --filter -o <sample>.count -m <memory.toGiga()-1>G reads`); technology mapping for 10XV1-4/DROPSEQ/SMARTSEQ mirrors upstream. Ext.args `--workflow … --filter` applied. Output `results/<aligner>/<sample>.count/`. |
| `STAR_GENOMEGENERATE` | `star_genomegenerate` | star 2.7.11b, samtools 1.21, gawk 5.1.0 | Same command: `samtools faidx` + gawk SAindexNbases heuristic from the `.fai` (14 cap), `--runMode genomeGenerate --genomeDir … --genomeFastaFiles … --sjdbGTFfile … --runThreadN … --genomeSAindexNbases … --limitGenomeGenerateRAM <memory-100000000>`. Output under `refs/star_index/`. |
| `STAR_GENOMEPARAMS_UPGRADE` | `star_genomeparams_upgrade` | gawk 5.3.1 | Same script: symlink the legacy index files, awk-rewrite `genomeParameters.txt` (versionGenome 20201 → 2.7.4a, append genomeType/Full + genomeTransformType/None + genomeTransformVCF/-), move to `refs/star_index_upgraded/`. Fires only when `star_index` is set and `star_index_legacy=true` (upstream `isStarIndexLegacy`). |
| `STAR_ALIGN` | `star_align` | star 2.7.10b | Same command: reads passed REVERSE first, `--readFilesCommand zcat --runDirPerm All_RWX --outWigType bedGraph --twopassMode Basic --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM <memory bytes>`, `--soloCBwhitelist` with the same `.gz`→uncompress handling (protocols without an upstream whitelist — dropseq/smartseq — get the literal `--soloCBwhitelist None`, STAR's required spelling for no whitelist; live-found: omitting the flag aborts with "--soloCBwhitelist is not defined"), `--soloType`/`--soloUMIlen` per protocol (10XV1/2→10, 10XV3/4→12, dropseq/smartseq→none), `--soloCellFilter CellRanger2.2 <expected_cells> 0.99 10` when set, `--soloFeatures <star_feature>` (+Velocyto publish rename). Solo.out tsv/mtx files gzipped in-place before publish, exactly like upstream. Index selection: upgraded legacy > user `star_index` > built. |
| `CELLRANGERARC_MKGTF` | `cellrangerarc_mkgtf` | cellranger-arc 2.0.2 | Same command as upstream (`cellranger-arc mkgtf` with the three biotype filters). Runs only when `build_cellranger_index=true`. |
| `CELLRANGERARC_MKREF` | `cellrangerarc_mkref` | cellranger-arc 2.0.2 | Same flow: auto-generated mkref config json (`organism: "refdata"`, `genome: ["<prefix>_reference"]`, `input_fasta`, `input_gtf`) or user `cellrangerarc_config`, then `cellranger-arc mkref --config=config --nthreads …`. Output at `refs/cellrangerarc_reference/` (the config's `genome` name; `cellrangerarc_reference` can point at an existing reference to skip building). |
| `CELLRANGERARC_COUNT` | `cellrangerarc_count` | cellranger-arc 2.0.2 | Same flow: fastqs staged under `fastqs/`, 2-row `lib.csv` (Gene Expression / Chromatin Accessibility), `cellranger-arc count --id=<sample> --libraries=… --reference=… --localcores … --localmem … [--expect-cells]`, outs tree relocated to `results/<aligner>/count/<sample>/outs/`. Deviation: the upstream samplesheet's `sample_type`/`fastq_barcode` columns are replaced by a fixed file-naming contract — see the sample-data requirements. |
| `MTX_TO_H5AD` | `mtx_to_h5ad_{raw,filtered,simpleaf,kallisto_raw,kallisto_filtered,star_raw,star_filtered}` | scanpy 1.10.2 / pandas / anndata | Same template scripts per aligner (`mtx_to_h5ad_cellranger.py` — read_10x_h5, also used for cellrangerarc exactly like upstream's `(input_aligner in ['cellranger','cellrangerarc','cellrangermulti']) ? 'cellranger' : input_aligner`; `mtx_to_h5ad_simpleaf.py`; `mtx_to_h5ad_kallisto.py` with standard/lamanno/nac branches; `mtx_to_h5ad_star.py` incl. the Velocyto layer code, dead upstream, kept verbatim), one rule per aligner×input_type. Raw/filtered gating mirrors the upstream channels: simpleaf emits only raw (upstream hardcodes `unfiltered-pl`); star/kallisto filtered conversions skip for protocols without a whitelist (dropseq/smartseq) — the upstream filtered dirs don't exist there. |
| `CELLBENDER_REMOVEBACKGROUND` | `cellbender_removebackground` | cellbender 0.3.2 | Same command `TMPDIR=. cellbender remove-background --cpu-threads … --estimator-multiple-cpu --input … --output <sample>.h5` (no `--cuda`: GPU profile is out of scope). Full output file set moved to `results/<aligner>/<sample>/cellbender_removebackground/`. Skipped for `cellrangerarc`, exactly like upstream. |
| `ANNDATA_BARCODES` | `anndata_barcodes` | anndata 0.11.4 / pandas | Same template script (barcode CSV → subset → write), same output name `<sample>_cellbender_filter_matrix.h5ad`. Skipped for `cellrangerarc` with the upstream subworkflow. |
| `CONCAT_H5AD` | `concat_h5ad_filtered`, `concat_h5ad_cellbender_filter`, `concat_h5ad_raw` | scanpy 1.10.2 | Same template script (`ad.concat(label="sample", merge="unique", index_unique="_")` + samplesheet join on `sample`). Upstream runs one process per input_type; the port has one rule per input_type. Gating mirrors the upstream channels: `filtered` skips for simpleaf (no filtered h5ads), star+dropseq and kallisto+dropseq (no filtered dirs), and smartseq (no whitelist); `raw` runs only when `skip_cellbender=true` or aligner=cellrangerarc (raw superseded by the CellBender-filtered h5ad otherwise). |
| `ANNDATAR_CONVERT` | `anndatar_convert_{filtered,cellbender_filter,raw}` + `anndatar_convert_combined_{…}` | anndataR 1.0.2, SeuratObject 5.5.0, SingleCellExperiment 1.32.0 | Same R template (read_h5ad → `as_Seurat()`/`as_SingleCellExperiment()` → saveRDS). Six rules: per sample and per combined h5ad, per input_type; type gating mirrors the concat rules. Upstream `dir.create(<sample>)` calls and versions.yml writing dropped (output dirs are pre-created by the engine; versions are recorded in `collect_versions`). |
| `softwareVersionsToYAML` + `collectFile` | `collect_versions` | — | Writes the same file `results/pipeline_info/nf_core_scrnaseq_software_mqc_versions.yml` consumed by MultiQC. Content is the port's pinned versions (upstream collates live tool versions from a channel topic, which has no oxo-flow equivalent); since containers are pinned, the recorded versions equal the executed ones. Only the active aligner's block is emitted, like the upstream channel topic. |
| `paramsSummaryMultiqc` + methods description | `workflow_summary`, `methods_description` | — | New default-ON rules producing the summary/methods MultiQC YAMLs from the copied-verbatim `assets/methods_description_template.yml` (the `${…}` placeholders are filled at render time; upstream fills them from the Nextflow workflow object, which has no oxo-flow equivalent). They run in the default config, so a single-sample default dry-run plan (`oxo-flow dry-run main.oxoflow --samples first:1`, as exercised by test/run.sh) shows 21 rules executing (19 baseline + these 2); with the two bundled samples the plan shows 29 running instances — documented new default behavior. |
| `MULTIQC` | `multiqc` | multiqc 1.34 | Same command (`multiqc --force [--title] --config <assets/multiqc_config.yml> .`) with inputs staged flat like the module's `stageAs '?/*'`; the input union covers the active aligner's web summaries/logs (FastQC + cellranger web_summary + simpleaf quants.h5ad + STAR Log.final.out). Default `assets/multiqc_config.yml` copied verbatim from upstream. |

**Not ported (with reasons):**

| Upstream branch | Reason |
|---|---|
| `aligner = cellrangermulti` (multiome VDJ/Ab-seq/CRO) | Structural: upstream feeds per-sample, per-modality fastq groups into `cellranger multi` via channel branching (`groupTuple` + EMPTY-file injection for missing modalities) and three index channels (GEX/VDJ/cellranger_multi_barcodes). A fixed rule input signature cannot express a variable number of modality input sets per sample — no oxo-flow analogue exists. |
| `PIPELINE_COMPLETION` (email/notification) | Structural: `workflow.onComplete`/`onError` hooks do not exist in the oxo-flow engine; the failure email would have to be sent by an external wrapper. |
| `skip_cellranger_renaming` (multi-lane samples) | One fastq pair per sample is supported; the staging rename hard-codes lane `L001`. |

**Other deliberate deviations** (documented per row above): FastQC is skipped
for `cellrangerarc` (five reads per sample cannot fit one static input
pattern; upstream runs it on all of them); the arc samplesheet columns are a
file-naming contract; `workflow_summary`/`methods_description` are new
default-ON rules; simpleaf/star/kallisto accept one explicit `whitelist` path
instead of upstream's automatic per-protocol mapping.

**Live-root-caused fixes** (engine 0.15.0, tx-ubuntu): tool-facing
threads/cores use `{effective_threads}` (rules declare 12/6 CPUs; a 4-core box
would oversubscribe); every container spec is quay.io-qualified (bare
`biocontainers/...` resolves to Docker Hub, not the pinned quay.io registry);
directory-moving rules `rm -rf` the engine-precreated output parent before
`mv` (the parent exists, so `mv` would nest the tree inside itself); the STAR
index nbases heuristic truncates with `int()` and a 14 cap (the 52kb fixture
genome rounded up to 7 where STAR requires 6 — "may cause seg-fault"); the
fixture GTF gives every gene two exons with an intron (single-exon
transcripts crash simpleaf's grangers intron pass: polars "invalid series
dtype: expected List, got null"); the fixture genome is padded to ~52kb (STAR
double-frees on the original 1.9 kb genome).

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
