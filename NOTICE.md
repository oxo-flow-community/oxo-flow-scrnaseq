oxo-flow-scrnaseq
Copyright (c) 2026 oxo-flow-community

This pipeline is a port of nf-core/scrnaseq
(https://github.com/nf-core/scrnaseq), version 4.2.0
(commit 3fc17b4f971a89e47c88337de71d0e777ffad8cc), authored by
The nf-core/scrnaseq team.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---------------------------------------------------------------------
Upstream license

This port is derived from nf-core/scrnaseq under the MIT license.
The upstream LICENSE must be included **verbatim** in this
repository: fetch it from the upstream repository at the ported
4.2.0 (commit 3fc17b4f971a89e47c88337de71d0e777ffad8cc) and place it
at LICENSE.upstream. (Apache-2.0 §4(d): attribution notices from the
Source form must be retained.)

The helper scripts in scripts/ are ports of files distributed with
nf-core/scrnaseq (MIT license, The nf-core/scrnaseq team):
  - filter_gtf_for_genes_in_genome.py  (bin/, from nf-core/rnaseq)
  - mtx_to_h5ad_cellranger.py          (modules/local/mtx_to_h5ad/templates/)
  - barcodes.py                        (modules/nf-core/anndata/barcodes/templates/)
  - concat_h5ad.py                     (modules/local/concat_h5ad/templates/)
  - anndatar_convert.R                 (modules/local/anndatar_convert/templates/)
---------------------------------------------------------------------
