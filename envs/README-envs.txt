The workflow in main.oxoflow uses the exact upstream container images (docker = the
nf-core module container strings). The *.yaml files here are conda alternatives for
local runs without docker. Rules running cellranger (cellranger_mkgtf, cellranger_mkref,
cellranger_count) are docker-only: the upstream nf-core modules explicitly refuse the
conda/mamba profiles, so no conda env is shipped for cellranger.

python-igraph and leidenalg (scanpy.yaml) are not pinned upstream
("conda-forge::python-igraph conda-forge::leidenalg"); they were resolved at port time
(2026-08-15) to the conda-forge releases current for scanpy 1.10.2.
