# HuEmProjection

R package for projecting query single-cell gene expression matrices onto the
human early embryo reference.

## Install in the Docker environment

From the repository root on the host:

```sh
docker run -it \
  -v /local/folder/HumanEarlyEmbryoRefProjection/:/home/docker/Install/project \
  --name eeptools \
  zhaocheng/shiny-tools-eeptools:v1 \
  /bin/bash
```

Inside the container:

```sh
mkdir -p /tmp/Rlib
cd /home/docker/Install
R CMD INSTALL -l /tmp/Rlib project/HuEm_stable_reference_projection_tool
Rscript -e '.libPaths(c("/tmp/Rlib", .libPaths())); source("project/HuEm_stable_reference_projection_tool/tests/test_function_loading.R")'
```

For this checkout, replace `/local/folder/HumanEarlyEmbryoRefProjection/` with:

```sh
/Users/cheng.zhao/Downloads/HumanEarlyEmbryoRefProjection/
```

## Basic use

```r
library(HuEmProjection)

# One-time setup after downloading the Release asset:
huem_set_model_dir("/path/to/model")

predict_out <- project_query_dataset(
  counts = "/path/to/query.counts.tsv.gz",
  cell_groups = "/path/to/cell.group.meta.tsv",
  run_milor = TRUE,
  cor_cutoff = 0.5
)

head(predict_out$full.anno)
```

Alternatively, download the model archive directly from a GitHub Release:

```r
huem_download_models(
  model_url = "https://github.com/<OWNER>/<REPO>/releases/download/v0.1.0/HuEmProjection-model-v0.1.0.tar.gz"
)
```

The release asset was generated outside the GitHub upload directory at:

```text
../HuEmProjection-model-v0.1.0.tar.gz
```

The counts file must be a tab-separated matrix with genes as rows, cells as
columns, and a `Gene` column containing gene symbols. The cell group file should
have two tab-separated columns without a header: cell barcode and group.

Demo data bundled with the package can be run with:

```r
predict_out <- project_demo_dataset(n = "n1000")
```

To test from the package source directory that the package namespace, exported
functions, and bundled data paths load correctly without running a full
projection:

```sh
Rscript tests/test_function_loading.R
```

To also force loading the external reference models:

```sh
HUEM_TEST_LOAD_REFERENCE=true Rscript tests/test_function_loading.R
```

## Useful files

- `R/package-api.R`: high-level user-facing API.
- `R/projection-core.R`: original projection functions adapted for package data paths.
- `inst/extdata/demo`: bundled demo datasets.
- `../HuEmProjection-model-v0.1.0.tar.gz`: upload this file to a GitHub Release, but do not commit it.
- `inst/scripts/projecting_query_dataset.R`: installed copy of the demo script.
