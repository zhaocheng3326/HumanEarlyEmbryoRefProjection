# HumanEarlyEmbryoRefProjection

Code, Docker instructions, and an R package for projecting query single-cell
gene expression datasets onto a stable human early embryo reference.

The user-facing R package is `HuEmProjection`, located in:

```text
HuEm_stable_reference_projection_tool/
```

Reference model files are distributed separately as a GitHub Release asset to
keep the Git repository small.

## Repository Structure

- `HuEm_stable_reference_projection_tool/`: R package source for `HuEmProjection`.

## Install the R Package

Install from GitHub:

```r
install.packages("remotes")

remotes::install_github(
  "zhaocheng3326/HumanEarlyEmbryoRefProjection",
  subdir = "HuEm_stable_reference_projection_tool"
)
```

Load the package:

```r
library(HuEmProjection)
```

## Download Reference Models

The model archive is attached to the `v0.1.0` release:

```r
huem_download_models(
  model_url = "https://github.com/zhaocheng3326/HumanEarlyEmbryoRefProjection/releases/download/v0.1.0/HuEmProjection-model-v0.1.0.tar.gz"
)
```

Alternatively, manually download and extract the release asset, then point the
package to the extracted `model` directory:

```r
huem_set_model_dir("/path/to/model")
```

## Basic Usage

Input counts should be a tab-separated gene-by-cell matrix with a `Gene` column
containing gene symbols. Optional cell group metadata should be a two-column
tab-separated file without a header: cell barcode and group.

```r
predict_out <- project_query_dataset(
  counts = "/path/to/query.counts.tsv.gz",
  cell_groups = "/path/to/cell.group.meta.tsv",
  run_milor = TRUE,
  cor_cutoff = 0.5
)

head(predict_out$full.anno)
```

Run the small bundled demo:

```r
predict_out <- project_demo_dataset(n = "small", run_milor = FALSE)
table(predict_out$full.anno$pred_EML)
```

## Docker Environment

A Docker image with the compatible R package environment is available:

```text
zhaocheng/shiny-tools-eeptools:v1
```

Run the container while preserving the image's internal packrat environment:

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

## Release Asset

Upload this local file to the GitHub release when making a new release:

```text
HuEmProjection-model-v0.1.0.tar.gz
```

Do not commit model archives or built package tarballs into git.

## Citation

If this repository is helpful, please cite:

Zhao et al., 2024. DOI: `10.1038/s41592-024-02493-2`

Original MNN paper:

Haghverdi et al., 2018. DOI: `10.1038/nbt.4091`

If using miloR aggregation, please also cite:

Dann et al., 2022. DOI: `10.1038/s41587-021-01033-z`

## License

This project is licensed under the GNU General Public License v3.0 or later.
