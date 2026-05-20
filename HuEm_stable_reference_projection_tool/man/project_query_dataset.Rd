\name{project_query_dataset}
\alias{project_query_dataset}
\alias{project_demo_dataset}
\alias{huem_read_counts}
\alias{huem_read_cell_groups}
\alias{huem_demo_file}
\alias{huem_download_models}
\alias{huem_model_dir}
\alias{huem_model_file}
\alias{huem_set_model_dir}
\alias{huem_load_reference}
\title{Project query cells onto the human early embryo reference}
\usage{
project_query_dataset(counts, cell_groups = NULL, run_milor = TRUE,
  cor_cutoff = 0.5, reference_names = c("SPH2016", "D3post",
  "Meistermann_2021", "CS7", "nBGuo", "Yan2013"))

project_demo_dataset(n = c("small", "n750", "n1000"), run_milor = TRUE,
  cor_cutoff = 0.5)
}
\arguments{
\item{counts}{A gene-by-cell count matrix/data frame, or a path to a tab-separated file with a \code{Gene} column.}
\item{cell_groups}{Optional data frame with \code{cell} and \code{group} columns, or a path to a two-column tab-separated file without a header.}
\item{run_milor}{Logical; whether to run Milo neighbourhood aggregation.}
\item{cor_cutoff}{Correlation cutoff used for filtering and annotation.}
\item{reference_names}{Reference datasets to use for MNN pairing.}
\item{n}{Bundled demo dataset size.}
}
\description{
Runs the complete offline projection workflow: optional Milo aggregation,
stable size-factor calculation, correlation to the reference, MNN pairing,
UMAP projection, SVM annotation, and pseudotime transfer.
}
\details{
Model files are not bundled in the GitHub source package. Configure them with
\code{huem_set_model_dir()}, the \code{HUEMPROJECTION_MODEL_DIR} environment
variable, or \code{huem_download_models()} after uploading the model archive to
a GitHub Release.
}
\value{
A list containing projected UMAP coordinates, annotation tables, MNN support,
raw metadata, and intermediate matrices.
}
