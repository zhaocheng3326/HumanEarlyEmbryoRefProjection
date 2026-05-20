\name{projection_core}
\alias{FunAdjust_shift_variance}
\alias{FunBigcorS}
\alias{FunCalCor}
\alias{FunCalCor_fast}
\alias{FunCalMNN_each}
\alias{FunCalSF}
\alias{FunCompute_correction_vectors}
\alias{FunMiloCal}
\alias{FunNWIN}
\alias{FunPredAnno}
\alias{FunProjCal}
\alias{FunRestricted_knn_cdx}
\alias{FunRestricted_mnn}
\alias{FunSmooth_gaussian_kernel}
\alias{FunUStage}
\alias{FunUS}
\title{Low-level projection helpers}
\description{
Low-level functions adapted from the original reference projection scripts.
They are exported for reproducibility and debugging. Most users should call
\code{\link{project_query_dataset}} instead.
}
\usage{
FunMiloCal(temp.counts, temp.counts.meta, temp.cal = TRUE)
FunCalSF(milo_out)
FunCalCor(sf_out)
FunNWIN(sf_out)
FunCalMNN_each(sf_out, ref_name)
FunProjCal(mnn.pairs.list, query.sce.ob, query.sce.cor.out, temp.max,
  D2_umap_model, Dmulti_umap_model, cor.cutoff)
FunPredAnno(predict_out, cor.cutoff = 0.5)
}
\arguments{
\item{temp.counts}{Gene-by-cell count matrix for query cells.}
\item{temp.counts.meta}{Cell metadata with at least \code{cell} and \code{group}.}
\item{temp.cal}{Logical; whether to run Milo neighbourhood aggregation.}
\item{milo_out}{Output from \code{FunMiloCal}.}
\item{sf_out}{Output from \code{FunCalSF}.}
\item{ref_name}{Reference dataset name.}
\item{mnn.pairs.list}{List of MNN pair tables, one entry per reference dataset.}
\item{query.sce.ob}{Query \code{SingleCellExperiment}.}
\item{query.sce.cor.out}{Query-reference correlation table.}
\item{temp.max}{Downsampling support table from \code{FunNWIN}.}
\item{D2_umap_model}{Two-dimensional uwot UMAP model.}
\item{Dmulti_umap_model}{Latent-space uwot UMAP model.}
\item{cor.cutoff}{Correlation cutoff used for filtering and annotation.}
\item{predict_out}{Projection result list from \code{FunProjCal}.}
}
\value{
Function-specific intermediate objects used by the projection workflow.
}
