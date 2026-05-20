# Compatibility shim for older examples that used:
#   source("main.function.R")
#
# The projection code is now provided by the HuEmProjection package. Install the
# package from this directory first with:
#   R CMD INSTALL .

if (!requireNamespace("HuEmProjection", quietly = TRUE)) {
  stop("HuEmProjection is not installed. Run `R CMD INSTALL .` from ",
       "HuEm_stable_reference_projection_tool first.", call. = FALSE)
}

suppressPackageStartupMessages(library(HuEmProjection))
