#!/usr/bin/env Rscript

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

pkg <- "HuEmProjection"

suppressPackageStartupMessages({
  library(pkg, character.only = TRUE)
})

expected_functions <- c(
  "FunAdjust_shift_variance",
  "FunBigcorS",
  "FunCalCor",
  "FunCalCor_fast",
  "FunCalMNN_each",
  "FunCalSF",
  "FunCompute_correction_vectors",
  "FunMiloCal",
  "FunNWIN",
  "FunPredAnno",
  "FunProjCal",
  "FunRestricted_knn_cdx",
  "FunRestricted_mnn",
  "FunSmooth_gaussian_kernel",
  "FunUStage",
  "FunUS",
  "huem_demo_file",
  "huem_download_models",
  "huem_load_reference",
  "huem_model_dir",
  "huem_model_file",
  "huem_read_cell_groups",
  "huem_read_counts",
  "huem_set_model_dir",
  "project_demo_dataset",
  "project_query_dataset"
)

exports <- getNamespaceExports(pkg)
missing_exports <- setdiff(expected_functions, exports)
if (length(missing_exports) > 0) {
  fail("Missing exported functions: ", paste(missing_exports, collapse = ", "))
}

not_functions <- expected_functions[
  !vapply(expected_functions, function(x) is.function(getExportedValue(pkg, x)), logical(1))
]
if (length(not_functions) > 0) {
  fail("Exported objects are not functions: ", paste(not_functions, collapse = ", "))
}

demo_file <- huem_demo_file("small.demo.counts.tsv.gz")
if (!file.exists(demo_file)) {
  fail("Cannot find bundled demo file: ", demo_file)
}

if (identical(tolower(Sys.getenv("HUEM_TEST_LOAD_REFERENCE")), "true")) {
  ref_env <- huem_load_reference(verbose = FALSE)
  required_reference_objects <- c("ref.umap", "ref_umap_nDim_model", "ref_svm_model")
  missing_reference_objects <- required_reference_objects[
    !vapply(required_reference_objects, exists, logical(1), envir = ref_env, inherits = FALSE)
  ]
  if (length(missing_reference_objects) > 0) {
    fail("Reference loading missing objects: ",
         paste(missing_reference_objects, collapse = ", "))
  }
} else {
  model_dir <- tryCatch(huem_model_dir(), error = function(e) NA_character_)
  if (is.na(model_dir)) {
    message("Model directory not configured; skipping reference model loading test")
  } else {
    model_file <- huem_model_file("model.preload.Rdata")
    if (!file.exists(model_file)) {
      fail("Cannot find model file: ", model_file)
    }
  }
}

cat("OK: HuEmProjection function loading test passed\n")
