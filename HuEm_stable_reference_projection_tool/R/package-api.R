huem_read_counts <- function(file, gene_column = "Gene") {
  data.table::fread(file, header = TRUE, sep = "\t") %>%
    tibble::column_to_rownames(gene_column)
}

huem_read_cell_groups <- function(file) {
  read.delim(file, header = FALSE, stringsAsFactors = FALSE, sep = "\t") %>%
    tibble::as_tibble() %>%
    stats::setNames(c("cell", "group"))
}

project_query_dataset <- function(counts,
                                  cell_groups = NULL,
                                  run_milor = TRUE,
                                  cor_cutoff = 0.5,
                                  reference_names = c(
                                    "SPH2016", "D3post", "Meistermann_2021",
                                    "CS7", "nBGuo", "Yan2013"
                                  )) {
  if (is.character(counts) && length(counts) == 1) {
    counts <- huem_read_counts(counts)
  }
  counts <- as.data.frame(counts)

  if (is.null(cell_groups)) {
    counts_meta <- tibble::tibble(cell = colnames(counts), group = "None")
  } else if (is.character(cell_groups) && length(cell_groups) == 1) {
    counts_meta <- huem_read_cell_groups(cell_groups)
  } else {
    counts_meta <- tibble::as_tibble(cell_groups)
  }

  if (!all(c("cell", "group") %in% colnames(counts_meta))) {
    stop("cell_groups must contain columns named 'cell' and 'group'.",
         call. = FALSE)
  }
  missing_cells <- setdiff(counts_meta$cell, colnames(counts))
  if (length(missing_cells) > 0) {
    stop("cell_groups contains cells that are absent from counts: ",
         paste(utils::head(missing_cells, 5), collapse = ", "),
         call. = FALSE)
  }

  counts_meta <- counts_meta %>%
    dplyr::mutate(pj = "query", EML = "query")

  milo_out <- FunMiloCal(counts, counts_meta, temp.cal = run_milor)
  sf_out <- FunCalSF(milo_out)
  sf_out$query.sce.cor.out <- FunCalCor(sf_out)
  sf_out$query.sce.cor.out.mean <- sf_out$query.sce.cor.out %>%
    tidyr::gather(ref_cell, cor, -query_cell) %>%
    dplyr::group_by(query_cell) %>%
    dplyr::top_n(20, cor) %>%
    dplyr::summarise(cor_top_mean = mean(cor), .groups = "drop")
  sf_out$NWIN <- FunNWIN(sf_out)

  mnn_pairs <- list()
  for (ref_name in reference_names) {
    message(ref_name)
    mnn_pairs[[ref_name]] <- FunCalMNN_each(sf_out, ref_name)
  }

  ref_env <- huem_load_reference()
  predict_out <- FunProjCal(
    mnn_pairs,
    query.sce.ob = sf_out$query.sce.ob,
    query.sce.cor.out = sf_out$query.sce.cor.out,
    temp.max = sf_out$NWIN$temp.max,
    D2_umap_model = ref_env$ref.umap,
    Dmulti_umap_model = ref_env$ref_umap_nDim_model,
    cor.cutoff = cor_cutoff
  )
  predict_out$HS <- milo_out$HS
  predict_out$raw.meta <- milo_out$raw.meta
  predict_out$query.sce.cor.out.mean <- sf_out$query.sce.cor.out.mean
  FunPredAnno(predict_out, cor.cutoff = cor_cutoff)
}

project_demo_dataset <- function(n = c("small", "n750", "n1000"),
                                 run_milor = TRUE,
                                 cor_cutoff = 0.5) {
  n <- match.arg(n)
  count_file <- switch(
    n,
    small = huem_demo_file("small.demo.counts.tsv.gz"),
    n750 = huem_demo_file("n750.demo.counts.tsv.gz"),
    n1000 = huem_demo_file("n1000.demo.counts.tsv.gz")
  )
  group_file <- switch(
    n,
    small = NULL,
    n750 = huem_demo_file("n750.cell.group.meta.tsv"),
    n1000 = huem_demo_file("n1000.cell.group.meta.tsv")
  )
  project_query_dataset(
    counts = count_file,
    cell_groups = group_file,
    run_milor = run_milor,
    cor_cutoff = cor_cutoff
  )
}
