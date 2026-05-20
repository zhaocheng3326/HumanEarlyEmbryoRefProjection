.huem_reference_env <- new.env(parent = emptyenv())
.huem_state <- new.env(parent = emptyenv())
.huem_state$reference_loaded <- FALSE

.huem_extdata_file <- function(...) {
  path <- system.file("extdata", ..., package = "HuEmProjection")
  if (!nzchar(path)) {
    stop("Cannot find package data file: ", file.path(...), call. = FALSE)
  }
  path
}

huem_demo_file <- function(filename) {
  .huem_extdata_file("demo", filename)
}

huem_set_model_dir <- function(path) {
  if (!dir.exists(path)) {
    stop("Model directory does not exist: ", path, call. = FALSE)
  }
  options(HuEmProjection.model_dir = normalizePath(path, mustWork = TRUE))
  invisible(getOption("HuEmProjection.model_dir"))
}

huem_model_dir <- function(must_exist = TRUE) {
  candidates <- c(
    getOption("HuEmProjection.model_dir", NULL),
    Sys.getenv("HUEMPROJECTION_MODEL_DIR", unset = NA_character_),
    system.file("extdata", "model", package = "HuEmProjection"),
    file.path(tools::R_user_dir("HuEmProjection", "data"), "model")
  )
  candidates <- candidates[nzchar(candidates) & !is.na(candidates)]
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) > 0) {
    return(normalizePath(existing[[1]], mustWork = TRUE))
  }
  if (!must_exist) {
    return(normalizePath(candidates[[length(candidates)]], mustWork = FALSE))
  }
  stop(
    "HuEmProjection model files were not found.\n",
    "Set a local model directory with `huem_set_model_dir('/path/to/model')`, ",
    "set the HUEMPROJECTION_MODEL_DIR environment variable, or download models ",
    "from the GitHub Release with `huem_download_models(model_url = '...')`.",
    call. = FALSE
  )
}

huem_model_file <- function(filename) {
  path <- file.path(huem_model_dir(), filename)
  if (!file.exists(path)) {
    stop("Cannot find model file: ", path, call. = FALSE)
  }
  path
}

huem_download_models <- function(model_url = getOption("HuEmProjection.model_url",
                                                      Sys.getenv("HUEMPROJECTION_MODEL_URL", "")),
                                 dest_dir = huem_model_dir(must_exist = FALSE),
                                 force = FALSE) {
  if (!nzchar(model_url)) {
    stop(
      "No model_url was provided. After uploading the model archive to a ",
      "GitHub Release, call `huem_download_models(model_url = '<release asset url>')` ",
      "or set option HuEmProjection.model_url / env HUEMPROJECTION_MODEL_URL.",
      call. = FALSE
    )
  }

  required_file <- file.path(dest_dir, "model.preload.Rdata")
  if (file.exists(required_file) && !force) {
    return(invisible(normalizePath(dest_dir, mustWork = TRUE)))
  }

  parent_dir <- dirname(dest_dir)
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  archive <- tempfile(fileext = ".tar.gz")
  utils::download.file(model_url, archive, mode = "wb", quiet = FALSE)
  if (dir.exists(dest_dir) && force) {
    unlink(dest_dir, recursive = TRUE)
  }
  utils::untar(archive, exdir = parent_dir)

  if (!file.exists(required_file)) {
    nested <- file.path(parent_dir, "model", "model.preload.Rdata")
    if (file.exists(nested) && normalizePath(dirname(nested)) != normalizePath(dest_dir, mustWork = FALSE)) {
      if (dir.exists(dest_dir)) {
        unlink(dest_dir, recursive = TRUE)
      }
      file.rename(file.path(parent_dir, "model"), dest_dir)
    }
  }
  if (!file.exists(required_file)) {
    stop("Downloaded archive did not create the expected model files in: ",
         dest_dir, call. = FALSE)
  }

  huem_set_model_dir(dest_dir)
  invisible(normalizePath(dest_dir, mustWork = TRUE))
}

huem_load_reference <- function(force = FALSE, verbose = FALSE) {
  if (.huem_state$reference_loaded && !force) {
    return(.huem_reference_env)
  }

  rm(list = ls(.huem_reference_env), envir = .huem_reference_env)
  .huem_reference_env$ref.umap <- uwot::load_uwot(
    file = huem_model_file("stable.ref.fastMNN.umap.model")
  )
  .huem_reference_env$ref_umap_nDim_model <- uwot::load_uwot(
    file = huem_model_file("stable.ref.fastMNN.umap.nDimN.model")
  )
  .huem_reference_env$ref_svm_model <- readRDS(
    huem_model_file("stable.ref.fastMNN.umap.nDim.svm.model.rds")
  )
  load(huem_model_file("model.preload.Rdata"),
       envir = .huem_reference_env,
       verbose = verbose)

  .huem_state$reference_loaded <- TRUE
  .huem_reference_env
}

.huem_bind_reference <- function(envir = parent.frame()) {
  list2env(as.list(huem_load_reference()), envir = envir)
  invisible(envir)
}

#' general function
FunRestricted_mnn <- function(left.data,  right.data, k=20) {
  # modified from https://rdrr.io/bioc/batchelor/src/R/MNN_tree.R#sym-.restricted_mnn
  k1 <- min(k,  nrow(left.data))
  k2 <- min(k, nrow(right.data))
  pairs <- findMutualNN(left.data, right.data, k1=k1, k2=k2)
  pairs$first <- pairs$first
  pairs$second <- pairs$second
  pairs
}

FunRestricted_knn_cdx <- function(left.data, right.data, k=20) {
  k <-  min(k,  nrow(left.data),nrow(right.data))
  Rindex.out <-  data.frame(left=rownames(left.data)) %>% tibble::as_tibble() %>% mutate(left=as.vector(left)) %>% tibble::rowid_to_column("Rindex")
  knns.out <- BiocNeighbors::queryKNN(left.data,right.data,k=k,get.index=TRUE,get.distance=FALSE)  %>% as.data.frame()
  rownames(knns.out) <- rownames(right.data)
  knns.out <- knns.out %>% tibble::rownames_to_column("right") %>% tibble::as_tibble() %>% gather(index,Rindex,-right) %>% left_join(Rindex.out ,by="Rindex") %>% select(-c(Rindex,index)) %>% select(left,right)
  return(knns.out)
}


FunBigcorS <- function (x, y = NULL,size = 2000, verbose = TRUE,...) {
  #quick way to calculated correlation , modified from post
  STR <- "Correlation"
  
  if (!is.null(y) & NROW(x) != NROW(y))
    stop("'x' and 'y' must have compatible dimensions!")
  NCOL <- ncol(x)
  if (!is.null(y))
    YCOL <- NCOL(y)
  REST <- NCOL%%size
  LARGE <- NCOL - REST
  NBLOCKS <- NCOL%/%size
  if (is.null(y))
    resMAT <- ff::ff(vmode = "double", dim = c(NCOL, NCOL))
  else resMAT <- ff::ff(vmode = "double", dim = c(NCOL, YCOL))
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0)
    GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  SPLIT <- split(1:NCOL, GROUP)
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  if (!is.null(y))
    COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
  timeINIT <- proc.time()
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (is.null(y)) {
      if (verbose)
        message("bigcor: ", sprintf("#%d: %s of Block %s and Block %s (%s x %s) ... ",
                                    i, STR, COMB[1], COMB[2], length(G1), length(G2)))
      RES <- cor(x[, G1], x[, G2], method="spearman")
      resMAT[G1, G2] <- RES
      resMAT[G2, G1] <- t(RES)
    }
    else {
      if (verbose)
        message("bigcor: ", sprintf("#%d: %s of Block %s and 'y' (%s x %s) ... ",
                                    i, STR, COMB[1], length(G1), YCOL))
      RES <- cor(x[, G1], y,method="spearman")
      resMAT[G1, ] <- RES
    }
    if (verbose) {
      timeNOW <- proc.time() - timeINIT
      message("bigcor: ", round(timeNOW[3], 2), " sec\n")
    }
    gc()
  }
  return(resMAT)
}
FunCompute_correction_vectors <- function(data1, data2, mnn1, mnn2, tdata2, sigma)  {     
  #modified from https://rdrr.io/github/LTLA/batchelor/src/R/mnnCorrect.R#sym-.compute_correction_vectors
  vect <- data1[mnn1,,drop=FALSE] - data2[mnn2,,drop=FALSE]
  averaged <- scuttle::sumCountsAcrossCells(t(vect), DataFrame(ID=mnn2), average=TRUE)
  cell.vect <- FunSmooth_gaussian_kernel(assay(averaged, withDimnames=FALSE), averaged$ID-1L, tdata2, sigma)
  t(cell.vect)
}
FunAdjust_shift_variance <- function(data1, data2, vect, sigma2=0.1) {
  # modified from https://rdrr.io/bioc/batchelor/src/R/RcppExports.R
  restrict1=seq_len(ncol(data1))-1
  restrict2=seq_len(ncol(data2))-1
  scaling <- .Call('_batchelor_adjust_shift_variance', PACKAGE = 'batchelor', data1, data2, vect, sigma2, restrict1, restrict2) %>% pmax(1)
  scaling * vect
}
FunSmooth_gaussian_kernel <-  function(averaged, index, mat, sigma2) {
  # modified from https://rdrr.io/bioc/batchelor/src/R/RcppExports.R
  .Call('_batchelor_smooth_gaussian_kernel', PACKAGE = 'batchelor', averaged, index, mat, sigma2)
}
FunUS = function(x) {
  x=unique(unlist(strsplit(x,":")))
  x=x[order(x)]
  return(data.frame(pred_EML=paste(x,collapse=":"),npred_EML=length(x)))
}

FunUStage = function(x) {
  x=unique(unlist(strsplit(x,":")))
  x=x[order(x)]
  return(data.frame(stage=paste(x,collapse=":")))
}

FunMiloCal <- function(temp.counts,temp.counts.meta, temp.cal=TRUE) {
  if (ncol(temp.counts) >2000) {
    temp.cal=TRUE
  }
  if (ncol(temp.counts)  < 200) {
    temp.cal=FALSE
  }
  if (temp.cal) {
    temp.list <- list()
    for (sp in unique(temp.counts.meta$group)) {
      print(sp)
      temp.M <- temp.counts.meta %>% filter(group==sp)
      if (nrow(temp.M) > 20) {
        data_sce <- SingleCellExperiment(list(counts=as.matrix(temp.counts[,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% scran::computeSumFactors()
        data_sce <- logNormCounts(data_sce)
  
        data_sce <- runPCA(data_sce, ncomponents=50)
        data_sce <- runUMAP(data_sce)
        
        set.seed(456)
        data_milo <- Milo(data_sce)
        #plotUMAP(data_milo)
        data_milo <- buildGraph(data_milo, k = 10, d = 50) ### default
        data_milo <- makeNhoods(data_milo, prop = 0.15, k = 21, d=30, refined = TRUE) ## 0.15 is better
        data_milo <- buildNhoodGraph(data_milo)
        
        data_graph <- nhoodGraph(data_milo)
        
        data_Nhood <- nhoods(data_milo)
        colnames(data_Nhood) <- colnames(data_milo)[as.numeric(colnames(data_Nhood))]
        rownames(data_Nhood) <- colnames(data_milo)
        
        data_Nhood <- data_Nhood %>% as.matrix() %>% as.data.frame() %>%  tibble::rownames_to_column("Scell") %>% gather("Hcell",n,-Scell) %>% filter(n > 0)  %>% tibble::as_tibble() %>% select(-n) %>% arrange(Hcell)
        data_Nhood <- data_Nhood %>% bind_rows(data_Nhood %>% select(Hcell) %>% unique() %>% mutate(Scell=Hcell)) %>% unique()
        data_edge <- as_edgelist(data_graph) %>% as.matrix() %>% as.data.frame() %>% setNames(c("nhood","nedge")) %>% mutate_all(as.numeric) %>% tibble::as_tibble()
        data_edge$nhCell <- colnames(data_milo)[data_edge$nhood]
        data_edge$neCell <- colnames(data_milo)[data_edge$nedge]
        data_edge <- data_edge %>% select(-c(nhood,nedge))
        
        #' merge the raw counts for the neighbourhood cells 
        temp.cell.list <- data_Nhood  %>% split(.,.$Hcell) %>% lapply(function(x) {x$Scell})
        temp.counts.merge <- list()
        for (n in names(temp.cell.list)) {
          temp.cells <- temp.cell.list[[n]]
          temp.counts.merge[[n]] <- as.data.frame(rowSums(temp.counts[,temp.cells,drop=F]) ) %>% setNames(n)
        }
        temp.counts.merge <- temp.counts.merge %>% do.call("cbind",.)
        
        sub.milo_out <- list()
        sub.milo_out$counts <- temp.counts.merge
        sub.milo_out$HS <- data_Nhood
        sub.milo_out$edge <- data_edge
        temp.list[[sp]] <-  sub.milo_out
      }
    }
    milo_out <- list()
    milo_out$counts <- temp.list %>% lapply(function(x) x$counts) %>% do.call("bind_cols",.)
    milo_out$HS  <- temp.list %>% lapply(function(x) x$HS) %>% do.call("bind_rows",.)
    milo_out$edge <- temp.list %>% lapply(function(x) x$edge) %>% do.call("bind_rows",.)
    milo_out$raw.meta <- temp.counts.meta
    milo_out$info <- "neighborhood calculation"
    print(milo_out$info)
  }else{
    milo_out <- list()
    milo_out$counts <- temp.counts
    milo_out$HS  <- temp.counts.meta %>% mutate(Scell=cell,Hcell=cell) %>% dplyr::select(Scell,Hcell)
    milo_out$edge <- milo_out$HS
    milo_out$raw.meta <- temp.counts.meta
    milo_out$info <- "passing to next step"
    print(milo_out$info)
  }
  return(milo_out)
}

FunCalSF<- function(milo_out) {
  .huem_bind_reference()
  
  set.seed(ref.seed)
  load(huem_model_file("ref.sce.ob.related.Rdata"), verbose = TRUE)
  ref.raw.data <- readRDS(huem_model_file("ref.raw.data.Rdata"))
  ref.sce.ob <- readRDS(huem_model_file("ref.sce.ob.rds"))
  
  temp.counts <- milo_out$counts
  temp.M <- milo_out$raw.meta %>% select(cell,group,EML,pj) %>% filter(cell %in% (colnames(temp.counts)))
  temp.counts <- temp.counts[rownames(temp.counts) %>% intersect(ref.sel.expG),temp.M$cell]
  sf_out <- list()
  
  
  
  #' align the reads name
  if (nrow(temp.counts) < length(ref.sel.expG)) {
    temp.counts.extra <- matrix(0,length(ref.sel.expG)-nrow(temp.counts),ncol(temp.counts)) %>% as.data.frame()
    rownames(temp.counts.extra) <- ref.sel.expG%>% setdiff(rownames(temp.counts))
    colnames(temp.counts.extra) <- colnames(temp.counts)
    temp.counts <- temp.counts %>% bind_rows(temp.counts.extra)
  }
  temp.counts <- temp.counts[ref.sel.expG,]
  
  if (nrow(temp.M) >50) {
    k_num=30
  }else{
    k_num=5
  }
  
  #' get stable sizeFactor
  #' modified code from https://rdrr.io/github/LTLA/batchelor/src/R/multiBatchNorm.R
  ref.sce.ob$query <- SingleCellExperiment(list(counts=as.matrix(temp.counts[,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% scran::computeSumFactors()
  query.index <- match("query",names(ref.sce.ob))
  
  ref.sce.stats$query <- list()
  #sf <-  sizeFactors(ref.sce.ob$query)
  sf <-   SingleCellExperiment(list(counts=as.matrix(temp.counts[,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% scran::computeSumFactors() %>% sizeFactors()
  sf <- sf/mean(sf)
  ave <- scuttle::calculateAverage(ref.sce.ob$query@assays@data$counts,  size.factors=sf )
  ref.sce.stats$query$sf <- sf
  ref.sce.stats$query$ave <- ave
  
  first <- ref.smallest
  second <- query.index 
  first.ave <- ref.sce.stats[[first]]$ave
  first.sum <- sum(first.ave)
  second.ave <- ref.sce.stats[[second]]$ave
  second.sum <- sum(second.ave)
  # Mimic calcAverage(cbind(first.ave, second.ave)).
  grand.mean <- (first.ave/first.sum + second.ave/second.sum)/2 * (first.sum + second.sum)/2
  min.mean <- 1
  keep <- grand.mean >= min.mean
  # Computing it twice for exactly equal results when order of batches is rearranged.
  kept.f <- first.ave[keep]
  kept.s <- second.ave[keep]
  curratio2 <- median(kept.f/kept.s)
  sf_out$rescale.out <- curratio2
  
  ref.rescaling.updated <- ref.rescaling
  ref.rescaling.updated[query.index] <- curratio2
  left.raw.data <- ref.raw.data
  query.sce.ob <- ref.sce.ob$query
  sizeFactors(query.sce.ob) <- ref.sce.stats$query$sf/ref.rescaling.updated[query.index]
  right.data <- ((scuttle::logNormCounts(query.sce.ob,center.size.factors=FALSE) %>% logcounts())[ref.mnnVGs,]  %>% cosineNorm( mode="all"))$matrix %>% t()
  

  
  sf_out$right.data <- right.data
  sf_out$query.sce.ob <- query.sce.ob
  sf_out$left.raw.data <- left.raw.data
  sf_out$k_num <- k_num
  return(sf_out)
}

FunCalCor_fast <- function(sf_out) {
  query.sce.cor=FunBigcorS(t(sf_out$right.data),t(sf_out$left.raw.data))
  query.sce.cor <- query.sce.cor[1:nrow(query.sce.cor),1:ncol(query.sce.cor)]
  rownames(query.sce.cor) <- rownames(sf_out$right.data)
  colnames(query.sce.cor) <-  rownames(sf_out$left.raw.data)
  query.sce.cor.out <- query.sce.cor %>% as.data.frame()%>% tibble::rownames_to_column("query_cell") #%>% gather(ref_cell,cor,-query_cell) %>% group_by()
  return(query.sce.cor.out)
}
FunCalCor <- function(sf_out) {
  query.sce.cor=cor(t(sf_out$right.data),t(sf_out$left.raw.data),method="spearman")
  query.sce.cor <- query.sce.cor[1:nrow(query.sce.cor),1:ncol(query.sce.cor)]
  rownames(query.sce.cor) <- rownames(sf_out$right.data)
  colnames(query.sce.cor) <-  rownames(sf_out$left.raw.data)
  query.sce.cor.out <- query.sce.cor %>% as.data.frame()%>% tibble::rownames_to_column("query_cell") #%>% gather(ref_cell,cor,-query_cell) %>% group_by()
  return(query.sce.cor.out)
}

FunNWIN <- function(sf_out) {
  right.data <- sf_out$right.data
  k_num <- sf_out$k_num
  
  
  win.size <- 200
  if (nrow( right.data) > win.size+50) {
    nwin <- round(nrow( right.data)/win.size+0.5)
    win.index <- rep(1:nwin,each=win.size)[1:nrow( right.data)]
    nrep <- 5
    set.seed(1);
    win.sel.index <- matrix(nrow=length(win.index),ncol=nrep)
    colnames(win.sel.index) <- paste0("rep",1:nrep)
    rownames(win.sel.index) <- rownames( right.data)
    win.sel.index[,"rep1"] <- sample(win.index)
    for (n in 2:nrep) {
      set.seed(n);
      win.sel.index[,paste0("rep",n)] <- sample(win.index)
    }
    temp.max <- win.sel.index %>% as.data.frame() %>% tibble::rownames_to_column("right") %>% gather(rep,nwin,-right)  %>% group_by(right) %>% summarise(nC=n())
  }else{
    nwin=1
    win.sel.index <- matrix(nrow=nrow( right.data),ncol=1,1)
    colnames(win.sel.index) <- paste0("rep",1:1)
    rownames(win.sel.index) <- rownames( right.data)
    temp.max <- data.frame(cell=rownames(sf_out$right.data))  %>% mutate(right=as.vector(cell),nC=1) %>% select(right,nC)
  }
  temp.ds.list <- list()
  for (nc in 1:ncol(win.sel.index)) {
    for (nw in (1:nwin)) {
      right.data.sel <- rownames(win.sel.index)[win.sel.index[,nc]==nw]
      if (length(right.data.sel) < win.size & nwin>1) {
        set.seed(100*nc+nw);right.data.extra.sel <- sample(rownames(win.sel.index) %>% setdiff(right.data.sel),win.size-length(right.data.sel))
        right.data.sel <- c(right.data.sel,right.data.extra.sel) ### fulfill the last win
      }
      temp.ds.list[[paste(nc,nw,sep=":")]] <- data.frame(right=right.data.sel,nc=nc,nw=nw,SID=paste(nc,nw,sep=":")) %>% tibble::as_tibble() %>% mutate_all(as.vector)
    }
  }
  temp.out <- list()
  
  
  temp.out$win.sel.index <- temp.ds.list %>% do.call("bind_rows",.)
  temp.out$temp.max <- temp.out$win.sel.index %>% group_by(right) %>% summarise(nC=n_distinct(SID)) %>%ungroup()
  temp.out$nwin  <- nwin
  temp.out$win.size  <- win.size
  return(temp.out)
}

FunCalMNN_each <- function(sf_out,ref_name) { 
  .huem_bind_reference()
  #' calculate the mnn pairs
  temp.ref.sce.ob <- readRDS(huem_model_file(paste0("ref.sce.ob.", ref_name, ".rds")))
  
  s <- ref_name
  k_num <- sf_out$k_num
  query.sce.ob <- sf_out$query.sce.ob
  temp.max <- sf_out$NWIN$temp.max
  win.sel.index <- sf_out$NWIN$win.sel.index 
  nwin <- sf_out$NWIN$nwin  
  win.size <- sf_out$NWIN$win.size
  
  mnn.pairs.list <- list()
  left.data.each <-  ((scuttle::logNormCounts(temp.ref.sce.ob,center.size.factors=FALSE) %>% logcounts())[ref.mnnVGs,]  %>% cosineNorm( mode="all"))$matrix %>% t()
  
  for (nc_nw in unique(win.sel.index$SID)) {
      right.data.sel <- win.sel.index %>% filter(SID==nc_nw) %>% pull(right)
      print(nc_nw)
      right.data.each <- ((scuttle::logNormCounts(query.sce.ob[,right.data.sel],center.size.factors=FALSE) %>% logcounts())[ref.mnnVGs,]  %>% cosineNorm( mode="all"))$matrix %>% t()
      mnn.sets <- FunRestricted_mnn(left.data.each, right.data.each,k=k_num)
      s1 <- mnn.sets$first
      s2 <- mnn.sets$second
      mnn.pairings.mdx <- DataFrame(left=s1, right=s2)
      mnn.pairs.list[[nc_nw]] <-  mnn.pairings.cdx <- DataFrame(left=rownames(left.data.each)[s1], right=rownames(right.data.each)[s2]) %>% tibble::as_tibble() %>% mutate(ref_pj=s,nc_nw=nc_nw)
  }
  return(mnn.pairs.list)
}

FunProjCal <- function(mnn.pairs.list,query.sce.ob,query.sce.cor.out,temp.max,D2_umap_model,Dmulti_umap_model,cor.cutoff) {
  .huem_bind_reference()
  load(huem_model_file("ref.mnn.correct.cal.Rdata"), verbose = TRUE)
  predict_out <- list()
  mnn.pairs.out <- mnn.pairs.list %>% lapply(function(x) {do.call("bind_rows",x)})  %>% do.call("bind_rows",.) %>% group_by(left,right) %>% summarise(nC=n()) %>% inner_join(temp.max ,by=c("right","nC"))%>% select(-nC)  %>% select(left,right) %>% rename(ref_cell=left,query_cell=right)   %>% inner_join(ref.meta %>% mutate(ref_cell=cell,ref_EML=cluster_EML,ref_pj=pj) %>% select(ref_cell,ref_EML,ref_pj),by="ref_cell")%>% select(ref_cell,query_cell,ref_EML,ref_pj) %>% unique() %>% left_join(query.sce.cor.out%>% gather(ref_cell,cor,-query_cell) %>% tibble::as_tibble(),by=c("ref_cell","query_cell")) %>% ungroup()

  #' the cells with conflict MNN pairs
  temp.check.cells <-  mnn.pairs.out %>% filter(ref_EML %in% c("TE","CTB","EVT","STB")) %>% pull(query_cell) %>% intersect(mnn.pairs.out %>% filter(ref_EML %in% c("Early_Amnion","Late_Amnion","Amnion" )) %>% pull(query_cell)) %>% unique()
  mnn.pairs.out.cor <- mnn.pairs.out  %>% group_by(query_cell,ref_EML) %>% top_n(20,cor) %>% summarise(cor_top_mean=mean(cor)) %>% ungroup()
  mnn.pairs.filter.out <-  mnn.pairs.out %>% filter(!query_cell %in% temp.check.cells) %>% left_join(mnn.pairs.out.cor,by=c("query_cell","ref_EML"))  %>% filter(cor > cor.cutoff | cor_top_mean > cor.cutoff  ) %>% select(ref_cell,query_cell) %>% unique()
  
  mnn.pairs.ref.EML.support <- mnn.pairs.out %>% ungroup() %>% select(-c(ref_pj,ref_cell)) %>% unique() %>% rename(sub_pred_EML=ref_EML)%>% mutate(sub_pred_EML=recode(sub_pred_EML,"2-4 cell"="Z4cell","Zygote"="Z4cell")) %>% left_join(ref.mergeCT,by="sub_pred_EML") %>% mutate(mergeCT=ifelse(is.na(mergeCT),sub_pred_EML,mergeCT)) %>% unique() #%>% filter(!query_cell %in% temp.check.cells)
  
  query.exp <-  (scuttle::logNormCounts(query.sce.ob,center.size.factors=FALSE) %>% logcounts())
  query.exp.l2 <- pmax(1e-8,  query.exp %>% cosineNorm(mode="l2norm",subset.row = ref.mnnVGs))
  query.exp <- scuttle::normalizeCounts(query.exp, size_factors=query.exp.l2, center_size_factors=FALSE, log=FALSE)[ref.mnnVGs,]
  query.exp.centered <- query.exp -ref.grand.centers
  w <- sqrt(ncol(query.exp.centered ))
  query.exp.scaled <- query.exp.centered/w
  
  query.svd.output <-  as.matrix(crossprod(query.exp.centered, ref.svd.out$u))
  set.seed(10)
   
  left.data <- ref.mnn.out
  right.data <- query.svd.output
  right.data <- batchelor:::.orthogonalize_other(right.data, NULL, ref.to.add )
  right.data <- batchelor:::.center_along_batch_vector(right.data, ref.overall.batch)
  trans.right <- t(right.data)
  mnn.sets  <- list()
  mnn.sets$first <- match(mnn.pairs.filter.out$ref_cell,rownames(left.data))
  mnn.sets$second <- match(mnn.pairs.filter.out$query_cell,rownames(right.data))
  s1 <- mnn.sets$first
  s2 <- mnn.sets$second
  correction.in <- FunCompute_correction_vectors(left.data, right.data, s1, s2, trans.right, sigma=0.1)
  correction.in <- FunAdjust_shift_variance(t(left.data),t(right.data),correction.in)
  right.data <- right.data + correction.in
  query.pca <- left.data %>%rbind(right.data) 
  
  set.seed(ref.seed)
  qr.umap <- umap_transform(query.pca[,1:ref.nPC],D2_umap_model) 
  
  temp.umap <- qr.umap  %>% as.data.frame() %>% setNames(c("UMAP_1","UMAP_2")) %>% tibble::rownames_to_column("cell") %>% inner_join(ref.meta %>% mutate(pj=paste0("ref_",pj)) %>% bind_rows(colData(query.sce.ob) %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% rename(devTime=group)) %>% select(cell,devTime,EML,pj),by="cell")
  
  #' rotation umap
  temp.umap <- temp.umap %>% rename(UMAP_2=UMAP_1, UMAP_1=UMAP_2)
  
  #' latent space
  set.seed(ref.seed)
  qr_latent <- umap_transform(query.pca[,1:ref.nPC],ref_umap_nDim_model) %>% as.data.frame()%>% setNames(paste("Dim",1:ref.svm.nDim,sep="_"))   %>% tibble::rownames_to_column("cell") %>% tibble::as_tibble()
  
  temp_svm_input <- qr_latent %>% filter(cell %in% rownames(query.svd.output)) %>% tibble::column_to_rownames("cell")
  
  predict_out$umap <- temp.umap
  predict_out$query.pca <- query.pca
  predict_out$mnn.pairs.ref.EML.support <- mnn.pairs.ref.EML.support
  predict_out$svm_input <- temp_svm_input
  predict_out$query.sce.cor.out.mean <- NULL
  #predict_out$cor.dis <- mnn.pairs.out  %>% ggplot()+geom_boxplot(mapping=aes(x=ref_EML,y=cor))+ theme_classic() + theme(axis.text.x=element_text(angle = 90))
  #predict_out$ref.mergeCT <- ref.mergeCT 
  return(predict_out)
}


FunPredAnno <- function(predict_out,cor.cutoff=0.5) {
  .huem_bind_reference()
  mnn.pairs.ref.EML.support <- predict_out$mnn.pairs.ref.EML.support
  query.sce.cor.out.mean <- predict_out$query.sce.cor.out.mean
  temp_svm_input <- predict_out$svm_input
  #ref.mergeCT <- predict_out$ref.mergeCT
  
  temp_pred <- list()
  for (temp_sel_EML in names(ref_svm_model)) {
    temp_pred[[temp_sel_EML]] <- predict(ref_svm_model[[temp_sel_EML]]$model, temp_svm_input, type = "prob") %>% tibble::as_tibble() %>% select(-other)%>% setNames(c("prop")) %>% mutate(query_cell=rownames(temp_svm_input),sub_pred_EML=paste0(temp_sel_EML,".prop"),ref_kappa=ref_svm_model[[temp_sel_EML]]$perf %>% pull(Kappa),ref_acc=ref_svm_model[[temp_sel_EML]]$perf %>% pull(Accuracy))
  }
  temp_pred_out <- temp_pred %>% do.call("bind_rows",.) %>% mutate(sub_pred_EML=recode(sub_pred_EML,"EightCell.prop"="8 cell.prop","AxMes.prop"="Axial Mes.prop")) %>% group_by(query_cell) %>% top_n(1,prop)%>% top_n(1,ref_kappa)%>% top_n(1,ref_acc) %>% rename(prediction_score_max=prop) %>% mutate(sub_pred_EML=gsub(".prop","",sub_pred_EML))%>% select(query_cell,sub_pred_EML,prediction_score_max) %>% mutate(pred_EML=recode(sub_pred_EML,"Late_Hypoblast"="Hypoblast","Early_Hypoblast"="Hypoblast","Early_EPI"="Epiblast","Late_EPI"="Epiblast","Early_Amnion"="Amnion","Late_Amnion"="Amnion"))%>% inner_join(temp_pred %>% do.call("bind_rows",.) %>% select(query_cell,sub_pred_EML,prop) %>% spread(sub_pred_EML,prop),by="query_cell")%>% ungroup()
  temp_pred_support_EML <- temp_pred_out %>% select(query_cell,sub_pred_EML) %>% left_join(ref.mergeCT,by="sub_pred_EML") %>% mutate(mergeCT=ifelse(is.na(mergeCT),sub_pred_EML,mergeCT))  %>% inner_join(mnn.pairs.ref.EML.support %>% rename(related_pred_EML=sub_pred_EML),by=c("query_cell","mergeCT")) %>% select(-mergeCT) %>% unique()
  temp_pred_out <- temp_pred_out %>% filter(sub_pred_EML %in% c(temp_pred_support_EML$related_pred_EML, temp_pred_support_EML$sub_pred_EML)) %>% bind_rows(temp_pred_out %>% filter(!sub_pred_EML %in% c(temp_pred_support_EML$related_pred_EML, temp_pred_support_EML$sub_pred_EML)) %>% mutate(sub_pred_EML="Ambiguous",pred_EML="Ambiguous"))
  
  
  temp_pred_out <- temp_pred_out %>% full_join( query.sce.cor.out.mean,by="query_cell") %>% mutate(pred_EML=ifelse(cor_top_mean < cor.cutoff, "low_cor", pred_EML))%>% mutate(sub_pred_EML=ifelse(cor_top_mean < cor.cutoff, "low_cor", sub_pred_EML)) %>% mutate(pred_EML=ifelse(prediction_score_max< 0.5 , "Ambiguous", pred_EML))%>% mutate(sub_pred_EML=ifelse(prediction_score_max< 0.5 , "Ambiguous", sub_pred_EML))
  
  
  # temp.umap %>% APPJ_UMAP(p)
  # temp_pred_out %>% ggplot()+geom_boxplot(mapping=aes(x=sub_pred_EML,y=prediction_score_max))+ theme(axis.text.x=element_text(angle = 90))
  # temp.umap %>% rows_update( temp_pred_out %>% tibble::as_tibble() %>% mutate(EML=sub_pred_EML,cell=query_cell) %>% select(cell,EML),by="cell") %>% APPJ_UMAP(p) 
  predict_out$anno <-  temp_pred_out %>% select(query_cell,pred_EML,sub_pred_EML,prediction_score_max)%>%  mutate(stage=ifelse(sub_pred_EML %in% stage.cate$Early,"Early","uncertained")) %>% mutate(stage=ifelse(sub_pred_EML %in% stage.cate$Late,"Late",stage)) 
  
  #' expand to full cells
  temp_pred_ac_out <- temp_pred_out%>% select(query_cell,pred_EML,sub_pred_EML,prediction_score_max)  %>% rename(Hcell=query_cell) %>% inner_join(predict_out$HS,by ="Hcell") %>% ungroup()%>% select(-Hcell) %>% rename(query_cell=Scell) %>% group_by(query_cell,pred_EML,sub_pred_EML) %>% summarise(n=n(),prediction_score_max=max(prediction_score_max)) %>% arrange(query_cell)
  
  temp1 <- temp_pred_ac_out %>% filter(!pred_EML %in% c("low_cor","Ambiguous")) %>% arrange(query_cell,pred_EML,sub_pred_EML) %>% group_by(query_cell) %>% top_n(1,n) %>% group_by(query_cell)  %>% top_n(1,prediction_score_max) %>% summarise(pred_EML=paste(unique(pred_EML),collapse=":"),sub_pred_EML=paste(sub_pred_EML,collapse=":"),prediction_score_max=unique(prediction_score_max))  %>%  mutate(stage=ifelse(sub_pred_EML %in% stage.cate$Early,"Early","uncertained")) %>% mutate(stage=ifelse(sub_pred_EML %in% stage.cate$Late,"Late",stage)) 
  temp2 <- temp_pred_ac_out %>% filter(!query_cell %in% temp1$query_cell) %>% filter(pred_EML %in% "Ambiguous") %>% select(query_cell,pred_EML,sub_pred_EML,prediction_score_max) %>% unique() %>% mutate(stage="uncertained")
  temp3 <- temp_pred_ac_out %>% filter(!query_cell %in% temp1$query_cell)%>% filter(!query_cell %in% temp2$query_cell) %>% filter(pred_EML %in% "low_cor") %>% select(query_cell,pred_EML,sub_pred_EML) %>% unique()%>% mutate(stage="uncertained",prediction_score_max=-1)
  
  predict_out$full.anno <- predict_out$raw.meta %>% mutate(query_cell=cell,pred_EML="nb_failed",sub_pred_EML="nb_failed",stage="unpredictable",prediction_score_max=-1) %>% select(query_cell,pred_EML,sub_pred_EML,stage,prediction_score_max)
  if (nrow(temp1) >0) {
    predict_out$full.anno <- predict_out$full.anno %>% rows_update(temp1,by="query_cell")
  }
  if (nrow(temp2) >0) {
    predict_out$full.anno <- predict_out$full.anno %>% rows_update(temp2,by="query_cell")
  }
  if (nrow(temp3) >0) {
    predict_out$full.anno <- predict_out$full.anno %>% rows_update(temp3,by="query_cell")
  }
  predict_out$full.anno <- predict_out$full.anno  %>% filter(!query_cell %in% predict_out$anno$query_cell ) %>% bind_rows(predict_out$anno )
  
  #' got the predicted pseudotime
  temp.anno <- predict_out$anno %>% select(query_cell,sub_pred_EML) 
  temp.pca.out <- predict_out$query.pca %>% as.data.frame() %>% setNames(paste0("PC_",1:50))%>% tibble::rownames_to_column("cell")
  temp.traj <- list()
  for (s in unique(temp.anno$sub_pred_EML)) {
    if (s %in% ref.pseudotime.sel$sub_rename_EML) {
      temp.query.cells <- temp.anno %>% filter(sub_pred_EML==s) %>% pull(query_cell)
      k=min(length(temp.query.cells),5)
      temp.ref.psdt <- ref.pseudotime.sel %>% filter(sub_rename_EML %in% c(s))
      
      temp.ref.umap.pos <-  temp.pca.out %>% filter(cell %in% temp.ref.psdt$cell)  %>% tibble::column_to_rownames("cell")
      temp.query.umap.pos <- temp.pca.out %>% filter(cell %in% temp.query.cells) %>% tibble::column_to_rownames("cell")
      temp.traj[[s]] <- FunRestricted_knn_cdx( temp.ref.umap.pos,temp.query.umap.pos,k=k)  %>% setNames(c("ref_cell","query_cell")) %>% inner_join(temp.ref.psdt %>% rename(ref_psdt=psdt,ref_cell=cell) %>% select(ref_cell,ref_psdt),by="ref_cell") %>% group_by(query_cell) %>% summarise(pred_psdt=mean(ref_psdt)) %>% inner_join(temp.anno,by="query_cell")
    }
  }
  temp.traj.out <- temp.traj %>% do.call("bind_rows",.)
  
  if (nrow(temp.traj.out)==0) {
    predict_out$anno <- predict_out$anno %>% mutate(pred_psdt=NA)
    predict_out$full.anno <- predict_out$full.anno %>% mutate(pred_psdt=NA)
    
  }else{
    temp.traj.psdt.ac.out <- predict_out$full.anno %>% select(query_cell,sub_pred_EML) %>% inner_join(predict_out$HS %>% rename(query_cell=Scell),by="query_cell") %>% inner_join(temp.traj.out  %>% rename(Hcell=query_cell),by=c("Hcell","sub_pred_EML")) %>% group_by(query_cell,sub_pred_EML) %>% summarise(pred_psdt=mean(pred_psdt)) %>% ungroup()
    
    predict_out$anno <- predict_out$anno %>% left_join(temp.traj.out,by=c("query_cell","sub_pred_EML"))
    predict_out$full.anno <- predict_out$full.anno %>% left_join(temp.traj.psdt.ac.out,by=c("query_cell","sub_pred_EML"))
  }
  
  return(predict_out)
}
