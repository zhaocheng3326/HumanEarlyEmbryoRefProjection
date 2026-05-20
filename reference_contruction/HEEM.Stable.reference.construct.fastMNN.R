#' ---
#' title: "human emrbyo reference construction"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---

# ##loading R library
#R4.0
rm(list=ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scran)
  library(batchelor)
  library(Seurat)
  library(SeuratWrappers)
  library(scuttle)
  library(SeuratDisk)
  library(ff)
  library(uwot)
})

source("src/project.setting.R") #### please modify this path as needed.
source("src/local.quick.fun.R")


suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)



#' loading data
load(paste0("tmp_data/gene.meta.Rdata"),verbose=T)

#' loading specific markers to identity the obsecu markers
#hm.AMvsTE.out <- readRDS(paste0("tmp_data/",TD,"/human.AmnionVsTE.sigMarker.rds")) %>% group_by(Type) %>% top_n(20,power)%>% select(gene,Type) %>% split(.,.$Type) %>% lapply(function(x){x$gene})


# introduce the marmoset Amnion and TE. 

model_file <- paste0(DIR,"/tmp_data/",TD,"/stable.ref.fastMNN.umap.model")
if (file.exists(paste0("tmp_data/",TD,"/stable.ref.fastMNN.Rdata"))) {
  load(paste0("tmp_data/",TD,"/stable.ref.fastMNN.Rdata"),verbose=T)
  ref.umap <- load_uwot(file =model_file)
}else{
  #' update reference annotation 
  meta.filter.updated <- readRDS(paste0("tmp_data/",TD,"/meta.filter.updated.rds"))
  meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("D3post","CS7","SPH2016","nBGuo","Meistermann_2021","Yan2013") ) %>% filter(cellType=="EM") %>%  mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))
  meta.filter <- meta.filter %>% rows_update(meta.filter.updated %>% filter(cell %in% meta.filter$cell) %>% select(cell,EML),by="cell")
  
  counts.filter <- (readRDS(paste0("tmp_data/",TD,"/counts.filter.rds")))[,meta.filter$cell]

  #' modify the cell names
  meta.filter <-meta.filter%>% mutate(cell=paste0("ref_",cell))
  colnames(counts.filter) <- paste0("ref_",colnames(counts.filter))
  
  expG.set <- list()
  for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) { 
    temp.cell <- meta.filter %>% filter(pj==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  
  sce.ob <- list()
  for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) { 
    print(b)
    temp.M <- meta.filter %>% filter(pj==b) 
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% scran::computeSumFactors()
    sce.ob[[b]] <- temp.sce
  }
  
  mBN.sce.ob <- multiBatchNorm(sce.ob$Yan2013,sce.ob$SPH2016,sce.ob$Meistermann_2021,sce.ob$nBGuo,sce.ob$D3post,sce.ob$CS7)
  names(mBN.sce.ob) <- c("Yan2013","SPH2016","Meistermann_2021","nBGuo","D3post","CS7")
  #pdf("tmp_data/temp.sf.pdf")
  mBN.sce.ob %>% lapply(function(x) {data.frame(cell=colnames(x),sf=sizeFactors(x)) %>% tbl_df() %>% return()})  %>% do.call("bind_rows",.) %>% inner_join(meta.filter %>% select(cell,pj)) %>% mutate(od=factor(pj,"Yan2013,Meistermann_2021,SPH2016,nBGuo,D3post,CS7" %>% strsplit(",") %>% unlist(),ordered = T)) %>% ggplot()+geom_violin(mapping=aes(x=od,y=sf),scale = "width",fill="royalblue")+theme_classic()
  #mBN.sce.ob <- multiBatchNorm(sce.ob$SPH2016,sce.ob$D3post,sce.ob$CS7)
  lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  
  temp.M <- meta.filter 
  temp.sel.expG <-rownames(lognormExp.mBN )
  
  sel.od <- "Yan2013,Meistermann_2021,SPH2016,nBGuo,D3post,CS7" %>% strsplit(",") %>% unlist()
  
  data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
  total.overall.VGs <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=5000)}) %>% lapply(VariableFeatures) %>% unlist()%>% unique()
  data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
  data.spt$Yan2013 <- FindVariableFeatures(data.spt$Yan2013,verbose=F,nfeatures=500) ### too few cells for this datasets
  #change the topVG to avoid the VG biased to pre-implantation embryos
  pre.VGs <- SelectIntegrationFeatures(object.list = data.spt[c("Yan2013","Meistermann_2021","nBGuo","SPH2016")], nfeatures = 2000)
  for (pj in c("Yan2013","Meistermann_2021","nBGuo","SPH2016")) {
    VariableFeatures(data.spt[[pj]]) <- pre.VGs
  }

  mnn.VGs.list <- list()
  mnn.VGs.list$VG3000 <- SelectIntegrationFeatures(object.list = data.spt[c("CS7","D3post","SPH2016")], nfeatures = 3000)
  mnn.VGs.list$VG3500 <- SelectIntegrationFeatures(object.list = data.spt[c("CS7","D3post","SPH2016")], nfeatures = 3500)
  mnn.VGs.list$VG4000 <- SelectIntegrationFeatures(object.list = data.spt[c("CS7","D3post","SPH2016")], nfeatures = 4000)
  mnn.VGs.list$VG4000 <- SelectIntegrationFeatures(object.list = data.spt[c("CS7","D3post","SPH2016")], nfeatures = 4000)
  mnn.VGs.list$VG.total <-  data.spt %>% lapply(VariableFeatures) %>% unlist() %>% unique() 
  mnn.VGs.list$total.overall.VGs <- total.overall.VGs

  mnn.VGs.list$VG.sel <- mnn.VGs.list$VG4000
  mnn.VGs <- mnn.VGs.list$VG.sel
  length(mnn.VGs)
  
  for (pj in names(data.spt)) {
    VariableFeatures(data.spt[[pj]]) <- mnn.VGs
  }
  
  #' release memory
  rm(data.merge)
  rm(counts.filter)
  rm(data.spt)
  
  mBN.sce.ob <- mBN.sce.ob[sel.od]
 
  
  mat.list <- list()
  for (idx in seq_along(mBN.sce.ob)) {
    temp <- (assay(mBN.sce.ob[[idx]], i="logcounts"))
    temp.l2 <- pmax(1e-8,  temp %>% cosineNorm(mode="l2norm",subset.row = mnn.VGs))
    mat.list[[idx]] <- scuttle::normalizeCounts(temp, size_factors=temp.l2, center_size_factors=FALSE, log=FALSE)[mnn.VGs,]
  }
  # modified from multibatchnorm, and fastMNN function
  grand.centers <- 0
  weights <- rep(1,length(mat.list)) #weight set to 1
  for (idx in seq_along(mat.list)) {
    centers <- rowMeans( mat.list[[idx]])
    grand.centers <- grand.centers + centers * weights[idx]
  }
  grand.centers <- grand.centers/sum(weights)
  centered <- lapply(mat.list, function(x) NULL)

  for (idx in seq_along(mat.list)) {
    current <- mat.list[[idx]] - grand.centers
    centered[[idx]] <- current
  }


  scaled <- vector("list", length(mat.list))
  for (idx in seq_along(centered)) {
    current <- centered[[idx]]
    w <- sqrt(ncol(current) / weights[idx])
    current <- current/w 
    scaled[[idx]] <- current
  }
  scaled <- do.call(cbind, scaled)
  centers <- grand.centers
  set.seed(10);svd.out <- BiocSingular::runSVD(scaled,k=0, nu=50,nv=0)
  #set.seed(10);svd.out <- BiocSingular::runSVD(scaled,k=0, nu=50,nv=0,BSPARAM = BSPARAM)
  output <- centered
  for (idx in seq_along(centered)) {
    output[[idx]] <- as.matrix(crossprod(centered[[idx]], svd.out$u))
  }
  names(output) <- names(mBN.sce.ob)
  
  left.data <- output[[sel.od[1]]]
  to.add <- NULL
  for (n in 2:length(output)) {
    right.data <- output[[sel.od[n]]]
    
    right.data <- batchelor:::.orthogonalize_other(right.data, NULL, to.add ) #.center_along_batch_vector iteration #Remove variation along the batch vector
    left.data <- batchelor:::.orthogonalize_other(left.data, NULL, to.add )
    
    mnn.sets <- FunRestricted_mnn(left.data, right.data)
    ave.out <- batchelor:::.average_correction(left.data, mnn.sets$first, right.data, mnn.sets$second)
    overall.batch <- colMeans(ave.out$averaged)
    to.add <- c(to.add,list(overall.batch))
    
    left.data <- batchelor:::.center_along_batch_vector(left.data, overall.batch)
    right.data <- batchelor:::.center_along_batch_vector(right.data, overall.batch)
    
    re.ave.out <- batchelor:::.average_correction(left.data, mnn.sets$first, right.data, mnn.sets$second)
    right.data <- batchelor:::.tricube_weighted_correction(right.data, re.ave.out$averaged, re.ave.out$second, k=20, ndist=3)
    left.data<- left.data %>% rbind(right.data)
  }
 
  ref.seed <- 10
  ref.nPC <- 50
  ref.mnn.out <- left.data
  set.seed(ref.seed)
  ref.umap <- umap(as.data.frame(ref.mnn.out[,c(1:ref.nPC)]), ret_model = TRUE, n_neighbors = 30L, metric = "cosine",min_dist = 0.3)
  temp.umap <- ref.umap$embedding %>% as.data.frame() %>% setNames(c("UMAP_1","UMAP_2")) %>% tibble::rownames_to_column("cell") %>% inner_join(meta.filter,by="cell")
  
  ref.grand.centers <- grand.centers
  ref.svd.out <- svd.out
  ref.mnn.out <- left.data
  ref.to.add <- to.add
  ref.ave.out <- ave.out
  ref.overall.batch <- overall.batch
  ref.sce.ob <- sce.ob
  ref.mnnVGs <- mnn.VGs
  ref.sel.expG <- sel.expG
  ref.sel.od <- sel.od
  ref.meta <- meta.filter
  ref.mnnVGs.list <- mnn.VGs.list
  
  #' for stable normalization
  ref.sce.stats <- list()
  for (n in names(ref.sce.ob)) {
    ref.sce.stats[[n]] <- list()
    sf <-  sizeFactors(ref.sce.ob[[n]])
    sf <- sf/mean(sf)
    ave <- scuttle::calculateAverage(ref.sce.ob[[n]]@assays@data$counts,  size.factors=sf )
    ref.sce.stats[[n]]$sf <- sf
    ref.sce.stats[[n]]$ave <- ave
  }
  ave.list <- lapply(ref.sce.stats,function(x){x$ave})
  sf.list <- lapply(ref.sce.stats,function(x){x$sf})
  min.mean <- 1
  nbatches <- length(ave.list)
  collected.ratios <- matrix(1, nbatches, nbatches)
  for (first in seq_len(nbatches-1L)) {
    first.ave <- ave.list[[first]]
    first.sum <- sum(first.ave)
    
    for (second in first + seq_len(nbatches - first)) {
      second.ave <- ave.list[[second]]
      second.sum <- sum(second.ave)
      
      # Mimic calcAverage(cbind(first.ave, second.ave)).
      grand.mean <- (first.ave/first.sum + second.ave/second.sum)/2 * (first.sum + second.sum)/2
      keep <- grand.mean >= min.mean
      
      # Computing it twice for exactly equal results when order of batches is rearranged.
      kept.f <- first.ave[keep]
      kept.s <- second.ave[keep]
      curratio1 <- median(kept.s/kept.f)
      curratio2 <- median(kept.f/kept.s)
      
      if (!is.finite(curratio1) || curratio1==0 || !is.finite(curratio2) || curratio2==0) {
        stop("median ratio of averages between batches is not finite")
      }
      
      collected.ratios[first,second] <- curratio1
      collected.ratios[second,first] <- curratio2
    }
  }
  # Rescaling to the lowest coverage batch.
  #ref.smallest <- which.min(apply(collected.ratios[-query.index,-query.index], 2, min, na.rm=TRUE))
  ref.smallest <- which.min(apply(collected.ratios, 2, min, na.rm=TRUE))
  ref.rescaling <- collected.ratios[,ref.smallest]
  
  
  #' calculate the raw data
  ref.raw.data <- list()
  for (x in seq_along(1:length(ref.rescaling))) {
    temp <-  ref.sce.ob[[x]] 
    sizeFactors(temp) <-  ref.sce.stats[[x]]$sf/ref.rescaling[x]
    
    ref.raw.data[[x]] <- ((scuttle::logNormCounts(temp,center.size.factors=FALSE) %>% logcounts())[ref.mnnVGs,]  %>% cosineNorm( mode="all" ))$matrix %>% t() 
  }
  ref.raw.data <- ref.raw.data %>% do.call("rbind",.)
  
  #ref.exclude.cells <- temp.umap %>% filter(pj=="nBGuo" & EML=="TE") %>% arrange(UMAP_2) %>% head(2) %>% tail(1) %>% pull(cell) 
  ref.exclude.cells <- temp.umap %>% filter(UMAP_1 > 0, UMAP_2< -2.5 ,EML=="TE", devTime %in% c("E6","E7"))  %>% pull(cell) 
  APHL_UMAP(temp.umap,ref.exclude.cells)
  #' saving object
  save(ref.sce.ob,ref.mnnVGs,ref.sel.expG,ref.sel.od,ref.meta,ref.mnn.out,ref.seed,ref.nPC,ref.sce.stats,ref.mnnVGs.list,ref.grand.centers,ref.svd.out,ref.to.add ,ref.ave.out,ref.overall.batch,ref.smallest,ref.rescaling,ref.raw.data,ref.exclude.cells,file=paste0("tmp_data/",TD,"/stable.ref.fastMNN.Rdata"))
  save_uwot(ref.umap, file = model_file)
  saveRDS(ref.sel.expG,paste0("tmp_data/",TD,"/stable.ref.fastMNN.ref.sel.expG.rds"))
  saveRDS(ref.mnnVGs,paste0("tmp_data/",TD,"/stable.ref.fastMNN.mnnVG.rds"))
  
  temp.umap <- ref.umap$embedding %>% as.data.frame() %>% setNames(c("UMAP_1","UMAP_2")) %>% tibble::rownames_to_column("cell") %>% inner_join(ref.meta,by="cell")
  saveRDS(temp.umap %>% tbl_df(),paste0("tmp_data/",TD,"/stable.ref.umap.rds"))
  
  #' try different parameters
  data.ob.umap.list <- list()
  for (nPC in c(15,20,25,30,35,40,50)) {
    set.seed(ref.seed)
    data.ob.umap.list[[paste0("PC",nPC)]] <-  (umap(as.data.frame(ref.mnn.out[,c(1:nPC)]), ret_model = TRUE, n_neighbors = 30L, metric = "cosine",min_dist = 0.3))$embedding %>% as.data.frame() %>% setNames(c("UMAP_1","UMAP_2")) %>% tibble::rownames_to_column("cell") %>% tbl_df()%>% mutate(nPC_sel=paste0("top:",nPC))
  }
  data.ob.umap.list <- data.ob.umap.list %>% do.call("bind_rows",.)
  saveRDS(data.ob.umap.list,paste0("tmp_data/",TD,"/stable.ref.umap.PC.para.list.rds"))
}


#' create the cluster information 
if (file.exists(paste0("tmp_data/",TD,"/stable.fastMNN.ref.cluster.rds"))) {
  ref.cluster <- readRDS(paste0("tmp_data/",TD,"/stable.fastMNN.ref.cluster.rds"))
}else{
  #ref.cluster<- paste0("C",Seurat:::RunModularityClusteringCpp(FindNeighbors(ref.mnn.out[,],verbose=F,k.param=15)$snn,modularity=1,resolution = 2,algorithm=1,nRandomStarts=10,nIterations=10,randomSeed=0,printOutput=T,edgefilename='')) 
  ref.cluster <- data.frame(cell=rownames(ref.mnn.out),SC = paste0("C",Seurat:::RunLeiden(FindNeighbors(ref.mnn.out[,],verbose=F,k.param=20)$snn,  method = "matrix",resolution.parameter = 2, random.seed=0,partition.type = "RBConfigurationVertexPartition"))) %>% tbl_df()
  saveRDS(ref.cluster,paste0("tmp_data/",TD,"/stable.fastMNN.ref.cluster.rds"))
}

