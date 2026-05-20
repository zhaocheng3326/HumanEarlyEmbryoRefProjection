#' ---
#' title: "traj for human reference"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library

rm(list=ls())
rewrite=FALSE

#' check whether in local computer
if (grepl("KI-",Sys.info()['nodename'])) {
  print("local computer")
  source("/Users/cheng.zhao/chzhao_bioinfo/PC/SnkM/SgCell.R")
  base_dir <- "/Users/cheng.zhao/Documents"
} else {
  print("On server")
  condaENV <- "/home/chenzh/miniconda3/envs/R4.2"  ####
  LBpath <- paste0(condaENV ,"/lib/R/library")
  .libPaths(LBpath)
  base_dir="/home/chenzh"
}



suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(tibble)
  library(scran)
  library(batchelor)
  #library(SeuratWrappers)
  library(slingshot)
  library(SCP)
})


# working directory
DIR <- "~/My_project/HumanEarlyEmbryoRefProjection" ### please modify this pathway as needed
knitr::opts_knit$set(root.dir=DIR)
setwd(DIR)

#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")
#source("src/figures.setting.R")

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter
#st_sfc <- sf::st_sfc


options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="Mar_2023"

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)


# slingshot object
stable.ref.umap <- readRDS(paste0("tmp_data/",TD,"/GS.stable.ref.umap.rds")) 
stable.ref.umap <- stable.ref.umap %>% rename(UMAP_2=UMAP_1, UMAP_1=UMAP_2) # %>% mutate(UMAP_2=-1*UMAP_2)


#' loading models
load("tmp_data/Mar_2023/stable.ref.fastMNN.Rdata",verbose=T)
rm(ref.sce.ob)
TF.genes <- read.delim("doc/TF_names_v_1.01.txt",stringsAsFactors = F,head=F)$V1


meta.filter <- stable.ref.umap %>% select(cell:rename_EML)
counts.filter <- (readRDS(paste0("tmp_data/",TD,"/counts.filter.rds")))[,meta.filter$cell]

sel.exp <- readRDS(paste0("tmp_data/",TD,"/lognormExp.only.human.EM.mBN.rds"))
TF.genes <- TF.genes %>% intersect(rownames(sel.exp))


#' cells with uncertained annotation
ref.further.exclude.cells <- readRDS(paste0("tmp_data/",TD,"/ref.further.exclude.cells.rds")) 

#' construct the overall object
temp.M <- stable.ref.umap
temp.M <- temp.M %>% filter(!cell %in% ref.further.exclude.cells )
table(duplicated(temp.M$cell))

temp.sel.expG<- rownames(sel.exp)
data.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
data.ob@assays$RNA@data <- as.matrix(sel.exp[temp.sel.expG,colnames(data.ob)])

mnn.embedding <-  ref.mnn.out[paste0("ref_",temp.M$cell),]
rownames(mnn.embedding) <- gsub("ref_","",rownames(mnn.embedding))
colnames(mnn.embedding) <-  paste0("mnn_",1:50)
umap.embedding <- temp.M  %>% select(cell,UMAP_1,UMAP_2) %>% tibble::column_to_rownames("cell") %>% as.matrix()
data.ob[["mnn"]] <- CreateDimReducObject(embeddings =mnn.embedding, key = "mnn_", assay = "RNA",global = TRUE)
data.ob[["UMAP"]] <- CreateDimReducObject(embeddings =umap.embedding , key = "UMAP_", assay = "RNA",global = TRUE)


#' general lineage branch
if (file.exists(paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.rds"))) {
  data.temp <- readRDS(paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.rds"))
  load(paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.Rdata"),verbose=T)
}else{
  traj.related.TF <- list()
  temp.sub.M <- temp.M %>% filter(sub_rename_EML %in% c("Zygote","2-4 cell","8 cell","Morula","Prelineage","ICM","Early_EPI","Late_EPI","Early_Hypoblast","Late_Hypoblast","TE","CTB")) # "PriS","Mesoderm","Amnion","AdvMes","STB","EVT"
  data.sub.temp <- subset(data.ob,cell=temp.sub.M$cell)
  data.sce <- as.SingleCellExperiment(data.sub.temp)
  
  data.sce <- slingshot(data.sce,  reducedDim = 'UMAP',clusterLabels = 'sub_rename_EML', start.clus = 'Zygote')
  data.crv <- SlingshotDataSet(data.sce)
  plot(reducedDims(as.SingleCellExperiment(data.ob))$UMAP,col="grey",pch=15,cex=0.5,xaxt="n",yaxt="n",bty="n",xlab = "", ylab = "", axes = FALSE)
  lines(data.crv, lwd=2, type = 'lineages', col = c("royalblue"))
  # lines(data.crv, lwd=2, col = c("royalblue"))
  
  
  data.pseudotime.na <-  slingPseudotime(data.crv, na=T) %>% as.data.frame() %>% setNames(c("EPI_traj","Hypoblast_traj","TE_traj"))%>% tibble::rownames_to_column("cell") %>% tbl_df() %>% inner_join(temp.M %>% select(cell,devTime,pj,rename_EML,sub_rename_EML),by="cell")
  data.pseudotime <-  slingPseudotime(data.crv, na=FALSE) %>% as.data.frame()  %>% setNames(c("EPI_traj","Hypoblast_traj","TE_traj"))%>% tibble::rownames_to_column("cell") %>% tbl_df() %>%  gather(traj,psdt,-cell) %>% inner_join(temp.M %>% select(cell,devTime,pj,rename_EML,sub_rename_EML),by="cell") # %>% gather(traj,psdt,
  #' scale TE, PE psedutime to make it comparable to EPI
  psdt.E4.EPItraj.max <-  data.pseudotime %>% filter(devTime=="E4" & traj=="EPI_traj") %>% pull(psdt) %>% max()
  psdt.E14.EPItraj.max <-  data.pseudotime %>% filter(devTime=="E14" & traj=="EPI_traj") %>% pull(psdt) %>% max()
  
  psdt.E4.PEtraj.max <-  data.pseudotime %>% filter(devTime=="E4" & traj=="Hypoblast_traj") %>% pull(psdt) %>% max()
  psdt.E14.PEtraj.max <-  data.pseudotime %>% filter(devTime=="E14" & traj=="Hypoblast_traj") %>% pull(psdt) %>% max()
  
  psdt.E4.TEtraj.max <-  data.pseudotime %>% filter(devTime=="E4" & traj=="TE_traj") %>% pull(psdt) %>% max()
  psdt.E14.TEtraj.max <-  data.pseudotime %>% filter(devTime=="E14" & traj=="TE_traj") %>% pull(psdt) %>% max()
  
  psdt.EPITE.scale.factor <- (psdt.E14.TEtraj.max-psdt.E4.TEtraj.max)/ (psdt.E14.EPItraj.max-psdt.E4.EPItraj.max) ### align psdt between EPI and TE
  psdt.EPIPE.scale.factor <- (psdt.E14.PEtraj.max-psdt.E4.PEtraj.max)/ (psdt.E14.EPItraj.max-psdt.E4.EPItraj.max) ### align psdt between EPI and TE
  #' a scaled psdt for cells belong to diff traj
  data.pseudotime.mod <- data.pseudotime %>% mutate(psdt=ifelse(traj=="TE_traj" & (!devTime %in% c("E1","E2","E3","E4")) & psdt > psdt.E4.TEtraj.max, (psdt-psdt.E4.TEtraj.max)/psdt.EPITE.scale.factor+psdt.E4.TEtraj.max,psdt )) %>% mutate(psdt=ifelse(traj=="Hypoblast_traj" & (!devTime %in% c("E1","E2","E3","E4")) & psdt > psdt.E4.PEtraj.max, (psdt-psdt.E4.PEtraj.max)/psdt.EPIPE.scale.factor+psdt.E4.PEtraj.max,psdt ))
  
  #' check
  data.pseudotime %>% split(.,.$traj) %>% lapply(function(x) {x%>% arrange(psdt) %>% tibble::rowid_to_column("n")  %>% select(cell,traj,n)}) %>% do.call("bind_rows",.) %>% inner_join(data.pseudotime.mod %>% split(.,.$traj) %>% lapply(function(x) {x%>% arrange(psdt) %>% tibble::rowid_to_column("n")  %>% select(cell,traj,n)}) %>% do.call("bind_rows",.),by=c("cell","traj")) %>% mutate(n=n.x-n.y) %>% pull(n) %>% max() ## same orders

  data.cellWeights <- slingCurveWeights(data.sce)
  
  traj_lineage_sel <- data.frame(sub_rename_EML=c("Morula","Prelineage","ICM","Early_EPI","Late_EPI"),traj="EPI_traj",sel_traj="EPI_main_traj")  %>% bind_rows(data.frame(sub_rename_EML=c("Morula","Prelineage","ICM","Early_Hypoblast","Late_Hypoblast"),traj="Hypoblast_traj",sel_traj="PE_main_traj"))%>% bind_rows( data.frame(sub_rename_EML=c("Morula","Prelineage","TE","CTB"),traj="TE_traj",sel_traj="TE_main_traj")) %>% tbl_df() %>% mutate_all(as.vector)
  
  data.pseudotime.sel <- data.pseudotime.mod %>% inner_join(traj_lineage_sel,by=c("traj","sub_rename_EML"),multiple = "all")
  #' modified from SCP tools https://github.com/zhanghao-njmu/SCP
  data.temp <- subset(data.ob,cell=data.pseudotime.sel$cell)
  data.temp@meta.data <- data.temp@meta.data %>% cbind((data.pseudotime.sel %>% select(cell,psdt,sel_traj) %>% spread(sel_traj,psdt) %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),])
  
  for (stj in unique(data.pseudotime.sel$sel_traj)) {
    temp.pseudotime <- data.pseudotime.sel %>% filter(sel_traj==stj)  %>% select(cell,psdt,sel_traj,rename_EML,devTime,pj) %>% arrange(psdt) 
    temp.pseudotime.vec <- temp.pseudotime %>% select(cell,psdt) %>% tibble::column_to_rownames("cell") %>% as.matrix()
    temp.sel.exp <- sel.exp[TF.genes,temp.pseudotime$cell]
    temp.sel.exp <- temp.sel.exp[rowSums(temp.sel.exp >0 ) >= 5,]
    gam_out <- list()
    Y_ordered <- temp.sel.exp %>% as.matrix()
    t_ordered <- temp.pseudotime.vec[,"psdt"]
    for (n in seq_len(nrow(Y_ordered))) {
      print(paste0(stj,n))
      feature_nm <- rownames(Y_ordered)[n]
      family_use <- "gaussian"
      sizefactror <- 1
      
      mod <- mgcv::gam(y ~ s(x, bs = "cs") + offset(rep(log(1),ncol(Y_ordered))),family = rep(family_use,ncol(Y_ordered)),data = data.frame(y = Y_ordered[feature_nm,,drop=T ], x = t_ordered))
      
      pre <- predict(mod, type = "link", se.fit = TRUE)
      upr <- pre$fit + (2 * pre$se.fit)
      lwr <- pre$fit - (2 * pre$se.fit)
      upr <- mod$family$linkinv(upr)
      lwr <- mod$family$linkinv(lwr)
      res <- summary(mod)
      fitted <- fitted(mod)
      pvalue <- res$s.table[[4]]
      dev.expl <- res$dev.expl
      r.sq <- res$r.sq
      fitted.values <- fitted * sizefactror
      upr.values <- upr * sizefactror
      lwr.values <- lwr * sizefactror
      exp_ncells <- sum(Y_ordered[feature_nm, ] > min(Y_ordered[feature_nm, ]), na.rm = TRUE)
      peaktime <- median(t_ordered[fitted.values > quantile(fitted.values, 0.99, na.rm = TRUE)])
      valleytime <- median(t_ordered[fitted.values < quantile(fitted.values, 0.01, na.rm = TRUE)])
      
      # ggplot(data = data.frame(
      #   x = t_ordered,
      #   raw = FetchData(data.ob, vars = feature_nm, slot = "data")[names(t_ordered), feature_nm, drop = TRUE],
      #   fitted = fitted.values,
      #   upr.values = upr.values,
      #   lwr.values = lwr.values
      # )) +
      #   geom_point(aes(x = x, y = raw), color = "black", size = 0.5) +
      #   geom_point(aes(x = x, y = fitted), color = "red", size = 0.5) +
      #   geom_path(aes(x = x, y = upr.values), color = "blue") +
      #   geom_path(aes(x = x, y = lwr.values), color = "green")
      
      # a <- data.frame(x = t_ordered, y = Y_ordered[feature_nm, ])
      # qplot(a$x, a$y)
      # length(unique(a$y) > 5)
      # a <- a[a$y > 0, ]
      # qplot(a$x, a$y)
      gam_out[[feature_nm]] <- list(
        features = feature_nm, exp_ncells = exp_ncells,
        r.sq = r.sq, dev.expl = dev.expl,
        peaktime = peaktime, valleytime = valleytime,
        pvalue = pvalue, fitted.values = fitted.values,
        upr.values = upr.values, lwr.values = lwr.values
      )
    }
    
    raw_matrix <- Y_ordered 
    fitted_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["fitted.values"]]))
    colnames(fitted_matrix) <- rownames(Y_ordered)
    fitted_matrix <- cbind(pseudotime = t_ordered, fitted_matrix)
    
    upr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["upr.values"]]))
    colnames(upr_matrix) <- rownames(Y_ordered)
    upr_matrix <- cbind(pseudotime = t_ordered, upr_matrix)
    
    lwr_matrix <- do.call(cbind, lapply(gam_out, function(x) x[["lwr.values"]]))
    colnames(lwr_matrix) <- rownames(Y_ordered)
    lwr_matrix <- cbind(pseudotime = t_ordered, lwr_matrix)
    
    DynamicFeatures <- as.data.frame(do.call(rbind.data.frame, lapply(gam_out, function(x) x[!names(x) %in% c("fitted.values", "upr.values", "lwr.values")]))) %>% filter(!is.na(peaktime))
    char_var <- c("features")
    numb_var <- colnames(DynamicFeatures)[!colnames(DynamicFeatures) %in% char_var]
    DynamicFeatures[, char_var] <- lapply(DynamicFeatures[, char_var, drop = FALSE], as.character)
    DynamicFeatures[, numb_var] <- lapply(DynamicFeatures[, numb_var, drop = FALSE], as.numeric)
    rownames(DynamicFeatures) <- DynamicFeatures[["features"]]
    DynamicFeatures[, "padjust"] <- p.adjust(DynamicFeatures[, "pvalue", drop = TRUE])
    raw_matrix <- raw_matrix[rownames(DynamicFeatures),]
    fitted_matrix <- fitted_matrix[,c("pseudotime",rownames(DynamicFeatures))]
    upr_matrix <- upr_matrix[,colnames(fitted_matrix)]
    lwr_matrix <- upr_matrix[,colnames(fitted_matrix)]
    
    res <- list(
      DynamicFeatures = DynamicFeatures,
      raw_matrix = raw_matrix,
      fitted_matrix = fitted_matrix,
      upr_matrix = upr_matrix,
      lwr_matrix = lwr_matrix,
      libsize = NULL,
      lineages = stj,
      family = family_use
    )
    data.temp@tools[[paste0("DynamicFeatures_", stj)]] <- res
  }
  
  psdt.genes <- list()
  
  # DynamicHeatmap for clustering
  ht1 <- DynamicHeatmap( srt = data.temp,lineages = c("EPI_main_traj"),n_split = 8,split_method = "kmeans-peaktime",cell_annotation = c("sub_rename_EML","devTime"),row_names_side="right",slot="data",min_expcells=10)#,features_label=c("DNMT1")
  temp.M <- data.pseudotime.sel %>% filter(traj=="EPI_traj") %>% arrange(psdt)
  temp.genes <- ht1$feature_metadata$features
  temp.cell.anno <- temp.M %>% select(cell,psdt,devTime,sub_rename_EML) %>% tibble::column_to_rownames("cell")
  temp.gene.anno <- ht1$feature_metadata %>% select(feature_split) %>% mutate(feature_split=recode(feature_split,"C2"="C1","C3"="C2","C4"="C2","C5"="C3","C6"="C4","C7"="C5","C7"="C6"))
  sel.exp[temp.genes,temp.M$cell] %>% FunPreheatmapNoLog() %>% pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_row=temp.gene.anno,annotation_col=temp.cell.anno,col=heat.col,gaps_row=(temp.gene.anno %>% split(.,.$feature_split) %>% lapply(nrow) %>% cumsum()))
  psdt.genes$EPI_traj <- data.temp@tools$DynamicFeatures_EPI_main_traj$DynamicFeatures %>% tbl_df() %>% rename(gene=features) %>% left_join(temp.gene.anno  %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% rowid_to_column("index") %>% rename(cluster=feature_split),by="gene") %>% mutate(traj="EPI_traj")
  
  ht1 <- DynamicHeatmap( srt = data.temp,lineages = c("PE_main_traj"),n_split = 8,split_method = "kmeans-peaktime",cell_annotation = c("sub_rename_EML","devTime"),row_names_side="right",slot="data",min_expcells=10)#,features_label=c("DNMT1")
  temp.M <- data.pseudotime.sel %>% filter(traj=="Hypoblast_traj") %>% arrange(psdt)
  temp.genes <- ht1$feature_metadata$features
  temp.cell.anno <- temp.M %>% select(cell,psdt,devTime,rename_EML) %>% tibble::column_to_rownames("cell")
  temp.gene.anno <- ht1$feature_metadata %>% select(feature_split) %>% mutate(feature_split=recode(feature_split,"C2"="C1","C3"="C1","C4"="C2","C5"="C3","C6"="C4","C7"="C5","C8"="C6"))
  sel.exp[temp.genes,temp.M$cell] %>% FunPreheatmapNoLog() %>% pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_row=temp.gene.anno,annotation_col=temp.cell.anno,col=heat.col,gaps_row=(temp.gene.anno %>% split(.,.$feature_split) %>% lapply(nrow) %>% cumsum()))
  psdt.genes$Hypoblast_traj <- data.temp@tools$DynamicFeatures_PE_main_traj$DynamicFeatures %>% tbl_df() %>% rename(gene=features) %>% left_join(temp.gene.anno  %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% rowid_to_column("index")%>% rename(cluster=feature_split),by="gene") %>% mutate(traj="Hypoblast_traj")
  
  ht1 <- DynamicHeatmap( srt = data.temp,lineages = c("TE_main_traj"),n_split = 7,split_method = "kmeans-peaktime",cell_annotation = c("sub_rename_EML","devTime"),row_names_side="right",slot="data",min_expcells=10)#,features_label=c("DNMT1")
  temp.M <- data.pseudotime.sel %>% filter(traj=="TE_traj") %>% arrange(psdt)
  temp.genes <- ht1$feature_metadata$features
  temp.cell.anno <- temp.M %>% select(cell,psdt,devTime,rename_EML) %>% tibble::column_to_rownames("cell")
  temp.gene.anno <- ht1$feature_metadata %>% select(feature_split) %>% mutate(feature_split=recode(feature_split,"C2"="C1","C3"="C2","C4"="C3","C5"="C4","C6"="C5","C7"="C6"))
  sel.exp[temp.genes,temp.M$cell] %>% FunPreheatmapNoLog() %>% pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_row=temp.gene.anno,annotation_col=temp.cell.anno,col=heat.col,gaps_row=(temp.gene.anno %>% split(.,.$feature_split) %>% lapply(nrow) %>% cumsum()))
  psdt.genes$TE_traj <- data.temp@tools$DynamicFeatures_TE_main_traj$DynamicFeatures %>% tbl_df() %>% rename(gene=features) %>% left_join(temp.gene.anno  %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% rowid_to_column("index")%>% rename(cluster=feature_split),by="gene") %>% mutate(traj="TE_traj")
  
  #' check two trajectories
  ht1 <- DynamicHeatmap( srt = data.temp,lineages = c("EPI_main_traj","TE_main_traj"),reverse_ht="TE_main_traj",n_split = 7,split_method = "kmeans-peaktime",cell_annotation = c("sub_rename_EML","devTime"),row_names_side="right",slot="data",min_expcells=10)#,features_label=c("DNMT1")
  temp.M <- data.pseudotime.sel %>% filter(traj %in% c("TE_traj","EPI_traj")) %>% arrange(psdt)
  temp.genes <- ht1$feature_metadata$features
  temp.cell.anno <- temp.M %>% select(cell,traj,psdt,devTime,sub_rename_EML) %>% mutate(SID=paste(traj,cell,sep=":")) # %>% tibble::column_to_rownames("cell")
  temp.gene.anno <- ht1$feature_metadata %>% tibble::rownames_to_column("gene") %>% tbl_df() %>%select(gene,feature_split,index) %>% left_join(psdt.genes$EPI_traj %>% mutate(EPI_traj_cluster=cluster) %>% filter(!is.na(EPI_traj_cluster)) %>% select(gene,EPI_traj_cluster),by="gene")%>% left_join(psdt.genes$TE_traj %>% mutate(TE_traj_cluster=cluster) %>% filter(!is.na(TE_traj_cluster)) %>% select(gene,TE_traj_cluster),by="gene") 
  sel.exp[temp.genes,temp.M %>% filter(traj=="EPI_traj") %>% arrange(desc(psdt)) %>% pull(cell)] %>% FunPreheatmapNoLog() %>% pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_row=(temp.gene.anno %>% select(-index)%>% tibble::column_to_rownames("gene")),annotation_col=(temp.cell.anno %>% filter(traj=="EPI_traj") %>% select(-SID) %>% tibble::column_to_rownames("cell")),col=heat.col,gaps_row=(temp.gene.anno %>% split(.,.$feature_split) %>% lapply(nrow) %>% cumsum()))
  sel.exp[temp.genes,temp.M %>% filter(traj=="TE_traj") %>% arrange(psdt) %>% pull(cell)] %>% FunPreheatmapNoLog() %>% pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_row=(temp.gene.anno %>% select(-index)%>% tibble::column_to_rownames("gene")),annotation_col=(temp.cell.anno %>% filter(traj=="TE_traj") %>% select(-SID) %>% tibble::column_to_rownames("cell")),col=heat.col,gaps_row=(temp.gene.anno %>% split(.,.$feature_split) %>% lapply(nrow) %>% cumsum()))
  psdt.genes$EPI_TE_traj <- temp.gene.anno  %>% mutate(traj="EPI_TE_traj")
  
  #' check two trajectories
  ht1 <- DynamicHeatmap( srt = data.temp,lineages = c("EPI_main_traj","PE_main_traj"),reverse_ht="PE_main_traj",n_split = 7,split_method = "kmeans-peaktime",cell_annotation = c("sub_rename_EML","devTime"),row_names_side="right",slot="data",min_expcells=10)#,features_label=c("DNMT1")
  temp.M <- data.pseudotime.sel %>% filter(traj %in% c("Hypoblast_traj","EPI_traj")) %>% arrange(psdt)
  temp.genes <- ht1$feature_metadata$features
  temp.cell.anno <- temp.M %>% select(cell,traj,psdt,devTime,sub_rename_EML) %>% mutate(SID=paste(traj,cell,sep=":")) # %>% tibble::column_to_rownames("cell")
  temp.gene.anno <- ht1$feature_metadata %>% tibble::rownames_to_column("gene") %>% tbl_df() %>%select(gene,feature_split,index) %>% left_join(psdt.genes$EPI_traj %>% mutate(EPI_traj_cluster=cluster) %>% filter(!is.na(EPI_traj_cluster)) %>% select(gene,EPI_traj_cluster),by="gene")%>% left_join(psdt.genes$Hypoblast_traj %>% mutate(PE_traj_cluster=cluster) %>% filter(!is.na(PE_traj_cluster)) %>% select(gene,PE_traj_cluster),by="gene") 
  sel.exp[temp.genes,temp.M %>% filter(traj=="EPI_traj") %>% arrange(desc(psdt)) %>% pull(cell)] %>% FunPreheatmapNoLog() %>% pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_row=(temp.gene.anno %>% select(-index)%>% tibble::column_to_rownames("gene")),annotation_col=(temp.cell.anno %>% filter(traj=="EPI_traj") %>% select(-SID) %>% tibble::column_to_rownames("cell")),col=heat.col,gaps_row=(temp.gene.anno %>% split(.,.$feature_split) %>% lapply(nrow) %>% cumsum()))
  sel.exp[temp.genes,temp.M %>% filter(traj=="Hypoblast_traj") %>% arrange(psdt) %>% pull(cell)] %>% FunPreheatmapNoLog() %>% pheatmap(cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,annotation_row=(temp.gene.anno %>% select(-index)%>% tibble::column_to_rownames("gene")),annotation_col=(temp.cell.anno %>% filter(traj=="Hypoblast_traj") %>% select(-SID) %>% tibble::column_to_rownames("cell")),col=heat.col,gaps_row=(temp.gene.anno %>% split(.,.$feature_split) %>% lapply(nrow) %>% cumsum()))
  psdt.genes$EPI_PE_traj <- temp.gene.anno  %>% mutate(traj="EPI_PE_traj")
  
  psdt.genes.out <- psdt.genes[c("EPI_traj","Hypoblast_traj","TE_traj")] %>% do.call("bind_rows",.)
  
  
  saveRDS(data.temp,file= paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.rds"))
  save(data.pseudotime.na,data.pseudotime.mod ,data.pseudotime.sel,data.pseudotime,data.sce,file=paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.Rdata"))
  saveRDS(data.pseudotime.sel, paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.psdt.rds"))
  saveRDS(data.pseudotime.mod,file=paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.psdt.mod.all.rds"))
  saveRDS(psdt.genes,paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.psdt.genes.rds"))
  saveRDS(data.crv,paste0("tmp_data/",TD,"/human.EM.main.traj.data.crv.rds"))
  #saveRDS(psdt.genes.out, paste0("tmp_data/",TD,"/human.EM.SCP.main.traj.psdt.genes.out.rds"))
}

