#' ---
#' title: "integrating three species together"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---


# ##loading R library

rm(list=ls())
rewrite=FALSE


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
  library(SeuratWrappers)
})

source("src/project.setting.R") #### please modify this path as needed.


## define the EM brief structure data
#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")


suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)



#' loading data
load(paste0("tmp_data/",TD,"/marmoset_counts_meta.Rdata"),verbose = T)
meta_marmo <- meta_marmo %>% filter(cellType=="EM")
VG_Marmoset <-rownames(readRDS(paste0("tmp_data/",TD,"/data.Marmoset.EM.rds")) @tools$RunFastMN@assays@data$reconstructed)
#' check the st of marmoset genes 
marmo.genes.anno %>% filter(ensembl_gene_id %in% c(count_marmo_genes %>% filter(source=="ENS_ID") %>% pull(GeneID))) %>% filter(hsapiens_homolog_orthology_type=="ortholog_one2one") %>% pull(hsapiens_homolog_associated_gene_name) %>% setdiff(count_marmo_genes$GeneID)
table(count_marmo_genes$source)
# we can tell except "MARCH9","MARC2","MARCH8","MARCH11","Y_RNA","Metazoa_SRP", the majority "ortholog_one2one" has been transfered to human gene names (decided to do the overlap directly)

meta_human <- readRDS(paste0("tmp_data/",TD,"/data.ob.umap.only.EM.rds")) %>% select(-c(UMAP_1,UMAP_2))
count_human <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))[,meta_human$cell]
VG_human <-VariableFeatures(readRDS(paste0("tmp_data/",TD,"/data.ob.only.EM.rds")))

load(paste0("tmp_data/",TD,"/cyno_counts_meta.update.Rdata"),verbose = T)
meta_cyno <- meta_cyno %>% filter(cellType=="EM")
VG_cyno <-rownames(readRDS(paste0("tmp_data/",TD,"/data.cyno.EM.rds")) @tools$RunFastMN@assays@data$reconstructed)
#' check the st of cyno genes 
#cyno.genes.anno %>% filter(ensembl_gene_id %in% c(count_cyno_genes %>% filter(source=="ENS_ID") %>% pull(GeneID))) %>% filter(hsapiens_homolog_orthology_type=="ortholog_one2one") %>% pull(hsapiens_homolog_associated_gene_name) %>% setdiff(count_cyno_genes$GeneID)
table(count_cyno_genes %>% filter(GeneName %in% rownames(count_cyno)) %>% select(source,GeneName) %>% unique() %>% group_by(GeneName) %>% summarise(source=paste(source,collapse=":"))%>% pull(source))
# gene name has been transfered to human gene names 


#' length of overlapped genes
intersect(rownames(count_human),rownames(count_marmo)) %>% intersect(rownames(count_cyno)) %>% length()

#' get the overlap of human marmoset
GI.list <- list()
GI.list$hu_mar_cyno <- list()
GI.list$hu_mar_cyno$ov <- intersect(rownames(count_human),rownames(count_marmo)) %>% intersect(rownames(count_cyno))
GI.list$hu_mar_cyno$human <- rownames(count_human)
GI.list$hu_mar_cyno$marmoset <- rownames(count_marmo)
GI.list$hu_mar_cyno$cyno <- rownames(count_cyno)
GI.list$hu_mar_cyno$all <- c(GI.list$hu_mar_cyno$human,GI.list$hu_mar_cyno$marmoset,GI.list$hu_mar_cyno$cyno) %>% unique()


metas.list <- list()
metas.list$hu_mar_cyno <- meta_human%>% mutate(species="human") %>% bind_rows(meta_marmo %>% mutate(species="marmoset"))%>% bind_rows(meta_cyno %>% mutate(species="cyno"))
table(duplicated(metas.list$hu_mar_cyno$cell))
table(metas.list$hu_mar_cyno$pj)

#' down sample the Cyno_Yang_2021 from 6980 to 2000
metas.list$hu_mar_cyno <- metas.list$hu_mar_cyno %>% filter(pj!="Cyno_Yang_2021") %>% bind_rows(metas.list$hu_mar_cyno %>% filter(pj=="Cyno_Yang_2021") %>% FunMaSF(2000))
table(metas.list$hu_mar_cyno$pj)

counts.list <- list()
counts.list$hu_mar_cyno <- (cbind(count_marmo[GI.list$hu_mar_cyno$ov ,],count_human[GI.list$hu_mar_cyno$ov ,]) %>% cbind(count_cyno[GI.list$hu_mar_cyno$ov ,]))[,metas.list$hu_mar_cyno $cell]
counts.list$hu_mar_cyno %>% dim()

counts <- counts.list$hu_mar_cyno 
meta <- metas.list$hu_mar_cyno

#' decrease mem
rm(counts.list)
rm(count_marmo)
rm(count_human)
rm(count_cyno)

expG.set <- list()
for (b in unique(meta$pj  %>% unique() %>% as.vector())) {
  temp.cell <- meta %>% filter(pj==b) %>% pull(cell)
  expG.set[[b]] <- rownames(counts )[rowSums(counts[,temp.cell] >=1) >=5]
}
sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()



if (rewrite) {
  print("rewrite")
  saveRDS(counts,paste0("tmp_data/",TD,"/cyno.human.marmoset.merge.counts.rds"))
  saveRDS(meta,paste0("tmp_data/",TD,"/cyno.human.marmoset.merge.meta.rds"))
}else{
  
  print("do not rewrite")
  
}


if (file.exists(paste0("tmp_data/",TD,"/human_marmoset_cyno.data.list.Rdata"))) {
  load(paste0("tmp_data/",TD,"/human_marmoset_cyno.data.list.Rdata"),verbose = T)
}else{
  sce.ob <- list()
  for (b in unique(meta$pj  %>% unique() %>% as.vector())) {
    print(b)
    temp.M <- meta %>% filter(pj==b)
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    sce.ob[[b]] <- temp.sce
  }
  temp.mBN.sce.ob <- multiBatchNorm(sce.ob$Yan2013,sce.ob$Meistermann_2021,sce.ob$SPH2016,sce.ob$nBGuo,sce.ob$D3post,sce.ob$CS7,sce.ob$Cyno_Ma_2019,sce.ob$Cyno_Nakamura_2016,sce.ob$Cyno_Yang_2021,sce.ob$Mamo_Bergmann_2022,sce.ob$Mamo_Boroviak_2018)
  lognormExp.3S.onlyhomo..mBN<- temp.mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  saveRDS(lognormExp.3S.onlyhomo..mBN,paste0("tmp_data/",TD,"/lognormExp.3S.onlyhomo..mBN.rds"))
  
  
  temp.M <-  meta %>% filter(pj %in% c("Yan2013","CS7","D3post","SPH2016","nBGuo","Meistermann_2021","Cyno_Ma_2019","Cyno_Nakamura_2016","Cyno_Yang_2021","Mamo_Bergmann_2022","Mamo_Boroviak_2018"))
  temp.sel.cell <- temp.M$cell
  data.merge <- CreateSeuratObject(counts[sel.expG,temp.sel.cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  data.spt <- SplitObject(data.merge, split.by = "species")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
 
  temp <- data.frame(gene=VG_human) %>% tibble::rowid_to_column("n") %>% tbl_df() %>% bind_rows(data.frame(gene=VG_cyno) %>% tibble::rowid_to_column("n") %>% tbl_df())%>% bind_rows(data.frame(gene=VG_Marmoset) %>% tibble::rowid_to_column("n") %>% tbl_df()) %>% filter(gene %in% GI.list$hu_mar_cyno$ov) %>% group_by(gene) %>% summarise(nS=n(),mS=mean(n),minS=min(n))
  top.sel.features <- temp %>% filter(nS==3) %>% arrange(mS) %>% bind_rows(temp %>% filter(nS==2) %>% arrange(mS)) %>% bind_rows(temp %>% filter(nS==1) %>% filter(gene %in% VG_human) %>% arrange(mS) %>% head(500)) %>% bind_rows(temp %>% filter(nS==1) %>% filter(gene %in% VG_Marmoset) %>% arrange(mS) %>% head(500))%>% bind_rows(temp %>% filter(nS==1) %>% filter(gene %in% VG_cyno) %>% arrange(mS) %>% head(500)) %>% unique() %>% pull(gene)
  
  
  data.list <- list()
  sce.norm.list <- list()
  #Yan2013,Meistermann_2021,SPH2016,nBGuo,D3post,CS7
  sce.norm.list$human <- multiBatchNorm(sce.ob$Yan2013,sce.ob$Meistermann_2021,sce.ob$SPH2016,sce.ob$nBGuo,sce.ob$D3post,sce.ob$CS7)
  names(sce.norm.list$human) <- c("Yan2013","Meistermann_2021","SPH2016","nBGuo","D3post","CS7")
 # sce.temp.M1.norm <-  mnnCorrect(sce.norm.list$human$Yan2013,sce.norm.list$human$Meistermann_2021,sce.norm.list$human$nBGuo,sce.norm.list$human$SPH2016,sce.norm.list$human$D3post, sce.norm.list$human$CS7, subset.row=top.sel.features)
  sce.temp.M1.norm <-  mnnCorrect(sce.norm.list$human$Yan2013,sce.norm.list$human$SPH2016,sce.norm.list$human$D3post,sce.norm.list$human$nBGuo,sce.norm.list$human$Meistermann_2021, sce.norm.list$human$CS7, subset.row=top.sel.features) 
  
  sce.norm.list$cyno <- multiBatchNorm(sce.ob$Cyno_Ma_2019,sce.ob$Cyno_Nakamura_2016,sce.ob$Cyno_Yang_2021)
  names(sce.norm.list$cyno) <- c("Cyno_Ma_2019","Cyno_Nakamura_2016","Cyno_Yang_2021")
  sce.temp.M2.norm  <-  mnnCorrect(sce.norm.list$cyno$Cyno_Ma_2019,sce.norm.list$cyno$Cyno_Nakamura_2016,sce.norm.list$cyno$Cyno_Yang_2021, subset.row=top.sel.features) 
  
  sce.norm.list$Marmoset <- multiBatchNorm(sce.ob$Mamo_Bergmann_2022,sce.ob$Mamo_Boroviak_2018)
  names(sce.norm.list$Marmoset) <- c("Mamo_Bergmann_2022","Mamo_Boroviak_2018")
  sce.temp.M3.norm  <-  mnnCorrect(sce.norm.list$Marmoset$Mamo_Bergmann_2022,sce.norm.list$Marmoset$Mamo_Boroviak_2018, subset.row=top.sel.features) 
  rm(sce.ob)
  
  temp.M1 <- meta %>% filter(species=="human")
  temp.M2 <- meta %>% filter(species=="cyno")
  temp.M3 <- meta %>% filter(species=="marmoset")
  data.list <- list()
  temp <- CreateSeuratObject(counts[top.sel.features,temp.M1$cell], meta.data =(temp.M1%>% tibble::column_to_rownames("cell")))  %>% NormalizeData(verbose = FALSE)%>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  temp@assays$RNA@data <- sce.temp.M1.norm@assays@data$corrected[,colnames(temp)]
  data.list$M1  <- temp %>% ScaleData(verbose=F)%>% RunPCA(verbose=F,npcs=20) %>% RunUMAP(dims=1:20,verbose=F) %>% FindNeighbors( dims = 1:20,verbose = FALSE) %>%  FindClusters(resolution = 0.1,verbose = FALSE)
  
  
  temp <- CreateSeuratObject(counts[top.sel.features,temp.M2$cell], meta.data =(temp.M2%>% tibble::column_to_rownames("cell")))  %>% NormalizeData(verbose = FALSE)%>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  temp@assays$RNA@data <- sce.temp.M2.norm@assays@data$corrected[rownames(temp),colnames(temp)]
  data.list$M2 <- temp   %>% ScaleData(verbose=F)%>% RunPCA(verbose=F,npcs=20) %>% RunUMAP(dims=1:20,verbose=F) %>% FindNeighbors( dims = 1:20,verbose = FALSE) %>%  FindClusters(resolution = 0.1,verbose = FALSE)
  
  temp <- CreateSeuratObject(counts[top.sel.features,temp.M3$cell], meta.data =(temp.M3%>% tibble::column_to_rownames("cell")))  %>% NormalizeData(verbose = FALSE)%>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  temp@assays$RNA@data <- sce.temp.M3.norm@assays@data$corrected[rownames(temp),colnames(temp)]
  data.list$M3 <- temp   %>% ScaleData(verbose=F)%>% RunPCA(verbose=F,npcs=20) %>% RunUMAP(dims=1:20,verbose=F) %>% FindNeighbors( dims = 1:20,verbose = FALSE) %>%  FindClusters(resolution = 0.1,verbose = FALSE)
  
  data.list <- data.list[c("M1","M2","M3")]
  save(data.list,top.sel.features,file=paste0("tmp_data/",TD,"/human_marmoset_cyno.data.list.Rdata"))
}

if (file.exists(paste0("tmp_data/",TD,"/human_marmoset_cyno.data.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/human_marmoset_cyno.data.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/human_marmoset_cyno.data.ob.umap.rds"))
}else{
  verbose=F;pc2 <- 20
  data.merge.features <-top.sel.features
  
  data.merge.anchors <- FindIntegrationAnchors(object.list = data.list[c("M1","M3","M2")],  anchor.features = data.merge.features, verbose = FALSE)
  data.merge <- IntegrateData(anchorset = data.merge.anchors, verbose = FALSE)
  DefaultAssay(object = data.merge) <- "integrated"
  data.ob <- data.merge %>% ScaleData(verbose = FALSE)%>% RunPCA( verbose = FALSE)  %>%  RunUMAP(dims = 1:pc2,verbose = FALSE)
  data.temp <- data.ob%>% FindNeighbors( dims = 1:pc2, verbose = FALSE,nn.method="annoy",annoy.metric="cosine") %>% FindClusters(resolution = 1, verbose =  FALSE)
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:LOC)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))%>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_25),by="cell")
  
  saveRDS(data.ob,paste0("tmp_data/",TD,"/human_marmoset_cyno.data.ob.rds"))
  saveRDS(data.ob.umap,paste0("tmp_data/",TD,"/human_marmoset_cyno.data.ob.umap.rds"))
  
  #' try the different parameters of PCs
  data.ob.umap.list <- list()
  for (nPC in c(15,20,25,30,35,40,50)) {
    data.temp <- data.ob %>% RunUMAP(1:nPC, verbose = FALSE) 
    data.ob.umap.list[[paste0("PC",nPC)]] <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% mutate(nPC_sel=paste0("top:",nPC))
  }
  data.ob.umap.list <- data.ob.umap.list %>% do.call("bind_rows",.)
  saveRDS(data.ob.umap.list,paste0("tmp_data/",TD,"/human_marmoset_cyno.data.ob.umap.PC.para.list.rds"))
}
