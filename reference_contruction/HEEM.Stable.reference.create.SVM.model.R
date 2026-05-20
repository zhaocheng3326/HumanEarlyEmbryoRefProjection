#' ---
#' title: "create the nDim UMAP model and SVM model for cell classification"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---

# ##loading R library
rm(list=ls())
rewrite=FALSE


suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  #library(scran)
  library(ff)
  #library(batchelor)
  library(Seurat)
  #library(SeuratWrappers)
  #library(scuttle)
  #library(SeuratDisk)
  library(uwot)
  library(caret)
  library(e1071)
})
source("src/project.setting.R") #### please modify this path as needed.


# model file
model_file <- paste0(DIR,"/tmp_data/",TD,"/stable.ref.fastMNN.umap.model")
model_db_file <- paste0(DIR,"/tmp_data/",TD,"/stable.ref.fastMNN.Rdata")
model_ref_file <- paste0(DIR,"/tmp_data/",TD,"/GS.stable.ref.umap.rds") ### with updated cell annotation
model_ndim_umap_file <- paste0(DIR,"/tmp_data/",TD,"/stable.ref.fastMNN.umap.nDimN.model")

#' loading model files
load(model_db_file,verbose=T)
ref.svm.nDim <- 20

# there is some random factors , careful to rewrite
if (rewrite) {
  ref_nDim_umap_model <- umap(as.data.frame(ref.mnn.out[,c(1:ref.nPC)]), ret_model = TRUE, n_neighbors = 30L, metric = "cosine",min_dist = 0.3,n_components=ref.svm.nDim)
  save_uwot(ref_nDim_umap_model, file = model_ndim_umap_file)
}

# there is some random factors , careful to rewrite
if (rewrite) {
  #' loading nDim UMAP model
  ref_nDim_umap_model <- load_uwot(file = model_ndim_umap_file)
  
  #' updated reference annotation
  ref.umap.updated.anno <- readRDS(model_ref_file) %>% mutate(rename_EML=ifelse(rename_EML %in% c("Amnion.Ecto","Late_Amnion"),"Amnion",rename_EML))%>% mutate(sub_rename_EML=ifelse(sub_rename_EML %in% c("Amnion.Ecto","Late_Amnion"),"Amnion",sub_rename_EML)) %>%  select(cell,rename_EML,sub_rename_EML) %>% mutate(cell=paste0("ref_",cell))
  
  ref.meta <- ref.meta %>% select(-rename_EML) %>% left_join(ref.umap.updated.anno,by="cell") %>% mutate(rename_EML=recode(rename_EML,"Amnion.Ecto"="Amnion")) %>% mutate(cluster_EML=ifelse(sub_rename_EML %in% c("Late_EPI","Early_EPI","Late_Hypoblast","Early_Hypoblast"),sub_rename_EML,rename_EML))%>% mutate(cluster_EML=ifelse(cluster_EML %in% c("Ambiguous","Unknown","EPI.PrE.INT"),"Ambiguous",cluster_EML)) %>% mutate(cluster_EML=ifelse(cluster_EML %in% c("Zygote","2-4 cell"),"Z4cell",cluster_EML)) %>% mutate(cluster_EML=recode(cluster_EML,"Axial Mes"="AxMes","8 cell"="EightCell"))
  
  
  #' cluster information
  ref.cluster <- readRDS(paste0("tmp_data/",TD,"/stable.fastMNN.ref.cluster.rds"))
  set.seed(ref.seed)
  
  #' stable prediction
  stable.ref.umap <- readRDS(paste0("tmp_data/",TD,"/GS.stable.ref.umap.rds")) %>% select(-rename_EML) %>% left_join(ref.meta %>% select(cell,rename_EML,cluster_EML) %>% mutate(cell=gsub("^ref_","",cell)),by="cell")
  #' rotation
  stable.ref.umap <-  stable.ref.umap %>% rename(UMAP_2=UMAP_1, UMAP_1=UMAP_2)
  
  ref.umap <- stable.ref.umap  %>% mutate(cell=paste0("ref_",cell))
  
  #' filter uncertained cells before training models
  ref.further.exclude.cells <- readRDS(paste0("tmp_data/",TD,"/ref.further.exclude.cells.rds")) %>% paste0("ref_",.) ## with uncertained lineage annotation
  
  ref.umap.filter <- ref.umap  %>% filter(!cell %in% ref.further.exclude.cells)
  
  #' training for SVM models
  set.seed(123)
  train_data <- ref_nDim_umap_model$embedding %>% as.data.frame() %>% setNames(paste("Dim",1:ref.svm.nDim,sep="_")) %>% tibble::rownames_to_column("cell") %>% inner_join(ref.umap.filter %>% select(cell,cluster_EML),by="cell") %>% filter(!cluster_EML %in% c("Ambiguous"))  
  
  svm_model <- list()
  for (sel_EML in unique(train_data$cluster_EML)) {
    print(paste("nDim",ref.svm.nDim,sel_EML,sep=":"))
    svm_model[[sel_EML]] <- FunSelEML_SVM_opt_model(train_data,ref.svm.nDim,sel_EML)
  }
  
  # saving object
  saveRDS(svm_model,paste0("tmp_data/",TD,"/stable.ref.fastMNN.umap.nDim.svm.model.rds"))
}
