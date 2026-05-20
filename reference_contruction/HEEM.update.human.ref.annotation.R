#' ---
#' title: "update human annotation"
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
  library(tidyr)
  library(pheatmap)
  library(tibble)
  #library(scran)
  #library(batchelor)
  #library(SeuratWrappers)
})

source("src/project.setting.R") #### please modify this path as needed.
source("src/local.quick.fun.R")


suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)


#' loading the ref cluster results
ref.cluster <- readRDS(paste0("tmp_data/",TD,"/stable.fastMNN.ref.cluster.rds"))

pj.dt <- data.frame(pj=c("CS7","Cyno_Ma_2019","Cyno_Nakamura_2016","Cyno_Yang_2021","D3post","Mamo_Bergmann_2022","Mamo_Boroviak_2018","Meistermann_2021","nBGuo","SPH2016","Yan2013"), pup=c("Tyser et al","Ma et al", "Nakamura et al","Yang et al","Xiang et al","Bergmann et al","Boroviak et al","Meistermann et al","Yanagida et al","Petropoulos et al","Yan et al"))

#' update AdvMes to ExE_Mes (loading update CS7 annotation)
#' update D3post from warmflash annotation
#' update SPH2016 using the Meistermann et al annotation
#' update SPH2016 ICM cells by Stirparo et al annotation and cluster information
meta.filter.updated <- readRDS(paste0("tmp_data/",TD,"/meta.filter.updated.rds")) %>% filter(!cellType %in% c("Blastoid","EPSC")) %>% select(cell,EML)  %>% filter(EML!="Oocyte")

#' SPH2016 ICM cells
SPH2016.ICM.cells <- readRDS(file=paste0("tmp_data/",TD,"/SPH2016.ICM.cells.rds")) %>% gsub("ref_","",.)
stable.ref.umap <- stable.ref.umap  %>% mutate(rename_EML=ifelse(cell %in% SPH2016.ICM.cells, "ICM",rename_EML))


#' stable reference UMAP
stable.ref.umap <- readRDS(paste0("tmp_data/",TD,"/stable.ref.umap.rds")) %>% mutate(cell=gsub("^ref_","",cell)) %>% select(-rename_EML) %>% rows_update(meta.filter.updated ,by=c("cell")) %>% left_join(rename.dt,by="EML") %>% mutate(rename_EML=ifelse(is.na(rename_EML),EML,rename_EML))%>% mutate(rename_EML=ifelse(devTime=="E3","8 cell",rename_EML))%>% mutate(rename_EML=ifelse(devTime=="E4","Morula",rename_EML))# %>% mutate(rename_EML=ifelse(pj=="CS7" & subCT=="PGC","PGC",rename_EML)) 

#' CS7 sub cluster of endoderm 
temp.CS7.M <- stable.ref.umap %>% filter(pj=="CS7" & EML=="Endoderm") %>% mutate(rename_EML=subCT) %>% mutate(rename_EML=recode(rename_EML,"DE(NP)"="DE","DE(P)"="DE","YS Endoderm"="YSE"))
stable.ref.umap <- stable.ref.umap %>% rows_update(temp.CS7.M %>% select(cell,rename_EML),by="cell") 
#' add pup
stable.ref.umap <- stable.ref.umap %>% inner_join(pj.dt,by="pj")

#' updated metas
saveRDS(stable.ref.umap,file=paste0("tmp_data/",TD,"/GS.stable.ref.umap.rds"))
