#' here is the setting for the whole project

#' check whether in local computer
if (grepl("KI-",Sys.info()['nodename'])) {
  print("local computer")
  #source("/Users/cheng.zhao/chzhao_bioinfo/PC/SnkM/SgCell.R")
  base_dir <- "/Users/cheng.zhao/Documents"
} else {
  print("On server")
  condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
  LBpath <- paste0(condaENV ,"/lib/R/library")
  .libPaths(LBpath)
  base_dir="/home/chenzh"
}

# working directory
DIR <- "~/My_project/HumanEarlyEmbryoRefProjection"
knitr::opts_knit$set(root.dir=DIR)
setwd(DIR)

#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter
#st_sfc <- sf::st_sfc


options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="Jan_2024"