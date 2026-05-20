
FunRestricted_mnn <- function(left.data,  right.data, k=20) {
  # modified from https://rdrr.io/bioc/batchelor/src/R/MNN_tree.R#sym-.restricted_mnn
  k1 <- min(k,  nrow(left.data))
  k2 <- min(k, nrow(right.data))
  pairs <- findMutualNN(left.data, right.data, k1=k1, k2=k2)
  pairs$first <- pairs$first
  pairs$second <- pairs$second
  pairs
}

FunStoufferP = function(x) {
  return(data.frame(merge_p_val=poolr::stouffer(x)$p))
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

FunSmooth_gaussian_kernel <-  function(averaged, index, mat, sigma2) {
  # modified from https://rdrr.io/bioc/batchelor/src/R/RcppExports.R
  .Call('_batchelor_smooth_gaussian_kernel', PACKAGE = 'batchelor', averaged, index, mat, sigma2)
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



FunMaSF <- function(temp,n) {
  set.seed(456)
  return(head(temp[sample(nrow(temp)),],n))
}

FunSelEML_SVM_opt_model <- function(train_data,nDim,sel_EML) {
  train_data_sub <- train_data %>% mutate(cluster_EML=ifelse(cluster_EML %in% sel_EML,cluster_EML,"other")) %>% mutate(cluster_EML=factor(cluster_EML)) %>% tibble::column_to_rownames("cell")
  train_ctrl <- trainControl(method="cv",number=5,classProbs =  TRUE)
  train_grid <- expand.grid(sigma = c(0.1, 0.5, 1), C = c(0.1, 1, 10))
  svm_test_model <- train(cluster_EML ~ ., data = train_data_sub, method = "svmRadial", trControl = train_ctrl, tuneGrid = train_grid,,metric="Kappa")
  #model_bestParameters <- svm_test_model$bestTune
  #finalModel <- svm(cluster_EML ~ ., data = train_data_sub, method = "svmRadial", cost = model_bestParameters$C, gamma = model_bestParameters$sigma)
  temp.out <- list()
  #temp.out$model <- svm_test_model$finalModel
  temp.out$model <- svm_test_model
  temp.out$para <- svm_test_model$bestTune
  temp.out$perf <- svm_test_model$results %>% tbl_df() %>% filter(sigma==svm_test_model$bestTune$sigma & C==svm_test_model$bestTune$C) %>% mutate(sel_EML=sel_EML,sel_nDim=nDim)
  return(temp.out)
}
