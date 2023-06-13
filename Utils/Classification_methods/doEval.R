suppressPackageStartupMessages({
  library(caret)
  library(doParallel)
  library(e1071)
  
  library(scran)
  library(bluster)
  library(scater)
})

# First method, KNN classifier
doKKnn_simple_union <- function(train.sce, train.cellTypes, test.sce, test.cellTypes, anchor = NULL){
  cat("do KkNN union classification \n")
  
  train.exprsMat <- as.data.frame(t(logcounts(train.sce)))
  train.cty <- droplevels(as.factor(train.cellTypes))
  test.exprsMat <- as.data.frame(t(logcounts(test.sce)))
  test.cty <- droplevels(as.factor(test.sce$cellTypes))

  train.tmp_cellType <- train.cty
  train.tmp_mat <- cbind(train.tmp_cellType, train.exprsMat)
  train.tmp_mat$train.tmp_cellType <- factor(train.tmp_mat$train.tmp_cellType)
  colnames(train.tmp_mat)[1] = "con"
  
  
  # predict
  test.tmp_cellType <- test.cty
  test.tmp_mat <- cbind(as.factor(test.tmp_cellType), test.exprsMat)
  colnames(test.tmp_mat)[1] <- "con"
  
  fits = kknn(con ~., train = train.tmp_mat, test = test.tmp_mat, k =7, distance = 2)
  
  # record
  res = confusionMatrix(fits$fitted.values, reference= test.tmp_mat$con)
  res.1 <- as.matrix(res, what="classes")
  res.2 <- as.matrix(res, what="overall")
  tt <- list(res.2, res.1)
  #}
  #names(tt) <- levels(test.cty)
  return(tt)
}

# 2nd method, SVM
######
doSvm_simple_union <- function(train.sce, train.cellTypes, test.sce, test.cellTypes,  anchor = NULL, worker=10){
  cat("do SVM union classification \n")
  
  train.exprsMat <- as.data.frame(t(logcounts(train.sce)))
  train.cty <- droplevels(as.factor(train.cellTypes))
  test.exprsMat <- as.data.frame(t(logcounts(test.sce)))
  test.cty <- droplevels(as.factor(test.sce$cellTypes))
  
  train.tmp_cellType <- train.cty
  train.tmp_mat <- cbind(train.tmp_cellType, train.exprsMat)
  train.tmp_mat$train.tmp_cellType <- factor(train.tmp_mat$train.tmp_cellType)
  colnames(train.tmp_mat)[1] = "con"
  
  
  # predict
  test.tmp_cellType <- test.cty
  test.tmp_mat <- cbind(as.factor(test.tmp_cellType), test.exprsMat)
  colnames(test.tmp_mat)[1] <- "con"
  
  classifier = svm(formula = con ~., data = train.tmp_mat, type = "C-classification", kernel = "linear", scale=FALSE)
  y_pred = predict(classifier, newdata = test.tmp_mat[-1])
  # record
  res = confusionMatrix(y_pred, reference =  test.tmp_mat$con)
  res.1 <- as.matrix(res, what="classes")
  res.2 <- as.matrix(res, what="overall")
  tt <- list(res.2, res.1)
  #}
  #names(tt) <- levels(test.cty)
  return(tt)
}
