# this file stores DE methods operate in R
# load packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(stringr)
})




# limma-voom, adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
doLimmavoom <- function(sce, cellTypes) { 
  # input must be , log-transformed data
  time.start = Sys.time()
  suppressPackageStartupMessages({library(limma)})
  suppressPackageStartupMessages({library(edgeR)})
  cat("do Limma-voom test \n")
  cat("input must be sce object with raw count")
  countMatrix = as.matrix(counts(sce))
  
  cty <- droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cty)) {
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 0))
    
    design <- stats::model.matrix(~tmp_celltype)
    
    dge <- edgeR::DGEList(countMatrix, group = tmp_celltype)
    dge <- edgeR::calcNormFactors(dge)
    vm <- limma::voom(dge, design = design, plot = FALSE)
    fit <- limma::lmFit(vm, design = design)
    
    #y <- methods::new("EList")
    #y$E <- exprsMat[keep, ]
    #y$E <- exprsMat
    #vm <- limma::voom(y, design = design)
    #fit <- limma::lmFit(vm, design = design)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
    if (!is.null(tt[[i]]$ID)) {
      tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
      rownames(tt[[i]]) <- tt[[i]]$ID
    }
    tt[[i]] = tt[[i]][order(tt[[i]][,"t"], decreasing = TRUE), ]
    
  }
  names(tt) <- levels(cty)
  time.end = Sys.time()
  time_elapse <- as.numeric(difftime(time.end, time.start, units = "secs"))
  to_return <- list(result = tt, time = time_elapse)
  return(to_return)
}


# method 11 scDD, adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
# and https://www.bioconductor.org/packages/release/bioc/vignettes/scDD/inst/doc/scDD.pdf
doScDD <- function(sce, cellTypes, worker=48, log2FCshre = 0){
  time.start = Sys.time()
  suppressPackageStartupMessages({
    library(scDD)
    library(BiocParallel)
  })
  cat("do scDD test \n")
  cat("input must be sce object with log-normalised count")
  
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  if(Sys.info()[1] == "Windows"){
    param=BiocParallel::SnowParam(workers=worker)
  } else {
    param=BiocParallel::MulticoreParam(workers=worker)
  }
  assayNames(sce) <- c("counts", "normcounts")
  
  cty = droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cty)){
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 2))
    names(tmp_celltype) <- colnames(sce)
    sce$condition <- tmp_celltype
    sce <- scDD(sce, prior_param = prior_param,
                condition = "condition", param = param,
                testZeroes = TRUE, categorize = FALSE)
    res <- results(sce)
    res_dd <- res
    #res_dd <- res[res[, "DDcategory"] %in% c("DE", "DP", "DM", "DB", "DZ"), ] # only return DD genes
    mat <- as.matrix(normcounts(sce))
    logFC <- (rowMeans(mat[,sce$condition == 1]) - rowMeans(mat[,sce$condition == 2]))
    logFC_pos_gene <- rownames(sce)[logFC > log2FCshre | logFC == Inf]
    res_dd <- res_dd[rownames(res_dd) %in% logFC_pos_gene,]
    res_dd <- res_dd[order(res_dd[,"combined.pvalue.adj"], decreasing = FALSE), ]
    
    tt[[i]] <- res_dd
    
  }
  names(tt) <- levels(cty)
  time.end = Sys.time()
  time_elapse <- as.numeric(difftime(time.end, time.start, units = "secs"))
  to_return <- list(result = tt, time = time_elapse)
  return(to_return)
}

doScDD_all <- function(sce, cellTypes, worker=48, log2FCshre = 0){
  # for reproducibility testing
  time.start = Sys.time()
  suppressPackageStartupMessages({
    library(scDD)
    library(BiocParallel)
  })
  cat("do scDD test \n")
  cat("input must be sce object with log-normalised count")
  
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  if(Sys.info()[1] == "Windows"){
    param=BiocParallel::SnowParam(workers=worker)
  } else {
    param=BiocParallel::MulticoreParam(workers=worker)
  }
  assayNames(sce) <- c("counts", "normcounts")
  
  cty = droplevels(as.factor(cellTypes))
  tt <- list()
  for (i in 1:nlevels(cty)){
    tmp_celltype <- (ifelse(cty == levels(cty)[i], 1, 2))
    names(tmp_celltype) <- colnames(sce)
    sce$condition <- tmp_celltype
    sce <- scDD(sce, prior_param = prior_param,
                condition = "condition", param = param,
                testZeroes = TRUE, categorize = FALSE)
    res <- results(sce)
    res_dd <- res
    #res_dd <- res[res[, "DDcategory"] %in% c("DE", "DP", "DM", "DB", "DZ"), ] # only return DD genes
    mat <- as.matrix(normcounts(sce))
    res_dd <- res_dd[order(res_dd[,"combined.pvalue.adj"], decreasing = FALSE), ]
    
    tt[[i]] <- res_dd
    
  }
  names(tt) <- levels(cty)
  time.end = Sys.time()
  time_elapse <- as.numeric(difftime(time.end, time.start, units = "secs"))
  to_return <- list(result = tt, time = time_elapse)
  return(to_return)
}


# adapt from https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
doDESeq2 <- function(sce, cellTypes) {
  # input must be raw count
  time.start = Sys.time()
  suppressPackageStartupMessages({
    library(DESeq2)
  })
  cat("do DEseq2 test \n")
  cat("input must be sce object with raw count")
  exprsMat = as.matrix(counts(sce)+1)
  cty <- droplevels(as.factor(cellTypes))
  
  tt <- list()
  for (i in 1:nlevels(cty)) {
    tmp_celltype <- factor(ifelse(cty == levels(cty)[i], "anchor", "others"), levels=c("others", "anchor"))
    dds <- DESeqDataSetFromMatrix(countData = exprsMat, 
                                  colData = data.frame(condition=tmp_celltype), 
                                  design = ~condition)
    dds <- DESeq(dds)
    tmp_res <- results(dds)
    tmp_res <- tmp_res[order(tmp_res[,"stat"], decreasing = TRUE), ]
    tt[[i]] <- tmp_res
  }
  names(tt) <- levels(cty)
  time.end = Sys.time()
  time_elapse <- as.numeric(difftime(time.end, time.start, units = "secs"))
  to_return <- list(result = tt, time = time_elapse)
  return(to_return)
}

doDESeq2_mc <- function(sce, cellTypes){ # equivilent implementation of previous function
  time.start = Sys.time()
  suppressPackageStartupMessages({library(DESeq2)})
  cat("do DEseq2 test \n")
  cat("input must be sce object with raw count")
  exprsMat = as.matrix(counts(sce)+1)
  cty <- droplevels(as.factor(cellTypes))
  time.end = Sys.time()
  #res_list = list()
  #for (i in 1:nlevels(cty)) {
  res_list_final <- list()
  res_list <- mclapply(c(1:nlevels(cty)), mc.cores = 10, function(i){
    tmp_time.start = Sys.time()
    tmp_celltype <- factor(ifelse(cty == levels(cty)[i], "anchor", "others"), levels=c("others", "anchor"))
    dds <- DESeqDataSetFromMatrix(countData = exprsMat, 
                                  colData = data.frame(condition=tmp_celltype), 
                                  design = ~condition)
    dds <- DESeq(dds)
    tmp_res <- results(dds)
    tmp_res <- tmp_res[order(tmp_res[,"stat"], decreasing = TRUE), ]
    tmp_time.end = Sys.time()
    tmp_time = as.numeric(difftime(tmp_time.end, tmp_time.start, units = "secs"))
    system(sprintf('echo "\n%s\n"', paste0("checkpoint 2", collapse="")))
    return(list(tmp_res, tmp_time))
  })
  
  
  res_list_final <- lapply(res_list, function(j){
    append(res_list_final, list(j[[1]]))
  })
  names(res_list_final) <- levels(cty)
  #}
  time_vector <- list()
  time_vector <- lapply(res_list, function(k){
    append(time_vector, k[[2]])
  })
  time_elapse = sum(unlist(time_vector)) + as.numeric(difftime(time.end, time.start, units = "secs"))
  #time.end = Sys.time()
  #time_elapse <- as.numeric(difftime(time.end, time.start, units = "secs"))
  to_return <- list(result = res_list_final, time = time_elapse)
  return(to_return)
}



























