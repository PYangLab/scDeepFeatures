#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(glue)
  source("../../Utils/utils.R")
  source("../../Utils/Feature_selection_methods/Traditional/doDE.R")
  source("../../Utils/Classification_methods/doEval.R")
  #print(doTtest)
})
 
method = "runScDD" ######
doFS = doScDD_all #####

######################### 
# Setting 4
n=1
nCells = 200
nGroups =10
j = 1
while(j<=10){
  cat(paste0("converting...nGroups=", nGroups,"...nCells=", paste0(nCells,collapse = ":"),"...iteration:", j,"\n"))
  # read subsampled data
  sce = readRDS(glue("../../Data/Simulation/Tabula_muris/Setting4/Simulation_setting4_seed_{n}_iteration{j}.rds"))
  rownames(sce) = str_replace_all(rownames(sce), "_", "-") 
  # split
  res.split = train_test_split(rna = sce, cty = sce$cellTypes, seedInput = 1234)
  train.rna = res.split$train.rna
  train.cty = res.split$train.cty
  test.rna = res.split$test.rna
  test.cty = res.split$test.cty
  
  
  ### filter lowly expressed genes, if it express less than 1% in within every cell types 
  # exprsMat <- logcounts(train.rna)
  # exprs_pct <- 0.01
  
  # label = train.rna$cellTypes
  # cts = unique(train.rna$cellTypes)
  
  # meanPct.list <- list()
  # for(i in 1:length(cts)) {
  #   idx <- which(label == cts[i])
  #   meanPct.list[[i]] <- (Matrix::rowSums(exprsMat[, idx, drop = FALSE] > 0)/sum(label == cts[i])) > exprs_pct
  # }
  # names(meanPct.list) <- cts
  # keep = rowSums(do.call(cbind, meanPct.list)) > 0 # remove any genes that are expressed lower than 5% in each and every cell type
  
  # train.rna = train.rna[keep, ]
  ###

  res <- doFS(train.rna, train.cty) #####
  dir.create(glue("../../Result/Fs_marker/Traditional/Setting4"),showWarnings = FALSE)
  dir.create(glue("../../Result/Fs_marker/Traditional/Setting4/{method}"),showWarnings = FALSE)
  saveRDS(res, file = glue("../../Result/Fs_marker/Traditional/Setting4/{method}/Res_{method}_Simulation_setting4_seed_{n}_iteration{j}_all.rds"))
  cat(glue("finished {j} datasets \n"))
  j=j+1
  }