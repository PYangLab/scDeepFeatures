#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(glue)
  library(stringr)
  source("../../Utils/utils.R")
  source("../../Utils/Feature_selection_methods/Traditional/doDE.R")
  source("../../Utils/Classification_methods/doEval.R")
  #print(doTtest)
})
 
method = "runDESeq2" ######
doFS = doDESeq2_mc #####


# ######################### 
# # Setting 2
n=1
j = 1
while(j<=10){
  for(nGroups in c(10,20)){#seq(5,20,5)
    for (nCells in seq(50,250,50)){
      cat(paste0(glue("doing feature selection {method}: ...nGroups="), nGroups, "...iteration:", j,"\n"))
      sce = readRDS(glue("../../Data/Simulation/Tabula_muris/Setting2/nGroups_{nGroups}/nCells_{nCells}/Simulation_setting2_nGroups_{nGroups}_nCells_{nCells}_iteration{j}.rds"))
      rownames(sce) = str_replace_all(rownames(sce), "_", "-")
      
      res.split <- train_test_split(rna = sce, cty = sce$cellTypes, seedInput = 1234)
      train.rna = res.split$train.rna
      train.cty = res.split$train.cty
      test.rna = res.split$test.rna
      test.cty = res.split$test.cty
      
      ## filter lowly expressed genes, if it express less than 1% in within every cell types 
      exprsMat <- logcounts(train.rna)
      exprs_pct <- 0.01
      
      label = train.rna$cellTypes
      cts = unique(train.rna$cellTypes)
      
      print(all(train.cty == label))
      meanPct.list <- list()
      for(i in 1:length(cts)) {
        idx <- which(label == cts[i])
        meanPct.list[[i]] <- (Matrix::rowSums(exprsMat[, idx, drop = FALSE] > 0)/sum(label == cts[i])) > exprs_pct
      }
      names(meanPct.list) <- cts
      keep = rowSums(do.call(cbind, meanPct.list)) > 0 # remove any genes that are expressed lower than 5% in each and every cell type
      
      train.rna = train.rna[keep, ]
      ###

      res <- doFS(train.rna, train.cty) #####
      dir.create(glue("../../Result/Fs_marker/Traditional/Setting2"),showWarnings = FALSE)
      dir.create(glue("../../Result/Fs_marker/Traditional/Setting2/{method}"),showWarnings = FALSE)
      saveRDS(res, file = glue("../../Result/Fs_marker/Traditional/Setting2/{method}/Res_{method}_Simulation_setting2_nGroups_{nGroups}_nCells_{nCells}_iteration{j}.rds"))
      cat(glue("finished {n} datasets \n"))
    }
    n=n+1
  }
  j=j+1
}
