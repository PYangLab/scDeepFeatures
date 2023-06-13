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
 
method = "runLimmavoom" ######
doFS = doLimmavoom #####

######################### 
# Setting 1
n=1 # seed number 
j=1 # repeat times
while(j<=10){
  for(nGroups in c(5, 10, 15, 20)){ # cell type number
    cat(paste0(glue("doing feature selection {method}: ...nGroups="), nGroups, "...iteration:", j,"\n"))
    # read simulated data
    sce = readRDS(glue("../../Data/Simulation/Tabula_muris/Setting1/Simulation_setting1_nCells_200_nGroups_{nGroups}__iteration{j}_seeds{n}.rds"))
    #rownames(sce) = str_replace_all(rownames(sce), "_", "-") 
    # split
    res.split = train_test_split(rna = sce, cty = sce$cellTypes, seedInput = 1234)
    train.rna = res.split$train.rna
    train.cty = res.split$train.cty
    test.rna = res.split$test.rna
    test.cty = res.split$test.cty
    
    ### filter lowly expressed genes, if it express less than 1% in within every cell types 
    exprsMat <- logcounts(train.rna)
    exprs_pct <- 0.01
    
    label = train.rna$cellTypes
    cts = unique(train.rna$cellTypes)
    
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
    dir.create(glue("../../Result/Fs_marker/Traditional/Setting1"),showWarnings = FALSE)
    dir.create(glue("../../Result/Fs_marker/Traditional/Setting1/{method}"),showWarnings = FALSE)
    saveRDS(res, file = glue("../../Result/Fs_marker/Traditional/Setting1/{method}/Res_{method}_Simulation_setting1_nCells_200_nGroups_{nGroups}__iteration{j}_seeds{n}.rds"))
    cat(glue("finished {n} datasets \n"))
    n=n+1
    
  }
  j=j+1
}
