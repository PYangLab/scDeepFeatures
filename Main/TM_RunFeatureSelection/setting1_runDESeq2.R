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
doFS = doDESeq2 #####

######################### 
# Setting 1
n=1 # seed number 
j=1 # repeat times
while(j<=10){
  for(nGroups in c(5, 10, 15, 20)){ # cell type number
    cat(paste0(glue("doing feature selection {method}: ...nGroups="), nGroups, "...iteration:", j,"\n"))
    # read simulated data
    sce = readRDS(glue("../../Data/Simulation/Tabula_muris/Setting1/Simulation_setting1_nCells_200_nGroups_{nGroups}__iteration{j}_seeds{n}.rds"))
    rownames(sce) = str_replace_all(rownames(sce), "_", "-") 
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

# ######################### 
# # Setting 2
# n=1
# for(nGroups in seq(5,20,5)){
#   for(nCells in seq(50,250,50)){
#     j=1
#     while(j<=10){
#       cat(paste0(glue("doing feature selection {method}: ...nGroups="), nGroups,"...nCells=", nCells,"...iteration:", j,"\n"))
#       # read simulated data
#       sce = readRDS(glue("./Simulation/Tabula_muris/Setting2/Simulation_setting2_nCells_{nCells}_nGroups_{nGroups}__iteration{j}_seeds{n}.rds"))
#       # split
#       res.split = train_test_split(rna = sce, cty = sce$cellTypes, seedInput = 1234)
#       train.rna = res.split$train.rna
#       train.cty = res.split$train.cty
#       test.rna = res.split$test.rna
#       test.cty = res.split$test.cty
#       
#       dir.create(glue("../../Result/Fs_marker/Traditional/Setting2"),showWarnings = FALSE)
#       dir.create(glue("../../Result/Fs_marker/Traditional/Setting2/{method}"),showWarnings = FALSE)
#       res <- doFS(train.rna, train.cty) #####
#       
#       saveRDS(glue("../../Result/Fs_marker/Traditional/Setting2/{method}/Res_{method}_Simulation_setting2_nCells_{nCells}_nGroups_{nGroups}__iteration{j}_seeds{n}.rds"))
#       n=n+1
#       j=j+1
#     }
#   }
#   
# }
# 
# ######################### 
# # Setting 3
# n=1
# 
# for(nGroups in c(10,20)){#seq(5,20,5)
#   imbratio = c(2,4,10)##
#   imbratio_each_num = 1##
#   cellratio=list(c(rep(200,nGroups/2),rep(100,nGroups/2)),c(rep(200,nGroups/2),rep(50,nGroups/2)),c(rep(200,nGroups/2),rep(20,nGroups/2)))
#   for(nCells in cellratio){
#     imbratio_each = imbratio[imbratio_each_num]##
#     j=1
#     while(j<=10){
#       cat(paste0(glue("doing feature selection {method}: ...nGroups="), nGroups,"...nCells=", paste0(nCells,collapse = ":"),"...iteration:", j,"\n"))
#       # read simulated data
#       sce = readRDS(glue("./Simulation/Tabula_muris/Setting3/Simulation_setting3_nCellsMajor_200_nGroups_{nGroups}_imbalanceRatio_1:{imbratio_each}__iteration{j}_seeds{n}.rds"))
#       # split
#       res.split = train_test_split(rna = sce, cty = sce$cellTypes, seedInput = 1234)
#       train.rna = res.split$train.rna
#       train.cty = res.split$train.cty
#       test.rna = res.split$test.rna
#       test.cty = res.split$test.cty
#       
#       dir.create(glue("../../Result/Fs_marker/Traditional/Setting3"),showWarnings = FALSE)
#       dir.create(glue("../../Result/Fs_marker/Traditional/Setting3/{method}"),showWarnings = FALSE)
#       res <- doFS(train.rna, train.cty) #####
#       
#       saveRDS(glue("../../Result/Fs_marker/Traditional/Setting3/{method}/Res_{method}_Simulation_setting3_nCellsMajor_200_nGroups_{nGroups}_imbalanceRatio_1:{imbratio_each}__iteration{j}_seeds{n}.rds"))
#       
#       n=n+1
#       j=j+1
#     }
#     imbratio_each_num = imbratio_each_num + 1 
#   }
#   
# }
# 
# ######################### 
# # Setting 4
# n=1
# for(nGroups in c(5,10)){
#   for(nCells in c(500)){
#     j=1
#     while(j<=10){
#       cat(paste0(glue("doing feature selection {method}: ...nGroups="), nGroups,"...nCells=", nCells,"...iteration:", j,"\n"))
#       # read simulated data
#       sce = readRDS(glue("./Simulation/Tabula_muris/Setting4/Simulation_setting4_nCells_{nCells}_nGroups_{nGroups}__iteration{j}_seeds{n}.rds"))
#       # split
#       res.split = train_test_split(rna = sce, cty = sce$cellTypes, seedInput = 1234)
#       train.rna = res.split$train.rna
#       train.cty = res.split$train.cty
#       test.rna = res.split$test.rna
#       test.cty = res.split$test.cty
#       
#       dir.create(glue("../../Result/Fs_marker/Traditional/Setting4"),showWarnings = FALSE)
#       dir.create(glue("../../Result/Fs_marker/Traditional/Setting4/{method}"),showWarnings = FALSE)
#       
#       res <- doFS(train.rna, train.cty) #####
#       saveRDS(glue("../../Result/Fs_marker/Traditional/Setting4/{method}/Res_{method}_Simulation_setting4_nCells_{nCells}_nGroups_{nGroups}__iteration{j}_seeds{n}.rds"))
#       
#       n=n+1
#       j=j+1
#     }
#   }
#   
# }
