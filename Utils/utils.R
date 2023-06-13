# load packages
suppressPackageStartupMessages({
  library(caret)
  library(rhdf5)
  library(HDF5Array)
})


# split train and test on 50% training set, 50% test set
# input "rna" can be sce, seurat, or matrix, "cty" is a vector of cell type
train_test_split <- function(rna, cty, seedInput = 1234){ 
  set.seed(seedInput)
  f <- createFolds(cty, 2)
  
  f_1_train = f[[1]]
  f_1_test = f[[2]]
  
  train.rna = rna[, f_1_train]
  train.cty = cty[f_1_train]
  
  test.rna = rna[, f_1_test]
  test.cty = cty[f_1_test]
  
  return(list(
    train.rna = train.rna,
    train.cty = train.cty,
    test.rna = test.rna,
    test.cty = test.cty
  ))
}



#' 
#' @param exprs_list a list of expression matrices
#' @param h5file_list a vector indictates the the h5 file names to output, corresponding to the expression matrix

write_h5_DL <- function(exprs_list, h5file_list) {
  
  if (length(unique(lapply(exprs_list, rownames))) != 1) {
    stop("rownames of exprs_list are not identical.")
  }
  
  for (i in seq_along(exprs_list)) {
    if (file.exists(h5file_list[i])) {
      warning("h5file exists! will rewrite it.")
      system(paste("rm", h5file_list[i]))
    }
    
    h5createFile(h5file_list[i])
    h5createGroup(h5file_list[i], "matrix")
    writeHDF5Array(t((exprs_list[[i]])), h5file_list[i], name = "matrix/data")
    h5write(rownames(exprs_list[[i]]), h5file_list[i], name = "matrix/features")
    h5write(colnames(exprs_list[[i]]), h5file_list[i], name = "matrix/barcodes")
    print(h5ls(h5file_list[i]))
    
  }
  
  
}


#' 
#' @param cellType_list a list of cell types
#' @param csv_list a vector indictates the the csv file names to output, corresponding to the cell type list

write_csv_DL <- function(cellType_list, csv_list) {
  
  for (i in seq_along(cellType_list)) {
    
    if (file.exists(csv_list[i])) {
      warning("csv_list exists! will rewrite it.")
      system(paste("rm", csv_list[i]))
    }
    
    names(cellType_list[[i]]) <- NULL
    write.csv(cellType_list[[i]], file = csv_list[i])
    
  }
}
