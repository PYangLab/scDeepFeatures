# scDeepFeatures: Deep learning-based feature selection for single-cell RNA sequencing data analysis

In this work, we explore the utility of various deep learning-based feature selection methods for scRNA-seq data analysis. We sample from Tabula Muris and Tabula Sapiens atlases to create scRNA-seq datasets with a range of data properties and evaluate the performance of traditional and deep learning-based feature selection methods for cell type classification, feature selection reproducibility and diversity, and computational time. Our study provides a reference for future development and application of deep learning-based feature selection methods for single-cell omics data analyses.


<img width=100% src="https://github.com/HaoHuang-USYD/DL_feature_selection/blob/main/img/main.png"/>

## Installation:
DL feature selection for scRNA-seq is developed using PyTorch 1.11.0 and Captum 0.5.0 and requires >=1 GPU to run. We recommend using conda enviroment to install and run the framework. We assume conda is installed. You can use the provided environment or install the environment by yourself accoring to your hardware settings. Note the following installation code snippets were tested on a Windows system (Win10 professional) with a NVIDIA GeForce RTX 3090 GPU. The installation process needs about 15 minutes.
 
### Installation using provided environment
Step 1: Create and activate the conda environment for matilda using our provided file
```
conda env create -f environment_DL_feature_selection.yaml
conda activate fsdl
```

Step 2:
Obtain DL_feature_selection by clonning the github repository:
```
git clone https://github.com/PYangLab/scDeepFeatures.git
```

## Preparing intput
The main function takes expression data (i.e. RNA) in `.h5` format and cell type labels in `.csv` format, with log-normalised count data. 
An example for creating .h5 file with Utils/utils.R from a singlecellexperiment object in the R environment is as below:
```
source("utils.R")

# singlecellexperiment object (sce) with cellTypes as cell type annotation
train.exprsmat = as.matrix(logcounts(sce))
write_h5_DL(exprs_list = list(rna = train.exprsmat),
                 h5file_list = c("data.h5")))
write_csv_DL(cellType_list =  list(rna = sce$cellTypes),
               csv_list = c("label.csv")))
```

### Example dataset
Without downloading all the datasets we sampled (see link provided below), one can use the example dataset saved in `Data/Example_dataset`
Training and testing on demo dataset will cost no more than 1 minute with GeForce RTX 3090 GPU.

## Running deep learning-based feature selection with the example dataset
### Training the MLP model and performing feature selection
```
cd Utils
cd Feature_selection_methods
cd Mlp

# training and perform feature selection for scRNA-seq 
# with the interest of feature selection, test data is the same as training data
python main.py --method [method] --train_data [path to training data]  --train_label [path to training label] --test_data [path to test data] --test_label [path to test label] --save_fs_eachcell [path to save feature selection results]
# Example run
python main.py --method DeepLift  --train_data ../../../Data/Example_dataset/data.h5 --train_label ../../../Data/Example_dataset/label.csv --test_data ../../../Data/Example_dataset/data.h5 --test_label ../../../Data/Example_dataset/label.csv --save_fs_eachcell ../../../Data/Example_dataset/
```
### Argument
+ `--method`: Method name of Deep learning-based feature selection (i.e. DeepLift, GradientShap, LRP, FeatureAblation, Occlusion, Lime)
+ `--train_data`: path to training scRNA-seq data (h5 file)
+ `--train_label`: path to training scRNA-seq label (csv file)
+ `--test_data`: same as train data
+ `--test_label`: same as train label
+ `--only_pos`: whether to select only up-regulated genes, default is "activated" , change to any other string to record the overall gene ranks
+ `--save_fs_eachcell`: path to save the feature selection results

Training and model config
+ `--batch_size`: Batch size (set as 64 by default)
+ `--epochs`: Number of epochs (set as 20 by default)
+ `--lr`: Learning rate.


Other config
+ `--seed`: The random seed for training.
+ `--augmentation`: Whether to augment simulated data.


## Feature selection 
Scripts for feature selection have been deposited in `Main/TM_RunFeatureSelection` and `Main/TS_RunFeatureSelection`, and `Utils/Feature_selection_methods/Mlp/`.

## Scripts for analysing feature selection results
Scripts for analysing the feature selection results have been deposited in `Main/TM_RunEvaluation`, `Main/TS_RunEvaluation` and `Main/Combine_Evaluation`.

## Citation
If you use this content, please cite:
```
Huang, H., Liu, C., Wagle, M.M. et al. Evaluation of deep learning-based feature selection for single-cell RNA sequencing data analysis. Genome Biol 24, 259 (2023). https://doi.org/10.1186/s13059-023-03100-x
```
View our publication here: 
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03100-x
