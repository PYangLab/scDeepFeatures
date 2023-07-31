import argparse
import os
import numpy as np
import math   
import h5py

import scipy

import torchvision.transforms as transforms
from torchvision.utils import save_image
 
from torch.utils.data import DataLoader
from torchvision import datasets
from torch.autograd import Variable
 
import torch.nn as nn
import torch.nn.functional as F
import torch
from tqdm import tqdm
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import random
import captum
from captum.attr import *
import time
from util import setup_seed, MyDataset, save_checkpoint,AverageMeter,accuracy,ToTensor, read_h5_data, read_fs_label, read_h5_feature, read_ct_name
from learn import Mlp_baseline, Mlp_baseline_LRP, train_model








#####################################################
setup_seed(1234)
cuda = True if torch.cuda.is_available() else False
device = torch.device("cuda:0") if cuda == True else torch.device("cpu") 
print("cuda is: " + str(cuda) + "\n" +
      "current device is: " + str(torch.cuda.get_device_name(0)))
FloatTensor = torch.cuda.FloatTensor if cuda else torch.FloatTensor
LongTensor = torch.cuda.LongTensor if cuda else torch.LongTensor










#####################################################
# add arguments
parser = argparse.ArgumentParser()
parser.add_argument('--method', type=str, default='IG', help='load generator model')
parser.add_argument('--seed', type=int, default=1234, help='random seed')
parser.add_argument('--seedInput', type=int, default=1234, help='which seed data we are running')
parser.add_argument('--train_data', type=str, default='../../../Data/Simulation/Tabula_muris_h5/Setting1/Train_Simulation_setting1_nCells_200_nGroups_5__iteration1_seeds1.h5', help='path of training data file')
parser.add_argument('--train_label', type=str, default='../../../Data/Simulation/Tabula_muris_h5/Setting1/Train_Cty_Simulation_setting1_nCells_200_nGroups_5__iteration1_seeds1.csv', help='path of training label file')
parser.add_argument('--test_data', type=str, default='../../../Data/Simulation/Tabula_muris_h5/Setting1/Test_Simulation_setting1_nCells_200_nGroups_5__iteration1_seeds1.h5', help='path of test data file')
parser.add_argument('--test_label', type=str, default='../../../Data/Simulation/Tabula_muris_h5/Setting1/Test_Cty_Simulation_setting1_nCells_200_nGroups_5__iteration1_seeds1.csv', help='path of test label file')
parser.add_argument('--n_epochs', type=int, default=20, help='num of training epochs')
parser.add_argument('--lr', type=float, default=2e-3, help='init learning rate')
parser.add_argument('--batch_size', type=int, default=64, help='batch size')
parser.add_argument('--marker_num', type=int, default=0, help='batch size') # add a tailor to this parameters
parser.add_argument('--only_pos', type=str, default="activated", help='select only positive marker')
parser.add_argument('--model_save_path', type=str, default='../../../Result/Fs_marker/Nn_test/', help='path for saving trained models') 
parser.add_argument('--save_fs_eachcell', type=str, default='../../../Result/Fs_marker/Nn_test/', help='path for saving fs gene name')
args = parser.parse_args()


train_data = read_h5_data(args.train_data)
train_label = read_fs_label(args.train_label)
test_data = read_h5_data(args.test_data)
test_label = read_fs_label(args.test_label)








features_name = read_h5_feature(args.train_data)
print("\n training dataset head of features: ", features_name[0:10])
print("\n training dataset feature dimension check: ", features_name.shape)

celltype_name = read_ct_name(args.train_label)
print("\n training set snapshot: ")
print(train_data)
print("\n training dataset celltype name and freq: ")
print(celltype_name)
print("\n training dataset shape is: ", train_data.shape,)


feature_num = train_data.shape[1]
rna_num_celltype = len(torch.unique(train_label))
print("\n training dataset feature num check 2: " + str(feature_num))
print("\n rna_num_celltype: " + str(rna_num_celltype))

assert(len(features_name) == feature_num)


if args.marker_num == 0:
    args.marker_num = feature_num
print("\n marker num for selection is: " + str(args.marker_num))



train_transformed_dataset = MyDataset(train_data,train_label)
train_dl = DataLoader(train_transformed_dataset, batch_size=args.batch_size, shuffle=True, num_workers=0,drop_last=True)
test_transformed_dataset = MyDataset(test_data,test_label)
valid_dl = DataLoader(test_transformed_dataset, batch_size=args.batch_size, shuffle=True, num_workers=0,drop_last=True)

#args.train_data[-68:]



print("========== now training the neural net ==========")
if args.method == "LRP" :
    model = Mlp_baseline_LRP(feature_num, rna_num_celltype)
else:
    model = Mlp_baseline(feature_num, rna_num_celltype)
print("\n usingi model: ")
print(model)
model.to(device)
train_model(model, train_dl, train_dl, lr=args.lr, n_epochs=args.n_epochs, FloatTensor=FloatTensor, LongTensor=LongTensor, classify_dim=rna_num_celltype, save_path=os.path.join(args.model_save_path, 'model_best.pth.tar'))   # last argument did not apply in main.py, save storage


print("========== now doing feature selection ==========")
if args.method == "DeepLift":
    deconv = DeepLift(model)
if args.method == "GradientShap": 
    deconv = GradientShap(model)
if args.method == "FeatureAblation":
    deconv = FeatureAblation(model)
if args.method == "Occlusion": 
    deconv = Occlusion(model)
if args.method == "Lime": 
    deconv = Lime(model)
if args.method == "LRP": 
    deconv = LRP(model)


print("\n run method: " + args.method)


rna_train_data_fs = train_data
rna_train_label_fs = train_label

rna_train_data_fs_topn_list = []
rna_train_label_fs_topn_list = []
rna_train_fs_index_list=[]
rna_train_fs_list=[]
all_feature = []


for i in range(rna_num_celltype):
    rna_train_index_fs= torch.where(rna_train_label_fs==i)
    rna_train_index_fs = [t.numpy() for t in rna_train_index_fs]
    rna_train_index_fs = np.array(rna_train_index_fs)
    rna_train_data_each_celltype_fs = rna_train_data_fs[rna_train_index_fs,:].reshape(-1,feature_num)
        
    attribution = torch.zeros(1,feature_num).to(device)
    full_attribution = torch.zeros(1,feature_num).to(device)
    print(rna_train_data_each_celltype_fs.size())
    if args.method == "FeatureAblation":
        attribution = attribution + deconv.attribute(rna_train_data_each_celltype_fs, target=i)
    elif args.method == "Occlusion":
        attribution = attribution + deconv.attribute(rna_train_data_each_celltype_fs, target=i, sliding_window_shapes=(3,))
    else:
        for j in tqdm(range(rna_train_data_each_celltype_fs.size(0)-1)):
            if args.method == "GradientShap":
                baselines = torch.zeros(1, rna_train_data_each_celltype_fs.size(1)).to(device)
                attribution = attribution + deconv.attribute(rna_train_data_each_celltype_fs[j:j+1,:], baselines, target=i)
            elif args.method == "Lime":
                tmp_attr = deconv.attribute(rna_train_data_each_celltype_fs[j:j+1,:], target=i, n_samples=500)
                attribution = attribution + tmp_attr
                full_attribution = torch.cat((full_attribution, tmp_attr), 0)
            else:
                attribution = attribution + deconv.attribute(rna_train_data_each_celltype_fs[j:j+1,:], target=i)
    print("--------------------------------the attribution shape is-------------------------------")
    print(attribution.shape)   

    attribution_mean = torch.sum((attribution),dim=0) # here is the score, mean and sum does not change the rank  
    attribution_mean_copy = attribution_mean # here is the copy for calculating name
    v,ind = attribution_mean_copy.sort(dim=0,descending=True)
    df = pd.DataFrame([  ind.cpu().numpy().astype("int64"), v.detach().cpu().numpy(), np.arange(features_name.shape[0]).astype("int64")]).T
    df = df.rename({0: "index", 1: "scores", 2: "rank"}, axis=1)
    sorted_df = df.sort_values(by=['index'], ascending=True)
    sorted_df = sorted_df.set_axis(features_name.tolist(), axis=0)
    # calculate logFC
    anchor_exp = rna_train_data_fs[rna_train_label_fs == i,]
    other_exp = rna_train_data_fs[rna_train_label_fs != i,]
    logFC = torch.mean(anchor_exp,0) - torch.mean(other_exp,0)
    sorted_df["logFC"] = logFC.cpu().numpy()
    # rank the df, by importance scores
    sorted_df_2 = sorted_df.sort_values(by=['scores'], ascending=False)
    if args.only_pos == "activated":
        sorted_df_2 = sorted_df_2[sorted_df_2["logFC"]>0]
    sorted_df_2 = sorted_df_2.iloc[0:args.marker_num,:]
    #print(sorted_df_2)
    
    print("\n feature selection of celltype:", i, "finish")
    if not os.path.exists(args.save_fs_eachcell):
        os.makedirs(args.save_fs_eachcell)
    sorted_df_2.to_csv(args.save_fs_eachcell+str(args.method)+"."+str(i) + "." + list(celltype_name.index.values)[i] + ".csv")
print("\n check result: \n")
print(sorted_df_2) #print to terminal for results checking
                                                                                           



























