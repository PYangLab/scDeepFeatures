import argparse
import os
import numpy as np
import math                                                                 

import torchvision.transforms as transforms
from torchvision.utils import save_image
 
from torch.utils.data import DataLoader
from torchvision import datasets
from torch.autograd import Variable
 
import torch.nn as nn
import torch.nn.functional as F
import torch
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import random
import captum
from captum.attr import *
import time
from tqdm import tqdm
from util import setup_seed, MyDataset, save_checkpoint,AverageMeter,accuracy,ToTensor, read_h5_data, read_fs_label


class Mlp_baseline_LRP(nn.Module):   
    def __init__(self, feature_num, rna_num_celltype):
        super(Mlp_baseline_LRP, self).__init__()
        self.feature_num = feature_num
        self.rna_num_celltype = rna_num_celltype

        self.layer1 = nn.Sequential(
            nn.Linear(self.feature_num, 1024), # 128
            #nn.BatchNorm1d(128),
            nn.ReLU(inplace=False),
            nn.Linear(1024, 512), # 128 128
            #nn.BatchNorm1d(128),
            nn.ReLU(inplace=False)
        )
        
        self.layer2 = nn.Sequential(
            nn.Linear(512, self.rna_num_celltype) # 128
        )

    # 前向传播
    def forward(self, x):
        x = x.view(x.size(0), -1)
        x = self.layer1(x)
        x = self.layer2(x)
        return x

class Mlp_baseline(nn.Module):   
    def __init__(self, feature_num, rna_num_celltype):
        super(Mlp_baseline, self).__init__()
        self.feature_num = feature_num
        self.rna_num_celltype = rna_num_celltype

        self.layer1 = nn.Sequential(
            nn.Linear(self.feature_num, 1024),
            #nn.BatchNorm1d(128),
            nn.LeakyReLU(inplace=False),
            nn.Linear(1024, 512),
            #nn.BatchNorm1d(128),
            nn.ReLU(inplace=False)
        )
        
        self.layer2 = nn.Sequential(
            nn.Linear(512, self.rna_num_celltype)
        )

    # 前向传播
    def forward(self, x):
        x = x.view(x.size(0), -1)
        x = self.layer1(x)
        x = self.layer2(x)
        return x


def train_model(model, train_dl, valid_dl, lr, n_epochs, FloatTensor, LongTensor, classify_dim=17, verbose=True,  save_path = ""):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.CrossEntropyLoss()
    
    best_top1_acc=  0
    best_each_celltype_top1 = []
    best_each_celltype_num=[]
    for i in range(classify_dim):
        best_each_celltype_top1.append(AverageMeter('Acc@1', ':6.2f'))
        best_each_celltype_num.append(0)
            
    for epoch in tqdm(range(1,n_epochs+1)):
        train_top1 = AverageMeter('Acc@1', ':6.2f')
        for i, batch_sample in enumerate(train_dl):
            model.train()
            train_data = batch_sample['data']
            train_label = batch_sample['label']

            # Configure input
            train_data = Variable(train_data.type(FloatTensor))
            train_label = Variable(train_label.type(LongTensor))

            #training
            optimizer.zero_grad()
            train_output = model(train_data)
            loss = criterion(train_output, train_label)
            loss.backward()
            optimizer.step()

            train_pred1,  = accuracy(train_output, train_label, topk=(1, ))
            train_top1.update(train_pred1[0], 1)
            
        valid_top1 = AverageMeter('Acc@1', ':6.2f')
        each_celltype_top1 = []
        each_celltype_num=[]
        for i in range(classify_dim):
            each_celltype_top1.append(AverageMeter('Acc@1', ':6.2f'))
            each_celltype_num.append(0)
               
        with torch.no_grad():
            for i, batch_sample in enumerate(valid_dl):
                test_data = batch_sample['data']
                test_label = batch_sample['label']
                model.eval()
                test_data = Variable(test_data.type(FloatTensor))
                test_label = Variable(test_label.type(LongTensor))

                test_output = model(test_data)
        
                valid_pred1,  = accuracy(test_output, test_label, topk=(1, ))
                valid_top1.update(valid_pred1[0], 1)
                
                for j in range(classify_dim):
                    if len(test_label[test_label==j])!=0:
                        pred1,  = accuracy(test_output[test_label==j,:], test_label[test_label==j], topk=(1, ))
                        each_celltype_top1[j].update(pred1[0],1)
                        each_celltype_num[j]=each_celltype_num[j] + len(test_label[test_label==j])
                                               
            if valid_top1.avg > best_top1_acc:

                best_top1_acc = valid_top1.avg
                for j in range(classify_dim):
                    best_each_celltype_top1[j] = each_celltype_top1[j].avg
                    best_each_celltype_num[j] = each_celltype_num[j]
 
            if epoch == n_epochs:
                print('Epoch : ',epoch, '\t', 'train prec:', train_top1.avg, 'valid_prec :', best_top1_acc)
                #torch.save(model.state_dict(), save_path) #commented to save disk space
                for j in range(classify_dim):
                    print('cell type : ',j, '\t', '\t',  'prec :', best_each_celltype_top1[j], 'number:',best_each_celltype_num[j])


