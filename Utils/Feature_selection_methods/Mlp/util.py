import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from torchvision import datasets
import torchvision.transforms as transforms
import random
import os
import umap
import anndata
import scanpy as sc
from torch.autograd import Variable
import h5py
import scipy

cuda = True if torch.cuda.is_available() else False
FloatTensor = torch.cuda.FloatTensor if cuda else torch.FloatTensor
LongTensor = torch.cuda.LongTensor if cuda else torch.LongTensor

def setup_seed(seed):
     torch.manual_seed(seed)
     torch.cuda.manual_seed_all(seed)
     np.random.seed(seed)
     random.seed(seed)
     torch.backends.cudnn.deterministic = True
     torch.backends.cudnn.benchmark = False
     os.environ['PYTHONHASHSEED']=str(seed)


class MyDataset(Dataset):
    def __init__(self, data, label):
        self.data = data
        self.label = label

    def __getitem__(self, index):#返回的是tensor
        img, target = self.data[index,:], self.label[index]
        sample = {'data': img, 'label': target}
        return sample

    def __len__(self):
        return len(self.data)
        

        
class ToTensor(object):
    def __call__(self, sample):
        data,label = sample['data'], sample['label']
        data = data.transpose((1, 0))
        return {'data': torch.from_numpy(data),
                'label': torch.from_numpy(label),
               }
               
class AverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self, name, fmt=':f'):
        self.name = name
        self.fmt = fmt
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

    def __str__(self):
        fmtstr = '{name} {val' + self.fmt + '} ({avg' + self.fmt + '})'
        return fmtstr.format(**self.__dict__)


def accuracy(output, target, topk=(1,)):
    """Computes the accuracy over the k top predictions for the specified values of k"""
    with torch.no_grad():
        maxk = max(topk)
        batch_size = target.size(0)

        _, pred = output.topk(maxk, 1, True, True)
        pred = pred.t()
        correct = pred.eq(target.view(1, -1).expand_as(pred))

        res = []
        for k in topk:
            correct_k = correct[:k].view(-1).float().sum(0, keepdim=True)
            res.append(correct_k.mul_(100.0 / batch_size))
            
        return res
    
def save_checkpoint(state, save):
    if not os.path.exists(save):
        os.makedirs(save)
    filename = os.path.join(save,'model_best.pth.tar')
    torch.save(state, filename)
    

def read_h5_data(data_path):
    data = h5py.File(data_path,"r")
    h5_data = data['matrix/data']
    sparse_data = scipy.sparse.csr_matrix(np.array(h5_data).transpose())
    data_fs = torch.from_numpy(np.array(sparse_data.todense()))
    data_fs = Variable(data_fs.type(FloatTensor))
    return data_fs

def read_h5_feature(data_path):
    data = h5py.File(data_path,"r")
    features = data['matrix/features'][:]
    features = features.astype(str)
    return features

    
def read_fs_label(label_path):
    label_fs = pd.read_csv(label_path,header=None,index_col=False)  #
    label_fs = pd.Categorical(label_fs.iloc[1:(label_fs.shape[0]),1]).codes
    label_fs = np.array(label_fs[:]).astype('int32')
    label_fs = torch.from_numpy(label_fs)#
    return label_fs

def read_ct_name(label_path):
    label_fs = pd.read_csv(label_path,header=None,index_col=False)  #
    ct_name = pd.DataFrame( pd.Categorical(label_fs.iloc[1:(label_fs.shape[0]),1]).value_counts())
    return ct_name


