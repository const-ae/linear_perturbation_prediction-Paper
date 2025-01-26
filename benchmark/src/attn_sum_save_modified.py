# -*- coding: utf-8 -*-
import os
import gc
import argparse
import json
import logging
import random
import math
import random
from functools import reduce
import numpy as np
import pandas as pd
from scipy import sparse
import scipy.io as sio
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, precision_recall_fscore_support, classification_report
import torch
from torch import nn
from torch.optim import Adam, SGD, AdamW
from torch.nn import functional as F
from torch.optim.lr_scheduler import StepLR, CosineAnnealingWarmRestarts, CyclicLR

from performer_pytorch import PerformerLM
import scanpy as sc
from utils import *

parser = argparse.ArgumentParser()
parser.add_argument("--bin_num", type=int, default=7, help='Number of bins.')
parser.add_argument("--gene_num", type=int, default=16906, help='Number of genes.')
parser.add_argument("--data_path", type=str, default='./data/data.h5ad', help='Path of data for generating the embeddings.')
parser.add_argument("--model_path", type=str, default='./model.pth', help='Path of model training on the data.')
parser.add_argument("--save_dir", type=str, default='./attention/', help='Directory of embeddings to save.')

args = parser.parse_args()

SEQ_LEN = args.gene_num + 1
CLASS = args.bin_num + 2

data_dir = args.data_path
model_dir = args.model_path
save_dir = args.save_dir

device = torch.device("cuda")
print('            =======  Config over  ======= \n')

data = sc.read_h5ad(data_dir)
methods = np.unique(data.obs['dataset'])
index_methods = data.obs['dataset']
index_labels = data.obs['celltype']
cellinds = list(set(index_labels.tolist()))
label_dict, label = np.unique(np.array(data.obs['celltype']), return_inverse=True)
data_counts = data.X

all_mtx = []
for cellind in cellinds:
    print(cellind)
    data_alpha = data_counts[index_labels == cellind]

    model = PerformerLM(
        num_tokens = CLASS,
        dim = 200,
        depth = 6,
        max_seq_len = SEQ_LEN,
        heads = 10,
        local_attn_heads = 0,
        g2v_position_emb = False
    )
    print(f'            =======  Model defined  ======= \n')

    ckpt = torch.load(model_dir)
    model.load_state_dict(ckpt['model_state_dict'], strict = False)
    model = model.to(device)
    print('            =======  Predict start  ======= \n')
    batch_size = data_alpha.shape[0]
    model.eval()
    with torch.no_grad():

        final_mtx = torch.zeros(batch_size, 200).to(device)
        for index in range(3):
            full_seq = data_alpha[index].toarray()[0]
            full_seq[full_seq > (CLASS - 2)] = CLASS - 2
            full_seq = torch.from_numpy(full_seq).long()
            full_seq = torch.cat((full_seq, torch.tensor([0]))).to(device)
            full_seq = full_seq.unsqueeze(0)
            x, _ = model(full_seq, output_attentions=True, return_encodings=True)
            emb_cls = x[:, -1, :]
            emb_avg = x.mean(dim=1)
            final_mtx[index,:] = emb_avg

    final_mtx =  final_mtx.mean(axis=0).detach().cpu().numpy()
    np.save(os.path.join(save_dir, 'full_attn_sum.npy'), final_mtx)
    print(f'            =======  Predict end  ======= \n')
