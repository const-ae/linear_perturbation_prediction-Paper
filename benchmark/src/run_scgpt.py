import argparse
from pathlib import Path
import tempfile
import shutil
import pickle
import json 
import time
import warnings
import copy
import os


# Import AnnData object with scanpy
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import torch
from torch import nn
from torch.nn import functional as F
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)
from torch_geometric.loader import DataLoader
import scgpt as scg
from scgpt.model import TransformerGenerator
from scgpt.loss import (
    masked_mse_loss,
    criterion_neg_log_bernoulli,
    masked_relative_error,
)
from scgpt.tokenizer import tokenize_batch, pad_batch, tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed, map_raw_id_to_vocab_id

import session_info

print("Hello from python")


parser = argparse.ArgumentParser(description='Run scGPT')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--patience', dest = 'patience', action = 'store', help = "How many epochs with increasing error are allowed", default = 1, type = int)
parser.add_argument('--epochs', dest = 'epochs', action = 'store', help = "How many epochs are run", default = 15, type = int)
parser.add_argument('--pool_size', dest = 'pool_size', action = 'store', help = "Number of cells used for the prediction", default = 100, type = int)
parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")
args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman",
#     "--test_train_config_id", "8443ed21d2ac4-f8716281f960b", "--working_dir",
#     "/scratch/ahlmanne/perturbation_prediction_benchmark", "--result_id", "0"])
print(args)

out_dir = args.working_dir + "/results/" + args.result_id

set_seed(args.seed)
# --------------------------------------------------------

if args.dataset_name == "norman_from_scfoundation":
  import sys
   # scfoundation uses a forked version of GEARS which
  sys.path.insert(0, "/g/huber/users/ahlmanne/projects/perturbation_prediction-benchmark/tmp/scfoundation/scfoundation_gears/")
  sys.path.append("/g/huber/users/ahlmanne/projects/perturbation_prediction-benchmark/tmp/scfoundation/model/")
  import gears.version
  assert gears.version.__version__ == '0.0.2'
else:
  import gears.version
  # scGPT is also built around the GEARS version 0.0.2
  assert gears.version.__version__ == '0.0.2'

print("GEARS version: " + str(gears.version.__version__))
print("GEARS location: " + str(gears.version.__file__))
from gears import PertData, GEARS
from gears.inference import compute_metrics, deeper_analysis, non_dropout_analysis
from gears.utils import create_cell_graph_dataset_for_prediction

# The configuration is based on the https://scgpt.readthedocs.io/en/v0.2.1/tutorial_perturbation.html
# Note that the values slightly differ from the values in https://github.com/bowang-lab/scGPT/blob/7301b51a72f5db321fccebb51bc4dd1380d99023/tutorials/Tutorial_Perturbation.ipynb

# settings for data prcocessing
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]
pad_value = 0  # for padding values
pert_pad_id = 2

n_hvg = 0  # number of highly variable genes
include_zero_gene = "all"  # include zero expr genes in training input, "all", "batch-wise", "row-wise", or False
max_seq_len = 1536

# settings for training
MLM = True  # whether to use masked language modeling, currently it is always on.
CLS = False  # celltype classification objective
CCE = False  # Contrastive cell embedding objective
MVC = False  # Masked value prediction for cell embedding
ECS = False  # Elastic cell similarity objective
cell_emb_style = "cls"
mvc_decoder_style = "inner product, detach"
amp = True
load_model = "/home/ahlmanne/huber/data/scgpt_models/scGPT_human"
load_param_prefixs = [
    "encoder",
    "value_encoder",
    "transformer_encoder",
]

# settings for optimizer
lr = 1e-4  # or 1e-4
batch_size = 64
eval_batch_size = 64
epochs = args.epochs
schedule_interval = 1
early_stop = 5

# settings for the model
embsize = 512  # embedding dimension
d_hid = 512  # dimension of the feedforward network model in nn.TransformerEncoder
nlayers = 12  # number of nn.TransformerEncoderLayer in nn.TransformerEncoder
nhead = 8  # number of heads in nn.MultiheadAttention
n_layers_cls = 3
dropout = 0.2  # dropout probability
use_fast_transformer = True  # whether to use fast transformer

# logging
log_interval = 100
device = "cuda"

# dataset and evaluation choices
pert_data_folder = Path("data/gears_pert_data/")
pert_data = PertData(pert_data_folder)
if args.dataset_name in ['norman', 'adamson', 'dixit']:
  pert_data.load(args.dataset_name)
else:
  pert_data.load(data_path = "data/gears_pert_data/" + args.dataset_name)


conds = pert_data.adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
gene_names = pert_data.adata.var["gene_name"].values.tolist() + ["ctrl"]
good_conds = np.array(conds)[[len(c) == 1 or (c[0] in gene_names and c[1] in gene_names) for c in [co.split("+") for co in conds]]]
pert_data.adata = pert_data.adata[[c in good_conds for c in pert_data.adata.obs["condition"]],:]

with open(args.working_dir + "/results/" + args.test_train_config_id) as json_file:
  set2conditions = json.load(json_file)

# Filter out problematic conditions
set2conditions["train"] = [c for c in set2conditions["train"] if c in good_conds]
set2conditions["test"] = [c for c in set2conditions["test"] if c in good_conds]
set2conditions["val"] = [c for c in set2conditions["val"] if c in good_conds]

print(set2conditions)
pert_data.set2conditions = set2conditions
pert_data.split = "custom"
pert_data.subgroup = None
pert_data.seed = 1
pert_data.train_gene_set_size = 0.75
pert_data.get_dataloader(batch_size = batch_size, test_batch_size = eval_batch_size)
logger = scg.logger 


# Load model
model_dir = Path(load_model)
model_config_file = model_dir / "args.json"
model_file = model_dir / "best_model.pt"
vocab_file = model_dir / "vocab.json"
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

pert_data.adata.var["id_in_vocab"] = [
    1 if gene in vocab else -1 for gene in pert_data.adata.var["gene_name"]
]
gene_ids_in_vocab = np.array(pert_data.adata.var["id_in_vocab"])
logger.info(
    f"match {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes "
    f"in vocabulary of size {len(vocab)}."
)
genes = pert_data.adata.var["gene_name"].tolist()

# model
with open(model_config_file, "r") as f:
    model_configs = json.load(f)
logger.info(
    f"Resume model from {model_file}, the model args will override the "
    f"config {model_config_file}."
)
embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]
n_layers_cls = model_configs["n_layers_cls"]

vocab.set_default_index(vocab["<pad>"])
gene_ids = np.array(
    [vocab[gene] if gene in vocab else vocab["<pad>"] for gene in genes], dtype=int
)
n_genes = len(genes)

ntokens = len(vocab)  # size of vocabulary
model = TransformerGenerator(
    ntokens,
    embsize,
    nhead,
    d_hid,
    nlayers,
    nlayers_cls=n_layers_cls,
    n_cls=1,
    vocab=vocab,
    dropout=dropout,
    pad_token=pad_token,
    pad_value=pad_value,
    pert_pad_id=pert_pad_id,
    do_mvc=MVC,
    cell_emb_style=cell_emb_style,
    mvc_decoder_style=mvc_decoder_style,
    use_fast_transformer=use_fast_transformer,
)
model_dict = model.state_dict()
pretrained_dict = torch.load(model_file)
pretrained_dict = {
    k: v
    for k, v in pretrained_dict.items()
    if any([k.startswith(prefix) for prefix in load_param_prefixs])
}
model_dict.update(pretrained_dict)
model.load_state_dict(model_dict)
model.to(device)

# Make scGPT compatible with GEARs Dataloader from version > 0.0.1
# The problem is that the new data loader provides a pert_idx which
# refers to the position in the gene2go vector.
# Code adapted from https://github.com/snap-stanford/GEARS/blob/e0c27b69f8a1a611d56c3f8c5f5a168cb2cde5f6/gears/pertdata.py#L128
if args.dataset_name == "norman_from_scfoundation":
  with open(os.path.join(pert_data.data_path, 'gene2go.pkl'), 'rb') as f:
          lookup_gene2go = pickle.load(f)
  pert_names = np.unique(list(lookup_gene2go.keys()))
else:
  path_ = os.path.join(pert_data.data_path, 'essential_all_data_pert_genes.pkl')
  with open(path_, 'rb') as f:
      essential_genes = pickle.load(f)
  with open(os.path.join(pert_data.data_path, 'gene2go_all.pkl'), 'rb') as f:
          lookup_gene2go = pickle.load(f)
  gene2go = {i: lookup_gene2go[i] for i in essential_genes if i in lookup_gene2go}
  pert_names = np.unique(list(gene2go.keys()))



criterion = masked_mse_loss
criterion_cls = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=lr)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, schedule_interval, gamma=0.9)
scaler = torch.cuda.amp.GradScaler(enabled=amp)


def train(model: nn.Module, train_loader: torch.utils.data.DataLoader) -> None:
    """
    Train the model for one epoch.
    """
    model.train()
    total_loss, total_mse = 0.0, 0.0
    start_time = time.time()

    num_batches = len(train_loader)
    for batch, batch_data in enumerate(train_loader):
        batch_size = len(batch_data.y)
        batch_data.to(device)
        x: torch.Tensor = batch_data.x  # (batch_size * n_genes, 2)
        # The scFoundation attaches the total counts to the array (pertdata.py#363)
        if args.dataset_name == "norman_from_scfoundation":
          ori_gene_values = x[:, 0].view(batch_size, n_genes+1)[:,0:n_genes]
        else:
          ori_gene_values = x[:, 0].view(batch_size, n_genes)
        # pert_flags = x[:, 1].long().view(batch_size, n_genes)
        pert_flags = torch.zeros((batch_size, n_genes), dtype = torch.long, device = device)
        for idx in range(batch_size):
            for pi in batch_data.pert_idx[idx]:
                if pi != -1:
                    pert_name = pert_names[pi]
                    gene_idx = np.where(pert_data.gene_names == pert_name)[0][0]
                    pert_flags[idx, gene_idx] = 1
        
        target_gene_values = batch_data.y  # (batch_size, n_genes)

        if include_zero_gene in ["all", "batch-wise"]:
            if include_zero_gene == "all":
                input_gene_ids = torch.arange(n_genes, device=device, dtype=torch.long)
            else:
                input_gene_ids = (
                    ori_gene_values.nonzero()[:, 1].flatten().unique().sort()[0]
                )
            # sample input_gene_id
            if len(input_gene_ids) > max_seq_len:
                input_gene_ids = torch.randperm(len(input_gene_ids), device=device)[
                    :max_seq_len
                ]
            input_values = ori_gene_values[:, input_gene_ids]
            input_pert_flags = pert_flags[:, input_gene_ids]
            target_values = target_gene_values[:, input_gene_ids]

            mapped_input_gene_ids = map_raw_id_to_vocab_id(input_gene_ids, gene_ids)
            mapped_input_gene_ids = mapped_input_gene_ids.repeat(batch_size, 1)

            # src_key_padding_mask = mapped_input_gene_ids.eq(vocab[pad_token])
            src_key_padding_mask = torch.zeros_like(
                input_values, dtype=torch.bool, device=device
            )

        with torch.cuda.amp.autocast(enabled=amp):
            output_dict = model(
                mapped_input_gene_ids,
                input_values,
                input_pert_flags,
                src_key_padding_mask=src_key_padding_mask,
                CLS=CLS,
                CCE=CCE,
                MVC=MVC,
                ECS=ECS,
            )
            output_values = output_dict["mlm_output"]

            masked_positions = torch.ones_like(
                input_values, dtype=torch.bool
            )  # Use all
            loss = loss_mse = criterion(output_values, target_values, masked_positions)

        model.zero_grad()
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("always")
            torch.nn.utils.clip_grad_norm_(
                model.parameters(),
                1.0,
                error_if_nonfinite=False if scaler.is_enabled() else True,
            )
            if len(w) > 0:
                logger.warning(
                    f"Found infinite gradient. This may be caused by the gradient "
                    f"scaler. The current scale is {scaler.get_scale()}. This warning "
                    "can be ignored if no longer occurs after autoscaling of the scaler."
                )
        scaler.step(optimizer)
        scaler.update()

        # torch.cuda.empty_cache()

        total_loss += loss.item()
        total_mse += loss_mse.item()
        if batch % log_interval == 0 and batch > 0:
            lr = scheduler.get_last_lr()[0]
            ms_per_batch = (time.time() - start_time) * 1000 / log_interval
            cur_loss = total_loss / log_interval
            cur_mse = total_mse / log_interval
            # ppl = math.exp(cur_loss)
            logger.info(
                f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
                f"lr {lr:05.4f} | ms/batch {ms_per_batch:5.2f} | "
                f"loss {cur_loss:5.2f} | mse {cur_mse:5.2f} |"
            )
            total_loss = 0
            total_mse = 0
            start_time = time.time()


def evaluate(model: nn.Module, val_loader: torch.utils.data.DataLoader) -> float:
    """
    Evaluate the model on the evaluation data.
    """
    model.eval()
    total_loss = 0.0
    total_error = 0.0

    with torch.no_grad():
        for batch, batch_data in enumerate(val_loader):
            batch_size = len(batch_data.y)
            batch_data.to(device)
            x: torch.Tensor = batch_data.x  # (batch_size * n_genes, 2)
            if args.dataset_name == "norman_from_scfoundation":
              ori_gene_values = x[:, 0].view(batch_size, n_genes+1)[:,0:n_genes]
            else:
              ori_gene_values = x[:, 0].view(batch_size, n_genes)
            # pert_flags = x[:, 1].long().view(batch_size, n_genes)
            pert_flags = torch.zeros((batch_size, n_genes), dtype = torch.long, device = device)
            for idx in range(batch_size):
                for pi in batch_data.pert_idx[idx]:
                    if pi != -1:
                        pert_name = pert_names[pi]
                        gene_idx = np.where(pert_data.gene_names == pert_name)[0][0]
                        pert_flags[idx, gene_idx] = 1
            target_gene_values = batch_data.y  # (batch_size, n_genes)

            if include_zero_gene in ["all", "batch-wise"]:
                if include_zero_gene == "all":
                    input_gene_ids = torch.arange(n_genes, device=device)
                else:  # when batch-wise
                    input_gene_ids = (
                        ori_gene_values.nonzero()[:, 1].flatten().unique().sort()[0]
                    )

                # sample input_gene_id
                if len(input_gene_ids) > max_seq_len:
                    input_gene_ids = torch.randperm(len(input_gene_ids), device=device)[
                        :max_seq_len
                    ]
                input_values = ori_gene_values[:, input_gene_ids]
                input_pert_flags = pert_flags[:, input_gene_ids]
                target_values = target_gene_values[:, input_gene_ids]

                mapped_input_gene_ids = map_raw_id_to_vocab_id(input_gene_ids, gene_ids)
                mapped_input_gene_ids = mapped_input_gene_ids.repeat(batch_size, 1)

                # src_key_padding_mask = mapped_input_gene_ids.eq(vocab[pad_token])
                src_key_padding_mask = torch.zeros_like(
                    input_values, dtype=torch.bool, device=input_values.device
                )
            with torch.cuda.amp.autocast(enabled=amp):
                output_dict = model(
                    mapped_input_gene_ids,
                    input_values,
                    input_pert_flags,
                    src_key_padding_mask=src_key_padding_mask,
                    CLS=CLS,
                    CCE=CCE,
                    MVC=MVC,
                    ECS=ECS,
                    do_sample=True,
                )
                output_values = output_dict["mlm_output"]

                masked_positions = torch.ones_like(
                    input_values, dtype=torch.bool, device=input_values.device
                )
                loss = criterion(output_values, target_values, masked_positions)
            total_loss += loss.item()
            total_error += masked_relative_error(
                output_values, target_values, masked_positions
            ).item()
    return total_loss / len(val_loader), total_error / len(val_loader)

def predict(
    model: TransformerGenerator, pert_list, pool_size = None
):
    adata = pert_data.adata
    ctrl_adata = adata[adata.obs["condition"] == "ctrl"]
    if pool_size is None:
        pool_size = len(ctrl_adata.obs)
    gene_list = pert_data.gene_names.values.tolist()
    for pert in pert_list:
        for i in pert:
            if i not in gene_list:
                raise ValueError(
                    "The gene is not in the perturbation graph. Please select from GEARS.gene_list!"
                )

    model.eval()
    device = next(model.parameters()).device
    with torch.no_grad():
        results_pred = {}
        for pert in pert_list:
            cell_graphs = create_cell_graph_dataset_for_prediction(
                pert, ctrl_adata, gene_list, device, num_samples=pool_size
            )
            loader = DataLoader(cell_graphs, batch_size=eval_batch_size, shuffle=False)
            preds = []
            for batch_data in loader:
                if args.dataset_name == "norman_from_scfoundation":
                  # Again some incompatibility around the pert_flags between 
                  # GEARS from scfoundation and GEARS version 0.2.0 as installed by flashattn-env
                  X = np.array(batch_data.x.view( len(batch_data.pert), n_genes+1)[:,0:n_genes].cpu())
                  pert_feats = np.zeros(X.shape)
                  pert_idx = [np.where(p == np.array(gene_list))[0][0] for p in pert]
                  for p in pert_idx:
                    pert_feats[:,int(np.abs(p))] = np.sign(p)
                  feature_mat = torch.Tensor(np.vstack([X.flatten(), pert_feats.flatten()])).T
                  batch_data.x = feature_mat
                pred_gene_values = model.pred_perturb(
                    batch_data, include_zero_gene, gene_ids=gene_ids, amp=amp
                )
                preds.append(pred_gene_values)
            preds = torch.cat(preds, dim=0)
            results_pred["_".join(pert)] = np.mean(preds.detach().cpu().numpy(), axis=0)

    return results_pred


best_val_loss = float("inf")
best_model = copy.deepcopy(model)
patience = args.patience

for epoch in range(0, epochs):
    epoch_start_time = time.time()
    train_loader = pert_data.dataloader["train_loader"]
    valid_loader = pert_data.dataloader["val_loader"]

    train(
        model,
        train_loader,
    )
    val_loss, val_mre = evaluate(
        model,
        valid_loader,
    )
    elapsed = time.time() - epoch_start_time
    logger.info("-" * 89)
    logger.info(
        f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
        f"valid loss/mse {val_loss:5.4f} |"
    )
    logger.info("-" * 89)

    if val_loss < best_val_loss:
        best_val_loss = val_loss
        best_model = copy.deepcopy(model)
        logger.info(f"Best model with score {best_val_loss:5.4f}")
        patience = 0
    else:
        patience += 1
        if patience >= early_stop:
            logger.info(f"Early stop at epoch {epoch}")
            break

    scheduler.step()

conds = pert_data.adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
split_conds = [x.split("+") for x in conds]
split_conds = [list(filter(lambda y: y != "ctrl", x)) for x in split_conds]

tmp_out_dir = tempfile.mkdtemp()
torch.save(best_model.state_dict(), f"{tmp_out_dir}/best_model.pt")

all_pred_vals = predict(best_model, split_conds, pool_size=args.pool_size)
all_pred_vals = {k: v.tolist() for k, v in all_pred_vals.items()}
with open(f"{tmp_out_dir}/all_predictions.json", 'w', encoding="utf8") as handle:
    json.dump(all_pred_vals, handle, indent = 4)
    
with open(f"{tmp_out_dir}/gene_names.json", 'w', encoding="utf8") as handle:
    json.dump(pert_data.adata.var["gene_name"].values.tolist(), handle, indent = 4)

# Make out_dir
shutil.move(tmp_out_dir, out_dir)



session_info.show()
print("Python done")
