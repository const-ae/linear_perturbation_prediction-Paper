
from pathlib import Path
import shutil
import anndata as ad
import numpy as np
import pandas as pd
import json
import pickle
import tempfile
import os
from scipy import sparse
import scanpy as sc
import argparse
import session_info

parser = argparse.ArgumentParser(description='Run scBERT')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--finetuning_epochs', dest = 'finetuning_epochs', action = 'store', help = "The number of epochs used for finetuning. scBERT has a default of 100", default = 100, type = int)
parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")

args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman_from_scfoundation",
#     "--finetuning_epochs", "1",
#     "--test_train_config_id", "ede703be5dd78-c8355c044d3ba", "--working_dir",
#     "/scratch/ahlmanne/perturbation_prediction_benchmark", "--result_id", "0"])
print(args)

out_dir = args.working_dir + "/results/" + args.result_id

np.random.seed(args.seed)
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
  assert gears.version.__version__ == '0.1.2'

from gears import PertData, GEARS
from gears.utils import filter_pert_in_go


pert_data_folder = Path("data/gears_pert_data/")
pert_data = PertData(pert_data_folder)
if args.dataset_name in ['norman', 'adamson', 'dixit']:
  pert_data.load(args.dataset_name)
else:
  pert_data.load(data_path = "data/gears_pert_data/" + args.dataset_name)

adata = pert_data.adata
with open(args.working_dir + "/results/" + args.test_train_config_id) as json_file:
  set2conditions = json.load(json_file)



adata.var = adata.var.set_index("gene_name", drop=False)

# Convert adata.X[ij] = log(c / (s / 10k) + 1) to my_counts
my_counts = np.expm1(adata.X) / 1e4
# The smallest non-zero element used to be a 1
one_equiv = np.array([np.min(my_counts[i,:].data) for i in range(my_counts.shape[0])])
my_counts = (my_counts / one_equiv[:,np.newaxis]).tocsr()

# The object must exactly match the shape and gene order of panglao
# NOTE: if a gene is not present in adata, it will be left blank
panglao = ad.read_h5ad('/home/ahlmanne/data/scbert/panglao_human.h5ad')
counts = sparse.lil_matrix((my_counts.shape[0],panglao.X.shape[1]),dtype=np.float32)
ref = panglao.var_names.values
obj = adata.var_names.values

# This is a more efficient way to copy over the content
matches = np.array([np.array([idx, np.where(obj == e)[0][0]]) if e in obj else np.array([idx, np.nan]) for idx, e in enumerate(ref)])
matches = matches[~np.isnan(matches[:,1])]
matches = matches.astype(int)
counts[:,matches[:,0]] = my_counts[:,matches[:,1]]

# TThe original inefficient method
# for i in range(len(ref)):
#     if ref[i] in obj:
#         loc = obj.index(ref[i])
#         counts[:,i] = my_counts[:,loc]

counts = counts.tocsr()
new_adata = ad.AnnData(counts, obs = adata.obs.copy(), var = panglao.var.copy())
sc.pp.filter_cells(new_adata, min_genes=200)
sc.pp.normalize_total(new_adata, target_sum=1e4)
sc.pp.log1p(new_adata, base=2) # NOTE: using the log2 differs from the original adata.X values

# Do finetuning 
finetuning_data = new_adata.copy()

# Subsample to speed-up training
sc.pp.subsample(finetuning_data, n_obs = min(finetuning_data.shape[0], 5000), random_state = 0)
finetuning_data.obs['celltype'] = finetuning_data.obs['condition']
# Filter out rare elements
for cond, cnt in zip(*np.unique(finetuning_data.obs['celltype'], return_counts = True)):
  if cnt <= 5:
    finetuning_data = finetuning_data[finetuning_data.obs['celltype'] != cond,:]
  
finetuning_data.obs["celltype"].cat.remove_unused_categories()
finetuning_dir = tempfile.mkdtemp()
finetuning_data.write_h5ad(finetuning_dir + "/finetuning_adata.h5ad")


def _find_free_port():
  # coppied from https://github.com/facebookresearch/detectron2/blob/main/detectron2/engine/launch.py#L15
  import socket
  
  sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
  # Binding to port 0 will cause the OS to find an available port for us
  sock.bind(("", 0))
  port = sock.getsockname()[1]
  sock.close()
  # NOTE: there is still a chance the port could be taken by other processes.
  return port

print("This calls a modified version of finetune.py (ie., setting 'strict=False' while loading)")
# I copied `./src/finetune_modified.py` to `/home/ahlmanne/prog/scBERT`
# 5 epochs are a lot less than the recommended 100 iterations
cmd = f"""
cd /home/ahlmanne/prog/scBERT && \
python \
  -m torch.distributed.launch --nproc_per_node=1 --master_port={_find_free_port()} \
  finetune_modified.py \
  --bin_num 5 \
  --epoch {args.finetuning_epochs} \
  --ckpt_dir {finetuning_dir}/ckpt_folder/ \
  --data_path {finetuning_dir}/finetuning_adata.h5ad \
  --model_path /home/ahlmanne/projects/perturbation_prediction-benchmark/data/panglao_pretrain.pth 
"""
print(f"Running:\n{cmd}")
os.system(cmd)


# These two values are needed to make the script happy
new_adata.obs['celltype'] =  "cellline"
new_adata.obs['dataset'] =  args.dataset_name

sub_adata = new_adata[new_adata.obs['condition'] == 'ctrl',:]
sub_adata = sub_adata[np.random.choice(range(sub_adata.shape[0]), size = 100), :]

conds = adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
embeddings = {}
for idx, genes in enumerate(conds):
  gene_list = list(set(genes.split("+")) - set(["ctrl"]))
  
  matches = np.isin(new_adata.var_names, gene_list)
  if genes == "ctrl":
    sub_adata_mod = sub_adata.copy()
  elif np.sum(matches) == len(gene_list):
    sub_adata_mod = sub_adata.copy()
    modified_expr = new_adata[new_adata.obs['condition'] == genes,matches].X.toarray()
    # This leaking some of the test data to the model. But it is still more elegant
    # than just setting the values to zero or some arbitrary value for overexpression
    sub_adata_mod.X[:,matches] = modified_expr[np.random.choice(range(modified_expr.shape[0]), sub_adata_mod.shape[0]),:]
  else:
    continue
  with  tempfile.TemporaryDirectory() as adata_dir:
    sub_adata_mod.write_h5ad(adata_dir + "/adata.h5ad")
    print("This calls a modified version of the attn_sum_save.py script")
    # I copied `./src/attn_sum_save_modified.py` to `/home/ahlmanne/prog/scBERT`
    cmd = f"""
    cd /home/ahlmanne/prog/scBERT && \
    python attn_sum_save_modified.py \
      --bin_num 5 \
      --gene_num {sub_adata.shape[1]} \
      --data_path {adata_dir}/adata.h5ad \
      --model_path {finetuning_dir}/ckpt_folder/finetune_best.pth \
      --save_dir "{adata_dir}/" 
    """
    print(f"Running:\n{cmd}")
    os.system(cmd)
    result = np.load(f"{adata_dir}/full_attn_sum.npy")
    embeddings[genes] = result 



# Using the embedding for each perturbation, predict expression!
def mean_gene_expr_for_condition(condition):
    sub_data = adata[adata.obs['condition'] == condition,:]
    return np.array(sub_data.X.mean(axis=0)).flatten()


training_df = (pd.DataFrame({"training": set2conditions.keys(),
              "condition": set2conditions.values()})
  .explode("condition"))

base_df = (pd.DataFrame({
 "condition": conds})
  .merge(training_df, how = "left", on = "condition", validate = "one_to_one")
  .assign(obs_mean = lambda x: [mean_gene_expr_for_condition(y) for y in x.condition])
  .assign(embedding = lambda x: [embeddings[y] if y in embeddings.keys() else np.nan for y in x.condition])
)

base_df = base_df.dropna()

train = base_df[base_df['training'] == "train"]
from sklearn.linear_model import Ridge
expr_predictor = Ridge(fit_intercept = True)
Y = np.vstack(train['obs_mean'])
emb_mat_train = np.vstack(train['embedding'])
expr_predictor.fit(X = emb_mat_train, y = Y)


emb_mat = np.vstack(base_df['embedding'])
res = expr_predictor.predict(X = emb_mat)
all_pred_vals = {k : res[i,:].tolist() for i,k in enumerate(base_df['condition'])}

tmp_out_dir = tempfile.mkdtemp()
with open(f"{tmp_out_dir}/all_predictions.json", 'w', encoding="utf8") as handle:
    json.dump(all_pred_vals, handle, indent = 4)
with open(f"{tmp_out_dir}/gene_names.json", 'w', encoding="utf8") as handle:
    json.dump(adata.var["gene_name"].values.tolist(), handle, indent = 4)

# Move results to out_dir
shutil.move(tmp_out_dir, out_dir)



session_info.show()
print("Python done")
