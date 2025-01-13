
from pathlib import Path
import shutil
import anndata as ad
import numpy as np
import pandas as pd
import torch
import json
import pickle
import tempfile
import os
import argparse
import session_info


parser = argparse.ArgumentParser(description='Run GeneFormer')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--model_type', dest = 'model_type', action = 'store', help = "The UCE model type ('4layers' or '33layers' ", default = "4layers")
parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")

args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman_from_scfoundation",
#     "--perturbation_type", "overexpress",  "--model_type", "4layers",
#     "--test_train_config_id", "ede703be5dd78-c8355c044d3ba", "--working_dir",
#     "/scratch/ahlmanne/perturbation_prediction_benchmark", "--result_id", "0"])
print(args)
print("GPU: " + str(torch.cuda.get_device_name()))

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



adata.var = adata.var.set_index("gene_name", drop=True)

# Convert adata.X[ij] = log(c / (s / 10k) + 1) to counts
counts = np.expm1(adata.X) / 1e4
# The smallest non-zero element used to be a 1
one_equiv = np.array([np.min(counts[i,:].data) for i in range(counts.shape[0])])
counts = (counts / one_equiv[:,np.newaxis]).tocsr()
new_adata = ad.AnnData(counts, obs = adata.obs, var = adata.var)

sub_adata = new_adata[new_adata.obs['condition'] == 'ctrl',:]
sub_adata = sub_adata[np.random.choice(range(sub_adata.shape[0]), 100),:]

conds = adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
embeddings = {}
for idx, genes in enumerate(conds):
  gene_list = list(set(genes.split("+")) - set(["ctrl"]))
  
  matches = np.isin(new_adata.var_names, gene_list)
  if genes == "ctrl":
    sub_adata_mod = sub_adata.copy()
  elif np.sum(matches) == len(gene_list):
    sub_adata_mod = sub_adata.copy()
    modified_expr = adata[adata.obs['condition'] == genes,matches].X.toarray()
    # This leaking some of the test data to the model. But it is still more elegant
    # than just setting the values to zero or some arbitrary value for overexpression
    sub_adata_mod.X[:,matches] = modified_expr[np.random.choice(range(modified_expr.shape[0]), sub_adata_mod.shape[0]),:]
  else:
    continue
  
  with  tempfile.TemporaryDirectory() as adata_dir:
    sub_adata_mod.write_h5ad(adata_dir + "/adata.h5ad")
    if args.model_type == "4layers":
      cmd = f"""
      cd /home/ahlmanne/prog/UCE && \
      python eval_single_anndata.py \
        --adata_path {adata_dir}/adata.h5ad \
        --dir "{adata_dir}/" \
        --species human \
        --model_loc /home/ahlmanne/data/universal_cell_embedding/4layer_model.torch \
        --batch_size 100 \
        --nlayers 4
      """
    elif args.model_type == "33layers":
      cmd = f"""
      cd /home/ahlmanne/prog/UCE && \
      python eval_single_anndata.py \
        --adata_path {adata_dir}/adata.h5ad \
        --dir "{adata_dir}/" \
        --species human \
        --model_loc /home/ahlmanne/data/universal_cell_embedding/33l_8ep_1024t_1280.torch \
        --batch_size 25 \
        --nlayers 33
      """
    else:
      raise ValueError(f"--model_type must be either '4layers' or '33layers', not '{args.model_type}'.")
    
    print(f"Running script for {genes} ({idx}):\n{cmd}")
    os.system(cmd)
    result = ad.read_h5ad(adata_dir + "/adata_uce_adata.h5ad")
    embeddings[genes] = result.obsm['X_uce'].mean(axis=0)



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
    json.dump(adata.var_names.values.tolist(), handle, indent = 4)

# Move results to out_dir
shutil.move(tmp_out_dir, out_dir)



session_info.show()
print("Python done")
