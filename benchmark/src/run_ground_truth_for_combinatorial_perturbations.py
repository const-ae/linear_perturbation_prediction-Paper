from gears import PertData, GEARS
from pathlib import Path
import numpy as np
import pandas as pd
import json
import tempfile
import shutil
import argparse
import session_info

parser = argparse.ArgumentParser(description='Collect ground truth for combinatorial data')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")

args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman",
#     "--test_train_config_id", "8443ed21d2ac4-f8716281f960b", "--working_dir",
#     "/scratch/ahlmanne/perturbation_prediction_benchmark", "--result_id", "0"])
print(args)

out_dir = args.working_dir + "/results/" + args.result_id

np.random.seed(args.seed)
# --------------------------------------------------------


pert_data_folder = Path("data/gears_pert_data/")
pert_data = PertData(pert_data_folder)
if args.dataset_name in ['norman', 'adamson', 'dixit']:
  pert_data.load(args.dataset_name)
else:
  pert_data.load(data_path = "data/gears_pert_data/" + args.dataset_name)
adata = pert_data.adata

def mean_gene_expr_for_condition(condition):
    sub_data = adata[adata.obs['condition'] == condition,:]
    return np.array(sub_data.X.mean(axis=0)).flatten()
  
def se_gene_expr_for_condition(condition):
    sub_data = adata[adata.obs['condition'] == condition,:]
    return sub_data.X.toarray().std(axis=0) / sub_data.shape[0]

conds = adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
split_conds = [x.split("+") for x in conds]
split_conds = [list(filter(lambda y: y != "ctrl", x)) for x in split_conds]

res_df = (pd.DataFrame({
 "condition": conds,
 "split_conditions": split_conds})
  .assign(obs_mean = lambda x: [mean_gene_expr_for_condition(y) for y in x.condition])
  .assign(obs_se = lambda x: [se_gene_expr_for_condition(y) for y in x.condition])
  .assign(n_cells = lambda x: [sum(adata.obs['condition'] == y) for y in x.condition])
)


all_pred_vals = dict(zip(res_df['condition'], res_df["obs_mean"]))
all_pred_vals = {k: v.tolist() for k, v in all_pred_vals.items()}

all_pred_se_vals = dict(zip(res_df['condition'], res_df["obs_se"]))
all_pred_se_vals = {k: v.tolist() for k, v in all_pred_se_vals.items()}

n_cells_vals = dict(zip(res_df['condition'], res_df["n_cells"]))

tmp_out_dir = tempfile.mkdtemp()
with open(f"{tmp_out_dir}/all_predictions.json", 'w', encoding="utf8") as handle:
    json.dump(all_pred_vals, handle, indent = 4)
with open(f"{tmp_out_dir}/all_predictions_se.json", 'w', encoding="utf8") as handle:
    json.dump(all_pred_se_vals, handle, indent = 4)
with open(f"{tmp_out_dir}/n_cells.json", 'w', encoding="utf8") as handle:
    json.dump(n_cells_vals, handle, indent = 4)
with open(f"{tmp_out_dir}/gene_names.json", 'w', encoding="utf8") as handle:
    json.dump(pert_data.adata.var["gene_name"].values.tolist(), handle, indent = 4)
    
# Move results to out_dir
shutil.move(tmp_out_dir, out_dir)


session_info.show()
print("Python done")
