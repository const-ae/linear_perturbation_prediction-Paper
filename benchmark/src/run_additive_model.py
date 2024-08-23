from gears import PertData, GEARS
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import json
import tempfile
import argparse
import session_info


parser = argparse.ArgumentParser(description='Run simple additive model7')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")

args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman",
#     "--test_train_config_id", "badee280e9f4d-f8716281f960b", "--working_dir",
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

with open(args.working_dir + "/results/" + args.test_train_config_id) as json_file:
  set2conditions = json.load(json_file)


def mean_gene_expr_for_condition(condition):
    sub_data = adata[adata.obs['condition'] == condition,:]
    return np.array(sub_data.X.mean(axis=0)).flatten()

def sum_single_predictions(my_conds):
    pred = baseline.copy()
    for c in my_conds:
        pred += single_lookup[c] - baseline
    return pred

conds = adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
split_conds = [x.split("+") for x in conds]
split_conds = [list(filter(lambda y: y != "ctrl", x)) for x in split_conds]

baseline = np.array(adata[adata.obs['condition'] == 'ctrl',:].X.mean(axis=0)).flatten()

training_df = (pd.DataFrame({"training": set2conditions.keys(),
              "condition": set2conditions.values()})
  .explode("condition"))

base_df = (pd.DataFrame({
 "condition": conds,
 "split_conditions": split_conds})
  .merge(training_df, how = "left", on = "condition")
)

single_expr_df = (base_df
  .assign(single = lambda x: [len(y) == 1 for y in x.split_conditions])
  .query('single')
  .assign(cond_name = lambda x: [y[0] for y in x.split_conditions])
  .assign(obs_mean = lambda x: [mean_gene_expr_for_condition(y) for y in x.condition])
)

single_lookup = dict(zip(single_expr_df['cond_name'], single_expr_df["obs_mean"]))

res_df = (base_df
  .assign(single = lambda x: [len(y) == 1 for y in x.split_conditions])
  .assign(simple_pred = lambda x: [sum_single_predictions(y) for y in x.split_conditions])
  # .assign(obs = lambda x: [mean_gene_expr_for_condition(y) for y in x.condition])
)


all_pred_vals = dict(zip(res_df['condition'], res_df["simple_pred"]))
all_pred_vals = {k: v.tolist() for k, v in all_pred_vals.items()}

tmp_out_dir = tempfile.mkdtemp()
with open(f"{tmp_out_dir}/all_predictions.json", 'w', encoding="utf8") as handle:
    json.dump(all_pred_vals, handle, indent = 4)
with open(f"{tmp_out_dir}/gene_names.json", 'w', encoding="utf8") as handle:
    json.dump(pert_data.adata.var["gene_name"].values.tolist(), handle, indent = 4)

# Move results to out_dir
shutil.move(tmp_out_dir, out_dir)


session_info.show()
print("Python done")
