
from pathlib import Path
import shutil
import numpy as np
import torch
import json
import tempfile
import argparse
import session_info

import cpa
import scanpy as sc
import anndata as ad

parser = argparse.ArgumentParser(description='Run scVI')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")

args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman_from_scfoundation",
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


# Convert adata.X[ij] = log(c / (s / 10k) + 1) to counts
counts = np.expm1(adata.X) / 1e4
counts = counts.tocsr()
# The smallest non-zero element used to be a 1
one_equiv = np.array([np.min(counts[i,:].data) for i in range(counts.shape[0])])
counts = (counts / one_equiv[:,np.newaxis]).tocsr()
new_adata = ad.AnnData(counts, obs = adata.obs, var = adata.var)

new_adata.obs['dose_val'] = ['1.0' if "ctrl" in x else "1.0+1.0" for x in new_adata.obs['condition']]
new_adata.obs['condition_harm'] = [x.replace("+ctrl", "").replace("ctrl+", "") for x in new_adata.obs['condition']]

cpa.CPA.setup_anndata(new_adata,
                    perturbation_key='condition_harm',
                    control_group='ctrl',
                    dosage_key='dose_val',
                    is_count_data=True,
                    max_comb_len=2
                   )
new_adata.obs['split'] = "train"
new_adata.obs['split'][new_adata.obs['condition'].isin(set2conditions['test'])] = "test"
new_adata.obs['split'][new_adata.obs['condition'].isin(set2conditions['val'])] = "ood"

model_params = {
    "n_latent": 32,
    "recon_loss": "nb",
    "doser_type": "linear",
    "n_hidden_encoder": 256,
    "n_layers_encoder": 4,
    "n_hidden_decoder": 256,
    "n_layers_decoder": 2,
    "use_batch_norm_encoder": True,
    "use_layer_norm_encoder": False,
    "use_batch_norm_decoder": False,
    "use_layer_norm_decoder": False,
    "dropout_rate_encoder": 0.2,
    "dropout_rate_decoder": 0.0,
    "variational": False,
    "seed": 8206,
}

trainer_params = {
    "n_epochs_kl_warmup": None,
    "n_epochs_adv_warmup": 50,
    "n_epochs_mixup_warmup": 10,
    "n_epochs_pretrain_ae": 10,
    "mixup_alpha": 0.1,
    "lr": 0.0001,
    "wd": 3.2170178270865573e-06,
    "adv_steps": 3,
    "reg_adv": 10.0,
    "pen_adv": 20.0,
    "adv_lr": 0.0001,
    "adv_wd": 7.051355554517135e-06,
    "n_layers_adv": 2,
    "n_hidden_adv": 128,
    "use_batch_norm_adv": True,
    "use_layer_norm_adv": False,
    "dropout_rate_adv": 0.3,
    "step_size_lr": 25,
    "do_clip_grad": False,
    "adv_loss": "cce",
    "gradient_clip_value": 5.0,
}

model = cpa.CPA(adata=new_adata,
                split_key='split',
                train_split='train',
                valid_split='test',
                test_split='ood',
                **model_params,
               )

with tempfile.TemporaryDirectory() as tmp_dir:
  model.train(max_epochs=2000,
              use_gpu=True,
              batch_size=2048,
              plan_kwargs=trainer_params,
              early_stopping_patience=5,
              check_val_every_n_epoch=5,
              save_path=tmp_dir,
             )


new_adata.layers['X_true'] = new_adata.X.copy()

# Override the newdata.X with random ctrl cells
ctrl_adata = new_adata[new_adata.obs['condition'] == 'ctrl'].copy()

new_adata.X = ctrl_adata.X[np.random.choice(ctrl_adata.n_obs, size=new_adata.n_obs, replace=True), :]


model.predict(new_adata, batch_size=2048)
print(f"CPA_pred max value {new_adata.obsm['CPA_pred'].max()}")
print(f"X_true max value {new_adata.layers['X_true'].max()}")
new_adata.layers['CPA_pred'] = new_adata.obsm['CPA_pred'].copy()

sc.pp.normalize_total(new_adata, target_sum=1e4, layer='CPA_pred')
sc.pp.log1p(new_adata, layer='CPA_pred')
sc.pp.normalize_total(new_adata, target_sum=1e4, layer='X_true')
sc.pp.log1p(new_adata, layer='X_true')

def mean_gene_expr_for_condition(condition):
    sub_data = new_adata[new_adata.obs['condition'] == condition,:]
    sub_data = sub_data.obsm['CPA_pred'].copy()
    sub_data = np.log1p(sub_data)
    return np.array(sub_data.mean(axis=0)).flatten()

conds = new_adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
all_pred_vals = {k: mean_gene_expr_for_condition(k).tolist() for k in conds}

tmp_out_dir = tempfile.mkdtemp()
with open(f"{tmp_out_dir}/all_predictions.json", 'w', encoding="utf8") as handle:
    json.dump(all_pred_vals, handle, indent = 4)
with open(f"{tmp_out_dir}/gene_names.json", 'w', encoding="utf8") as handle:
    json.dump(pert_data.adata.var["gene_name"].values.tolist(), handle, indent = 4)

# Move results to out_dir
shutil.move(tmp_out_dir, out_dir)



session_info.show()
print("Python done")
