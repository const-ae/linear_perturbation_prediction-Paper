from pathlib import Path
import shutil
import numpy as np
import json
import tempfile
import argparse
import session_info


import sys
 # scfoundation uses a forked version of GEARS which
sys.path.insert(0, "/g/huber/users/ahlmanne/projects/perturbation_prediction-benchmark/tmp/scfoundation/scfoundation_gears/")
sys.path.append("/g/huber/users/ahlmanne/projects/perturbation_prediction-benchmark/tmp/scfoundation/model/")
import gears.version
assert gears.version.__version__ == '0.0.2'
from gears import PertData, GEARS
from gears.utils import filter_pert_in_go

parser = argparse.ArgumentParser(description='Run scfoundation with gears')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--epochs', dest = 'epochs', action = 'store', help = "How many epochs are run", default = 15, type = int)
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

with open(args.working_dir + "/results/" + args.test_train_config_id) as json_file:
  set2conditions = json.load(json_file)

print(set2conditions)
set2conditions["train"] = list(set(set2conditions["train"]).difference(set(["IER5L+ctrl", "LYL1+IER5L"])))
set2conditions["test"] = list(set(set2conditions["test"]).difference(set(["IER5L+ctrl", "LYL1+IER5L"])))
set2conditions["val"] = list(set(set2conditions["val"]).difference(set(["IER5L+ctrl", "LYL1+IER5L"])))
print(set2conditions)

pert_data.set2conditions = set2conditions
pert_data.split = "custom"
pert_data.subgroup = None
pert_data.seed = 1
pert_data.train_gene_set_size = 0.75

# These are based on https://github.com/biomap-research/scFoundation/blob/69b0710660aded7e071d4f67e9f3ac03b096e587/GEARS/run_sh/run_singlecell_maeautobin-0.1B-res0-norman.sh
batch_size=6
accumulation_steps=5
test_batch_size=6
hidden_size=512
model_type="maeautobin"
bin_set="autobin_resolution_append"
# This file was downloaded separately from https://hopebio2020-my.sharepoint.com/:f:/g/personal/dongsheng_biomap_com/Eh22AX78_AVDv6k6v4TZDikBXt33gaWXaz27U9b1SldgbA
singlecell_model_path="/home/ahlmanne/huber/data/scfoundation_model/models.ckpt"
finetune_method="frozen"
train_gene_set_size=0.75
mode="v1"
highres=0 # 0
lr=0.01 #1e-3

pert_data.get_dataloader(batch_size = batch_size, test_batch_size = test_batch_size)
gears_model = GEARS(pert_data, device = 'cuda')
gears_model.model_initialize(hidden_size = hidden_size, 
                             model_type = model_type,
                             bin_set=bin_set,
                             load_path=singlecell_model_path,
                             finetune_method=finetune_method,
                             accumulation_steps=accumulation_steps,
                             mode=mode,
                             highres=highres)

temp_dir = tempfile.TemporaryDirectory()
gears_model.train(epochs = args.epochs, result_dir=temp_dir.name,lr=lr)
temp_dir.cleanup()


tmp_out_dir = tempfile.mkdtemp()
gears_model.save_model(f'{tmp_out_dir}/gears_model')

conds = pert_data.adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
split_conds = [x.split("+") for x in conds]
split_conds = [list(filter(lambda y: y != "ctrl", x)) for x in split_conds]

all_pred_vals = gears_model.predict(split_conds)
all_pred_vals = {k: v.tolist() for k, v in all_pred_vals.items()}
with open(f"{tmp_out_dir}/all_predictions.json", 'w', encoding="utf8") as handle:
    json.dump(all_pred_vals, handle, indent = 4)
with open(f"{tmp_out_dir}/gene_names.json", 'w', encoding="utf8") as handle:
    json.dump(pert_data.adata.var["gene_name"].values.tolist(), handle, indent = 4)

# Move results to out_dir
shutil.move(tmp_out_dir, out_dir)



session_info.show()
print("Python done")
