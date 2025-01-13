
from pathlib import Path
import shutil
import anndata as ad
import numpy as np
import pandas as pd
import torch
import json
import pickle
import tempfile
import argparse
import session_info

parser = argparse.ArgumentParser(description='Run GeneFormer')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--perturbation_type', dest = 'perturbation_type', action = 'store', required = True, help = "The Geneformer perturbation type ('delete' or 'overexpress'). For more information see: https://geneformer.readthedocs.io/en/latest/geneformer.in_silico_perturber.html")
parser.add_argument('--ncells_training', dest = 'ncells_training', action = 'store', default = 300 ,type = int)
parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")

args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman_from_scfoundation",
#     "--perturbation_type", "overexpress",
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
# The smallest non-zero element used to be a 1
one_equiv = np.array([np.min(counts[i,:].data) for i in range(counts.shape[0])])
counts = (counts / one_equiv[:,np.newaxis]).tocsr()
new_adata = ad.AnnData(counts, obs = adata.obs, var = adata.var)
new_adata.obs['n_counts'] = np.round(new_adata.X.sum(axis=1))

with open("/home/ahlmanne/prog/Geneformer/geneformer/ensembl_mapping_dict_gc95M.pkl", "rb") as f:
  gene2ensembl = pickle.load(f)

if args.dataset_name in ['adamson', 'norman', 'replogle_k562_essential', 'replogle_rpe1_essential']:
  # Already has ENSEMBL IDs
  new_adata.var['ensembl_id'] = new_adata.var.index
else:
  # Add ENSEMBL IDs 
  gene_info = pd.DataFrame.from_dict(gene2ensembl, orient = "index", columns = ["ensembl_id"])
  gene_info = gene_info.drop_duplicates(subset = "ensembl_id", keep = "first")
  if args.dataset_name == "dixit":
    new_adata.var = new_adata.var.set_index('gene_name', drop=False)
  
  new_adata.var = new_adata.var.join(gene_info, how = "left", validate = "one_to_one")
  new_adata = new_adata[:, new_adata.var.dropna().index]

print(f"Matched {new_adata.shape[1]} rows out of {adata.shape[1]}")

ctrl_adata = new_adata[new_adata.obs['condition'] == "ctrl",:].copy()
new_adata = new_adata[new_adata.obs['condition'].isin(set2conditions['train']), :]

# Geneformer
from geneformer import TranscriptomeTokenizer
from geneformer import Classifier

from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import perturber_utils as pu
GENEFORMER_MODEL_LOCATION = "/home/ahlmanne/prog/Geneformer/gf-12L-95M-i4096"


import tempfile
from datetime import datetime

trained_model_dir = tempfile.mkdtemp()

# Tokenization
def tokenize_adata(adata, output_dir = None):
  if not output_dir:
    output_dir = tempfile.mkdtemp()
  with  tempfile.TemporaryDirectory() as adata_dir:
    adata.write_h5ad(adata_dir + "/adata.h5ad")
    tk = TranscriptomeTokenizer({"condition": "condition"})
    tk.tokenize_data(data_directory = adata_dir, 
                     output_directory = output_dir, 
                     output_prefix = "tokenized", 
                     file_format="h5ad")
  return f"{output_dir}/tokenized.dataset"

ctrl_adata_tok = tokenize_adata(ctrl_adata)
train_adata_tok = tokenize_adata(new_adata)

# Classification (i.e., pretraining)
cc = Classifier(classifier="cell",
                cell_state_dict = {"state_key": "condition", "states": "all"},
                freeze_layers = 2,
                num_crossval_splits = 1,
                max_ncells_per_class = args.ncells_training,
                forward_batch_size=200,
                nproc=1)
cc.prepare_data(input_data_file = train_adata_tok,
                output_directory=trained_model_dir,
                output_prefix="classification",
                split_id_dict=None)
all_metrics = cc.validate(model_directory=GENEFORMER_MODEL_LOCATION,
                          prepared_input_data_file=f"{trained_model_dir}/classification_labeled_train.dataset",
                          id_class_dict_file=f"{trained_model_dir}/classification_id_class_dict.pkl",
                          output_directory=trained_model_dir,
                          output_prefix="classification_result")

today_date = datetime.now().strftime('%y%m%d')
trained_model_location = f"{trained_model_dir}/{today_date}_geneformer_cellClassifier_classification_result/ksplit1"


# In Silico Perturbation

# isp.perturb_data internally calculates the embedding after perturbation, but there is no way to extract it
# So I copied the code from https://huggingface.co/ctheodoris/Geneformer/blob/main/geneformer/in_silico_perturber.py#L421
# and https://huggingface.co/ctheodoris/Geneformer/blob/main/geneformer/in_silico_perturber.py#L761
# to export the embedding
def get_perturbed_embedding(isp: InSilicoPerturber, model, input_data):
    from geneformer.emb_extractor import get_embs

    isp.max_len = pu.get_model_input_size(model)
    layer_to_quant = pu.quant_layers(model) + isp.emb_layer

    ### filter input data ###
    # general filtering of input data based on filter_data argument


    def make_group_perturbation_batch(example):
        example_input_ids = example["input_ids"]
        example["tokens_to_perturb"] = isp.tokens_to_perturb
        indices_to_perturb = [
            example_input_ids.index(token) if token in example_input_ids else None
            for token in isp.tokens_to_perturb
        ]
        indices_to_perturb = [
            item for item in indices_to_perturb if item is not None
        ]
        if len(indices_to_perturb) > 0:
            example["perturb_index"] = indices_to_perturb
        else:
            # -100 indicates tokens to overexpress are not present in rank value encoding
            example["perturb_index"] = [-100]
        if isp.perturb_type == "delete":
            example = pu.delete_indices(example)
        elif isp.perturb_type == "overexpress":
            example = pu.overexpress_tokens(
                example, isp.max_len, isp.special_token
            )
            example["n_overflow"] = pu.calc_n_overflow(
                isp.max_len,
                example["length"],
                isp.tokens_to_perturb,
                indices_to_perturb,
            )
        return example


    perturbed_data = input_data.map(
        make_group_perturbation_batch, num_proc=isp.nproc
    ) 

    if isp.perturb_type == "overexpress":
        input_data = input_data.add_column(
            "n_overflow", perturbed_data["n_overflow"]
        )
        input_data = input_data.map(
            pu.truncate_by_n_overflow_special, num_proc=isp.nproc
        )


    perturbation_cls_emb = get_embs(
      model,
      perturbed_data,
      "cls",
      layer_to_quant,
      isp.pad_token_id,
      isp.forward_batch_size,
      token_gene_dict=isp.token_gene_dict,
      summary_stat=None,
      silent=True,
    )
    return perturbation_cls_emb.cpu().detach().numpy()

model = pu.load_model(
    "CellClassifier", num_classes = len(new_adata.obs['condition'].values.categories),
    model_directory = trained_model_location, mode="eval"
)
ctrl_dataset = pu.load_and_filter(
    None, 1, ctrl_adata_tok
)
start_isp = InSilicoPerturber(perturb_type=args.perturbation_type,
                      model_type="CellClassifier",
                      emb_mode="cls",
                      max_ncells=150,
                      emb_layer=0,
                      forward_batch_size=100,
                      nproc=1) 
ctrl_dataset = start_isp.apply_additional_filters(ctrl_dataset)

conds = adata.obs["condition"].cat.remove_unused_categories().cat.categories.tolist()
embeddings = {}
for idx, genes in enumerate(conds):
  gene_list = list(set(genes.split("+")) - set(["ctrl"]))
  mappable_genes = [gene2ensembl[x] for x in gene_list if x in gene2ensembl.keys()]
  if genes == "ctrl":
    embeddings[genes] = np.zeros(512)
  elif len(mappable_genes) == len(gene_list):
    isp = InSilicoPerturber(perturb_type=args.perturbation_type,
                          genes_to_perturb= mappable_genes,
                          model_type="CellClassifier",
                          emb_mode="cls",
                          max_ncells=150,
                          emb_layer=0,
                          forward_batch_size=100,
                          nproc=1) 
    local_emb = get_perturbed_embedding(isp, model = model, input_data = ctrl_dataset)
    embeddings[genes] = local_emb.mean(axis=0)
  else:
    pass

 
 
 
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
