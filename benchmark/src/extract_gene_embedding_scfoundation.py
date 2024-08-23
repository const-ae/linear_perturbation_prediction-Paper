from pathlib import Path
import numpy as np
import pandas as pd
import tempfile
import argparse
import session_info

import torch
import anndata as ad

parser = argparse.ArgumentParser(description='Extract gene embedding from scFoundation')

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")
args = parser.parse_args()
# args = parser.parse_args(["--working_dir", "/scratch/ahlmanne/perturbation_prediction_benchmark", 
#                           "--result_id", "0"])
print(args)

# Both files were downloaded from the scFoundation's Github
demo_adata = ad.read_h5ad("/home/ahlmanne/huber/data/scfoundation_model/demo.h5ad")
singlecell_model_path="/home/ahlmanne/huber/data/scfoundation_model/models.ckpt"
ckp = torch.load(singlecell_model_path)
gene_pos_emb = ckp['gene']['state_dict']['model.pos_emb.weight'].cpu().numpy()
gene_names = demo_adata.var.gene_name.tolist()
gene_names = gene_names + ["log10TotalCount1", "log10TotalCount2", "<pad>"]


df = pd.DataFrame(data = gene_pos_emb.transpose(), columns = gene_names)

outfile = args.working_dir + "/results/" + args.result_id
df.to_csv(outfile, sep = "\t", index = False)



session_info.show()
print("Python done")
