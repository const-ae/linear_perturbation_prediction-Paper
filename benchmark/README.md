# Overview

This folder contains the Python and R scripts to reproduce the results presented in our paper.

The `submission` folder contains a file called `run_perturbation_benchmark.R` that launches a series of slurm jobs which run `GEARS`, `scFoundation`, `scGPT`, and our linear models on the three single perturbation and one double perturbation datasets.

The jobs are defined using my own workflow manager (creatively named [MyWorkflowManager](https://github.com/const-ae/MyWorkflowManager)). Each job is defined by the executed script and the command line arguments used to call it. In addition, each job is run in a separate environment with reproducible dependencies. 

# Prepare environments

The first time you run `R` (make sure to use version 4.4.1) in the `benchmark` folder, the `.Rprofile` file is loaded and calls `source("renv/activate.R")` which automatically prepares [`renv`](https://rstudio.github.io/renv/articles/renv.html). This process can take a few minutes. Then call 
```
renv::restore()
```
to install all dependencies.

For the Python environments, I used `conda` from [Miniforge](https://conda-forge.org/download/), which you can recreate with
```
 conda env create --file conda_environments/<env>-environment.yml
```
where you replace `<env>` with the name of the environment.

# Reproduce all results

To reproduce all results, open the `run_perturbation_benchmark.R`, adjust the working directory (line 1) and the folder for the intermediate results (line 5). Then, source the `make_single_perturbation_jobs` and `make_double_perturbation_jobs` functions and call
```r
double_pert_jobs <- make_double_perturbation_jobs(datasets = c("norman_from_scfoundation"), seeds = 1:5)
run_job(double_pert_jobs, priority = "normal")

single_pert_jobs <- make_single_perturbation_jobs(datasets = c("adamson", 'replogle_k562_essential', 'replogle_rpe1_essential'), seeds = 1:2)
run_job(single_pert_jobs, priority = "normal")
```

This will launch approximately 100 slurm jobs, which will finish within a few hours to a few days. When all jobs complete, the results will be in `result_file_path(single_pert_jobs)` and `result_file_path(double_pert_jobs)`. These are the input files used in the notebooks in the parent folder.

# Reproduce some results

If you don't want to reproduce all results, you can also directly call the scripts in `src` directly with the correct command line arguments. I will illustrate the process with the Adamson data. 

First, we will will setup some folders:

```shell
# This could be anywhere
mkdir working_dir
mkdir working_dir/results
mkdir -p data/gears_pert_data
```

Next, make sure to load the necessary conda environment:
```shell
conda activate gears_env2
```

Now we will call `src/prepare_perturbation_data` to download the data and split all perturbations into a training test, and validation set:
```shell
python3 src/prepare_perturbation_data.py \
  --dataset_name adamson \
  --seed 1 \
  --working_dir /tmp/working_dir \
  --result_id test_train_split
```

In the next step we will use the output (`/tmp/working_dir/results/test_train_split`) to run the linear model (`src/run_linear_pretrained_model.R`):
```shell
Rscript --no-restore src/run_linear_pretrained_model.R \
    --dataset_name adamson \
    --test_train_config_id test_train_split \
    --pca_dim 10 \
    --gene_embedding training_data \
    --pert_embedding training_data \
    --working_dir /tmp/working_dir \
    --result_id linear_prediction
```

Similarly, we can now run scFoundation (`src/run_scfoundation.py`), GEARS (`src/run_gears.py`), or scGPT (`src/run_scgpt.py`). We can also calculate the ground truth by calling `src/ground_truth_combinatorial_prediction`. Just remember to load the right conda environment each time before executing the script.

```shell
conda activate gears_env2
python3 src/run_ground_truth_for_combinatorial_perturbations.py \
    --dataset_name adamson \
    --test_train_config_id test_train_split \
    --working_dir /tmp/working_dir \
    --result_id ground_truth

conda activate flashattn_env    
python3 src/run_scgpt.py \
    --dataset_name adamson \
    --test_train_config_id test_train_split \
    --working_dir /tmp/working_dir \
    --result_id ground_truth    

conda activate gears_env2
python3 src/run_gears.py \
    --dataset_name adamson \
    --test_train_config_id test_train_split \
    --working_dir working_dir \
    --result_id ground_truth
```

Each script produces a JSON file with the predictions of the gene expression for each perturbation and one JSON file listing the genes (as the order might differ). You can load the results with R and plot them:

```r
library(tidyverse)

test_train <- rjson::fromJSON(file = "/tmp/working_dir/results/test_train_split") %>%
  enframe(name = "train", value = "perturbation") %>%
  unnest(perturbation)

load_results <- function(path){
  preds <- rjson::fromJSON(file = file.path(path, "all_predictions.json"))
  names <- rjson::fromJSON(file = file.path(path, "gene_names.json"))
  tibble(perturbation = names(preds), prediction = unname(preds), gene = list(names)) %>%
    unnest(c(prediction, gene)) %>%
    mutate(perturbation = map_chr(str_split(perturbation, pattern = "[+_]", n = 2), \(x) {
      tmp <- if(all(x == "ctrl" | x == "")) "ctrl" else if(length(x) == 2) x else c(x, "ctrl")
      paste0(tmp, collapse = "+")
    }))
}

ground_truth <- load_results("/tmp/working_dir/results/ground_truth/") %>%
  rename(truth = prediction)
preds <- bind_rows(linear = load_results("/tmp/working_dir/results/linear_results/"),
                   gears = load_results("/tmp/working_dir/results/gears_results/"), 
                   scGPT = load_results("/tmp/working_dir/results/sgpt_results/"), 
                   .id = "method")

highest_expr_genes <- ground_truth %>%
  filter(perturbation == "ctrl") %>%
  slice_max(truth, n = 1000) %>%
  select(gene, baseline = truth)

res <- inner_join(preds, ground_truth, by = c("gene", "perturbation")) %>%
  inner_join(highest_expr_genes, by = "gene") %>%
  summarize(dist = sqrt(sum((truth - prediction)^2)), 
            pearson_delta = cor(prediction - baseline, truth - baseline), 
            .by = c("perturbation", "method")) 

res %>%
  inner_join(test_train, by = "perturbation") %>%
  filter( train != "train") %>%
  ggplot(aes(x = method, y = dist)) +
    ggbeeswarm::geom_quasirandom()
```



