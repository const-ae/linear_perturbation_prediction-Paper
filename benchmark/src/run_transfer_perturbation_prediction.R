library(SingleCellExperiment)
library(tidyverse)
Sys.setenv("BASILISK_EXTERNAL_CONDA"="/g/easybuild/x86_64/Rocky/8/haswell/software/Miniforge3/24.1.2-0")

pa <- argparser::arg_parser("Run transfer perturbation prediction")
pa <- argparser::add_argument(pa, "--dataset_name", type = "character", help = "The name of the dataset") 
pa <- argparser::add_argument(pa, "--test_train_config_id", type = "character", help = "The ID of the test/train/holdout run") 
pa <- argparser::add_argument(pa, "--reference_data", type = "character", help = "The name of the reference dataset")
pa <- argparser::add_argument(pa, "--pca_dim", type = "integer", default = 30, nargs = 1, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--ridge_penalty", type = "numeric", default = 0.1, nargs = 1, help = "The ridge penalty")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --dataset_name adamson
#                             --test_train_config_id fe90cfa1ab324-aab07fc0fac08
#                             --reference_data replogle_k562_essential
#                             --pca_dim 4
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

# Load data
folder <- "data/gears_pert_data"
sce <- zellkonverter::readH5AD(file.path(folder, pa$dataset_name, "perturb_processed.h5ad"))
ref <- zellkonverter::readH5AD(file.path(folder, pa$reference_data, "perturb_processed.h5ad"))
set2condition <- rjson::fromJSON(file = file.path(pa$working_dir, "results", pa$test_train_config_id))

# Only keep valid condtions
sce <- sce[,sce$condition %in% unlist(set2condition)]

# Clean up the colData(sce) a bit
sce$condition <- droplevels(sce$condition)
sce$clean_condition <- stringr::str_remove(sce$condition, "\\+ctrl")
training_df <- tibble(training = names(set2condition), condition = set2condition) %>%
  unnest(condition)
colData(sce) <- colData(sce) %>%
  as_tibble() %>%
  tidylog::left_join(training_df, by = "condition") %>%
  DataFrame()

ref$condition <- droplevels(ref$condition)
ref$clean_condition <- stringr::str_remove(ref$condition, "\\+ctrl")

baseline <- MatrixGenerics::rowMeans2(assay(sce, "X")[,sce$condition == "ctrl",drop=FALSE])
ref_baseline <- MatrixGenerics::rowMeans2(assay(ref, "X")[,ref$condition == "ctrl",drop=FALSE])

# Pseudobulk everything
psce <- glmGamPoi::pseudobulk(sce, group_by = vars(condition, clean_condition, training))
ref_psce <- glmGamPoi::pseudobulk(ref, group_by = vars(condition, clean_condition))

# Compare perturbation and gene names
ref_conds <- unique(ref_psce$clean_condition)
conds <- unique(psce$clean_condition)
# conds <- setdiff(conds, "ctrl")
print("Number of perturbations that are reference: ")
print(table(conds %in% ref_conds))

psce$in_ref <- psce$clean_condition %in% ref_conds

# Subtract baseline
assay(psce, "centered_mat", withDimnames = FALSE) <- as.matrix(assay(psce, "X")) - baseline
assay(ref_psce, "centered_mat", withDimnames = FALSE) <- as.matrix(assay(ref_psce, "X")) - ref_baseline

# Work on training data!
train_data <- psce[,psce$training == "train"]

# Relate the embeddings with ridge regression
train_idx_matches <- match(train_data$clean_condition,ref_psce$clean_condition)
train_pert_match <- seq_len(ncol(train_data))[! is.na(train_idx_matches)]
train_ref_match <- na.omit(train_idx_matches)

idx_matches <- match(psce$clean_condition,ref_psce$clean_condition)
pert_match <- seq_len(ncol(psce))[! is.na(idx_matches)]
ref_match <- na.omit(idx_matches)

# Build embeddings
pert_pca <- lemur:::pca(assay(train_data, "centered_mat"), n = pa$pca_dim)
ref_pca <- lemur:::pca(assay(ref_psce, "centered_mat"), n = pa$pca_dim)

X_train <- cbind(1, t(ref_pca$embedding[,train_ref_match,drop=FALSE]))
beta <- lemur:::ridge_regression(pert_pca$embedding[,train_pert_match,drop=FALSE], X_train,
                                 ridge_penalty = c(0, rep(pa$ridge_penalty, times = min(ncol(X_train)-1, pa$pca_dim))))

pred <- matrix(NA, nrow = nrow(psce), ncol = ncol(psce))
pred[,pert_match] <- pert_pca$coordsystem %*% beta %*% rbind(1, ref_pca$embedding[,ref_match,drop=FALSE]) + pert_pca$offset + baseline
rownames(pred) <- rownames(psce)

# Print some nice summary statistics
bind_cols(tibble(cond = psce$clean_condition, training = psce$training),
          as_tibble(t(assay(psce, "centered_mat"))) %>% rename_with(~ paste0("obs-", .x)),
          as_tibble(t(pred))  %>% rename_with(~ paste0("pred-", .x)))  %>%
  pivot_longer(starts_with(c("pred-", "obs-")), names_sep = "-", names_to = c(".value", "gene")) %>%
  summarize(r2 = cor(obs, pred), .by = c(cond, training)) %>%
  group_by(training) %>%
  skimr::skim(r2)

# Store output
colnames(pred) <- psce$condition
tmp <- as.list(as.data.frame(pred))
tmp_out_dir <- tempdir()
write_lines(rjson::toJSON(tmp), file.path(tmp_out_dir, "all_predictions.json"))
write_lines(rjson::toJSON(rownames(pred)), file.path(tmp_out_dir, "gene_names.json"))
file.rename(tmp_out_dir,out_dir)


#### Session Info
sessionInfo()


