library(SingleCellExperiment)
library(tidyverse)
Sys.setenv("BASILISK_EXTERNAL_CONDA"="/g/easybuild/x86_64/Rocky/8/haswell/software/Miniforge3/24.1.2-0")

pa <- argparser::arg_parser("Run linear pretrained model for double perturbations")
pa <- argparser::add_argument(pa, "--dataset_name", type = "character", help = "The name of the dataset") 
pa <- argparser::add_argument(pa, "--test_train_config_id", type = "character", help = "The ID of the test/train/holdout run") 
pa <- argparser::add_argument(pa, "--pca_dim", type = "integer", default = 10, nargs = 1, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--ridge_penalty", type = "numeric", default = 0.1, nargs = 1, help = "The ridge penalty")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")

pa <- argparser::add_argument(pa, "--gene_embedding", type = "character", default = "training_data", help = "The path to a tsv with a gene embedding or the string 'training_data'")
pa <- argparser::add_argument(pa, "--pert_embedding", type = "character", default = "training_data", help = "The path to a tsv with a perturbation embedding or the string 'training_data'")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --dataset_name norman_from_scfoundation
#                             --test_train_config_id ede703be5dd78-c8355c044d3ba
#                             --pca_dim 4
#                             --gene_embedding training_data
#                             --pert_embedding training_data
#                             --ridge_penalty 0.1
#                             --working_dir /scratch/ahlmanne/perturbation_prediction_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

solve_y_axb <- function(Y, A = NULL, B = NULL, A_ridge = 0.01, B_ridge = 0.01){
  stopifnot(is.matrix(Y) || is(Y, "Matrix"))
  stopifnot(is.null(A) || is.matrix(A) || is(A, "Matrix"))
  stopifnot(is.null(B) || is.matrix(B) || is(B, "Matrix"))
  
  center <- rowMeans(Y)
  Y <- Y - center
  
  if(! is.null(A) && ! is.null(B)){
    stopifnot(nrow(Y) == nrow(A))
    stopifnot(ncol(Y) == ncol(B))
    # fit <- lm.fit(kronecker(t(B), A), as.vector(Y))
    tmp <- as.matrix(Matrix::solve(t(A) %*% A + Matrix::Diagonal(ncol(A)) * A_ridge) %*% t(A) %*% Y %*% t(B) %*% Matrix::solve(B %*% t(B) + Matrix::Diagonal(nrow(B)) * B_ridge))
  }else if(is.null(B)){
    fit <- lm.fit(A, Y)
    tmp <- as.matrix(Matrix::solve(t(A) %*% A + Matrix::Diagonal(ncol(A)) * A_ridge) %*% t(A) %*% Y)
  }else if(is.null(A)){
    fit <- lm.fit(t(B), t(Y))
    tmp <- as.matrix(Y %*% t(B) %*% Matrix::solve(B %*% t(B) + Matrix::Diagonal(nrow(B)) * B_ridge))
  }else{
    stop("Either A or B must be non-null")
  }
  tmp[is.na(tmp)] <- 0
  list(K = tmp, center = center)
}


# Load data
folder <- "data/gears_pert_data"
sce <- zellkonverter::readH5AD(file.path(folder, pa$dataset_name, "perturb_processed.h5ad"))
set2condition <- rjson::fromJSON(file = file.path(pa$working_dir, "results", pa$test_train_config_id))
if(! "ctrl" %in% set2condition$train){
  set2condition$train <- c(set2condition$train, "ctrl")
}

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

gene_names <- rowData(sce)[["gene_name"]]
rownames(sce) <- gene_names

baseline <- MatrixGenerics::rowMeans2(assay(sce, "X")[,sce$condition == "ctrl",drop=FALSE])

# Pseudobulk everything
psce <- glmGamPoi::pseudobulk(sce, group_by = vars(condition, clean_condition, training))
assay(psce, "change") <- assay(psce, "X") - baseline


train_data <- psce[,psce$training == "train"]


# Get embeddings
gene_emb <- if(pa$gene_embedding == "training_data"){
  pca <- irlba::prcomp_irlba(as.matrix(assay(train_data, "X")), n = pa$pca_dim)
  rownames(pca$x) <- rownames(train_data)
  pca$x
}else if(pa$gene_embedding == "identity"){
  tmp <- Matrix::Diagonal(n = nrow(psce))
  rownames(tmp) <- rownames(psce)
  tmp
}else if(pa$gene_embedding == "zero"){
  tmp <- Matrix::Diagonal(n = nrow(psce)) - Matrix::Diagonal(n = nrow(psce))
  rownames(tmp) <- rownames(psce)
  tmp
}else if(pa$gene_embedding == "random"){
  tmp <- matrix(rnorm(nrow(psce) * pa$pca_dim), nrow = nrow(psce), ncol = pa$pca_dim)
  rownames(tmp) <- rownames(psce)
  tmp
}else{
  t(as.matrix(data.table::fread(pa$gene_embedding)))
}

pert_emb <- if(pa$pert_embedding == "training_data"){
  pca <- irlba::prcomp_irlba(as.matrix(assay(train_data, "X")), n = pa$pca_dim)
  rownames(pca$x) <- rownames(train_data)
  t(pca$x)
}else if(pa$pert_embedding == "identity"){
  tmp <- Matrix::Diagonal(n = ncol(psce))
  colnames(tmp) <- psce$clean_condition
  tmp
}else if(pa$pert_embedding == "zero"){
  tmp <- Matrix::Diagonal(n = ncol(psce)) - Matrix::Diagonal(n = ncol(psce))
  colnames(tmp) <- psce$clean_condition
  tmp
}else if(pa$pert_embedding == "random"){
  tmp <- matrix(rnorm(ncol(psce) * pa$pca_dim), nrow = pa$pca_dim, ncol = ncol(psce))
  colnames(tmp) <- psce$clean_condition
  tmp
}else{
  as.matrix(data.table::fread(pa$pert_embedding))
}
if(! "ctrl" %in% colnames(pert_emb)){
  pert_emb <- cbind(pert_emb, ctrl = rep(0, nrow(pert_emb)))
}

cond1 <- str_remove(train_data$condition, "\\+.+")
cond2 <- str_remove(train_data$condition, ".+\\+")

pert_matches1 <- match(cond1, colnames(pert_emb))
pert_matches2 <- match(cond2, colnames(pert_emb))
stopifnot(length(pert_matches2) == length(pert_matches1))
if(length(pert_matches1) <= 1){
  stop("Too few matches between clean_conditions and pert_embedding")
}
gene_matches <- match(rownames(gene_emb), rownames(train_data))
if(sum(! is.na(gene_matches)) <= 1){
  stop("Too few matches between gene names and gene_embedding")
}


gene_emb_sub <- gene_emb[! is.na(gene_matches),,drop=FALSE]
pert_emb_training <- pert_emb[,na.omit(pert_matches1),drop=FALSE] + pert_emb[,na.omit(pert_matches1),drop=FALSE]
Y <-  assay(train_data, "change")[na.omit(gene_matches), ! is.na(pert_matches1) & ! is.na(pert_matches2),drop=FALSE]
coefs <- solve_y_axb(Y = Y, A = gene_emb_sub, B = pert_emb_training,
                     A_ridge = pa$ridge_penalty, B_ridge = pa$ridge_penalty)


cond1_all <- str_remove(psce$condition, "\\+.+")
cond2_all <- str_remove(psce$condition, ".+\\+")
pert_matches_all1 <- match(cond1_all, colnames(pert_emb))
pert_matches_all2 <- match(cond2_all, colnames(pert_emb))
pert_emb_all <- pert_emb[,na.omit(pert_matches_all1),drop=FALSE] + pert_emb[,na.omit(pert_matches_all2),drop=FALSE]
colnames(psce) <- psce$clean_condition

baseline <- baseline[na.omit(gene_matches)]

pred <- as.matrix(gene_emb_sub %*% coefs$K %*% pert_emb_all + coefs$center + baseline)


rownames(pred) <- rownames(psce)[na.omit(gene_matches)]
colnames(pred) <- psce$clean_condition

# Print some nice summary statistics
bind_cols(tibble(cond = psce$clean_condition, training = psce$training), 
          as_tibble(t(assay(psce, "change"))) %>% rename_with(~ paste0("obs-+-", .x)),
          as_tibble(t(pred - baseline))  %>% rename_with(~ paste0("pred-+-", .x)))  %>%
  pivot_longer(starts_with(c("pred-+-", "obs-+-")), names_sep = "-\\+-", names_to = c(".value", "gene")) %>%
  tidylog::drop_na() %>%
  mutate( is_training = training == "train") %>%
  summarize(r2 = cor(obs, pred),
            l2 = sqrt(sum((obs - pred)^2)),
            .by = c(cond, is_training)) %>%
  group_by(is_training) %>%
  skimr::skim(r2, l2)

# Store output
tmp <- as.list(as.data.frame(pred))
tmp_out_dir <- file.path(tempdir(), "prediction_storage")
dir.create(tmp_out_dir)
write_lines(rjson::toJSON(tmp), file.path(tmp_out_dir, "all_predictions.json"))
write_lines(rjson::toJSON(rownames(pred)), file.path(tmp_out_dir, "gene_names.json"))
file.rename(tmp_out_dir,out_dir)


#### Session Info
sessionInfo()


