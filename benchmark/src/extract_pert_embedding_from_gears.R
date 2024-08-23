library(tidyverse)

pa <- argparser::arg_parser("Extract perturbation embedding from GEARS' GO similarity")
pa <- argparser::add_argument(pa, "--pca_dim", type = "integer", default = 10, nargs = 1, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --pca_dim 10
#                             --working_dir /scratch/ahlmanne/perturbation_prediction_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

# -----------------
outfile <- file.path(pa$working_dir, "results/", pa$result_id)


# Get GO similarities
df_jaccard <- read_csv("data/gears_pert_data/go_essential_all/go_essential_all.csv")
genes <- unique(c(df_jaccard$source, df_jaccard$target))

sp_mat <- df_jaccard %>%
  mutate(source = factor(source, levels = genes), target = factor(target, levels = genes)) %>%
  (\(x) Matrix::sparseMatrix(i = as.integer(x$source), j = as.integer(x$target), x = x$importance,
                             dims = c(length(genes), length(genes)), dimnames = list(genes, genes)))
gr <- igraph::graph_from_adjacency_matrix(sp_mat, mode = "undirected", weighted = TRUE)
emb <- igraph::embed_adjacency_matrix(gr, no = 2, which = "lm", scaled = TRUE)

pert_emb <- t(emb$X)
colnames(pert_emb) <- genes


write_tsv(as_tibble(pert_emb), outfile)

