setwd("/g/huber/users/ahlmanne/projects/perturbation_prediction-benchmark")
library(tidyverse)
library(MyWorkflowManager)
source("submission/wrap_scripts.R")
init("/scratch/ahlmanne/perturbation_prediction_benchmark")


make_single_perturbation_jobs <- function(datasets = c('adamson', 'replogle_k562_essential', 'replogle_rpe1_essential'),
                                   seeds = 1L){
  jobs <- lapply(datasets, \(dataset){
    scgpt_emb <- scgpt_extract_embedding()
    scfoundation_emb <- scfoundation_extract_embedding()
    gears_go_emb <- gears_extract_pert_embedding(list(pca_dim = 10L))
    replogle_k562_emb <- pca_extract_pert_embedding(list(dataset_name = "replogle_k562_essential"))
    replogle_rpe1_emb <- pca_extract_pert_embedding(list(dataset_name = "replogle_rpe1_essential"))
    
    inner_jobs <- lapply(seeds, \(se){
      config_job <- prepare_perturbation_data(list(dataset = dataset,  seed = se))
      default_params <- list(dataset_name = dataset, test_train_config_id = config_job$result_id, seed = se)
      
      
      pert_jobs <- list(
        scgpt = scgpt_combinatorial_prediction(default_params, dep_jobs = list(config_job), duration = "30:00:00"),
        gears = gears_combinatorial_prediction(default_params, dep_jobs = list(config_job)),
        ground_truth = ground_truth_combinatorial_prediction(default_params, dep_jobs = list(config_job)),
        lpm_selftrained = linear_pretrained_model_prediction(c(list(gene_embedding = "training_data", pert_embedding = "training_data"), default_params), dep_jobs = list(config_job)),
        lpm_randomPertEmb = linear_pretrained_model_prediction(c(list(gene_embedding = "training_data", pert_embedding = "random"), default_params), dep_jobs = list(config_job)),
        lpm_randomGeneEmb = linear_pretrained_model_prediction(c(list(gene_embedding = "random", pert_embedding = "training_data"), default_params), dep_jobs = list(config_job)),
        lpm_scgptGeneEmb = linear_pretrained_model_prediction(c(list(gene_embedding = result_file_path(scgpt_emb), pert_embedding = "training_data"), default_params), dep_jobs = list(config_job)),
        lpm_scFoundationGeneEmb = linear_pretrained_model_prediction(c(list(gene_embedding = result_file_path(scfoundation_emb), pert_embedding = "training_data"), default_params), dep_jobs = list(config_job)),
        lpm_gearsPertEmb = linear_pretrained_model_prediction(c(list(gene_embedding = "training_data", pert_embedding = result_file_path(gears_go_emb)), default_params), dep_jobs = list(config_job)),
        lpm_k562PertEmb = linear_pretrained_model_prediction(c(list(gene_embedding = "training_data", pert_embedding = result_file_path(replogle_k562_emb)), default_params), dep_jobs = list(config_job)),
        lpm_rpe1PertEmb = linear_pretrained_model_prediction(c(list(gene_embedding = "training_data", pert_embedding = result_file_path(replogle_rpe1_emb)), default_params), dep_jobs = list(config_job))
      )
      pert_jobs
    })
    names(inner_jobs) <- seeds
    purrr::list_flatten(inner_jobs, name_spec = "{outer}-{inner}")
  })
  names(jobs) <- datasets
  jobs <- purrr::list_flatten(jobs, name_spec = "{outer}-{inner}")
  
  params <- list(
    job_ids = map_chr(jobs, "result_id"),
    names = names(jobs)
  )
  collect_perturbation_predictions(params, dep_jobs = jobs)
}

make_double_perturbation_jobs <- function(datasets = c('norman_from_scfoundation'),
                                   seeds = 1L){
  jobs <- lapply(datasets, \(dataset){
    inner_jobs <- lapply(seeds, \(se){
      config_job <- prepare_perturbation_data(list(dataset = dataset,  seed = se))
      default_params <- list(dataset_name = dataset, test_train_config_id = config_job$result_id, seed = se)
      mem <- if(dataset == "norman_from_scfoundation") "80GB" else "40GB"
      long_dur <- if(dataset == "norman_from_scfoundation") "6-00:00:00" else "10:00:00"
      pert_jobs <- list(
        scgpt = scgpt_combinatorial_prediction(default_params, dep_jobs = list(config_job), memory = mem, duration = long_dur),
        gears = gears_combinatorial_prediction(default_params, dep_jobs = list(config_job), memory = mem),
        additive_model = additive_model_combinatorial_prediction(default_params, dep_jobs = list(config_job)),
        ground_truth = ground_truth_combinatorial_prediction(default_params, dep_jobs = list(config_job), memory = mem),
        scfoundation = scfoundation_combinatorial_prediction(c(list(epochs = 5),default_params),, dep_jobs = list(config_job), memory = mem)
      )
      pert_jobs
    })
    names(inner_jobs) <- seeds
    purrr::list_flatten(inner_jobs, name_spec = "{outer}-{inner}")
  })
  names(jobs) <- datasets
  jobs <- purrr::list_flatten(jobs, name_spec = "{outer}-{inner}")
  
  params <- list(
    job_ids = map_chr(jobs, "result_id"),
    names = names(jobs)
  )
  collect_perturbation_predictions(params, dep_jobs = jobs)
}

# Launch double perturbation jobs
double_pert_jobs <- make_double_perturbation_jobs(datasets = c("norman_from_scfoundation"), seeds = 1:5)
stat <- map_chr(double_pert_jobs$dependencies, job_status); table(stat)
write_rds(double_pert_jobs, "tmp/double_perturbation_jobs.RDS")
run_job(double_pert_jobs, priority = "normal")
file.copy(file.path(result_file_path(double_pert_jobs), "predictions.RDS"), to = "output/double_perturbation_results_predictions.RDS", overwrite = TRUE)
file.copy(file.path(result_file_path(double_pert_jobs), "parameters.RDS"), to = "output/double_perturbation_results_parameters.RDS", overwrite = TRUE)

# Launch single perturbation jobs
single_pert_jobs <- make_single_perturbation_jobs(datasets = c("adamson", 'replogle_k562_essential', 'replogle_rpe1_essential'), seeds = 1:2)
stat <- map_chr(single_pert_jobs$dependencies, job_status); table(stat)
write_rds(single_pert_jobs, "tmp/single_perturbation_jobs.RDS")
run_job(single_pert_jobs, priority = "normal")
file.copy(file.path(result_file_path(single_pert_jobs), "predictions.RDS"), to = "output/single_perturbation_results_predictions.RDS", overwrite = TRUE)
file.copy(file.path(result_file_path(single_pert_jobs), "parameters.RDS"), to = "output/single_perturbation_results_parameters.RDS", overwrite = TRUE)

stats_df <- bind_rows(
  enframe(single_pert_jobs$dependencies, name = "name", value = "job"),
  enframe(double_pert_jobs$dependencies, name = "name", value = "job"),
) %>%
  mutate(stats = map(job, \(j) read_delim(stats_file_path(j), col_names = c("metric", "value"), delim = " ", col_types = "cd"))) %>%
  unnest(stats) %>%
  dplyr::select(-job)
write_tsv(stats_df, "output/single_perturbation_jobs_stats.tsv")




