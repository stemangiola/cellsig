source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")
args <- commandArgs(trailingOnly = TRUE)

is_hierarchy = args[1]
contrast = args[2]
ranking_method = args[3] 
ranking_stat = args[4]
bayes = args[5] %>% as.integer()
selection = args[6] 
optimisation = args[7] 
output_file = args[8]

.is_hierarchy = is_hierarchy %>% when(
  (.)=="hierarchical" ~ TRUE,
  (.)=="non_hierarchical" ~ FALSE
)

.contrast_method = contrast %>% when(
  (.)=="mean_contrast" ~ mean_contrast,
  (.)=="pairwise_contrast" ~ pairwise_contrast
)

.ranking_method = ranking_method %>% when(
  (.)=="edgR" ~ rank_edgR_quasi_likelihood,
  (.)=="bayes" ~ rank_bayes
)

.rank_stat = ranking_stat %>% when(
  (.)=="_" ~ NULL,
  ~ (.)
)

.bayes  =
  if(bayes==1L){
    readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/benchmark_results_t_helper/counts_bayes_imputed_t_helper_tree.rds")
  } else {NULL}

readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/benchmark_results_t_helper/counts_imputed_t_helper_tree.rds") %>%
  dplyr::rename(symbol = feature) %>% 
  main(.sample = sample, 
       .symbol = symbol,
       .count = NULL,
       .cell_type = cell_type,
       .is_hierarchy=.is_hierarchy,
       .contrast_method=.contrast_method, 
       .ranking_method=.ranking_method, 
       .rank_stat=.rank_stat, 
       .bayes=.bayes, 
       .selection_method=selection, .kmax=60, .discard_number=2000, .reduction_method = "PCA", .dims=2,
       .optimisation_method=optimisation, .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
       .is_complete = TRUE) %>% 
  saveRDS(file = output_file, compress = "xz")
