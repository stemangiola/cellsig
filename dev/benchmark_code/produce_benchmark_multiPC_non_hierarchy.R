source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
is_hierarchy = args[1]
contrast = args[2]
rank = args[3] 
rank_stat = args[4]
selection = args[5] 
optimisation = args[6]
dims = args[7] %>% as.integer
bayes = args[8] %>% as.integer
output_file = args[9]


.is_hierarchy = is_hierarchy %>% when(
  (.)=="hierarchical" ~ TRUE,
  (.)=="non_hierarchical" ~ FALSE
)

.contrast_method = contrast %>% when(
  (.)=="mean_contrast" ~ mean_contrast,
  (.)=="pairwise_contrast" ~ pairwise_contrast
)

.ranking_method = rank %>% when(
  (.)=="edgR" ~ rank_edgR_quasi_likelihood,
  (.)=="edgR_robust" ~ rank_edgR_robust_likelihood_ratio,
  (.)=="bayes" ~ rank_bayes
)

.rank_stat = rank_stat %>% when(
  (.)=="_" ~ NULL,
  ~ (.)
)

.bayes  =
  if(bayes==1L){
    readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_bayes_imputed_non_hierarchy.rds")
  } else {NULL}


readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_imputed_non_hierarchy.rds") %>% 
  
  main(.sample=sample, .symbol=symbol, .count = NULL, .cell_type = cell_type,
       .is_hierarchy=.is_hierarchy, 
       .contrast_method=.contrast_method, 
       .ranking_method=.ranking_method, 
       .rank_stat=.rank_stat, 
       .bayes=.bayes, 
       .tree = NULL,
       .selection_method=selection, .kmax=60, .discard_number=2000, .reduction_method = "PCA",
       .dims=dims,
       .optimisation_method = optimisation, .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
       .is_complete = TRUE) %>% 
  
  saveRDS(file = output_file, compress = "xz")