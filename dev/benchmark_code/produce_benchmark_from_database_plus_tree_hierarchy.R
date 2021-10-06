source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
is_hierarchy = args[1]
contrast_name = args[2]
rank_name = args[3] 
rank_stat = args[4]
bayes = args[5] %>% as.integer()
.selection_method = args[6] 
.optimisation_method = args[7] 
output_file = args[8]


.is_hierarchy = is_hierarchy %>% when(
  (.)=="hierarchical" ~ TRUE,
  (.)=="non_hierarchical" ~ FALSE
)

.contrast_method = contrast_name %>% when(
  (.)=="mean_contrast" ~ mean_contrast,
  (.)=="pairwise_contrast" ~ pairwise_contrast
)

.ranking_method = rank_name %>% when(
  (.)=="edgR" ~ rank_edgR_quasi_likelihood,
  (.)=="edgR_robust" ~ rank_edgR_robust_likelihood_ratio,
  (.)=="bayes" ~ rank_bayes
)

.rank_stat = rank_stat %>% when(
  (.)=="_" ~ NULL,
  ~ (.)
)

.bayes  =
  if(bayes==1){
    readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/cellsig_theoretical_transcript_abundance_distribution.rds") %>% 
      # take the first row of the duplicated data so that each gene in a cell type has only one set of quantiles
      nest(data = - c(symbol, cell_type)) %>% 
      mutate(data = map(data, ~.x[1, ])) %>% 
      # ensure all cell types have the same set of genes
      add_count(symbol) %>%
      filter(n == max(n)) %>%
      unnest(data)
    } else {NULL}


readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_hierarchy.rds") %>% 
  main(.is_hierarchy=.is_hierarchy,
       .contrast_method = .contrast_method, .ranking_method = .ranking_method, .rank_stat = .rank_stat, .bayes = .bayes,
       .selection_method = .selection_method, .kmax = 60, .discard_number = 2000, .reduction_method = "PCA",
       .optimisation_method= .optimisation_method) %>% 
  saveRDS(file = glue("{output_file}", compress = "xz"))