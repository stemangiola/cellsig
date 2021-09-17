source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")

preprocess <- function(.transcript, .level) {
  
  # load data
  .transcript %>%
    
    dplyr::rename("symbol" = ".feature", sample = ".sample") %>% 
    # tidybulk(sample, symbol, count_scaled) %>% for imputed counts data
    # tidybulk(.sample = sample, .transcript = symbol, .abundance = count, .abundance_scaled = counts_scaled) %>%
    
    # filter for cells at the level of interest. .level == level_1
    filter(is.na(!!as.symbol(.level))==FALSE) %>%
    
    nest(data = -c(symbol, !!as.symbol(.level))) %>%
    
    # for a cell type some samples may miss genes in other samples: so for the same cell type genes may have different number of samples
    mutate(n_samples_per_gene = map_int(
      data,
      ~ .x %>%
        distinct(sample) %>%
        nrow)) %>%
    
    unnest(data) %>%
    nest(data = -!!as.symbol(.level)) %>%
    
    mutate(n_samples= map_int(
      data,
      ~ .x %>%
        distinct(sample) %>%
        nrow
    )) %>%
    
    unnest(data) %>% 
    
    mutate(ratio_imputed_samples = 1 - n_samples_per_gene / n_samples) %>% 
    
    # filter(ratio_imputed_samples < 0.2) %>% 
  
    # nest by ancestor
    nest(data = - !!as.symbol(pre(.level)))
  
}

counts_imputed_hierarchy <- scale_input_counts(counts_imputed, .is_hierarchy = TRUE)

counts_bayes_imputed_hierarchy <- 
  counts_bayes_imputed %>% 
  select(-level) %>% 
  scale_input_counts(.is_hierarchy = TRUE)

ranked_PW_L4 <- counts_imputed_L4 %>% 
  do_ranking(.ranking_method = rank_edgR_quasi_likelihood, 
             .contrast_method = pairwise_contrast,
             .rank_stat = "logFC")

naive_PW_L4 <- ranked_PW_L4 %>% 
  naive_selection(1)

# preprocessing step
# -> calculate the frequency of imputation for each celltype/gene


# Ranking step, e.g. epithelial vs all
# rank epithelial
# left_join frequency_of_imputation based on target cell type (epithelial) and gene
# Filter(frequency_of_imputation < 0.2)
# proceed to optimisation

# x <- counts_imputed %>% 
#   
#   dplyr::rename("symbol" = ".feature", sample = ".sample") %>% 
#   
#   # tidybulk(sample, symbol, count_scaled) %>% for imputed counts data
#   tidybulk(.sample = sample, .transcript = symbol, .abundance = count, .abundance_scaled = counts_scaled) %>%
#   
#   # filter for cells at the level of interest. .level == level_1
#   filter(is.na(level_3)==FALSE) %>% 
#   
#   # Level 3
#   # Do it just ones per level
#   nest(data = -c(symbol, level_3)) %>%
#   
#   # for a cell type some samples may miss genes in other samples: so for the same cell type genes may have different number of samples
#   mutate(n_samples_per_gene = map_int(
#     data,
#     ~ .x %>%
#       distinct(sample) %>%
#       nrow)) %>%
#   
#   unnest(data) %>%
#   nest(data = -level_3) %>%
#   
#   mutate(n_samples= map_int(
#     data,
#     ~ .x %>%
#       distinct(sample) %>%
#       nrow
#   )) %>%
#   
#   unnest(data) 
#   
# x %>% 
#   
#   mutate(ratio_imputed_samples = 1 - n_samples_per_gene / n_samples) %>% 
#   
#   filter(ratio_imputed_samples < 0.2)%>% nrow()
# 
# # DO manu times for each cell type in your benchmark. Think about where to add this filter
# nk_primed %>%
#  filter(ratio_imputed < .20)





