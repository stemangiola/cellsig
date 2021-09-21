source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")

preprocess <- function(.transcript, .level) {
  
  # load data
  .transcript %>%
    
    dplyr::rename("symbol" = ".feature", sample = ".sample") %>% 
    # tidybulk(sample, symbol, count_scaled) %>% for imputed counts data
    tidybulk(.sample = sample, .transcript = symbol, .abundance = count, .abundance_scaled = counts_scaled) %>%
    
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

readRDS("dev/intermediate_data/counts_imputed.rds") %>% 
  scale_input_counts(.is_hierarchy = TRUE) %>% 
  saveRDS("dev/intermediate_data/counts_imputed_hierarchy.rds", compress = "xz")

readRDS("dev/intermediate_data/counts_imputed.rds") %>% 
  scale_input_counts(.is_hierarchy = FALSE) %>% 
  saveRDS("dev/intermediate_data/counts_imputed_non_hierarchy.rds", compress = "xz")

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

    # nest by ancestor
    nest(data = - !!as.symbol(pre(.level)))
  
}

readRDS("dev/intermediate_data/counts_bayes_imputed.rds") %>% 
  select(-level) %>% 
  scale_input_counts(.is_hierarchy = TRUE) %>% 
  saveRDS("dev/intermediate_data/counts_bayes_imputed_hierarchy.rds", compress = "xz")

readRDS("dev/intermediate_data/counts_bayes_imputed.rds") %>% 
  select(-level) %>% 
  scale_input_counts(.is_hierarchy = FALSE) %>% 
  saveRDS("dev/intermediate_data/counts_bayes_imputed_non_hierarchy.rds", compress = "xz")

# ranked_MC_L4 <- tt_L4 %>%
#   do_ranking(.sample = sample, .symbol=symbol, 
#     .ranking_method = rank_edgR_quasi_likelihood, 
#              .contrast_method = mean_contrast,
#              .rank_stat = "logFC")
# 
# naive_PW_L4 <- ranked_MC_L4 %>%
#   do_selection(.sample = sample, .symbol = symbol, .selection_method = "silhouette",
#                .reduction_method = "PCA", .discard_number = 50, .dims = 2)
  # naive_selection(1)

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





