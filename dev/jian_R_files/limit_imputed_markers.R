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
    
    mutate(ratio_imputed_samples = 1 - n_samples_per_gene / n_samples_in_tot) %>% 
    
    filter(ratio_imputed_samples < 0.2) %>% 
  
    # nest by ancestor
    nest(data = - !!as.symbol(pre(.level)))
  
}



x <- counts_imputed %>% 
  
  dplyr::rename("symbol" = ".feature", sample = ".sample") %>% 
  
  # tidybulk(sample, symbol, count_scaled) %>% for imputed counts data
  tidybulk(.sample = sample, .transcript = symbol, .abundance = count, .abundance_scaled = counts_scaled) %>%
  
  # filter for cells at the level of interest. .level == level_1
  filter(is.na(level_3)==FALSE) %>% 
  
  # Level 3
  # Do it just ones per level
  nest(data = -c(symbol, level_3)) %>%
  
  # for a cell type some samples may miss genes in other samples: so for the same cell type genes may have different number of samples
  mutate(n_samples_per_gene = map_int(
    data,
    ~ .x %>%
      distinct(sample) %>%
      nrow)) %>%
  
  unnest(data) %>%
  nest(data = -level_3) %>%
  
  mutate(n_samples= map_int(
    data,
    ~ .x %>%
      distinct(sample) %>%
      nrow
  )) %>%
  
  unnest(data) 
  
x %>% 
  
  mutate(ratio_imputed_samples = 1 - n_samples_per_gene / n_samples) %>% 
  
  filter(ratio_imputed_samples < 0.2)%>% nrow()

# DO manu times for each cell type in your benchmark. Think about where to add this filter
nk_primed %>%
 filter(ratio_imputed < .20)





