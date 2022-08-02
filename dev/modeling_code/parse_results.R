library(tidyverse)
library(cellsig)
library(tidybulk)
library(tidySummarizedExperiment)
library(magrittr)
library(furrr)
library(data.tree)
library(yaml)
library(glue)
plan(multisession, workers = 10)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
directory_in = args[1]
tree_in_path = args[2]

counts_bayes = 
  dir(directory_in, pattern = "result.rds", full.names = T) %>%
  future_map_dfr(~ {
    x = readRDS(.x) %>% mutate(file = .x)
    attr(x, "fit") = NULL
    attr(x, "rng") = NULL
    x
  }) 

counts_bayes %>% saveRDS(glue("{directory_in}/counts_bayes.rds"), compress = "xz") 

# If no tree provided, just une one level
tree_in = 
  tree_in_path %>% 
  when(
    is.na(.) ~ from_dataframe_to_one_level_tree(counts_bayes, cell_type),
    ~ read_yaml(.) %>% as.Node
  )


counts_bayes_parsed = 
  counts_bayes %>%
  
  #unite( "sample", c(cell_type, level), remove = FALSE) %>%
  tree_and_signatures_to_database(
    tree_in, 
    ., 
    cell_type, 
    cell_type, 
    .feature,  
    `50%`
  )  

 counts_bayes_parsed %>% saveRDS(glue("{directory_in}/counts_bayes_parsed.rds"), compress = "xz") 

counts_bayes_imputed = 
  counts_bayes_parsed %>%

  # Convert to SE
  as_SummarizedExperiment(cell_type, .feature, c(`10%`, `50%`,  `90%`, log_mean, log_sd) ) %>%
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts
  when("level_4" %in% colnames(.) ~ impute_missing_abundance(.,~ level_4, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd), ~ (.))) %>%
  when("level_3" %in% colnames(.) ~ impute_missing_abundance(.,~ level_3, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd), ~ (.))) %>%
  when("level_2" %in% colnames(.) ~ impute_missing_abundance(.,~ level_2, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd), ~ (.))) %>%
  when("level_1" %in% colnames(.) ~ impute_missing_abundance(.,~ level_1, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd), ~ (.))) %>% 
  
  # Convert back to tibble
  as_tibble() %>%
  
  # Merge the imputed column
  when("level_1" %in% colnames(.) ~ mutate(., .imputed = .imputed), ~ (.)) %>%
  when("level_2" %in% colnames(.) ~ mutate(., .imputed = .imputed | .imputed.1), ~ (.)) %>%
  when("level_3" %in% colnames(.) ~ mutate(., .imputed = .imputed | .imputed.2), ~ (.)) %>%
  when("level_4" %in% colnames(.) ~ mutate(., .imputed = .imputed | .imputed.3), ~ (.)) %>% 
  
  dplyr::select(-contains(  ".imputed." )) 

 counts_bayes_imputed %>% saveRDS(glue("{directory_in}/counts_bayes_parsed_imputed.rds"), compress = "xz") 
