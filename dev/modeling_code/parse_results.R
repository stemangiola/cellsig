library(tidyverse)
library(cellsig)
library(tidybulk)
library(tidySummarizedExperiment)
library(magrittr)

# level_4_patition_801_input_result$rng %>% 
#   rstan::summary("Y_gen", c(0.1, 0.9)) %$%
#   summary %>% 
#   as_tibble() %>% 
#   rowid_to_column(var = ".feature_idx") 
# %>% 
#   left_join(
#     df %>% 
#       mutate(database_for_cell_type = factor(database_for_cell_type), .feature = factor(.feature)) %>% 
#       select(database_for_cell_type, cell_type, .feature, count_scaled, count, level_4, multiplier) %>% 
#       mutate(.feature_idx = as.integer(.feature))
#   ) 

library(furrr)
plan(multisession, workers = 10)
local_dir = "."

counts_bayes = 
  dir(sprintf("%s/dev/modeling_results/", local_dir), pattern = "result.rds", full.names = T) %>%
  grep("result.rds", ., value = T, fixed = TRUE) %>%
  future_map_dfr(~ {
    x = readRDS(.x) %>% mutate(file = .x)
    attr(x, "fit") = NULL
    attr(x, "rng") = NULL
  }) %>%
  
  unite( "sample", c(cell_type, level), remove = FALSE) %>%
  tree_and_signatures_to_database(tree, ., sample, cell_type, .feature,  `50%`)  

counts_bayes %>% saveRDS("dev/counts_bayes.rds", compress = "xz")

counts_bayes %>%

  # Convert to SE
  as_SummarizedExperiment(sample, .feature, c(`10%`, `50%`,  `90%`, log_mean, log_sd) ) %>%
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts
  impute_missing_abundance(~ level_4, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd)) %>%
  impute_missing_abundance(~ level_3, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd)) %>%
  impute_missing_abundance(~ level_2, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd)) %>%
  impute_missing_abundance(~ level_1, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd)) %>% 
  
  # Convert back to tibble
  as_tibble() %>%
  
  # Merge the imputed column
  mutate(.imputed = .imputed | .imputed.1 | .imputed.2 | .imputed.3  ) %>%
  select(-c(  .imputed.1 , .imputed.2 ,.imputed.3 )) %>%

  saveRDS("dev/counts_bayes_imputed.rds", compress = "xz")
