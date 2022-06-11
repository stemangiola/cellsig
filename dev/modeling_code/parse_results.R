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

library(data.tree)
library(yaml)


# counts_bayes_OLD <- readRDS("~/counts_bayes.rds")
# counts_bayes_new = counts_bayes %>% bind_rows(counts_bayes_OLD %>% anti_join(counts_bayes, by = c(".feature", "cell_type")))
directory = glue("{local_dir}/dev/benchmark_database_crossvalidation/")
prefix = "kamran_crossvalidation_dataset_1"


counts_bayes = 
  dir(directory, pattern = "result.rds", full.names = T) %>%
  future_map_dfr(~ {
    x = readRDS(.x) %>% mutate(file = .x)
    attr(x, "fit") = NULL
    attr(x, "rng") = NULL
    x
  }) 

job::job({ counts_bayes %>% saveRDS(glue("{directory}/counts_bayes_{prefix}.rds"), compress = "xz") })

counts_bayes_parsed = 
  counts_bayes %>%
  
  #unite( "sample", c(cell_type, level), remove = FALSE) %>%
  tree_and_signatures_to_database(
    read_yaml("dev/tree_kamran.yaml") %>% as.Node, 
    ., 
    cell_type, 
    cell_type, 
    .feature,  
    `50%`
  )  

job::job({ counts_bayes_parsed %>% saveRDS(glue("{directory}/counts_bayes_{prefix}_parsed.rds"), compress = "xz") })

counts_bayes_imputed = 
  counts_bayes_parsed %>%

  # Convert to SE
  as_SummarizedExperiment(cell_type, .feature, c(`10%`, `50%`,  `90%`, log_mean, log_sd) ) %>%
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts
  impute_missing_abundance(~ level_4, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd)) %>%
  impute_missing_abundance(~ level_3, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd)) %>%
  impute_missing_abundance(~ level_2, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd)) %>%
  impute_missing_abundance(~ level_1, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd)) %>% 
  
  # Convert back to tibble
  as_tibble() %>%
  
  # Merge the imputed column
  mutate(.imputed = .imputed | .imputed.1 | .imputed.2 | .imputed.3  ) %>%
  dplyr::select(-c(  .imputed.1 , .imputed.2 ,.imputed.3 )) 

job::job({ counts_bayes_imputed %>% saveRDS(glue("{directory}/counts_bayes_{prefix}_parsed_imputed.rds"), compress = "xz") })
