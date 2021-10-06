library(tidyverse)
library(cellsig)
library(tidybulk)
library(tidySummarizedExperiment)

dir(sprintf("%s/dev/modeling_results/", local_dir), pattern = "result.rds", full.names = T) %>%
  grep("result.rds", ., value = T, fixed = TRUE) %>%
  map_dfr(~ readRDS(.x)) %>%
  
  unite( "sample", c(cell_type, level), remove = FALSE) %>%
  tree_and_signatures_to_database(tree, ., sample, cell_type, .feature,  `50%`)  %>%

  # Convert to SE
  as_SummarizedExperiment(sample, .feature, c(`2.5%`, `25%`,  `50%`,  `75%`, `97.5%`) ) %>%
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts
  impute_missing_abundance(~ level_4, .abundance = c(`2.5%`, `25%`,  `50%`,  `75%`, `97.5%`)) %>%
  impute_missing_abundance(~ level_3, .abundance = c(`2.5%`, `25%`,  `50%`,  `75%`, `97.5%`)) %>%
  impute_missing_abundance(~ level_2, .abundance = c(`2.5%`, `25%`,  `50%`,  `75%`, `97.5%`)) %>%
  impute_missing_abundance(~ level_1, .abundance = c(`2.5%`, `25%`,  `50%`,  `75%`, `97.5%`)) %>% 
  
  # Convert back to tibble
  as_tibble() %>%
  
  # Merge the imputed column
  mutate(.imputed = .imputed | .imputed.1 | .imputed.2 | .imputed.3  ) %>%
  select(-c(  .imputed.1 , .imputed.2 ,.imputed.3 )) %>%

  saveRDS("dev/counts_bayes_imputed.rds", compress = "xz")
