library(tidyverse)
library(magrittr)
library(tidybulk)
library(tidySummarizedExperiment)

options("tidybulk_do_validate"= FALSE) 

readRDS("dev/counts.rds") %>% 
  
  
  
  # Convert to SE
  as_SummarizedExperiment(.sample, .feature, c(count, count_scaled)) %>%
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts
  impute_missing_abundance(~ cell_type, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_5, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_4, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_3, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_2, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_1, .abundance = c(count, count_scaled)) %>% 
  
  # Convert back to tibble
  as_tibble() %>%
  
  # Merge the imputed column
  mutate(.imputed = .imputed | .imputed.1 | .imputed.2 | .imputed.3 |.imputed.4 |.imputed.5  ) %>%
  select(-c( .imputed.1 , .imputed.2 , .imputed.3 ,.imputed.4 ,.imputed.5  )) %>%
  
  # Save
  saveRDS("dev/counts_imputed.rds", compress = "xz")
