library(tidyverse)
library(magrittr)
library(tidybulk)
library(tidySummarizedExperiment)

readRDS("dev/counts.rds") %>% 
  
  
  
  # Convert to SE
  as_SummarizedExperiment(.sample, .feature, count) %>%
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts
  impute_missing_abundance(~ cell_type, suffix="") %>%
  impute_missing_abundance(~ level_5, suffix="") %>%
  impute_missing_abundance(~ level_4, suffix="") %>%
  impute_missing_abundance(~ level_3, suffix="") %>%
  impute_missing_abundance(~ level_2, suffix="") %>%
  impute_missing_abundance(~ level_1, suffix="") %>% 
  
  # Convert back to tibble
  as_tibble() %>%
  
  # Merge the imputed column
  mutate(.imputed = .imputed.x | .imputed.y | .imputed.x.x |.imputed.y.y |.imputed.x.x.x| .imputed.y.y.y ) %>%
  select(-c( .imputed.x , .imputed.y , .imputed.x.x ,.imputed.y.y ,.imputed.x.x.x, .imputed.y.y.y)) %>%
  
  # Save
  saveRDS("dev/counts_imputed.rds", compress = "xz")
