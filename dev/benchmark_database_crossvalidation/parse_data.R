# parse_data

library(here)
library(tidyverse)
library(data.tree)
library(yaml)
library(cellsig)
library(tidybulk)

here("dev/benchmark_database_crossvalidation/training_data_1.rds" ) %>%
  readRDS() %>% 
  tree_and_signatures_to_database(
    read_yaml("dev/tree_kamran.yaml") %>% as.Node, 
    ., 
    sample, 
    cell_type, 
    symbol,  
    count
  ) %>%
  rename(
    .sample = sample,
    .feature = symbol
  ) %>%
  identify_abundant(.sample, .feature, count) %>%
  scale_abundance(.sample, .feature, count) %>%
  select(-count_scaled) %>%
  saveRDS(here("dev/benchmark_database_crossvalidation/training_data_1_parsed.rds")) 
