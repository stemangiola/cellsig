# parse_data

library(here)
library(tidyverse)
library(data.tree)
library(yaml)
library(cellsig)

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
  rename(.feature = symbol) %>%
  saveRDS(here("dev/benchmark_database_crossvalidation/training_data_1_parsed.rds")) 
