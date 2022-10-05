README
================

Create a hierarchical signature data frame from a tree and and signature
database

``` r
library(tidyverse)
library(cellsig)
library(tidybulk)
data("tree")

counts_first_db_raw = readRDS("dev/counts_first_db_raw.rds")
load("dev/counts_second_db_raw.rda")
counts_second_db_raw = new_data %>% select(sample, symbol, count, dataset, cell_type)

signatures = 
  counts_first_db_raw %>%
  select(-cell_type_original) %>%
  bind_rows(counts_second_db_raw %>% rename(database = dataset)) %>% 
  select(sample, cell_type, symbol, count)

counts = 
  tree_and_signatures_to_database(
    tree,
    signatures,
    sample,
    cell_type,
    symbol, 
    count
  )

counts
```
