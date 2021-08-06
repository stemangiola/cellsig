library(tidyverse)
library(cellsig)
library(tidybulk)
library(data.tree)
data("tree")

tree_and_signatures_to_database = function(tree, signatures, .sample, .cell_type, .symbol, .count){
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  
  signatures %>%
    
    # Add tree info
    left_join(
      tree %>%
        data.tree::Clone() %>%
        ToDataFrameTypeColFull(fill=NA) ,
      by = quo_name(.cell_type)
    ) %>%
    filter(level_1 %>% is.na %>% `!`) %>%
    
    # Reduce size
    mutate_if(is.character, as.factor) %>% 
    droplevels %>% 
    mutate(!!.count := !!.count %>% as.integer) %>%
    
    # Filter only symbol existing
    filter(!!.symbol %>% is.na %>% `!`) %>%
    
    # Aggregate
    aggregate_duplicates(!!.sample, !!.symbol, !!.count) %>% 
    select(-one_of("merged_transcripts"))
}

# Create a hierarchical signature data frame from a tree and and signature database
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