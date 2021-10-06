# adapt input data cell types to the subset tree of interest
library(tidyverse)
library(cellsig)
library(tidybulk)
library(data.tree)

adapt_tree <- function(.data, .tree, .subtree=NULL){
  
  # if no .subtree of interest is given, the original dataframe will be returned
  if (is.null(.subtree)){.data}
  
  # if subtree of interest is given, cell types in data will be modified accordingly
  else{
    
    # subtree <- tree_subset(.tree = .tree, .node = .node)
    
    .tree %>%
      ToDataFrameTypeColFull(fill=NA) %>% 
      
      # find ancestors of each cell type in the full tree
      unite("ancestors", contains("level"), sep=".", na.rm = TRUE) %>% 
      mutate(ancestors = map2(
        ancestors, cell_type,
        ~.x %>% 
          paste("Tissue", ., sep=".") %>% 
          str_split("\\.") %>% 
          unlist() %>% 
          .[. != .y] %>% 
          # arrange the ancestors from most recent to the furthest
          rev()
      )) %>% 
      
      # get all the cell types in the subtree of interest
      mutate(cell_type_subtree = .subtree$Get("level") %>% names() %>% list()) %>% 
      
      # if the cell type in the full tree is present in the subtree, leave it as it is;
      # otherwise find its most recent ancestor present in the subtree
      mutate(cell_type_to_be = pmap_chr(
        list(cell_type, ancestors, cell_type_subtree),
        ~ if (! ..1 %in% ..3){
          ..2[which(..2 %in% ..3)[1]]
        } else {..1}
      )) %>% 
      select(cell_type, cell_type_to_be) %>% 
      
      left_join(.data, ., by = "cell_type") %>% 
      
      # remove the old cell type
      select(-cell_type) %>% 
      
      # use cell type to be as the new cell type
      rename("cell_type" = "cell_type_to_be")
    
  }
  
}

tree_subset <- function(.tree, .node) {
  
  tree_subset_tibble(.tree = .tree, .node = .node) %>% 
    
    # reverse the row order of tibble to make the root node on top
    # .[rev(1: nrow(.)), ] %>%
    
    update_tree_edge_matrix() %>% 
    
    as.phylo() %>% 
    
    as.Node(replaceUnderscores=FALSE)
}

tree_subset_tibble <- function(.tree, .node){
  
  # Args:
  # .tree: a Node R6 object from data.tree
  # .node: a character specifying the cell type of interest
  
  tree_tbl <- .tree %>% as.phylo %>% as_tibble
  
  if (.node == "Tissue") { return( tree_tbl %>% filter(label == "Tissue") )
    
  } else {
    
    return(sibling(tree_tbl, .node) %>% 
             bind_rows( 
               tree_subset_tibble(.tree, parent(tree_tbl, .node) %>% .$label),
               .)
    )
  }
  
}

update_tree_edge_matrix <- function(.tree_tbl) {
  
  # Args:
  # .tree_tbl: a tree_tibble object from tidytree
  
  tree_phylo <- .tree_tbl %>% as.phylo
  
  .tree_tbl %>% 
    
    mutate(node = map_int(label, 
                          ~ if(.x %in% tree_phylo$tip.label){
                            which(.x == tree_phylo$tip.label)
                          } else {
                            which(.x == tree_phylo$node.label) + Ntip(tree_phylo)
                          }
    )) %>% 
    
    mutate(parent = rep(parent %>% unique() %>% rank(), 
                        times = table(parent)) + Ntip(tree_phylo)
    )
  
}

# Load cell differentiation tree
data("tree")

t_helper_tree <- tree_subset(tree, "t_helper_h1")

# Database #1
counts_first_db_raw = readRDS("dev/raw_data/counts_first_db_raw.rds")
counts_first_db_raw <- counts_first_db_raw %>% 
  select(sample, symbol, count, database, cell_type)

# Database #2
counts_second_db_raw <- readRDS("dev/raw_data/counts_second_db_raw.rds")
counts_second_db_raw <- counts_second_db_raw %>% 
  select(sample, symbol, count, database=dataset, cell_type)

# Database #3
counts_third_db_raw <- readRDS("dev/raw_data/counts_third_db_raw.rds")
counts_third_db_raw <- counts_third_db_raw %>% 
  select(sample, symbol, count, database=dataset, cell_type)

counts =
  
  # Merge dataset
  counts_first_db_raw %>%
  bind_rows(counts_second_db_raw) %>%
  bind_rows(counts_third_db_raw) %>%
  select(sample, cell_type, symbol, count, database)

rm(counts_first_db_raw, counts_second_db_raw, counts_third_db_raw)

counts %>%  
  adapt_tree(tree, t_helper_tree)
  

# this shows joining the tree dataframe won't remove cell types that are not present in the tree
# counts_t_helper_tree %>% 
#   # Add tree info
#   left_join(
#     t_helper_tree %>%
#       ToDataFrameTypeColFull(fill=NA),
#     by = "cell_type"
#   )

# Add tree structure to raw data frame
# tree_and_signatures_to_database = function(tree, signatures, .sample, .cell_type, .symbol, .count){
#   .sample = enquo(.sample)
#   .cell_type = enquo(.cell_type)
#   .symbol = enquo(.symbol)
#   .count = enquo(.count)
#   
#   signatures %>%
#     
#     # Add tree info
#     left_join(
#       tree %>%
#         data.tree::Clone() %>%
#         ToDataFrameTypeColFull(fill=NA) ,
#       by = quo_name(.cell_type)
#     ) %>%
#     filter(level_1 %>% is.na %>% `!`) %>%
#     
#     # Reduce size
#     mutate_if(is.character, as.factor) %>% 
#     droplevels %>% 
#     mutate(!!.count := !!.count %>% as.integer) %>%
#     
#     # Filter only symbol existing
#     filter(!!.symbol %>% is.na %>% `!`) %>%
#     
#     # Aggregate
#     aggregate_duplicates(!!.sample, !!.symbol, !!.count) %>% 
#     select(-one_of("merged_transcripts"))
# }
