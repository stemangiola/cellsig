library(ape)
library(tidytree)
library(data.tree)
library(tidyverse)

data("tree")

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




# this block of code returns a tibble ready to be joined with raw counts data
# tree_tbl <- tree_subset_tibble(tree, "t_CD4_memory") %>% 
#   filter(label != "Tissue") %>% 
#   mutate(parent = factor(parent, labels = glue("level_{1: n_distinct(parent)}") )) %>% 
#   pivot_wider(names_from = parent, values_from = label) %>% 
#   select(-c(node, branch.length)) %>% 
#   unite("cell_type", remove = FALSE, na.rm = TRUE)

# Modified function do not use!
# tree_and_signatures_to_database = function(tree, .interest, signatures, .sample, .cell_type, .symbol, .count){
#   .sample = enquo(.sample)
#   .cell_type = enquo(.cell_type)
#   .symbol = enquo(.symbol)
#   .count = enquo(.count)
#   
#   signatures %>%
#     
#     # Add tree info
#     left_join(
#       tree_subset(tree, .interest) %>% 
#         filter(label != "Tissue") %>% 
#         mutate(parent = factor(parent, labels = glue("level_{1: n_distinct(parent)}") )) %>% 
#         pivot_wider(names_from = parent, values_from = label) %>% 
#         select(-c(node, branch.length)) %>% 
#         unite("cell_type", remove = FALSE, na.rm = TRUE),
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