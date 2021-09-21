# devtools::install_github("stemangiola/tidybulk@for-jian", force = TRUE)
# compare efficiency of methods on subset tree
library(tidyverse)
library(cellsig)
library(tidybulk)
library(tidySummarizedExperiment)
library(data.tree)
library(tidytree)
library(ape)

adapt_tree <- function(.data, .tree, .subtree=NULL){
  
  # if no .node of interest is given, the original dataframe will be returned
  if (is.null(.subtree)){.data}
  
  # if node of interest is given, a subtree is generated and cell types in data will be modified accordingly
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
  
  adapt_tree(tree, t_helper_tree) %>% 
  
  # Parse into hierarchical dataset
  tree_and_signatures_to_database(t_helper_tree, ., sample, cell_type, symbol, count)  %>%
  
  # Remove redundant samples
  remove_redundancy(sample, symbol, count, correlation_threshold = 0.999, top = 500, method = "correlation") %>%
  droplevels() %>% 
  
  # Eliminate suspicious samples
  filter(!grepl("GSM3722278|GSM3722276|GSM3722277", sample)) %>%
  
  # eliminate genes that are not in all cell types level 1
  nest(data = -c(level_1, symbol)) %>%
  add_count( symbol) %>%
  filter(n==4) %>%
  select(-n) %>%
  unnest(data) %>%
  
  # Convert to SE
  as_SummarizedExperiment(sample, symbol, count) %>% 
  
  # Scale with first degree imputation. 
  # This because there are no common genes to all samples
  impute_missing_abundance(~ cell_type, suffix="") %>% 
  identify_abundant() %>%
  scale_abundance() %>%
  filter(!.imputed) %>% 
  
  select(-.imputed) %>% 
  
  # Just needed for the old version
  # select(-one_of("exposure_rate")) %>%
  
  # Calculate exposure for Bayes model
  mutate(exposure_rate = -log(multiplier)) %>% 
  
  # Save
  saveRDS("dev/counts_t_helper_tree.rds", compress = "xz")


# Perform imputation
readRDS("dev/counts_t_helper_tree.rds") %>% 

  # Convert to SE
  as_SummarizedExperiment(sample, feature, count_scaled) %>%
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts
  impute_missing_abundance(~ cell_type, suffix="") %>%
  impute_missing_abundance(~ level_5, suffix="") %>%
  impute_missing_abundance(~ level_4, suffix="") %>%
  impute_missing_abundance(~ level_3, suffix="") %>%
  impute_missing_abundance(~ level_2, suffix="") %>%
  impute_missing_abundance(~ level_1, suffix="") %>% 
  
  # Convert back to tibble
  as_tibble() %>% 
  
  mutate(.imputed = if_any(contains("imputed"), ~ .x != 0)) %>% 
  
  # Merge the imputed column
  # mutate(.imputed = .imputed | .imputed.1 | .imputed.2 |.imputed.3 |.imputed.4| .imputed.5 )
  # mutate(.imputed = .imputed.x | .imputed.y | .imputed.x.x |.imputed.y.y |.imputed.x.x.x| .imputed.y.y.y ) %>%
  
  select(-matches("imputed\\.\\d")) %>%
  
  # Save
  saveRDS("dev/intermediate_data/counts_imputed_t_helper_tree.rds", compress = "xz")

# create slurm workflow file

input_directory <- "dev/benchmark_code/"
output_directory <- "dev/benchmark_results_t_helper/"
tab <- "\t"


tibble(is_hierarchy = c("hierarchical", "non_hierarchical"),
       contrast = c("mean_contrast", "mean_contrast"),
       ranking_method = c("bayes", "edgR"),
       ranking_stat = c("_", "PValue"),
       bayes = c(1, 0),
       selection = c("silhouette", "silhouette"),
       optimisation = c("curvature", "penalty")
       ) %>% 
  slice(1) %>% 
  unite("arguments", everything(), sep=" ", remove=FALSE) %>% 
  unite("output_file", -c(arguments, bayes), remove=FALSE) %>% 
  mutate(output_file = glue("{output_directory}{output_file}.rds")) %>% 
  mutate(R_script = glue("{input_directory}produce_benchmark_t_helper.R") ) %>% 

  mutate(command = glue("{output_file}:\n{tab}Rscript {R_script} {arguments} {output_file}")) %>%
  
  pull(command) %>% 
  purrr::prepend("CATEGORY=yes_no_hierarchy\nMEMORY=30000\nCORES=2\nWALL_TIME=86400") %>% 
  # 
  # mutate(SLURM_command = glue::glue("sbatch Rscript {R_command}")) %>% 
  # pull(SLURM_command) %>%
  write_lines("./dev/benchmark_code/benchmark_t_helper.makeflow")

# use the signature output for evaluation
hierarchical_mean_contrast_bayes___silhouette_curvature <- 
  readRDS("dev/benchmark_results/hierarchical_mean_contrast_bayes___silhouette_curvature.rds") %>% 
  mutate(stream = "hierarchical_mean_contrast_bayes_silhouette_curvature", .before=1)

hierarchical_mean_contrast_bayes___silhouette_curvature_t_helper <- 
  readRDS("dev/benchmark_results_t_helper/hierarchical_mean_contrast_bayes___silhouette_curvature.rds") %>% 
  mutate(stream = "hierarchical_mean_contrast_bayes_silhouette_curvature_t_helper", .before=1)

cibersortx <- 
  read_delim("dev/jian_R_files/cibersortx/CIBERSORTx_Job21_phenoclass_1.CIBERSORTx_Job21_reference_1.bm.K999.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  pull(NAME) %>% 
  list() %>% 
  tibble(stream = "cibersortx", signature = .)

x <- 
  bind_rows(
    hierarchical_mean_contrast_bayes___silhouette_curvature_t_helper,
    hierarchical_mean_contrast_bayes___silhouette_curvature
    ) %>% 
  nest(signature = - stream) %>% 
  mutate(signature = map(signature, ~ .x %>% pull(signature) %>% unlist %>% unique)) %>% 
  bind_rows(cibersortx)
  
rm(hierarchical_mean_contrast_bayes___silhouette_curvature,
   hierarchical_mean_contrast_bayes___silhouette_curvature_t_helper,
   cibersortx)

# Evaluation
counts_imputed_non_hierarchy <- 
  readRDS("dev/intermediate_data/counts_imputed_non_hierarchy.rds")

y <- x %>% 
  
  # silhouette evaluation
  mutate(silhouette = map(
    signature, 
    ~ silhouette_evaluation(
      .signature = .x,
      .reduction_method = "PCA",
      .non_hierarchical_counts = counts_imputed_non_hierarchy,
      .sample = sample,
      .symbol = symbol)
  )) %>% 

  mutate(avg_silhouette = map_dbl(silhouette, ~ mean(.x$sil_width))) %>%
  
  mutate(silhouette = map(
    silhouette, 
    ~ .x %>% 
      group_by(cell_type) %>% 
      summarise(cluster_silhouette = mean(sil_width), cluster_size = n()) %>% 
      distinct() %>% 
      ungroup()
  )) 

y %>% 
  
  # deconvolution evaluation
  # for each mixture, combine with the signatures from all methods
  expand_grid(mix100, .) %>% 
  
  mutate(deconvolution = map2(
    signature, mix, 
    ~ deconvolution_evaluation(
      .signature = .x, 
      .mix=.y, 
      .non_hierarchical_counts = counts_imputed_non_hierarchy,
      .sample = sample,
      .symbol = symbol)
  )) %>% 
  
  # mse by method
  mutate(MSE = map_dbl(
    deconvolution,
    ~ mean((.x$estimated_proportion - .x$proportion)^2)
  )) %>% 
  nest(data=-stream) %>% 
  # mutate(median_MSE_over_mixes = map_dbl(data, ~ median(.x$MSE))) %>%
  mutate(mean_MSE_over_mixes = map_dbl(data, ~ mean(.x$MSE))) %>%
  unnest(data) %>% 
  
  # mse by cell type
  unnest(deconvolution) %>% 
  mutate(squared_error = (estimated_proportion - proportion)^2) %>% 
  nest(data = -c(stream, cell_type)) %>% 
  mutate(mean_MSE_for_cell_type = map_dbl(data, ~ mean(.x$squared_error))) %>% 
  unnest(data) %>% 
  
  select(-c(signature, mixture_ID, mix, replicate, estimated_proportion, proportion, squared_error))
  
  
y %>% 
  unnest(silhouette) %>% 
  filter(cell_type %in% c("t_helper_h1", "t_helper_h2"," t_helper_h17", "t_reg")) %>% 
  group_by(cell_type) %>% 
  arrange(desc(cluster_silhouette), .by_group = TRUE) %>% 
  ungroup




