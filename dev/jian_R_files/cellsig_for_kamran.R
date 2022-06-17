# devtools::install_github("stemangiola/nanny@convert-to-S3", force = TRUE)
# devtools::install_github("stemangiola/tidybulk@dev", force = TRUE)
# devtools::install_github("stemangiola/tidybulk@for-jian", force = TRUE)
# devtools::install_github("stemangiola/cellsig@dev", force = TRUE)

library(yaml)
library(tidytext)
library(data.tree)
library(tidytree)
library(ape)
library(glue)
library(rlang)
library(factoextra)
library(stringr)
library(scales)
library(KernSmooth)
library(splus2R)
library(data.tree)
library(cluster)
library(tidyverse)
library(tidybulk)
library(cellsig)
library(patchwork)
library(tidySummarizedExperiment)
library(treemap)





# to run cell signature
# input dataframe must contain these columns: sample | symbol (i.e. gene) | count | cell type
# tree is the cell differentiation tree

# load the input gene expression database for cell types in comparison
#counts = readRDS("dev/jian_R_files/raw_counts_kam.rds")

# load the cell differentiation tree as a data.tree object
#kamran_tree = read_yaml("dev/jian_R_files/kamran_tree.yaml") %>% 
#  as.Node


# get cell signature
#counts_imputed_kamran %>% 
#  main(.sample = sample, .symbol = symbol, .count = count, .cell_type = cell_type,
#       .is_hierarchy=TRUE,
#       .tree = kamran_tree,
#       .contrast_method = pairwise_contrast, .ranking_method = rank_edgR_quasi_likelihood, .rank_stat="PValue",
#       .selection_method = "silhouette", .reduction_method = "PCA", .dims=4, .discard_number = 1000,
#       .optimisation_method = "penalty",
#       .is_complete = TRUE) %>% 
  
#  saveRDS("dev/jian_R_files/cellsignature_kamran.rds", compress = "xz")

