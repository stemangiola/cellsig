library(R.utils)
library(tidyverse)
library(tidybulk)
library(ggplot2)
library(plotly)
library(stringr)
library(cluster)
library(factoextra)

load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")

tt_simple <- readRDS("dev/intermediate_data/tt_simple.rds")

counts_imputed <- 
  readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_imputed.rds") %>% 
  rename(symbol = feature)

new_tree <- read_yaml("dev/jian_R_files/new_tree.yaml") %>% as.Node

# Functions ======================================================

produce_cibersortx_bulk_rnaseq_input <- function(.expression_df, .transcript, .sample, .cell_type, .count, 
                                                 .dir, .tree=NULL, .suffix=NULL){
  
  # Args: 
  # .expression_df is a tibble with shape: transcript | sample | cell_type | count
  # .transcript: tell the function which column in the data represents gene symbol
  # .sample: tell the function which column in the data represents sample
  # .cell_type: tell the function which column in the data represents cell_type
  # .count: tell the function which column in the data represents count
  # .tree: optional argument. if .tree = NULL, cell types will be sourced from those present in the .expression_df column
  #       if provided, the cell types used for cibersortx will be the leaves of the tree
  # .dir: the path to the directory where the two output files should be saved
  # .suffix: optional argument which adds suffix to the output filename: reference.txt and phenoclass.txt
  
  .transcript = enquo(.transcript)
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)
  
  .dir <- ifelse(grepl("\\/$", .dir), .dir, paste0(.dir, "/"))
  
  ## ref_names is produced to create header for reference file
  ref_names <- .expression_df %>% 
    {
      if (is.null(.tree)){
        (.)
      } else {
        (.) %>% 
          filter(!!.cell_type %in% as.phylo(.tree)$tip.label)
      }
      
    } %>% 
    mutate(!!.sample := str_replace_all(!!.sample, "\\.", "_")) %>%
    unite(cell_sample, c(!!.cell_type, !!.sample), sep = ".") %>% 
    pull(cell_sample) %>% 
    unique()
  
  
  # create reference file
  .expression_df %>% 
    
    # filter for leaf cell types in the tree
    {
      if (is.null(.tree)){
        (.)
      } else {
        (.) %>% 
          filter(!!.cell_type %in% as.phylo(.tree)$tip.label)
      }
      
    } %>% 
    select(!!.sample, !!.cell_type, !!.count, !!.transcript) %>% 
    mutate(!!.sample := str_replace_all(!!.sample, "\\.", "_")) %>%
    pivot_wider(names_from = c(!!.cell_type, !!.sample), values_from = !!.count, names_sep=".") %>% 
    `names<-`(ref_names %>% 
                str_extract(".*(?=\\.)") %>% 
                prepend(quo_name(.transcript))
    ) %>% 
    # save files as txt files (tab separated files)
    write_tsv(glue("{.dir}reference{.suffix}.txt"))
  
  # create phenoclass file
  tibble(cell_type = 
           {
             if (is.null(.tree)){
               .expression_df %>% 
                 distinct(!!.cell_type) %>% 
                 pull %>% 
                 as.character
             } else {
               as.phylo(.tree)$tip.label
             }
             
           }
         
  ) %>% 
    bind_cols(
      tibble(ref_names, value=1L) %>% 
        pivot_wider(names_from = ref_names, values_from = value)
    ) %>% 
    pivot_longer(-cell_type, names_to="ref_names", values_to="membership") %>% 
    
    # assign membership according to cibersortx tutorial 6: 
    # "1" indicates membership of the reference sample to the class as defined in that row,
    # "2" indicates the class that the sample will be compared against
    mutate(membership = if_else(str_detect(ref_names, cell_type), 1L, 2L)) %>% 
    
    pivot_wider(names_from = ref_names, values_from = membership) %>% 
    
    # set col_names to FALSE following cibersortx format for phenoclass
    write_tsv(glue("{.dir}phenoclass{.suffix}.txt"), col_names = FALSE)
  
}

# No hierarchy CIBERSORTx files creation=============================

## create reference file

reference <- counts_imputed %>% 
  filter(cell_type %in% as.phylo(new_tree)$tip.label) %>% 
  select(sample, cell_type, count_scaled, symbol) %>% 
  mutate(sample = str_replace_all(sample, "\\.", "_")) %>%
  pivot_wider(names_from = c(cell_type, sample), values_from = count_scaled, names_sep=".")

## vector ref_names is produced below yet used here to create header for reference file
ref_names <- counts_imputed %>% 
  filter(cell_type %in% as.phylo(new_tree)$tip.label) %>% 
  mutate(sample = str_replace_all(sample, "\\.", "_")) %>%
  unite(cell_sample, c(cell_type, sample), sep = ".") %>% 
  pull(cell_sample) %>% 
  unique()

header <- ref_names %>% 
  str_extract(".*(?=\\.)") %>% 
  insert(ats=1, "symbol")

names(reference) <- header

reference

## create phenoclass file
cell_types <- as.phylo(new_tree)$tip.label

ref_tibble <- tibble(ref_names, value=1:length(ref_names)) %>% 
  pivot_wider(names_from = ref_names, values_from = value)

phenoclass <- tibble(cell_types) %>% 
  bind_cols(ref_tibble) %>% 
  pivot_longer(-cell_types, names_to="ref_names", values_to="membership") %>% 
  mutate(membership = if_else(str_detect(ref_names, cell_types), 1, 2)) %>% 
  pivot_wider(names_from = ref_names, values_from = membership)


# save files as txt files (tab separated files)
write_tsv(reference, "./dev/jian_R_files/cibersortx/reference.txt")

write_tsv(phenoclass, "./dev/jian_R_files/cibersortx/phenoclass.txt", col_names = FALSE)

# Analysis on the signature matrix generated by CIBERSORTx========================
CIBERSORTx_Job20_signature_bm_K999 <- 
  read_delim("cibersortx_NH/CIBERSORTx_Job20_signature.bm.K999.txt", "\t", 
             escape_double = FALSE, trim_ws = TRUE)

cibersort_signature <- CIBERSORTx_Job20_signature_bm_K999 %>% 
  pull(NAME)


CIBERSORTx_Job21_phenoclass_1_CIBERSORTx_Job21_reference_1_bm_K999 <- 
  read_delim("dev/jian_R_files/cibersortx/CIBERSORTx_Job21_phenoclass_1.CIBERSORTx_Job21_reference_1.bm.K999.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE)

cibersort_signature <- CIBERSORTx_Job21_phenoclass_1_CIBERSORTx_Job21_reference_1_bm_K999 %>% 
  pull(NAME)

# OBSOLETE =============
cibersort_silhouette <- tt_non_hierarchy %>% 
  unnest(tt) %>% 
  unnest(data) %>% 
  filter(symbol %in% cibersort_signature) %>% 
  sil_func("root", METHOD) %>% 
  select(signature, reduced_dimensions, silhouette, real_size) %>% 
  mutate(method = "cibersortx", .before = signature)

cibersort_silhouette

ciber_pca <- ciber_sil %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  # stat_ellipse(type = 't') +
  ggtitle("ciber_pca") +
  theme(
    plot.title = element_text(hjust=0.5) )

ciber_pca

ggsave("ciber_pca.png", ciber_pca)

