# devtools::github_install("stemangiola/nanny@convert-to-S3")
# devtools::github_install("stemangiola/tidybulk@dev")

library(tidyverse)
library(plotly)
library(nanny)
library(ggrepel)
library(GGally)
library(tidyHeatmap)
library(furrr)
plan(multiprocess, workers=5)

# To be loaded after all libraries
library(tidybulk)


tt <- 
  
  # Load dataset
  cellsig::counts %>%
  tidybulk(sample, symbol, count) %>%

  # Group by level because otherwise samples are duplicated
  nest(data = -level) %>%
  
  # Redefine factors inside each level
  mutate(data = future_map(data, ~ droplevels(.x))) %>%
  
  # Fill missing data. There are many genes that
  # are not shared by the majority of samples
  mutate(data = future_map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
  
  # Scale for future PCA plotting
  mutate(data = future_map(data, ~ scale_abundance(.x)))


tt_naive <-  
  tt %>%

  # Scale and reduce dimensions
  mutate(data = future_map(
    data,
    ~ .x %>%
      scale_abundance() %>%
      reduce_dimensions(method="MDS") %>%
      reduce_dimensions(method="PCA") %>%
      reduce_dimensions(method="tSNE")
  )) %>%

  # Cluster
  mutate(data = future_map(data,~ cluster_elements(.x, method="SNN")))


counts_imm_epi_naive <-
  tt %>%

  # Investigate one level
  filter(level==1) %>%

  # investigate one cell type pair
  mutate(comparison_data = map(
    data,
    ~ .x %>%
      filter(cell_type %in% c("immune_cell", "epithelial")) %>%
      mutate(cell_type = as.character(cell_type) )
    )) %>%

  #test. We run on the two populations but we select data for all populations
  mutate(markers = map(
    comparison_data,
    ~ .x %>%
      test_differential_abundance( ~ cell_type, action="only")
    )) %>%

  # Add marker info to original data
  mutate(data = map2(data, markers, ~ left_join(.x, .y))) %>%
  select(-comparison_data, - markers) %>%
  unnest(data) %>%

  # Nest
  nest_subset(data = -symbol) %>%

  # Select markers
  filter(FDR < 0.05 & abs(logFC) > 2) %>%
  filter(logCPM > mean(logCPM)) %>%
  arrange(logFC %>% desc()) %>%
  slice(1:10) %>%

  unnest(data)

######################################################################################################

# cell_sig() extracts top10 differentially expressed genes from different cell type pairs

cell_sig <- function(input, pair) {
  
  input %>%
    
    # Investigate one level
    filter(level==1) %>%
    
    # investigate one cell type pair
    mutate(comparison_data = map(
      data,
      ~ .x %>%
        filter(cell_type %in% pair) %>%
        mutate(cell_type = as.character(cell_type) )
    )) %>%
    
    #test. We run on the two populations but we select data for all populations
    mutate(markers = map(
      comparison_data,
      ~ .x %>%
        test_differential_abundance( ~ cell_type, sample, symbol, count, action="only")
    )) %>%
    
    # Add marker info to original data
    mutate(data = map2(data, markers, ~ left_join(.x, .y))) %>%
    select(-comparison_data, - markers) %>%
    unnest(data) %>%
    
    # Nest
    nest_subset(data = -symbol) %>%
    
    # Select markers
    filter(FDR < 0.05 & abs(logFC) > 2) %>%
    filter(logCPM > mean(logCPM)) %>%
    arrange(logFC %>% desc()) %>%
    slice(1:10) %>%
    unnest(data)
}

# extract top10 differentially expressed genes from 6 cell type pairs in level 1
counts_imm_epi_de <- cell_sig(tt, c("immune_cell", "epithelial"))
counts_imm_endo_de <- cell_sig(tt, c("immune_cell", "endothelial"))
counts_imm_fib_de <- cell_sig(tt, c("immune_cell", "fibroblast"))
counts_epi_endo_de <- cell_sig(tt, c("epithelial", "endothelial"))
counts_epi_fib_de <- cell_sig(tt, c("epithelial", "fibroblast"))
counts_endo_fib_de <- cell_sig(tt, c("endothelial", "fibroblast"))

# collage all the differentially expressed genes between cell types
sig_level1 <- bind_rows(counts_imm_epi_de, counts_imm_endo_de, counts_imm_fib_de, 
               counts_epi_endo_de, counts_epi_fib_de, counts_endo_fib_de) %>%
  tidybulk(sample, symbol, count) %>% 
  select(cell_type, sample, symbol, count) %>%
  
  # remove duplicate genes that arise during pairwise comparison and reduce dimensions
  distinct() %>%
#  scale_abundance() %>% 
  reduce_dimensions(sample, symbol, count, method = "PCA")

# plot PCA 
sig_level1 %>% 
  pivot_sample(sample) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + 
  geom_point() +
  theme_bw()

