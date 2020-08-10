# devtools::install_github("stemangiola/nanny@convert-to-S3")
# devtools::install_github("stemangiola/tidybulk@dev")

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

load("data/counts.rda")

tt <- 
  
  # Load dataset
  counts %>%
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


counts_endo_epi_naive <-
  tt %>%

  # Investigate one level
  filter(level==1) %>%

  # investigate one cell type pair
  mutate(comparison_data = map(
    data,
    ~ .x %>%
      filter(cell_type %in% c("endothelial", "epithelial")) %>%
      mutate(cell_type = as.character(cell_type) )
    )) %>%

  #test. We run on the two populations but we select data for all populations
  mutate(markers = map(
    comparison_data,
    ~ .x %>%
      
      # Differential transcription
      test_differential_abundance(
        ~ 0 + cell_type, 
        .contrasts = c(
          "cell_typeendothelial - cell_typeepithelial", 
          "cell_typeepithelial - cell_typeendothelial"
        ),
        action="only"
      ) 
    
  )) %>%
  
  # Select rank from each contrast
  mutate(markers = map(
    markers, ~ .x %>%
      
      # Group by contrast. Comparisons both ways.
      # This results in for example
      #
      # A tibble: 1 x 2
      # `cell_typeendothelial - cell_typeepithelial` `cell_typeepithelial - cell_typeendothelial`
      # <list>                                       <list>                        
      # 
      pivot_longer(
        cols = contains("___"),
        names_to = c("stats", "contrast"), 
        values_to = ".value", 
        names_sep="___"
      ) %>% 
      
      # Markers selection
      nest(stat_df = -contrast) %>%
      
      # Reshape inside each contrast
      mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>%
      
      # Rank
      mutate(stat_df = map(stat_df, ~.x %>%
                             filter(FDR < 0.05 & abs(logFC) > 2) %>%
                             filter(logCPM > mean(logCPM)) %>%
                             arrange(logFC %>% desc()) %>%
                             slice(1:10)        
      )) %>%
      unnest(stat_df)
  )) %>%

  # Add marker info to original data
  mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
  select(-comparison_data, - data) %>%
  unnest(markers) 


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
    filter(FDR < 0.05 & logFC > 2) %>%
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

# permutation
counts_epi_imm_de <- cell_sig(tt, c("epithelial", "immune_cell"))
counts_endo_imm_de <- cell_sig(tt, c("endothelial", "immune_cell"))
counts_fib_imm_de <- cell_sig(tt, c("fibroblast", "immune_cell"))
counts_endo_epi_de <- cell_sig(tt, c("endothelial", "epithelial"))
counts_fib_epi_de <- cell_sig(tt, c("fibroblast", "epithelial"))
counts_fib_endo_de <- cell_sig(tt, c("fibroblast", "endothelial"))

# collage all the differentially expressed genes between cell types
sig_level1 <- bind_rows(counts_imm_epi_de, counts_imm_endo_de, counts_imm_fib_de, 
               counts_epi_endo_de, counts_epi_fib_de, counts_endo_fib_de,
               counts_epi_imm_de, counts_endo_imm_de, counts_fib_imm_de,
               counts_endo_epi_de, counts_fib_epi_de, counts_fib_endo_de) %>%
  tidybulk(sample, symbol, count_scaled) %>% 
  select(sample, symbol, count_scaled, cell_type) %>%
  
  # remove duplicate genes that arise during pairwise comparison and reduce dimensions
  distinct() %>%
  # group_by(cell_type) %>% 
  reduce_dimensions(sample, symbol, count_scaled,  method = "PCA")



# plot PCA 
sig_level1 %>% 
  pivot_sample(sample) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + 
  geom_point() +
  theme_bw()


# Compare marker gene counts distribution
sig_level1 %>% 
  ggplot(aes(x = cell_type, y = log(count_scaled + 1), colour = cell_type)) +
  geom_boxplot() +
  geom_jitter(aes(alpha=0.01)) +
  facet_wrap(~symbol) +
  theme(axis.text.x = element_text(angle = 45))

# marker gene boxplot between endothelial and epithelial cells
counts_epi_endo_de %>% 
  ggplot(aes(x = symbol, y = log(count_scaled + 1), colour = cell_type)) +
  geom_boxplot() +
  geom_point(aes(alpha = 0.01)) +
  theme(axis.text.x = element_text(angle = 45))
