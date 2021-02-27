# Collection of all the signature genes selected upon silhouette score

library(tidyverse)
library(ggplot2)
library(plotly)
library(future)
library(furrr)
library(tidybulk)
library(cluster)
library(proxy)
library(factoextra)
library(stringr)


# Load data
load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")

# Functions================================================

## 1 preprocess data

### 1.1 string manipulation that converts level of interest (e.g "level_5") to its ancestor level (e.g "level_4")
pre <- function(.level) {
  .level %>% 
    str_split("_") %>% 
    {as.numeric(.[[1]][2])-1} %>% 
    paste("level", ., sep = "_")
}

### 1.2 preprocess

preprocess <- function(.data, .level) {
  
  # load data
  .data %>%
    
    tidybulk(sample, symbol, count) %>%
    
    # filter for the cell types of interest for gene marker selection
    filter(is.na(!!as.symbol(.level))==F) %>%
    
    # Imputation of missing data within each level_5
    # impute_missing_abundance(~ !!as.symbol(.level)) %>%
    
    # Group by ancestor
    nest(data = - !!as.symbol(pre(.level))) %>%
    
    # Eliminate genes that are present in some but all cell types
    # (can be still present in a subset of samples from each cell-type)
    mutate(data = map(
      data,
      ~ .x %>%
        nest(data = -c(symbol, !!as.symbol(.level))) %>%
        add_count(symbol) %>%
        filter(n == max(n)) %>%
        unnest(data)
    )) %>%
    
    # Imputation of missing data within each level_5
    mutate(data = map(data, ~ .x %>% impute_missing_abundance(~ !!as.symbol(.level)))) %>%
    
    # scale count for further analysis
    mutate(data=map(data, ~ .x %>%
                      identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
                      scale_abundance()
    ))
}

## 2 contrast functions

### 2.1 pairwise comparisons
get_contrasts_from_df = function(.data, .level){
  
  .data %>% 
    distinct(!!as.symbol(.level)) %>% 
    mutate(!!as.symbol(.level) := paste0(.level, !!as.symbol(.level))) %>% 
    
    # Permute
    mutate(cell_type2 := !!as.symbol(.level)) %>% 
    expand(!!as.symbol(.level), cell_type2) %>% 
    filter(!!as.symbol(.level) != cell_type2) %>% 
    
    # Create contrasts
    mutate(contrast = sprintf("%s - %s", !!as.symbol(.level), cell_type2)) %>%
    pull(contrast)
  
}

### 2.2 create a contrast vector for limma::makeContrasts() or tidybulk::test_differential abundance()

mean_contrast <- function(.data, .level){
  
  # find all cell types
  cell_types <- .data %>% 
    distinct(!!as.symbol(.level)) %>% 
    pull() %>% 
    as.vector()
  
  # format cell_types with prefix
  cell_types <- paste0(.level, cell_types)
  
  # initialise a vector called contrasts
  contrasts <- 1: length(cell_types)
  
  # create all contrasts and store them in contrasts
  for(i in 1:length(cell_types) ){
    background = paste(cell_types[-i], collapse = "+")
    divisor = length(cell_types[-i])
    contrasts[i] <- sprintf("%s-(%s)/%s", cell_types[i], background, divisor)
  }
  
  return(contrasts)
}

## 3 marker ranking & selection

select_markers_for_each_contrast = function(.markers, .sig_size){
  .markers %>%
    
    # Group by contrast. Comparisons both ways.
    pivot_longer(
      cols = contains("___"),
      names_to = c("stats", "contrast"), 
      values_to = ".value", 
      names_sep="___"
    ) %>% 
    
    # Markers selection within each pair of contrast
    nest(stat_df = -contrast) %>%
    
    # Reshape inside each contrast
    mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>%
    
    # Rank
    mutate(stat_df = map(stat_df, ~.x %>%
                           filter(FDR < 0.05 & logFC > 2) %>%
                           filter(logCPM > mean(logCPM)) %>%
                           arrange(logFC %>% desc()) %>%
                           slice(1: .sig_size)
                         
    )) %>%
    
    unnest(stat_df)
}

## 4 marker collection for each contrast

### pairwise contrast method
contrast_PW <- function(.tt, .level){
  .tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + !!as.symbol(.level), 
                                    .contrasts = get_contrasts_from_df(.x, .level),
                                    action="only") 
    ))
}

# mean contrast method
contrast_MC <- function(.tt, .level){
  .tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + !!as.symbol(.level), 
                                    .contrasts = mean_contrast(.x, .level),
                                    action="only") 
    ))
}

### select signature genes and processing
sig_select <- function(.contrast, .level, .sig_size) {
  .contrast %>% 
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, .sig_size))) %>%
    
    # Add original data info to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(!!as.symbol(pre(.level)), markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, .level, "") %>% str_replace(.level, ""))
}

## 2.5 Silhouette function
sil_func0 <- function(.markers, .method){
  .markers %>% 
    nest(rdim = - level_0) %>% 
    mutate(rdim = map(rdim, ~ .x %>% 
                        distinct(sample, symbol, count_scaled, cell_type))) %>% 
    mutate(rdim = map(rdim, ~ .x %>%
                        reduce_dimensions(sample, symbol, count_scaled,
                                          method = .method,
                                          action = "add",
                                          transform = log1p,
                                          # check_duplicates is for Rtsne method
                                          check_duplicates = FALSE) %>% 
                        
                        # save symbols for calculating real_size while reducing replicated rows resulted from symbol
                        nest(data_symbol = c(symbol, count_scaled))
    )) %>%
    
    # calculate the dissimilarity matrix with PC values
    mutate(distance = map(rdim, ~ .x %>%
                            select(contains(str_sub(.method, end = -2L))) %>%
                            factoextra::get_dist(method = "euclidean")
    )) %>%
    
    # calculate silhouette score
    mutate(sil = map2(rdim, distance, 
                      ~ silhouette(as.numeric(as.factor(`$`(.x, cell_type))), .y)
    )) %>% 
    mutate(sil = map(sil, ~ .x %>% summary())) %>%
    mutate(sil = map(sil, ~ .x %>% `$`(avg.width) ))%>% 
    mutate(sil = unlist(sil)) %>% 
    
    # obtain the actual number of signature genes
    mutate(real_size = map_int(rdim, ~ .x$data_symbol %>% 
                                 map_int(~ n_distinct(.x$symbol)) %>% 
                                 unlist() %>% 
                                 unique() ))
  
}


# Preprocess===============================================

tt_simple <- readRDS("dev/intermediate_data/tt_simple.rds")

contrast_all <- 
  # tibble(level = 1:5) %>% 
  # mutate(level = paste("level", level, sep = "_")) %>% 
  # 
  # # preprocess data
  # mutate(tt = map(level, ~ counts %>% 
  #                   mutate(level_0 = "cell") %>% 
  #                   preprocess(.x))) %>% 
  tt_all %>% 
  
  # generate contrast by pairwise comparison
  mutate(contrast_PW = map2(tt, level, ~ contrast_PW(.x, .y) )) %>% 
  
  # generate contrast by mean contrast method
  mutate(contrast_MC = map2(tt, level, ~ contrast_MC(.x, .y) ))

# Signature genes =======================================================

## Identify optimal signature sizes for each cell type (of the ancestor nodes) for PCA method
opPCA_sig_PW <- c(44, 17, 60, 5, 1, 4, 15, 1, 1, 6, 0, 3, 4, 1, 3)
opPCA_sig_MC <- c(60, 41, 29, 31, 1, 4, 15, 1, 1, 6, 0, 3, 4, 1, 2)

## create a collection of all signature genes from each node at each level
sig_collect <- contrast_all %>% 
  unnest(tt) %>% 
  mutate(opPCA_sig_PW = opPCA_sig_PW) %>% 
  mutate(opPCA_sig_MC = opPCA_sig_MC) %>% 
  mutate(markers_PW = pmap(list(contrast_PW, level, opPCA_sig_PW), ~ sig_select(..1, ..2, ..3))) %>% 
  mutate(markers_MC = pmap(list(contrast_MC, level, opPCA_sig_MC), ~ sig_select(..1, ..2, ..3))) %>% 
  mutate(markers_PW = map(markers_PW, ~.x %>%
                                pull(symbol) %>% 
                                unique())) %>% 
  mutate(markers_MC = map(markers_MC, ~.x %>%
                                pull(symbol) %>% 
                                unique())) %>% 
  select(-c(level_0, level_1, level_2, level_3, level_4, data))

# Hierarchy + Pairwise ==============================================

# Combine signature genes from all levels into a single vector
all_sig_PW <- sig_collect %>% 
  pull(markers_PW) %>% 
  unlist() %>% 
  unique()

hierarchy_PW_sil <- tt_simple %>% 
  # filter signature genes selected by H+PW method from the non-hierarchical data
  mutate(data = map(data, ~ .x %>% filter(symbol %in% all_sig_PW))) %>% 
  unnest(data) %>% 
  sil_func0("PCA")

hierarchy_PW_sil

hierarchy_PW_pca <- hierarchy_PW_sil %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  # stat_ellipse(type = 't') +
  ggtitle("H_PW_pca") +
  theme_bw()

ggsave("hierarchy_PW_pca.png", hierarchy_PW_pca)

# Hierarchy + Mean Contrast ===================================

# Combine signature genes from all levels into a single vector
all_sig_MC <- sig_collect %>% 
  pull(markers_MC) %>%
  unlist() %>% 
  unique()

hierarchy_MC_sil <- tt_simple %>% 
  # filter signature genes selected by H+PW method from the non-hierarchical data
  mutate(data = map(data, ~ .x %>% filter(symbol %in% all_sig_MC))) %>% 
  unnest(data) %>% 
  sil_func0("PCA")

hierarchy_MC_sil

hierarchy_MC_pca <- hierarchy_MC_sil %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  # stat_ellipse(type = 't') +
  ggtitle("H_MC_pca") +
  theme_bw()

ggsave("hierarchy_MC_pca.png", hierarchy_MC_pca)
