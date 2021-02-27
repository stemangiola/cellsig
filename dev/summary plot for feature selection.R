
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

# All functions

## Functions for hierarchy analysis

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


## 6 Silhouette score calculation for nodes in a level using .method of dimension reduction
sil_func <- function(.markers, .level, .method){
  .markers %>%
    nest(rdim = - !!as.symbol(pre(.level))) %>%
    mutate(rdim = map(rdim, ~ .x %>%
                       distinct(sample, symbol, count_scaled, !!as.symbol(.level)))) %>%
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
                      ~ silhouette(as.numeric(as.factor(`$`(.x, !!as.symbol(.level)))), .y)
    )) %>%
    mutate(sil = map(sil, ~ .x %>% summary())) %>%
    mutate(sil = map(sil, ~ .x %>% 
                       `$`(avg.width) ))%>% 
    mutate(sil = unlist(sil)) %>% 
    
    # obtain the actual number of signature genes
    mutate(real_size = map_int(rdim, ~ .x$data_symbol %>% 
                                 map_int(~ n_distinct(.x$symbol)) %>% 
                                 unlist() %>% 
                                 unique() ))
  
}

sil_tb <- function(.contrast, .level, .sig_size, .method) {
  tibble(sig_size = 1: .sig_size) %>% 
    
    # select signature genes for each sig_size at each level
    mutate(sil_df = map(sig_size, ~ sig_select(.contrast, .level, .x))) %>% 
    
    # calculate silhouette score for each ancestor cell type under sil_df
    mutate(sil_df = map(sil_df, ~ sil_func(.x, .level, .method)))
}


# Functions for no hierarchy analysis

## 2.1

mean_contrast0 <- function(.data){
  
  prefix <- "cell_type"
  
  # find all cell types
  cell_types <- .data %>% 
    distinct(cell_type) %>% 
    pull() %>% 
    as.vector()
  
  # format cell_types with prefix
  cell_types <- paste0(prefix, cell_types)
  
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

get_contrasts_from_df0 = function(.data){
  
  .data %>% 
    distinct(cell_type) %>% 
    mutate(cell_type = paste0("cell_type",cell_type)) %>% 
    
    # Permute
    mutate(cell_type2 = cell_type) %>% 
    expand(cell_type, cell_type2) %>% 
    filter(cell_type != cell_type2) %>% 
    
    # Create contrasts
    mutate(contrast = sprintf("%s - %s", cell_type, cell_type2)) %>%
    pull(contrast)
  
}

## 2.3 marker collection for each contrast

# mean contrast method
contrast_MC0 <- function(.tt){
  .tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = mean_contrast0(.x),
                                    action="only")  ))
}

# pairwise contrast method
contrast_PW0 <- function(.tt){
  .tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = get_contrasts_from_df0(.x),
                                    action="only")  ))
}

# marker selection & processing
sig_select0 <- function(.contrast, .sig_size) {
  .contrast %>% 
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, .sig_size))) %>%
    
    # Add original data info to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(level_0, markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
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

sil_tb0 <- function(.contrast, .sig_size, .method) {
  tibble(sig_size = 1: .sig_size) %>% 
    
    # select signature genes for each sig_size at each level
    mutate(sil_df = map(sig_size, ~ sig_select0(.contrast, .x))) %>% 
    
    # calculate silhouette score for each ancestor cell type under sil_df
    mutate(sil_df = map(sil_df, ~ sil_func0(.x, .method)))
}

# Hierarchical analysis

## Preprocess data at ALL LEVELS EXCEPT FOR no hierarchy========================

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

saveRDS(contrast_all, "contrast_all.rds")

contrast_all <- readRDS("dev/intermediate_data/contrast_all.rds")

# Hierarchy + Pairwise ============

sig_size <- 20

H_PW <- contrast_all %>%
  # slice(3) %>%
  
  # calculate silhouette score for all signature sizes at each level of interest
  mutate(sil = map2(contrast_PW, level, ~ sil_tb(.x, .y, sig_size, "tSNE"))) %>%
  
  # select what is needed to plot sil score for each ancestor node(type) changing over sig_size
  mutate(sil = map(sil, ~ .x %>% 
                     unnest(sil_df) %>% 
                     mutate(ancestor_type = select(., contains("level_")) %>% 
                              as_vector() ) %>% 
                     select(sig_size, real_size, sil, ancestor_type)
                            )) %>% 
  select(sil) %>% 
  unnest(sil) %>% 
  mutate(analysis = "PW+H")


# Hierarchy + Mean contrast ============

sig_size <- 20

H_MC <- contrast_all %>%
  # slice(3) %>% 
  mutate(sil = map2(contrast_MC, level, ~ sil_tb(.x, .y, sig_size, "tSNE"))) %>% 
  mutate(sil = map(sil, ~ .x %>% 
                     unnest(sil_df) %>% 
                     mutate(ancestor_type = select(., contains("level_")) %>% 
                              as_vector() ) %>% 
                     select(sig_size, real_size, sil, ancestor_type)
  )) %>% 
  select(sil) %>% 
  unnest(sil) %>% 
  mutate(analysis = "MC+H")


# No hierarchy Analysis

## Preprocess for non-hierarchical analysis ====================================
tt_simple <- counts %>%
  
  mutate(level_0 = "cell") %>% 
  
  tidybulk(sample, symbol, count) %>%
  
  # filter for the cell types of interest for gene marker selection
  filter(is.na(cell_type)==F) %>%
  
  # Group by ancestor
  nest(data = - level_0) %>%
  
  # Eliminate genes that are present in some but all cell types
  # (can be still present in a subset of samples from each cell-type)
  mutate(data = map(
    data,
    ~ .x %>%
      nest(data = -c(symbol, cell_type)) %>%
      add_count(symbol) %>%
      filter(n == max(n)) %>%
      unnest(data)
  )) %>%
  
  # Imputation of missing data within each level_5
  mutate(data = map(data, ~ .x %>% impute_missing_abundance(~ cell_type))) %>%
  
  # scale count for further analysis
  mutate(data=map(data, ~ .x %>%
                    identify_abundant(factor_of_interest = cell_type) %>%
                    scale_abundance()
  ))

saveRDS(tt_simple, "tt_simple.rds")

tt_simple <- readRDS("dev/intermediate_data/tt_simple.rds")

contrast_NH_MC <- tt_simple %>% 
  contrast_MC0()

saveRDS(contrast_NH_MC, "contrast_NH_MC.rds")

contrast_NH_MC <- readRDS("dev/intermediate_data/contrast_NH.rds")

contrast_NH_PW <- tt_simple %>% 
  contrast_PW0()

saveRDS(contrast_NH_PW, "contrast_NH_PW.rds")

contrast_NH_PW <- readRDS("dev/intermediate_data/contrast_NH_pairwise.rds")

## No hierarchy + Mean contrast ===============
sig_size <- 20

NH_MC <- 
  sil_tb0(contrast_NH_MC, sig_size, "tSNE") %>% 
  unnest(sil_df) %>% 
  mutate(ancestor_type = level_0) %>% 
  select(sig_size, real_size, sil, ancestor_type) %>% 
  mutate(analysis = "MC+NH")

## No hierarchy + Pairwise comparison ===============
sig_size <- 20

NH_PW <- 
  sil_tb0(contrast_NH_PW, sig_size, "tSNE") %>% 
  unnest(sil_df) %>% 
  mutate(ancestor_type = level_0) %>% 
  select(sig_size, real_size, sil, ancestor_type) %>% 
  mutate(analysis = "PW+NH")


# Combine all sil_score from all three types of analysis

final_PCA <- bind_rows(H_PW, H_MC, NH_MC, NH_PW)

summary_plot_PCA <- final_PCA %>% 
  ggplot(aes(sig_size, sil, 
             group=interaction(ancestor_type, analysis), 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5))

summary_plot_PCA

ggsave("summary_plot_PCA.png", summary_plot_PCA)

final_tSNE <- bind_rows(H_PW, H_MC, NH_MC, NH_PW)

summary_plot_tSNE <- final_tSNE %>% 
  ggplot(aes(sig_size, sil, 
             group=interaction(ancestor_type, analysis), 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5))

summary_plot_tSNE

ggsave("summary_plot_tSNE.png", summary_plot_tSNE)

# Below is a summary plot comparing silhouette score for four selection methods + CIBERSORT
# It requires code from these files: cibersortx.R, hierarchy.R, no_hierarchy.R====================

sil_all_methods <- bind_rows(ciber_sil, hierarchy_PW_sil, hierarchy_MC_sil, NH_PW_sil, NH_MC_sil) %>% 
  mutate(method = c("CIBERSORTx", "H + PW", "H + MC", "NH + PW", "NH + MC")) %>% 
  select(method, sil, real_size) %>% 
  arrange(desc(sil))

sil_all_methods_plot <- sil_all_methods %>% 
  ggplot(aes(method, sil, fill = method)) +
  geom_col() +
  geom_text(aes(label = signif(sil, 3)), vjust = 1.5)

ggsave("sil_all_methods_plot.png", sil_all_methods_plot)

