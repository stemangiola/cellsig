
library(tidyverse)
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

## 1 preprocess data

### 1.1 string manipulation that converts level of interest (e.g "level_5") to its ancestor level (e.g "level_4")
pre <- function(.level) {
  .level %>% 
    str_split("_") %>% 
    {as.numeric(.[[1]][2])-1} %>% 
    paste("level", ., sep = "_")
}

### 1.2 preprocess

preprocess <- function(.data, LEVEL) {
  
  # load data
  .data %>%
    
    tidybulk(sample, symbol, count) %>%
    
    # filter for the cell types of interest for gene marker selection
    filter(is.na(!!as.symbol(LEVEL))==F) %>%
    
    # Imputation of missing data within each level_5
    # impute_missing_abundance(~ !!as.symbol(LEVEL)) %>%
    
    # Group by ancestor
    nest(data = - !!as.symbol(pre(LEVEL))) %>%
    
    # Eliminate genes that are present in some but all cell types
    # (can be still present in a subset of samples from each cell-type)
    mutate(data = map(
      data,
      ~ .x %>%
        nest(data = -c(symbol, !!as.symbol(LEVEL))) %>%
        add_count(symbol) %>%
        filter(n == max(n)) %>%
        unnest(data)
    )) %>%
    
    # Imputation of missing data within each level_5
    mutate(data = map(data, ~ .x %>% impute_missing_abundance(~ !!as.symbol(LEVEL)))) %>%
    
    # scale count for further analysis
    mutate(data=map(data, ~ .x %>%
                      identify_abundant(factor_of_interest = !!as.symbol(LEVEL)) %>%
                      scale_abundance()
    ))
}

## 2 contrast functions

### 2.1 pairwise comparisons
get_contrasts_from_df = function(.data, LEVEL){
  
  .data %>% 
    distinct(!!as.symbol(LEVEL)) %>% 
    mutate(!!as.symbol(LEVEL) := paste0(LEVEL, !!as.symbol(LEVEL))) %>% 
    
    # Permute
    mutate(cell_type2 := !!as.symbol(LEVEL)) %>% 
    expand(!!as.symbol(LEVEL), cell_type2) %>% 
    filter(!!as.symbol(LEVEL) != cell_type2) %>% 
    
    # Create contrasts
    mutate(contrast = sprintf("%s - %s", !!as.symbol(LEVEL), cell_type2)) %>%
    pull(contrast)
  
}

### 2.2 create a contrast vector for limma::makeContrasts() or tidybulk::test_differential abundance()

make_contrasts <- function(.data, LEVEL){
  
  # find all cell types
  cell_types <- .data %>% 
    distinct(!!as.symbol(LEVEL)) %>% 
    pull() %>% 
    as.vector()
  
  # format cell_types with prefix
  cell_types <- paste0(LEVEL, cell_types)
  
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

select_markers_for_each_contrast = function(.markers, sig_size){
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
                           slice(1:sig_size)
                         
    )) %>%
    
    unnest(stat_df)
}

## 4 marker collection for each contrast

contrast <- function(tt, LEVEL, sig_size){
  tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + !!as.symbol(LEVEL), 
                                    .contrasts = get_contrasts_from_df(.x, LEVEL),
                                    action="only") 
    )) %>%
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, sig_size))) %>%
    
    # Add original data info to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(!!as.symbol(pre(LEVEL)), markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, LEVEL, "") %>% str_replace(LEVEL, ""))
}


## 5 calculate the area of confidence ellipses and the sum of their areas

ellip_area <- function(.markers, LEVEL){
  
  # reduce dimension
  PCA <- .markers %>%
    distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)) %>% 
    reduce_dimensions(sample, symbol, count_scaled, method = "PCA", action = "add", transform = log1p)
  
  # number of unique symbols for each cell type
  real_size <- PCA %>% 
    nest(data=-!!as.symbol(LEVEL)) %>% 
    mutate(real_size=map_int(data, ~ n_distinct(.x$symbol)))
  
  area <- PCA %>%   
    # remove non-numerical data to form a numerical data frame
    select(!!as.symbol(LEVEL), PC1, PC2) %>%
    
    # normalize principle component values
    mutate(across(c("PC1", "PC2"), ~ .x %>% scale())) %>% 
    
    # nest by cell_type so as to calculate ellipse area for each cell type
    nest(PC = - !!as.symbol(LEVEL)) %>% 
    
    # obtain covariance matrix for each cell type
    mutate(cov = map(PC, ~ cov(.x))) %>% 
    
    # calculate the eigenvalues for the covariance matrix of each cell type
    mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% 
    
    # transformation
    mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>%
    
    # below is the actual area for each ellipse
    mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% 
    
    # collect size of each cluster as factors for weights
    mutate(cluster_size = map_int(PC, ~ nrow(.x))) %>%
    
    # weight each area by the inverse of its cluster size
    mutate(weighted_area = map2_dbl(area, cluster_size, ~ .x / .y))
  
  left_join(area, real_size)
  
}

## 5.2 Ellipse area calculation for a series of sig_sizes
area_df_func <- function(.area_df, LEVEL){
  .area_df %>% 
    nest(markers = - !!as.symbol(pre(LEVEL))) %>% 
    mutate(ellip = map(markers, ~ .x %>% ellip_area(LEVEL)))
}

## 6 Silhouette score calculation for a series of sig_sizes
sil_func <- function(.sil_df, LEVEL){
  .sil_df %>%
    nest(pca = - !!as.symbol(pre(LEVEL))) %>%
    mutate(pca = map(pca, ~ .x %>%
                       distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>%
    mutate(pca = map(pca, ~ .x %>%
                       reduce_dimensions(sample, symbol, count_scaled,
                                         method = "PCA",
                                         action = "add",
                                         transform = log1p)
    )) %>%
    
    # calculate the dissimilarity matrix with PC values
    mutate(distance = map(pca, ~ .x %>%
                            select(contains("PC")) %>%
                            dist()
    )) %>%
    
    # calculate silhouette score
    mutate(sil = map2(pca, distance,
                      ~ silhouette(as.numeric(as.factor(`$`(.x, !!as.symbol(LEVEL)))), .y)
    )) %>%
    mutate(sil_info = map(sil, ~ .x %>% summary())) %>%
    mutate(sil_score = map(sil_info, ~ .x %>% `$`(avg.width))) %>%
    mutate(sil_score = unlist(sil_score)) %>%
    mutate(real_size=map_int(pca, ~ n_distinct(.x$symbol) ))
  
}


sil_score_at_level <- function(tt, level, sig_size) {
  tibble(sig_size = 1:sig_size) %>% 
    mutate(sil_df = map(sig_size, ~ contrast(tt, level, .x))) %>% 
    mutate(sil_df = map(sil_df, ~ .x %>% sil_func(level)))
}

# preprocess data at all levels except for no hierarchy

tt_all <- tibble(level = 1:5) %>% 
  mutate(level = paste("level", level, sep = "_")) %>% 
  mutate(tt = map(level, ~ counts %>% 
                    mutate(level_0 = "cell") %>% 
                    preprocess(.x)))

saveRDS(tt_all, "tt_all.rds")

# Pairwise + hierarchy

sig_size <- 20

pairwise <- tt_all %>%
  slice(3:5) %>% 
  mutate(sil = map2(tt, level, ~ sil_score_at_level(.x, .y, sig_size))) %>% 
  mutate(sil = map(sil, ~ .x %>% 
                     unnest(sil_df) %>% 
                     mutate(ancestor_type = select(., contains("level_")) %>% 
                              as_vector() ) %>% 
                     select(sig_size, real_size, sil_score, ancestor_type)
                            )) %>% 
  select(sil) %>% 
  unnest(sil) %>% 
  mutate(analysis = "PW+H")

# Mean contrast + hierarchy

contrast <- function(tt, LEVEL, sig_size){
  tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + !!as.symbol(LEVEL), 
                                    .contrasts = make_contrasts(.x, LEVEL),
                                    action="only") 
    )) %>%
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, sig_size))) %>%
    
    # Add original data info to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(!!as.symbol(pre(LEVEL)), markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, LEVEL, "") %>% str_replace(LEVEL, ""))
}

sig_size <- 20

mean_contrast <- tt_all %>%
  slice(3:5) %>% 
  mutate(sil = map2(tt, level, ~ sil_score_at_level(.x, .y, sig_size))) %>% 
  mutate(sil = map(sil, ~ .x %>% 
                     unnest(sil_df) %>% 
                     mutate(ancestor_type = select(., contains("level_")) %>% 
                              as_vector() ) %>% 
                     select(sig_size, real_size, sil_score, ancestor_type)
  )) %>% 
  select(sil) %>% 
  unnest(sil) %>% 
  mutate(analysis = "MC+H")

# Mean contrast + No hierarchy

# Functions====================================================

make_contrasts0 <- function(.data, cell_type){
  
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


## 2.3 marker collection for each contrast

contrast0 <- function(tt, cell_type, sig_size){
  tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = make_contrasts0(.x, cell_type),
                                    action="only") 
    )) %>%
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, sig_size))) %>%
    
    # Add original data info to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(level_0, markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
}


## 2.5 Silhouette function
sil_func0 <- function(.sil_df, cell_type){
  .sil_df %>% 
    nest(pca = - level_0) %>% 
    mutate(pca = map(pca, ~ .x %>% 
                       distinct(sample, symbol, count_scaled, cell_type))) %>% 
    mutate(pca = map(pca, ~ .x %>% 
                       reduce_dimensions(sample, symbol, count_scaled,
                                         method = "PCA",
                                         action = "add",
                                         transform = log1p)
    )) %>% 
    
    # calculate the dissimilarity matrix with PC values
    mutate(distance = map(pca, ~ .x %>% 
                            select(contains("PC")) %>% 
                            dist()
    )) %>% 
    
    # calculate silhouette score
    mutate(sil = map2(pca, distance, 
                      ~ silhouette(as.numeric(as.factor(`$`(.x, cell_type))), .y)
    )) %>% 
    mutate(sil_info = map(sil, ~ .x %>% summary())) %>% 
    mutate(sil_score = map(sil_info, ~ .x %>% `$`(avg.width))) %>% 
    mutate(sil_score = unlist(sil_score)) %>% 
    mutate(real_size=map_int(pca, ~ n_distinct(.x$symbol) ))
  
}

sil_score_no_hierarchy <- function(tt, cell_type, sig_size) {
  tibble(sig_size = 1:sig_size) %>% 
    mutate(sil_df = map(sig_size, ~ contrast0(tt, cell_type, .x))) %>% 
    mutate(sil_df = map(sil_df, ~ .x %>% sil_func0(cell_type)))
}

# preprocess===============

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

# =================
sig_size <- 2

no_hierarchy <- 
  tibble(sig_size = 1:sig_size) %>% 
  mutate(sil_df = map(sig_size, ~ contrast0(tt_simple, cell_type, .x))) %>% 
  mutate(sil_df = map(sil_df, ~.x %>% sil_func0(cell_type)))




# Commbine all sil_score from all three types of analysis

final <- bind_rows(pairwise, mean_contrast)

final %>% 
  ggplot(aes(sig_size, sil_score, 
             group=interaction(ancestor_type, analysis), 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5))


# Test (To be deleted)

pairwise %>% 
  ggplot(aes(sig_size, sil_score, group=ancestor_type, color=analysis)) +
  geom_line() +
  geom_point()

mean_contrast %>% 
  ggplot(aes(sig_size, sil_score, group=ancestor_type, color=analysis)) +
  geom_line() +
  geom_point()
