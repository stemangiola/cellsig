devtools::install_github("stemangiola/nanny@convert-to-S3", force = TRUE)
devtools::install_github("stemangiola/tidybulk@dev", force = TRUE)


library(tidyverse)
library(plotly)
library(nanny)
library(ggrepel)
library(GGally)
library(tidyHeatmap)
library(future)
library(furrr)
plan(multiprocess, workers=5)

# To be loaded after all libraries
library(tidybulk)
library(cluster)
library(proxy)
library(factoextra)
library(stringr)

# OLD Functions for data of old format===============================================================================
## 1 preprocess data


preprocess <- function(.data) {
  # Load dataset
  .data %>%
    
    # Remove entries with NA in gene symbol
    filter(symbol %>% is.na %>% `!`) %>% 
    
    tidybulk(sample, symbol, count) %>%
    
    # Group by level because otherwise samples are duplicated
    nest(data = -level) %>%
    # filter(level ==3) %>%
    
    # Redefine factors inside each level
    mutate(data = map(data, ~ droplevels(.x))) %>%
    
    # Remove redundancy
    mutate(data = map(data, ~aggregate_duplicates(.x))) %>% 
    
    # Fill missing data. There are many genes that are not shared by the majority of samples
    mutate(data = map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
    
    
    # Scale for future PCA plotting
    mutate(data = map(
      data, ~ .x %>% 
        identify_abundant(factor_of_interest = cell_type) %>%
        scale_abundance() 
    ))
  
}



## 2 contrast functions

### 2.1 pairwise comparisons
get_contrasts_from_df = function(.data){
  .data %>% 
    
    distinct(cell_type) %>% 
    
    # Permute
    mutate(cell_type2 = cell_type) %>% 
    expand(cell_type, cell_type2) %>% 
    filter(cell_type != cell_type2) %>% 
    
    # Create contrasts
    mutate(contrast = sprintf("cell_type%s - cell_type%s", cell_type, cell_type2)) %>%
    pull(contrast)
  
}

### 2.2 create a contrast vector for limma::makeContrasts() or tidybulk::test_differential abundance()

make_contrasts <- function(tt, LEVEL){
  
  prefix = "cell_type"
  
  # find all cell types
  cell_types <- tt %>% 
    filter(level == LEVEL) %>% 
    unnest(data) %>% 
    distinct(cell_type) %>% 
    pull() 
  
  # format cell_types with prefix
  cell_types <- paste(prefix, cell_types, sep="")
  
  # initialise a vector called contrasts
  contrasts <- 1: length(cell_types)
  
  # create all contrasts and store them in contrasts
  for(i in 1: length(cell_types) ){
    background = paste(cell_types[-i], collapse = "+")
    divisor = length(cell_types[-i])
    contrasts[i] <- sprintf("%s-(%s)/%s", cell_types[i], background, divisor)
  }
  
  return(contrasts)
}


### 2.3 contrast with no average background

make_contrasts2 <- function(tt, LEVEL){
  
  prefix = "cell_type"
  
  # find all cell types
  cell_types <- tt %>% 
    filter(level == LEVEL) %>% 
    unnest(data) %>% 
    distinct(cell_type) %>% 
    pull() 
  
  # format cell_types with prefix
  cell_types <- paste(prefix, cell_types, sep="")
  
  # initialise a vector called contrasts
  contrasts <- 1: length(cell_types)
  
  # create all contrasts and store them in contrasts
  for(i in 1: length(cell_types) ){
    background = paste(cell_types[-i], collapse = "+")
    contrasts[i] <- sprintf("%s-(%s)", cell_types[i], background)
  }
  
  return(contrasts)
}

## 3 marker selection

select_markers_for_each_contrast = function(.data, sig_size){
  .data %>%
    
    # Group by contrast. Comparisons both ways.
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
                           filter(FDR < 0.05 & logFC > 2) %>%
                           filter(logCPM > mean(logCPM)) %>%
                           arrange(logFC %>% desc()) %>%
                           slice(1:sig_size)
                         
    )) %>%
    
    unnest(stat_df)
}


## 4

contrast <- function(tt, LEVEL, sig_size){
  tt %>%
    
    # Investigate one level
    filter(level == LEVEL) %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = get_contrasts_from_df(.x),
                                    action="only") 
    )) %>%
    
    # Select rank from each contrast
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, sig_size))) %>%
    
    # Add original data info to markers
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
}


## 5 calculate the area of confidence ellipses and the sum of their areas
PCA_level5

ellip_area <- function(all_contrasts){
  # reduce dimension
  PCA <- all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="add", transform = log1p)
  
  real_size <- PCA %>% 
    nest(data = -cell_type) %>% 
    mutate(real_size = map_int(data, ~ n_distinct(.x$symbol)))
  
  area <- PCA %>%   
    # remove non-numerical data to form a numerical data frame
    select(cell_type, PC1, PC2) %>%
    
    # normalize principle component values
    mutate(across(c("PC1", "PC2"), ~ .x %>% scale())) %>% 
    
    # nest by cell_type so as to calculate ellipse area for each cell type
    nest(PC = -cell_type) %>% 
    
    # obtain covariance matrix for each cell type
    mutate(cov = map(PC, ~ cov(.x))) %>% 
    
    # calculate the eigenvalues for the covariance matrix of each cell type
    mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% 
    
    # transformation
    mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>% # unnest(area)
    
    # below is the actual area for each ellipse
    mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% 
    
    # collect size of each cluster as factors for weights
    mutate(cluster_size = map_int(PC, ~ nrow(.x))) %>%
    
    # weight each area by the inverse of its cluster size
    mutate(weighted_area = map2_dbl(area, cluster_size, ~ .x / .y))
  
  left_join(area, real_size)
  
}


# NEW Functions for new data format=================================================================================================

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


### 2.3 contrast with no average background

make_contrasts2 <- function(.data, LEVEL){
  
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
  for(i in 1: length(cell_types) ){
    background = paste(cell_types[-i], collapse = "+")
    contrasts[i] <- sprintf("%s-(%s)", cell_types[i], background)
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

## 5.2 
area_df_func <- function(.area_df, LEVEL){
  .area_df %>% 
    nest(markers = - !!as.symbol(pre(LEVEL))) %>% 
    mutate(ellip = map(markers, ~ .x %>% ellip_area(LEVEL)))
}

## 6 Silhoette function
sil_func <- function(.sil_df, LEVEL){
  .sil_df %>% 
    nest(markers = - !!as.symbol(pre(LEVEL))) %>% 
    mutate(pca = map(markers, ~ .x %>% 
                       distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>% 
    mutate(pca = map(pca, ~ .x %>% 
                       reduce_dimensions(sample, symbol, count_scaled,
                                         method = "PCA",
                                         action = "add",
                                         transform = log1p)
    )) %>% 
    # mutate(pca_norm = pca) %>% 
    # mutate(pca_norm = map(pca_norm, ~ .x %>% 
    #                         mutate(across(c("PC1", "PC2"), scale) ))) %>%
    
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