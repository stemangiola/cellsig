devtools::install_github("stemangiola/nanny@convert-to-S3", force = TRUE)
devtools::install_github("stemangiola/tidybulk@dev", force = TRUE)


library(tidyverse)
library(plotly)
library(future)
library(furrr)
library(tidybulk)
library(cluster)
library(proxy)
library(factoextra)
library(stringr)
plan(multiprocess, workers=5)


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






# NEW Functions for new data format ======================================================

# Functions for Hierarchical Analysis

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
    contrasts[i] <- sprintf("%s - (%s)/%s", cell_types[i], background, divisor)
  }
  
  return(contrasts)
}

### 2.3 contrast with no average background (OBSOLETE)

mean_contrast2 <- function(.data, .level){
  
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
  for(i in 1: length(cell_types) ){
    background = paste(cell_types[-i], collapse = "+")
    contrasts[i] <- sprintf("%s-(%s)", cell_types[i], background)
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

### mean contrast method
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

## 5 calculate the area of confidence ellipses and the sum of their areas

## 5.1 calculate the area of confidence ellipses and the sum of their areas
ellipse <- function(.rdim, .level, .method) {
  .rdim %>% 
    
    # remove non-numerical data to form a numerical data frame
    select(!!as.symbol(.level), 
           contains(str_sub(.method, end=-2L))) %>% 
    
    # normalize principle component values
    mutate(across(contains(str_sub(.method, end=-2L)), scale)) %>% 
    
    # nest by cell_type so as to calculate ellipse area for each cell type
    nest(dims = - !!as.symbol(.level)) %>% 
    
    # obtain covariance matrix for each cell type
    mutate(cov = map(dims, ~ cov(.x))) %>% 
    
    # calculate the eigenvalues for the covariance matrix of each cell type
    mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% 
    
    # transformation
    mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>%
    
    # below is the actual area for each ellipse
    mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% 
    
    # collect size of each cluster as factors for weights
    mutate(cluster_size = map_int(dims, ~ nrow(.x))) %>%
    
    # weight each area by the inverse of its cluster size
    mutate(weighted_area = map2_dbl(area, cluster_size, ~ .x / .y))
}

## 5.2 Ellipse area calculation 
ellip_func <- function(.markers, .level, .method){
  .markers %>% 
    
    # nest by ancestor cell types
    nest(rdim = - !!as.symbol(pre(.level))) %>%
    
    # reduce dimension
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
    
    mutate(real_size = map_int(rdim, ~ .x$data_symbol %>% 
                                 map_int(~ n_distinct(.x$symbol)) %>% 
                                 unlist() %>% 
                                 unique() )) %>% 
    
    mutate(area_df = map(rdim, ~ ellipse(.x, .level, .method) )
    )
  
}

## 5.3 Scale serialised ellip_func() output (a tibble called ellip_tb) for plotting
ellip_scale <- function(.ellip_tb, .level) {
  .ellip_tb %>% 
    unnest(ellip) %>%
    unnest(area_df) %>%
    
    # nest by ancestor cell type to rescale area for all sig_sizes
    nest(cell_data = - !!as.symbol(pre(LEVEL))) %>%
    mutate(cell_data = map(cell_data, ~ .x %>% 
                             mutate(rescaled_area = area %>% 
                                      scale(center = F))
    )) %>%
    
    # nest by ancestor cell type to summarise areas for each real_size/sig_size
    mutate(plot_data = map(cell_data, ~ .x %>%
                             # sum all areas for each real_size for an ancestor node
                             group_by(real_size) %>%
                             summarise(sig_size,
                                       stded_sum=sum(area, na.rm = T),
                                       wted_sum = sum(weighted_area, na.rm = T),
                                       rescaled_sum= sum(rescaled_area, na.rm = T)) %>%
                             # remove duplicate rows
                             distinct(real_size, sig_size, stded_sum, wted_sum, rescaled_sum) %>% 
                             pivot_longer(ends_with("sum"), names_to='area_type', values_to="area_value")
    ))
}

## 6 Silhoette function
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

### for summary plot
sil_tb <- function(.contrast, .level, .sig_size, .method) {
  tibble(sig_size = 1: .sig_size) %>% 
    
    # select signature genes for each sig_size at each level
    mutate(sil_df = map(sig_size, ~ sig_select(.contrast, .level, .x))) %>% 
    
    # calculate silhouette score for each ancestor cell type under sil_df
    mutate(sil_df = map(sil_df, ~ sil_func(.x, .level, .method)))
}

# Functions for Non-hierarchical Analysis
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
    contrasts[i] <- sprintf("%s - (%s)/%s", cell_types[i], background, divisor)
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

## 2.4 Ellipse area

## 2.4.1 calculate the area of confidence ellipses and the sum of their areas
ellipse0 <- function(.rdim, .method) {
  .rdim %>% 
    
    # remove non-numerical data to form a numerical data frame
    select(cell_type, 
           contains(str_sub(.method, end=-2L))) %>% 
    
    # normalize principle component values
    mutate(across(contains(str_sub(.method, end=-2L)), scale)) %>% 
    
    # nest by cell_type so as to calculate ellipse area for each cell type
    nest(dims = - cell_type) %>% 
    
    # obtain covariance matrix for each cell type
    mutate(cov = map(dims, ~ cov(.x))) %>% 
    
    # calculate the eigenvalues for the covariance matrix of each cell type
    mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% 
    
    # transformation
    mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>%
    
    # below is the actual area for each ellipse
    mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% 
    
    # collect size of each cluster as factors for weights
    mutate(cluster_size = map_int(dims, ~ nrow(.x))) %>%
    
    # weight each area by the inverse of its cluster size
    mutate(weighted_area = map2_dbl(area, cluster_size, ~ .x / .y))
}

## 2.4.2 Ellipse area calculation 
ellip_func0 <- function(.markers, .method){
  .markers %>% 
    
    # nest by ancestor cell types
    nest(rdim = - level_0) %>%
    
    # reduce dimension
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
    
    mutate(real_size = map_int(rdim, ~ .x$data_symbol %>% 
                                 map_int(~ n_distinct(.x$symbol)) %>% 
                                 unlist() %>% 
                                 unique() )) %>% 
    
    mutate(area_df = map(rdim, ~ ellipse0(.x, .method) ))
  
}

## 2.4.3 Scale serialised ellip_func() output (a tibble called ellip_tb) for plotting
ellip_scale0 <- function(.ellip_tb) {
  .ellip_tb %>% 
    unnest(ellip) %>%
    unnest(area_df) %>%
    
    # nest by ancestor cell type to rescale area for all sig_sizes
    nest(cell_data = - level_0) %>%
    mutate(cell_data = map(cell_data, ~ .x %>% 
                             mutate(rescaled_area = area %>% 
                                      scale(center = F))
    )) %>%
    
    # nest by ancestor cell type to summarise areas for each real_size/sig_size
    mutate(plot_data = map(cell_data, ~ .x %>%
                             # sum all areas for each real_size for an ancestor node
                             group_by(real_size) %>%
                             summarise(sig_size,
                                       stded_sum=sum(area, na.rm = T),
                                       wted_sum = sum(weighted_area, na.rm = T),
                                       rescaled_sum= sum(rescaled_area, na.rm = T)) %>%
                             # remove duplicate rows
                             distinct(real_size, sig_size, stded_sum, wted_sum, rescaled_sum) %>% 
                             pivot_longer(ends_with("sum"), names_to='area_type', values_to="area_value")
    ))
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

### for summary plot
sil_tb0 <- function(.contrast, .sig_size, .method) {
  tibble(sig_size = 1: .sig_size) %>% 
    
    # select signature genes for each sig_size at each level
    mutate(sil_df = map(sig_size, ~ sig_select0(.contrast, .x))) %>% 
    
    # calculate silhouette score for each ancestor cell type under sil_df
    mutate(sil_df = map(sil_df, ~ sil_func0(.x, .method)))
}

# Functions for Shiny App
format_name <- function(.method) {
  paste("markers", .method, sep = "_")
}

cell_sig_select <- function(.markers) {
  .markers %>% 
    # obtain cell types in a node and nest by it to extract signatures for all cell types
    mutate(cell_type = str_extract(contrast_pretty, "([a-z]|\\_)+(?=\\s)")) %>% 
    nest(signature = - cell_type) %>% 
    mutate(signature = map(signature, ~ .x %>% 
                             pull(symbol) %>% 
                             unique()
    ))
}