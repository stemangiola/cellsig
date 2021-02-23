# No hierarchy + mean contrast 

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

# To be loaded after all libraries
library(tidybulk)
library(cluster)
library(proxy)
library(factoextra)
library(stringr)

# 1 Load data=============================================================================

load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")


# 2 Functions=============================================================================

## 2.1

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


## 2.2 marker ranking & selection

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

## 2.4.1 calculate the area of confidence ellipses and the sum of their areas

ellip_area0 <- function(.markers, cell_type){
  
  # reduce dimension
  PCA <- .markers %>%
    distinct(sample, symbol, count_scaled, cell_type) %>% 
    reduce_dimensions(sample, symbol, count_scaled, method = "PCA", action = "add", transform = log1p)
  
  # number of unique symbols for each cell type
  real_size <- PCA %>% 
    nest(data=- cell_type) %>% 
    mutate(real_size=map_int(data, ~ n_distinct(.x$symbol)))
  
  area <- PCA %>%   
    # remove non-numerical data to form a numerical data frame
    select(cell_type, PC1, PC2) %>%
    
    # normalize principle component values
    mutate(across(c("PC1", "PC2"), ~ .x %>% scale())) %>% 
    
    # nest by cell_type so as to calculate ellipse area for each cell type
    nest(PC = - cell_type) %>% 
    
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

## 2.4.2 
area_df_func0 <- function(.area_df, cell_type){
  .area_df %>% 
    nest(markers = - level_0) %>% 
    mutate(ellip = map(markers, ~ .x %>% ellip_area(cell_type)))
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
                      ~ silhouette(as.numeric(as.factor(`$`(.x, cell_type))), .y)
    )) %>% 
    mutate(sil_info = map(sil, ~ .x %>% summary())) %>% 
    mutate(sil_score = map(sil_info, ~ .x %>% `$`(avg.width))) %>% 
    mutate(sil_score = unlist(sil_score)) %>% 
    mutate(real_size=map_int(pca, ~ n_distinct(.x$symbol) ))
  
}

# 3 No hierarchy + mean contrast Analysis=============================================================================

## 3.1 preprocess data==================================================================

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

## 3.2 Ellipse Area Analysis

### 3.2.1 Single sig_size calculation
sig_size <- 20

tt_simple %>% 
  contrast0(cell_type, sig_size) %>% 
  nest(markers = - level_0) %>% 
  mutate(ellip = map(markers, ~ .x %>% ellip_area0(cell_type)))

### 3.2.2 Serial sig_size calculation

sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
area_tb0 <- 
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(area_df = map(sig_size, ~ contrast0(tt_simple, cell_type, .x))) %>% 
  mutate(area_df = map(area_df, ~ area_df_func0(.x, cell_type)))


# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
area_data0 <- area_tb %>% 
  
  unnest(area_df) %>% 
  unnest(ellip) %>% 
  
  # nest by ancestor cell type to rescale area for all sig_sizes
  nest(cell_data = - level_0) %>% 
  mutate(cell_data = map(cell_data, ~ .x %>% mutate(rescaled_area = area %>% scale(center = F)))) %>% 
  unnest(cell_data) %>%
  
  # nest by ancestor cell type to summarise areas for each real_size/sig_size
  nest(cell_data = - level_0) %>% 
  mutate(plot_data = map(cell_data, ~ .x %>% 
                           group_by(real_size) %>%
                           summarise(sig_size,
                                     stded_sum=sum(area, na.rm = T), 
                                     wted_sum = sum(weighted_area, na.rm = T), 
                                     rescaled_sum= sum(rescaled_area, na.rm = T)) %>% 
                           pivot_longer(ends_with("sum"), names_to='area_type', values_to="area_value")
  ))

cell_elli0 <- area_data0 %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

sig_size <- 4


## 3.3 Silhouette Analysis

### 3.3.1 Single sig_size calculation

sig_size <- 11

tt_simple %>%
  
  # select markers and nest by ancestor cell type
  contrast0(cell_type, sig_size) %>%
  nest(markers = - level_0) %>% 
  
  # reduce dimensions
  mutate(pca = map(markers, ~ .x %>% 
                     distinct(sample, symbol, count_scaled, cell_type))) %>% 
  mutate(pca = map(pca, ~ .x %>% 
                     reduce_dimensions(
                       sample, symbol, count_scaled, 
                       method = "PCA", 
                       action="add", 
                       transform = log1p))) %>% 
  
  # calculate the dissimilarity matrix with PC values
  mutate(distance = map(pca, ~ .x %>% 
                          select(contains("PC")) %>% 
                          dist()
  )) %>% 
  
  # calculate silouette score
  mutate(sil = map2(pca, distance, 
                    ~ silhouette(as.numeric(as.factor(`$`(.x, cell_type))), .y)
  )) %>% 
  mutate(sil_info = map(sil, ~ .x %>% summary())) %>% 
  mutate(sil_score = map(sil_info, ~ .x %>% `$`(avg.width)))

## 3.3.2 serial sig_size calculation

sig_size <- 20

sil_tb0 <- 
  tibble(sig_size = 1:sig_size) %>% 
  # slice(1) %>% 
  mutate(sil_df = map(sig_size, ~ contrast0(tt_simple, cell_type, .x))) %>% 
  mutate(sil_df = map(sil_df, ~.x %>% sil_func0(cell_type)))


sil_data0 <- sil_tb %>% 
  unnest(sil_df) %>% 
  nest(plot_data = - level_0)

cell_sil0 <- sil_data0 %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

sig_size <- 4


# 4 Plot

PCA_simple <- tt_simple %>%
  contrast0(cell_type, sig_size) %>% 
  nest(markers = - level_0) %>% 
  mutate(pca = map(markers, ~ .x %>% 
                     distinct(sample, symbol, count_scaled, cell_type))) %>% 
  mutate(pca = map(pca, ~ .x %>% 
                     reduce_dimensions(
                       sample, symbol, count_scaled, 
                       method = "PCA", 
                       action="add", 
                       transform = log1p)))


sig_size <- 4

PCA0_cell <- PCA_simple %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()
