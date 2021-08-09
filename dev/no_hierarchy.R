# No hierarchy + mean contrast

devtools::install_github("stemangiola/nanny@convert-to-S3", force = TRUE)
devtools::install_github("stemangiola/tidybulk@dev", force = TRUE)


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

# 1 Load data=============================================================================

# tt_all <- readRDS("dev/intermediate_data/tt_all.rds")

tt_simple <- readRDS("dev/intermediate_data/tt_simple.rds")

# 2 Functions=============================================================================

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

## 2.2 marker ranking & selection

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

## 2.4.1 calculate the area of confidence ellipses and the sum of their areas

## 5.1 calculate the area of confidence ellipses and the sum of their areas
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

## 5.2 Ellipse area calculation 
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

## 5.3 Scale serialised ellip_func() output (a tibble called ellip_tb) for plotting
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

# 3 No hierarchy + mean contrast Analysis

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

tt_simple <- readRDS("dev/intermediate_data/tt_simple.rds")

contrast_NH_MC <- tt_simple %>% contrast_MC0()

saveRDS(contrast_NH_MC, "contrast_NH_MC.rds")

contrast_NH_MC <- readRDS("dev/intermediate_data/contrast_NH.rds")

## 3.2 Ellipse Area Analysis

### 3.2.1 Single sig_size calculation
sig_size <- 20

contrast_NH_MC %>%
  sig_select0(sig_size) %>%
  nest(markers = - level_0) %>%
  mutate(ellip = map(markers, ~ ellip_func0(.x, "PCA")))

### 3.2.2 Serial sig_size calculation

sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
ellip_tb0 <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(ellip = map(sig_size, ~ sig_select0(contrast_NH_MC, .x))) %>%
  mutate(ellip = map(ellip, ~ ellip_func0(.x, "PCA")))

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
ellip_data0 <- ellip_tb0 %>% ellip_scale0()

cell_elli0 <- ellip_data0 %>%
  pluck("plot_data", 1) %>%
  ggplot(aes(real_size, area_value, colour=area_type)) +
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

cell_elli0

ggsave("cell_elli0.png", cell_elli0)

# Plot
sig_size <- 8

PCA0_cell_elli <- ellip_tb0 %>%
  pluck("ellip", 8) %>% 
  pluck("rdim", 1) %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA0_cell_elli

## 3.3 Silhouette Analysis

### 3.3.1 Single sig_size calculation

sig_size <- 46

NH_MC_sil <- contrast_NH_MC %>%
  sig_select0(sig_size) %>%
  sil_func0("PCA")


NH_MC_pca <- NH_MC_sil %>%
  pluck("rdim", 1) %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  # stat_ellipse(type = 't') +
  ggtitle("NH_MC_pca") +
  theme_bw()

NH_MC_pca

ggsave("NH_MC_pca.png", NH_MC_pca)

## 3.3.2 Serial sig_size calculation

sig_size <- 60

sil_tb0 <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ sig_select0(contrast_NH_MC, .x))) %>%
  mutate(sil_df = map(sil_df, ~ sil_func0(.x, "PCA")))

sil_data0 <- sil_tb0 %>%
  unnest(sil_df) %>%
  nest(plot_data = - level_0)

cell_sil0 <- sil_data0 %>%
  pluck("plot_data", 1) %>%
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

cell_sil0

# find out which sig_size gives the maximum silhouette score
sil_data0 %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  which.max()

# find out what the maximum silhouette score is
sil_data0 %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  max()

# Plot
sig_size <- 16

PCA0_cell_sil <- sil_tb0 %>%
  pluck("sil_df", 16) %>% 
  pluck("rdim", 1) %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA0_cell_sil

# 4 No hierarchy + pairwise comparison Analysis

## 4.1 preprocess data==================================================================

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

tt_simple <- readRDS("dev/intermediate_data/tt_simple.rds")

contrast_NH_PW <- tt_simple %>% 
  contrast_PW0()

saveRDS(contrast_NH_PW, "contrast_NH_PW.rds")

contrast_NH_PW <- readRDS("dev/intermediate_data/contrast_NH_pairwise.rds")

## 4.2 Ellipse Area Analysis

### 4.2.1 Single sig_size calculation
sig_size <- 20

contrast_NH_PW %>%
  sig_select0(sig_size) %>%
  nest(markers = - level_0) %>%
  mutate(ellip = map(markers, ~ ellip_func0(.x, "PCA")))

### 4.2.2 Serial sig_size calculation

sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
ellip_tb0 <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(ellip = map(sig_size, ~ sig_select0(contrast_NH_PW, .x))) %>%
  mutate(ellip = map(ellip, ~ ellip_func0(.x, "PCA")))

saveRDS(ellip_tb0, "ellip_tb0.rds")

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
ellip_data0 <- ellip_tb0 %>% ellip_scale0()

cell_elli0 <- ellip_data0 %>%
  pluck("plot_data", 1) %>%
  ggplot(aes(real_size, area_value, colour=area_type)) +
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

cell_elli0

# Plot

sig_size <- 3

PCA0_cell_elli <- ellip_tb0 %>%
  pluck("ellip", 3) %>% 
  pluck("rdim", 1) %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA0_cell_elli

## 4.3 Silhouette Analysis

### 4.3.1 Single sig_size calculation

sig_size <- 3

NH_PW_sil <- contrast_NH_PW %>%
  # select markers and nest by ancestor cell type
  sig_select0(sig_size) %>%
  sil_func0("PCA")

NH_PW_pca <- NH_PW_sil %>%
  pluck("rdim", 1) %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  # stat_ellipse(type = 't') +
  ggtitle("NH_PW_pca") +
  theme_bw()

NH_PW_pca

ggsave("NH_PW_pca.png", NH_PW_pca)

## 4.3.2 Serial sig_size calculation

sig_size <- 30

sil_tb0 <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ sig_select0(contrast_NH_PW, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func0("PCA")))


sil_data0 <- sil_tb0 %>%
  unnest(sil_df) %>%
  nest(plot_data = - level_0)


cell_sil0 <- sil_data0 %>%
  pluck("plot_data", 1) %>%
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

cell_sil0

sil_data0 %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  which.max()

sil_data0 %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  max()

# Plot

sig_size <- 3

PCA0_cell_sil <- sil_tb0 %>%
  pluck("sil_df", 3) %>% 
  pluck("rdim", 1) %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA0_cell_sil
