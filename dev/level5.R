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

# Load data=================================================================================================

# load("dev/database/orig_counts.rda")

load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")

# Literal parameters========================================================================================

LEVEL = "level_5"

# Functions=================================================================================================

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

# Analysis=========================================================================================================

# 1 Setup data frame & preprocessing==============================================================================

tt <-
  
  # load data
  counts %>%
  
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

# 2 No hierachy + mean contrast on level 5================================================================================

# 3 No selection of markers========================================================================================
tt_naive <-  
  tt %>%
  
  # Scale and reduce dimensions
  mutate(data = map(
    data,
    ~ .x %>%
      reduce_dimensions(method="PCA") %>% 
      reduce_dimensions(method="MDS") %>%
      reduce_dimensions(method="tSNE")
  )) 
# %>% 
# # Cluster
# mutate(data = map(data, ~ cluster_elements(.x, method="SNN")))

PCA_naive_tCD4_memory <- tt_naive %>% 
  pluck("data", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_tCD8_memory <- tt_naive %>% 
  pluck("data", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_t_helper <- tt_naive %>% 
  pluck("data", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

# 4 Ellipse area==================================================================================================

## 4.1 Single sig_size calculation
sig_size <- 20

tt %>% 
  contrast(LEVEL, sig_size) %>% 
  nest(markers = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(ellip = map(markers, ~ .x %>% ellip_area(LEVEL)))

## 4.2 Serial sig_size calculation

sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
area_tb <- 
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(area_df = map(sig_size, ~ contrast(tt, LEVEL, .x))) %>% 
  mutate(area_df = map(area_df, ~ area_df_func(.x, LEVEL)))
  

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
area_data <- area_tb %>% 
  
  unnest(area_df) %>% 
  unnest(ellip) %>% 
  
  # nest by ancestor cell type to rescale area for all sig_sizes
  nest(cell_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(cell_data = map(cell_data, ~ .x %>% mutate(rescaled_area = area %>% scale(center = F)))) %>% 
  unnest(cell_data) %>%
 
  # nest by ancestor cell type to summarise areas for each real_size/sig_size
  nest(cell_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(cell_data, ~ .x %>% 
                           group_by(real_size) %>%
                           summarise(sig_size,
                                     stded_sum=sum(area, na.rm = T), 
                                     wted_sum = sum(weighted_area, na.rm = T), 
                                     rescaled_sum= sum(rescaled_area, na.rm = T)) %>% 
                           pivot_longer(ends_with("sum"), names_to='area_type', values_to="area_value")
                           ))
  
tCD4_memory_elli <- area_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

sig_size <- 4

tCD8_memory_elli <- area_data %>% 
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

sig_size <- 5

t_helper_elli <- area_data %>% 
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

sig_size <- 3

# %>% 
# group_by %>% 
# make_the_sum %>% 
# select_only_columns_of_interest %>%
# pivot_longer(areas, names_to area_type, area_value) %>% 
# ggplot(aes(real_size, area_value, color = area_type)) + geom_line()

area_sum <- 1: sig_size

for (i in 1: sig_size) {
  sig_size <- i
  print(sig_size)
  area_sum[i] <- tt %>% 
    contrast(LEVEL, sig_size) %>% 
    ellip_area()
}

tb <- tibble(size = 1: sig_size, total_area = area_sum)
ggplot(tb, aes(size, total_area)) +
  geom_line() +
  geom_point()

which.min(area_sum)

min(area_sum)

# when using get_contrasts_from_df(), sig_size = 44 gives the minimum total_area: 8.840;

# when using contrast1, sig_size = 8  gives the elbow point in all three scaling methods. 

# the PCA plots 336 most variable genes, which is when sig_size = 56.
# second minimum when sig_size <- 28, gives 0.204)

# when using contrast2, sig_size = 39  gives the minimum total_area: 8.014

sig_size <- 4
sig_size <- 8
sig_size <- 58
sig_size <- 39

# 5 Silhouette analysis=============================================================================

## 5.1 Single sig_size calculation

sig_size <- 11

tt %>%
  
  # select markers and nest by ancestor cell type
  contrast(LEVEL, sig_size) %>%
  nest(markers = - !!as.symbol(pre(LEVEL))) %>% 
 
   # reduce dimensions
  mutate(pca = map(markers, ~ .x %>% 
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>% 
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
                    ~ silhouette(as.numeric(as.factor(`$`(.x, !!as.symbol(LEVEL)))), .y)
                    )) %>% 
  mutate(sil_info = map(sil, ~ .x %>% summary())) %>% 
  mutate(sil_score = map(sil_info, ~ .x %>% `$`(avg.width)))

## 5.2 serial sig_size calculation

sig_size <- 20

sil_tb <- 
  tibble(sig_size = 1:sig_size) %>% 
  # slice(1) %>% 
  mutate(sil_df = map(sig_size, ~ contrast(tt, LEVEL, .x))) %>% 
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL)))
 

sil_data <- sil_tb %>% 
  unnest(sil_df) %>% 
  nest(plot_data = - !!as.symbol(pre(LEVEL)))
                           
tCD4_memory_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

sig_size <- 4

tCD8_memory_sil <- sil_data %>%
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

sig_size <- 1


t_helper_sil <- sil_data %>%
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

sig_size <- 3

for (i in 1:sig_size){
  sig_size = i
  print(sig_size)
  all_contrasts <- contrast(tt, LEVEL, sig_size)
  
  PCA_level5 <- all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="add", transform = log1p)
  
  real_size <- PCA_level5 %>% 
    nest(data = -cell_type) %>% 
    mutate(real_size = map_int(data, ~ n_distinct(.x$symbol)))
  
  distance <- PCA_level5 %>% 
    select(contains("PC")) %>% 
    dist()
  
  sil <- silhouette(as.numeric(as.factor(PCA_level5$cell_type)), distance)
  sil_info <- summary(sil)
  sil_scores[i] <- sil_info$avg.width
}

tb2 <- tibble(size = 1:sig_size, sil_scores = sil_scores)
ggplot(tb2, aes(size, sil_scores)) +
  geom_line() +
  geom_point()

max(sil_scores)
which.max(sil_scores)

# when using get_contrasts_from_df(), sig_size = 48 gives maximum silhouette score: 0.349; 
# when using contrast1, sig_size = 3 gives maximum silhouette score:  -0.09251488. (!bad)
# when using contrast2, sig_size = 43 gives maximum silhouette score: 0.358

sig_size <- 48
sig_size <- 3
sig_size <- 43


# 6 Plots==============================================================================================
## 6.1 Plot Markers
all_contrasts %>% 
  nest(markers= - !!as.symbol(pre(LEVEL))) %>%
  pluck("markers", 3) %>% 
  tidyr::nest(data = -c(contrast, symbol)) %>%
  slice(1:20) %>%
  unnest(data) %>%
  ggplot(aes(!!as.symbol(LEVEL), count_scaled + 1, color=contrast_pretty)) + 
  geom_point(size = 0.5) + 
  facet_wrap(~contrast_pretty + symbol, scale="free_y") + 
  scale_y_log10() +
  theme_bw() +
  theme(
    text = element_text(size=6), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) 


## 6.2 Plot PCA

PCA_level5 <- tt %>%
  contrast(LEVEL, sig_size) %>% 
  nest(markers= - !!as.symbol(pre(LEVEL))) %>% 
  mutate(pca = map(markers, ~ .x %>% 
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>% 
  mutate(pca = map(pca, ~ .x %>% 
                     reduce_dimensions(
                       sample, symbol, count_scaled, 
                       method = "PCA", 
                       action="add", 
                       transform = log1p)))
  

# plot PCA for PC values normalised within clusters
# group_by(cell_type) %>% 
# mutate(across(c("PC1", "PC2"), ~ (.x - mean(.x))/sd(.x) ))

# plot PCA for PC values normalised globally
# mutate(across(c("PC1", "PC2"), ~ (.x - mean(.x))/sd(.x) )) %>% 

# plot PCA for PC values normalised first globally then within clusters
# group_by(cell_type) %>%
# mutate(across(c("PC1", "PC2"), ~ (.x - mean(.x))/sd(.x) ))

sig_size <- 4

PCA5_tCD4_memory <- PCA_level5 %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

sig_size <- 5

PCA5_tCD8_memory <- PCA_level5 %>% 
  pluck("pca", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

sig_size <- 3

PCA5_t_helper <- PCA_level5 %>% 
  pluck("pca", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD4_memory %>% ggplotly(tooltip = c("label", LEVEL))
PCA5_tCD8_memory %>% ggplotly(tooltip = c("label", LEVEL))
PCA5_t_helper %>% ggplotly(tooltip = c("label", LEVEL))

sig_size <- 4

PCA5_tCD4_memory_sil <- PCA_level5 %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

sig_size <- 1

PCA5_tCD8_memory_sil <- PCA_level5 %>% 
  pluck("pca", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

sig_size <- 3

PCA5_t_helper_sil <- PCA_level5 %>% 
  pluck("pca", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()


# 7 NEW =================

```{r library, message = FALSE}
# Load libraries
library(tidyverse)
library(plotly)
library(future)
library(furrr)
library(tidybulk)
library(cluster)
library(proxy)
library(factoextra)
library(stringr)
```


```{r data}
# Load data
tt_all <- readRDS("intermediate_data/tt_all.rds")

tt_L5 <- tt_all %>% 
  pluck("tt", 5)
```

```{r literal}
# select level of interest
LEVEL = "level_5"
```


```{r functions, message=F}
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

## 5 
##5.1 calculate the area of confidence ellipses and the sum of their areas
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

## 6 Silhouette score calculation for a series of sig_sizes
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
```


# 1 Pre-analysis

## 1.1 Preprocess data
```{r preprocess, message=F, results = FALSE, warning = FALSE}
# 1 Setup data frame & preprocessing

# tt_L5 <- preprocess(counts, LEVEL)

```

```{r cell_types, message=F, warning=F}
# View cell types on level_1
tt_L5 %>% 
  unnest(data) %>% 
  select(!!as.symbol(LEVEL)) %>% 
  distinct()

```

## 1.2 Cluster without marker selection for contrast

```{r naive, results = FALSE, warning = FALSE, message=F}
# 2 No selection of markers
tt_naive_L5 <-
  tt_L5 %>%
  
  # Scale and reduce dimensions
  mutate(data = map(
    data,
    ~ .x %>%
      reduce_dimensions(method="PCA")
    # reduce_dimensions(method="MDS") %>%
    # reduce_dimensions(method="tSNE")
  ))

```

```{r PCA_naive_tCD4_memory, message=F}
PCA_naive_tCD4_memory <- tt_naive_L5 %>% 
  pluck("data", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_tCD4_memory
```

```{r PCA_naive_tCD8_memory, message=F}
PCA_naive_tCD8_memory <- tt_naive_L5 %>% 
  pluck("data", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_tCD8_memory
```

```{r PCA_naive_t_helper, message=F}
PCA_naive_t_helper <- tt_naive_L5 %>% 
  pluck("data", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_t_helper
```

# 2 Hierarchy + Pairwise Analysis

## 2.0 Generate pairwise contrast

```{r pairwise_contrast, message=F, warning = FALSE}
contrast_PW_L5 <- tt_L5 %>% 
  contrast_PW(LEVEL)
```

## 2.1 Ellipse Area Analysis

```{r ellipse, message=F, results = FALSE, warning = FALSE,}
sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
ellip_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(ellip = map(sig_size, ~ sig_select(contrast_PW_L5, LEVEL, .x))) %>%
  mutate(ellip = map(ellip, ~ ellip_func(.x, LEVEL, "PCA") ))

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
ellip_data <- ellip_tb %>% ellip_scale(LEVEL)

```

### 2.1.1 Signature size selection & PCA plot for t_CD4_memory cell

Signature size selection plot for t_CD4_memory cell:
  
  ```{r tCD4_memory_elli, message=F}
tCD4_memory_elli <- ellip_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

tCD4_memory_elli
```

The elbow point indicates optimal sig_size is $4$. PCA at sig_size = $4$: 
  
  ```{r PCA5_tCD4_memory, message=F, message=F, warning = FALSE}
sig_size <- 4

PCA5_tCD4_memory <- ellip_tb %>% 
  pluck("ellip", 4) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD4_memory
```

### 2.1.2 Signature size selection & PCA plot for t_CD8_memory cell

Signature size selection plot for t_CD8_memory cell:
  
  ```{r tCD8_memory, message=F}
tCD8_memory <- ellip_data %>% 
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

tCD8_memory
```

The elbow point indicates optimal sig_size is $5$. PCA at sig_size = $5$: 
  
  ```{r PCA5_tCD8_memory, message=F, warning = FALSE}
sig_size <- 5

PCA5_tCD8_memory <- ellip_tb %>% 
  pluck("ellip", 5) %>% 
  pluck("rdim", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD8_memory

```


### 2.1.3 Signature size selection & PCA plot for t_helper cell

Signature size selection plot for t_helper cell:
  
  ```{r t_helper, message=F}
t_helper <- ellip_data %>% 
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

t_helper

```

The elbow point indicates optimal sig_size is $3$. PCA at sig_size = $3$: 
  
  ```{r PCA5_t_helper, message=F, warning = FALSE}
sig_size <- 3

PCA5_t_helper <- ellip_tb %>% 
  pluck("ellip", 3) %>% 
  pluck("rdim", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_t_helper
```

## 2.2 Silhouette Analysis

```{r silhouette, message=F, results = FALSE, warning = FALSE}
sig_size <- 60

sil_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ sig_select(contrast_PW_L5, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL, "PCA")))

sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL)))

```

### 2.2.1 Signature size selection & PCA plot for t_CD4_memory cell

Signature size selection plot for t_CD4_memory cell:
  
  ```{r tCD4_memory_sil, message=F}
tCD4_memory_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

tCD4_memory_sil
```

```{r max_sig1, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $4$. PCA at sig_size = $4$:
  
  ```{r PCA5_tCD4_memory_sil, message=F, warning = FALSE}
sig_size <- 4

PCA5_tCD4_memory_sil <- sil_tb %>% 
  pluck("sil_df", 4) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD4_memory_sil

```

### 2.2.2 Signature size selection & PCA plot for t_CD8_memory cell

Signature size selection plot for t_CD8_memory cell:
  
  ```{r tCD8_memory_sil, message=F}
tCD8_memory_sil <- sil_data %>%
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

tCD8_memory_sil
```

```{r max_sig2, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 2) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 2) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $1$. PCA at sig_size = $1$:
  
  ```{r PCA5_tCD8_memory_sil, message=F, warning = FALSE}
sig_size <- 1

PCA5_tCD8_memory_sil <- sil_tb %>% 
  pluck("sil_df", 1) %>% 
  pluck("rdim", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD8_memory_sil
```

### 2.2.3 Signature size selection & PCA plot for t_helper cell

Signature size selection plot for t_helper cell:
  
  ```{r t_helper_sil, message=F}
t_helper_sil <- sil_data %>%
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

t_helper_sil
```

```{r max_sig3, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 3) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 3) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $3$. PCA at sig_size = $3$:
  
  ```{r PCA5_t_helper_sil, message=F, warning = FALSE}
sig_size <- 3

PCA5_t_helper_sil <- sil_tb %>% 
  pluck("sil_df", 3) %>% 
  pluck("rdim", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_t_helper_sil
```

# 3 Hierarchy + mean contrast analysis

## 3.0 Generate mean contrast

```{r mean_contrast, message=F, warning=F}
contrast_MC_L5 <- tt_L5 %>% 
  contrast_MC(LEVEL)
```

## 3.1 Ellipse Area Analysis

```{r ellipse_3.1, warning = FALSE, message=F}
sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
ellip_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(ellip = map(sig_size, ~ sig_select(contrast_MC_L5, LEVEL, .x))) %>%
  mutate(ellip = map(ellip, ~ ellip_func(.x, LEVEL, "PCA")))

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
ellip_data <- ellip_tb %>% ellip_scale(LEVEL)
```

### 3.1.1 Signature size selection & PCA plot for t_CD4_memory cell

Signature size selection plot for t_CD4_memory cell:
  
  ```{r tCD4_memory_elli_3.1.1, message=F}
tCD4_memory_elli <- ellip_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

tCD4_memory_elli
```

The elbow point indicates optimal sig_size is $4$. PCA at sig_size = $4$: 
  
  ```{r PCA5_tCD4_memory_3.1.1, warning = FALSE, message=F}
sig_size <- 4

PCA5_tCD4_memory <- ellip_tb %>% 
  pluck("ellip", 4) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD4_memory
```

### 3.1.2 Signature size selection & PCA plot for t_CD8_memory cell

Signature size selection plot for t_CD8_memory cell:
  
  ```{r tCD8_memory_3.1.2, warning = FALSE, message=F}
tCD8_memory <- ellip_data %>% 
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

tCD8_memory
```

The elbow point indicates optimal sig_size is $5$. PCA at sig_size = $5$: 
  
  ```{r PCA5_tCD8_memory_3.1.2, warning = FALSE, message=F}
sig_size <- 5

PCA5_tCD8_memory <- ellip_tb %>% 
  pluck("ellip", 5) %>% 
  pluck("rdim", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD8_memory
```

### 3.1.3 Signature size selection & PCA plot for t_helper cell

Signature size selection plot for t_helper cell:
  
  ```{r t_helper_3.1.3, warning = FALSE, message=F}
t_helper <- ellip_data %>% 
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

t_helper
```

The elbow point indicates optimal sig_size is $3$. PCA at sig_size = $3$: 
  
  ```{r PCA5_t_helper_3.1.3, warning = FALSE, message=F}
sig_size <- 3

PCA5_t_helper <- ellip_tb %>% 
  pluck("ellip", 3) %>% 
  pluck("rdim", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_t_helper
```

## 3.2 Silhouette Analysis

```{r silhouette_3.2, warning = FALSE, message=F}
sig_size <- 60

sil_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ sig_select(contrast_MC_L5, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL, "PCA")))

sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL)))
```


### 3.2.1 Signature size selection & PCA plot for t_CD4_memory cell

Signature size selection plot for t_CD4_memory cell:
  
  ```{r tCD4_memory_sil_3.2.1, warning = FALSE, message=F}
tCD4_memory_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

tCD4_memory_sil
```

```{r max_sig_3.2.1, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $4$. PCA at sig_size = $4$:
  
  ```{r PCA5_tCD4_memory_sil_3.2.1, warning = FALSE, message=F}
sig_size <- 4

PCA5_tCD4_memory_sil <- sil_tb %>% 
  pluck("sil_df", 4) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD4_memory_sil
```

### 3.2.2 Signature size selection & PCA plot for t_CD8_memory cell

Signature size selection plot for t_CD8_memory cell:
  
  ```{r tCD8_memory_sil_3.2.2, warning = FALSE, message=F}
tCD8_memory_sil <- sil_data %>%
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

tCD8_memory_sil
```

```{r max_sig_3.2.2, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 2) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 2) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $1$. PCA at sig_size = $1$:
  
  ```{r PCA5_tCD8_memory_sil_3.2.2, warning = FALSE, message=F}
sig_size <- 1

PCA5_tCD8_memory_sil <- sil_tb %>% 
  pluck("sil_df", 1) %>% 
  pluck("rdim", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_tCD8_memory_sil
```

### 3.2.3 Signature size selection & PCA plot for t_helper cell

Signature size selection plot for t_helper cell:
  
  ```{r t_helper_sil_3.2.3, warning = FALSE, message=F}
t_helper_sil <- sil_data %>%
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

t_helper_sil
```

```{r max_sig_3.2.3, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 3) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 3) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $2$. PCA at sig_size = $2$:
  
  ```{r PCA5_t_helper_sil_3.2.3, warning = FALSE, message=F}
sig_size <- 2

PCA5_t_helper_sil <- sil_tb %>% 
  pluck("sil_df", 2) %>% 
  pluck("rdim", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA5_t_helper_sil
```