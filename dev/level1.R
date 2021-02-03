#==========================================================================================
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
library(limma)

# Load data========================================================================================================

load("data/counts.rda")
new_counts <- counts

load("dev/orig_counts.rda")
orig_counts <- counts

rm(counts)

# Literal parameters===============================================================================================

LEVEL = 1

SIG_SIZE

# Functions=========================================================================================================

## 1 preprocess data

preprocess <- function(counts){
  
  # convert to tidybulk object
  counts %>%
    tidybulk(sample, symbol, count) %>%
    
    # Group by level because otherwise samples are duplicated
    nest(data = -level) %>%
    
    # Redefine factors inside each level
    mutate(data = map(data, ~ droplevels(.x))) %>%
    
    # Fill missing data. There are many genes that
    # are not shared by the majority of samples
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

contrast <- function(tt, LEVEL){
  tt %>%
    
    # Investigate one level
    filter(level == LEVEL) %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = make_contrasts(tt, LEVEL),
                                    action="only") 
    )) %>%
    
    # Select rank from each contrast
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, sig_size))) %>%
    
    # Select markers from each cell type
    # mutate(markers = map(markers, ~ select_markers_for_each_cell_type(.x, sig_size))) %>% 
    
    # Add original data info to markers
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
}


## 5 calculate the area of confidence ellipses and the sum of their areas

ellip_area <- function(all_contrasts){
  # reduce dimension
  all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p) %>% 
    
    # remove non-numerical data to form a numerical dataframe
    select(-sample) %>% 
    
    # normalise principle component values
    mutate(across(c("PC1", "PC2"), ~ (.x - mean(.x))/sd(.x) )) %>% 
    
    # nest by cell_type so as to calculate ellipse area for each cell type
    nest(PC = -cell_type) %>% # pluck("PC", 1) %>% select(contains('PC')) %>% cov()
    
    # # alternatively normalise PC values within clusters to standardise the variability within clusters contributed by sample size
    # mutate(PC = map_depth(PC, 2, ~ (.x - mean(.x))/sd(.x) )) %>%
    # mutate(PC = map(PC, ~ as_tibble(.x))) %>%
    
    # obtain covariance matrix for each cell type
    mutate(cov = map(PC, ~ cov(.x))) %>% # pluck("cov", 1)
    
    # calculate the eigenvalues for the covariance matrix of each cell type
    mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% # pluck("eigval", 1)
    
    # transformation
    mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>% # unnest(area)
    
    # below is the actual area for each ellipse
    mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% # unnest(area)
    
    # # calculate the mean area score
    # pull(area) %>% 
    # 
    # mean()
    
    # collect size of each cluster as factors for weights
    mutate(c_size = map_int(PC, ~ nrow(.x))) %>%
    
    # weight each area by the inverse of its cluster size
    mutate(weighted_area = map2_dbl(area, c_size, ~ .x / .y)) %>%
    
    # calculated the weighted sum of ellipse areas for clusters
    pull(weighted_area) %>%
    
    sum()
}

# For debugging you have 2 options
## call debugonce(get_constrasts_from_df), run the function that calls this function, which is constrast(tt), 
## when get_constrasts_from_df is called, you can just run any line to see the output
## insert "browser()" inside the function, load the function, run the function

# Analysis===============================================================================================

# 1 Setup data frame, pre-processing data================================================================

tt <- preprocess(orig_counts)

# 2 Inspect cell types on level 1========================================================================

tt %>% filter(level==LEVEL) %>% unnest(data) %>% select(cell_type) %>% distinct()

tt %>% filter(level==LEVEL) %>% unnest(data)

# 3 No selection of markers==================================================================================
tt_naive <-  
  tt %>%
  
  # Scale and reduce dimensions
  mutate(data = map(
    data,
    ~ .x %>%
      identify_abundant(factor_of_interest = cell_type) %>% 
      scale_abundance() %>%
      reduce_dimensions(method="PCA") %>% 
      reduce_dimensions(method="MDS") %>%
      reduce_dimensions(method="tSNE")
  )) %>%
  
  # Cluster
  mutate(data = map(data, ~ cluster_elements(.x, method="SNN")))

# 4 Ellipse area==================================================================================================

# 4.1 Single sig_size calculation

sig_size <- 60

tt %>% 
  contrast(., LEVEL) %>% 
  ellip_area()

# 4.2 Serial sig_size calculation

sig_size <- 60

area_sum <- 1: sig_size

for (i in 1: sig_size) {
  sig_size = i
  print(sig_size)
  area_sum[i] <- tt %>% 
    contrast(., LEVEL) %>% 
    ellip_area()
}

tb <- tibble(size = 1: sig_size, total_area = area_sum)
ggplot(tb, aes(size, total_area)) +
  geom_line() +
  geom_point()

which.min(area_sum)

# when using get_contrasts_from_df(), sig_size = 57 gives the minimum total_area: 4.364;
# when using contrast1, sig_size = 60 gives the minimum total_area: 6.005;
# when using contrast2, sig_size = 34 gives the minimum total_area: 5.570 (NaN for sig_size = 1 !)

min(area_sum)

sig_size <- 57
sig_size <- 60
sig_size <- 34

# 5 Silhouette analysis===============================================================================

## 5.1 Single sig_size calculation
sig_size <- 50

all_contrasts <- contrast(tt, LEVEL)

PCA_level1 <- all_contrasts %>% 
  distinct(sample, symbol, count_scaled, cell_type) %>%
  reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p)

# calculate the dissimilarity matrix
distance <- PCA_level1 %>%
  select(contains("PC")) %>% 
  dist()

# calculate silhouette score for each data point
sil <- silhouette(as.numeric(PCA_level1$cell_type), distance)
sil_info <- summary(sil)

# find the average silhouette width (or score) for all data points
avg_score <- sil_info$avg.width
avg_score

clus_size <- sil_info$clus.sizes
clus_avg_widths <- sil_info$clus.avg.widths
sum(clus_avg_widths/clus_size)

# plot all silhouette scores and store the data in sil_data
sil_data <- fviz_silhouette(sil)
# visualise
sil_data


# identify samples that have negative silhouette values
neg_sil_index <- which(sil[, "sil_width"] < 0)
neg_sil_index
sil[neg_sil_index, , drop = FALSE]

## 5.2 Serial sig_size calculation

sig_size <- 60

sil_scores <- 1:sig_size

for (i in 1: sig_size){
  sig_size = i
  print(sig_size)
  all_contrasts <- contrast(tt, LEVEL)
  
  PCA_level1 <- all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p)
  
  distance <- PCA_level1 %>% 
    select(contains("PC")) %>% 
    dist()
  
  sil <- silhouette(as.numeric(PCA_level1$cell_type), distance)
  sil_info <- summary(sil)
  sil_scores[i] <- sil_info$avg.width
}

tb2 <- tibble(size = 1 :sig_size, sil_scores = sil_scores )
ggplot(tb2, aes(size, sil_scores)) +
  geom_line() +
  geom_point()

max(sil_scores)
which.max(sil_scores) 

# when using get_contrasts_from_df(), sig_size = 30 gives maximum silhouette score: 0.729; 
# when using contrast1, sig_size = 50 gives maximum silhouette score: 0.706
# when using contrast2, sig_size = 31 gives maximum silhouette score: 0.678

sig_size <- 30
sig_size <- 50
sig_size <- 31

# 6 Plots==============================================================================================

## 6.1 Plot Markers

all_contrasts %>% 
  ggplot(aes(cell_type, count_scaled + 1, color=contrast_pretty)) + 
  geom_point(size = 0.5) + 
  facet_wrap(~contrast_pretty + symbol, scale="free_y") + 
  scale_y_log10() +
  theme_bw() +
  theme(
    text = element_text(size=6), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) 


## 6.2 Plot PCA

all_contrasts <- contrast(tt, LEVEL)


PCA_level1 <- all_contrasts %>%
  distinct(sample, symbol, count_scaled, cell_type) %>%
  reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p) # %>% 
  
  # plot PCA for PC values normalised within clusters
  # group_by(cell_type) %>% 
  # mutate(across(c("PC1", "PC2"), ~ (.x - mean(.x))/sd(.x) ))
  
  # plot PCA for PC values normalised globally
  # mutate(across(c("PC1", "PC2"), ~ (.x - mean(.x))/sd(.x) )) %>% 
  
  # plot PCA for PC values normalised first globally then within clusters
  # group_by(cell_type) %>%
  # mutate(across(c("PC1", "PC2"), ~ (.x - mean(.x))/sd(.x) ))

PCA_level1

PCA1 <- PCA_level1 %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA1

PCA1 %>% ggplotly(tooltip = c("label", "cell_type"))


# 7 new ====================================================================================================


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
load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")
```

```{r literal}
# select level of interest
LEVEL = "level_1"
```



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
    
    # create an ancestor node for cell types on level_1
    mutate(level_0 = "cell") %>% 
    
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



# 1 Pre-analysis

## 1.1 Preprocess data
```{r preprocess, message=F, warning = FALSE}
# 1 Setup data frame & preprocessing

tt <- preprocess(counts, LEVEL)

```

```{r cell_types, message=F, warning=F}

# View cell types on ancestor level
tt %>% 
  unnest(data) %>% 
  select(!!as.symbol(pre(LEVEL))) %>% 
  distinct()

```

## 1.2 Cluster without marker selection for contrast

```{r naive, results = FALSE, warning = FALSE, message=F}
# 2 No selection of markers
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

```

```{r PCA_cell, message=F}
PCA_cell <- tt_naive %>% 
  pluck("data", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_cell
```

# 2 Hierarchy + Pairwise Analysis

## 2.1 Ellipse Area Analysis


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



### 2.1.1 Signature size selection & PCA plot for cell

# Signature size selection plot for cell:
  

cell_elli <- area_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

cell_elli



The elbow point indicates optimal sig_size is $7$. PCA at sig_size = $7$: 
  

sig_size <- 7

PCA_level1 <- tt %>% 
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

PCA1_cell_elli <- PCA_level1 %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA1_cell_elli


## 2.2 Silhouette Analysis


sig_size <- 20

sil_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ contrast(tt, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL)))


sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL)))




### 2.2.1 Signature size selection & PCA plot for cell

Signature size selection plot for cell:
  
  ```{r cell_sil, message=F}
cell_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

cell_sil
```


The peak is reached when sig_size = $4$. PCA at sig_size = $4$:
  
  ```{r PCA5_tCD4_memory_sil, message=F, warning = FALSE}
sig_size <- 14

PCA_level1 <- tt %>% 
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

PCA1_cell_sil <- PCA_level1 %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA1_cell_sil

```

# 3 Hierarchy + mean contrast analysis

```{r mean_contrast, message=F}
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
```

## 3.1 Ellipse Area Analysis

```{r ellipse_3.1, warning = FALSE, message=F}
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

```

### 3.1.1 Signature size selection & PCA plot for cell

Signature size selection plot for t_CD4_memory cell:
  
  ```{r cell_elli_3.1.1, message=F}
cell_elli <- area_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

cell_elli
```


The elbow point indicates optimal sig_size is $4$. PCA at sig_size = $4$: 
  
  ```{r PCA5_tCD4_memory_3.1.1, warning = FALSE, message=F}
sig_size <- 4

PCA_level1 <- tt %>% 
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

PCA1_cell_elli <- PCA_level5 %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA1_cell_elli
```

## 3.2 Silhouette Analysis

```{r silhouette_3.2, warning = FALSE, message=F}
sig_size <- 20

sil_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ contrast(tt, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL)))


sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL)))

```


### 3.2.1 Signature size selection & PCA plot for cell

Signature size selection plot for t_CD4_memory cell:
  
  ```{r cell_sil_3.2.1, warning = FALSE, message=F}
cell_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

cell_sil
```

The peak is reached when sig_size = $4$. PCA at sig_size = $4$:
  
  ```{r PCA1_cell_sil_3.2.1, warning = FALSE, message=F}
sig_size <- 14

PCA_level1 <- tt %>% 
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

PCA1_cell_sil <- PCA_level5 %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA1_cell_sil

```



