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

# Load data=================================================================================================

load("data/counts.rda")
new_counts <- counts

load("orig_counts.rda")
orig_counts <- counts

rm(counts)

# Literal parameters========================================================================================

LEVEL = 2

SIG_SIZE

# Functions=================================================================================================

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

contrast <- function(tt, LEVEL, sig_size){
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
    
    # Add original data info to markers
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
}


## 5 calculate the area of confidence ellipses and the sum of their areas
PCA_level1

ellip_area <- function(all_contrasts){
  # reduce dimension
  PCA <- all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="add", transform = log1p)
  
  real_size <- PCA %>% 
    nest(data = -cell_type) %>% 
    mutate(real_size = map_int(data, ~ n_distinct(.x$symbol)))
  
  area <- PCA %>%   
    # remove non-numerical data to form a numerical dataframe
    select(cell_type, PC1, PC2) %>%
    
    # normalise principle component values
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

# For debugging you have 2 options
## call debugonce(get_constrasts_from_df), run the function that calls this function, which is constrast(tt), 
## when get_constrasts_from_df is called, you can just run any line to see the output
## insert "browser()" inside the function, load the function, run the function

# Analysis=========================================================================================================

# 1 Setup data frame & preprocessing==============================================================================
tt <- preprocess(orig_counts)


# 2 Inspect cell types on level 2=================================================================================
tt %>% filter(level==LEVEL) %>% unnest(data) %>% select(cell_type) %>% distinct()

tt %>% filter(level==LEVEL) %>% unnest(data)

# check how many genes there are for each cell type
tt %>% filter(level == LEVEL) %>% unnest(data) %>% group_by(cell_type) %>% select(symbol) %>% count()

sig_size <- 60

all_contrasts <- contrast(tt, LEVEL)

all_contrasts %>% group_by(cell_type) %>% select(symbol) %>% count()

all_contrasts %>% group_by(contrast) %>% select(symbol) %>% count()

all_contrasts %>% group_by(contrast, cell_type) %>% select(symbol) %>% count()

# 3 No selection of markers========================================================================================
tt_naive <-  
  tt %>%
  
  # Scale and reduce dimensions
  mutate(data = map(
    data,
    ~ .x %>%
      # identify_abundant(factor_of_interest = cell_type) %>% 
      # scale_abundance() %>%
      reduce_dimensions(method="PCA") %>% 
      reduce_dimensions(method="MDS") %>%
      reduce_dimensions(method="tSNE")
  )) 
  # %>% 
  # # Cluster
  # mutate(data = map(data, ~ cluster_elements(.x, method="SNN")))

PCA_naive <- tt_naive %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive

# 4 Ellipse area==================================================================================================

## 4.1 Single sig_size calculation
sig_size <- 60

tt %>% 
  contrast(., LEVEL) %>% 
  ellip_area()

## 4.2 Serial sig_size calculation

sig_size <- 60

# create a tibble that stores the confidence ellipse area output for each signature size
area_tb <- tibble(sig_size = 1:sig_size) %>%
  # slice(1:3) %>% 
  mutate(area_df = map(
    sig_size, ~ contrast(tt, LEVEL, .x) %>% ellip_area()
  ))

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
area_tb %>% 
  unnest(area_df) %>% 
  nest(cell_data = -cell_type) %>% 
  mutate(cell_data = map(cell_data, ~ .x %>% mutate(rescaled_area = area %>% scale(center = F)))) %>% 
  unnest(cell_data) %>%
  # select(real_size) %>% 
  # max()
  group_by(real_size) %>% 
  summarise(stded_sum=sum(area, na.rm = T), 
            wted_sum = sum(weighted_area, na.rm = T), 
            rescaled_sum= sum(rescaled_area, na.rm = T)) %>% 
  pivot_longer(c(stded_sum, wted_sum, rescaled_sum), names_to='area_type', values_to="area_value") %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")


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
  contrast(., LEVEL) %>% 
  ellip_area()
}

tb <- tibble(size = 1: sig_size, total_area = area_sum)
ggplot(tb, aes(size, total_area)) +
  geom_line() +
  geom_point()

which.min(area_sum)

min(area_sum)

# when using get_contrasts_from_df(), sig_size = 44 gives the minimum total_area: 8.840;

# when using contrast1, sig_size = 27  gives the minimum total_area: 8.713; 
# (when using contrast1, weighted sum of area, sig_size = 58 gives the minimum:0.2029, but 58 out of range? 
# the PCA plots 336 most variable genes, which is when sig_size = 56.
# second minimum when sig_size <- 28, gives 0.204)

# when using contrast2, sig_size = 39  gives the minimum total_area: 8.014

sig_size <- 44
sig_size <- 27
sig_size <- 58
sig_size <- 39

# 5 Silhouette analysis=============================================================================

## 5.1 Single sig_size calculation
sig_size <- 11

all_contrasts <- contrast(tt, LEVEL, sig_size)

PCA_level2 <- all_contrasts %>% 
  distinct(sample, symbol, count_scaled, cell_type) %>%
  reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p)

# calculate the dissimilarity matrix
distance <- PCA_level2 %>%
  select(contains("PC")) %>% 
  dist()

# calculate silhouette score for each data point
sil <- silhouette(as.numeric(PCA_level2$cell_type), distance)
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

## 5.2 serial sig_size calculation

sig_size <- 60

sil_scores <- 1:sig_size

for (i in 1: sig_size){
  sig_size = i
  print(sig_size)
  all_contrasts <- contrast(tt, LEVEL, sig_size)
  
  PCA_level2 <- all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p)
  
  distance <- PCA_level2 %>% 
    select(contains("PC")) %>% 
    dist()
  
  sil <- silhouette(as.numeric(PCA_level2$cell_type), distance)
  sil_info <- summary(sil)
  sil_scores[i] <- sil_info$avg.width
}

tb2 <- tibble(size = 1 :sig_size, sil_scores = sil_scores )
ggplot(tb2, aes(size, sil_scores)) +
  geom_line() +
  geom_point()

max(sil_scores)
which.max(sil_scores) 

# when using get_contrasts_from_df(), sig_size = 48 gives maximum silhouette score: 0.349; 
# when using contrast1, sig_size = 5 gives maximum silhouette score: 0.335 (try sig_size=7, 58) 
# when using contrast2, sig_size = 43 gives maximum silhouette score: 0.358

sig_size <- 48
sig_size <- 5
sig_size <- 43


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
 
all_contrasts <- contrast(tt, LEVEL, sig_size)


PCA_level2 <- all_contrasts %>%
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

PCA_level2

PCA2 <- PCA_level2 %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA2

PCA2 %>% ggplotly(tooltip = c("label", "cell_type"))


# 7 Unused Materials====================================================================================================
nest(data = -cell_type) %>%
  mutate(x_center = map_dbl(data, ~ mean(.x$PC1) )) %>%
  mutate(y_center = map_dbl(data, ~ mean(.x$PC2) )) %>%
  
  
  unnest(data)


sil <- function(.data) {
  PC <- select(.data %>% starts_with("PC"))
  distance <- dist(PC)
  silhouette(.data, distance)
}


cluster_elements(sample, symbol, count_scaled, method = "kmeans", centers = 4) %>% 
  ggplot() +
  geom_point(aes(x= PC1, y= PC2, color=cell_type)) 

mutate(a_dis = map_dbl(data, ~ sqrt(.x$PC1 - .x$x_center)^2 + (.x$PC2 - .x$y_center)^2))

# 3
select_markers_for_each_cell_type <- function(.data, sig_size) {
  .data %>%
    # split the contrast column into two so as to group all the contrasts into 4 groups based on cell_type
    separate(col = contrast, into = c("contrast1", "contrast2"), sep = "-", remove = F) %>% 
    
    # remove white spaces before names, make contrast easy to operate on
    mutate(across(c("contrast1", "contrast2"), ~ trimws(.x))) %>% 
    
    group_by(contrast) %>%
    arrange( logFC %>% desc) %>%
    slice(1:sig_size) 
  %>%
    
    
    # group by cell_type, which is in the contrast column
    # nest(stat_df = -contrast) %>%
    group_by(contrast1) %>% 
    dplyr::distinct(symbol, .keep_all = T) %>% 
    group_by(contrast2) %>% 
    arrange(logFC %>% desc) %>% 
    slice(1:sig_size) %>% 
    
    
    # mutate(dup = map(stat_df, ~ .x %>% duplicated(symbol))) %>% 
    # 
    # mutate(stat_df = map(stat_df, ~ x. %>% subset(dup == T))) %>% 
    
    
    # remove any duplicate genes within one cell_type that from different contrasts
    mutate(stat_df = map(stat_df, ~ .x %>% dplyr::distinct(symbol, .keep_all = T))) %>% 
    
    # mutate(stat_df = map(stat_df, ~ .x %>% slice(1:sig_size))) %>% 
    
    mutate(stat_df = map(stat_df, ~ .x %>% nest(stat_df2 = -contrast2))) %>% 
    
    
    unnest(stat_df)
}

# # alternatively normalise PC values within clusters to standardise the variability within clusters contributed by sample size
# mutate(PC = map_depth(PC, 2, ~ (.x - mean(.x))/sd(.x) )) %>%
# mutate(PC = map(PC, ~ as_tibble(.x))) %>%
# ...
# # calculate the mean area score
# pull(area) %>% 
# mean()

