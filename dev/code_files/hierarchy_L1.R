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
tt_all <- readRDS("dev/intermediate_data/tt_all.rds")

tt_L1 <- tt_all %>% 
  pluck("tt", 1)

L1 <- tt_L1 %>% 
  select("data") %>% 
  unnest(data) %>% 
  group_by(level_1, symbol) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(level_1) %>% 
  summarise(min=min(n))

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








# 7 NEW ===================================================================================
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
library(scales)
```


```{r data, message=FALSE, warning=FALSE}
# Load data
tt_all <- readRDS("intermediate_data/tt_all.rds")

tt_L1 <- tt_all %>% 
  pluck("tt", 1)
```

```{r literal}
# select level of interest
LEVEL = "level_1"

METHOD = "PCA"
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
    
    # Differential transcription gives all the statistics for each gene in each contrast
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

# ====================================================================

# calculates silhoutte score for each new marker selected for each contrast
sil_for_each_contrast <- function(.sig_select_data, .level, .method) {
  
  .sig_select_data %>% 
    
    # nest by pairwise contrast instead of ancestor level before
    nest(markers = - contrast_pretty) %>% 
    
    # calculate silhouette score for PCA of cell types resolved by the selected 1 markers
    mutate(sil_pair = map(markers, ~ sil_func(.x, .level, .method))) %>% 
    
    # make the silhouette score explicit
    mutate(sil_pair = map_dbl(sil_pair, ~ .x$sil)) %>% 
    
    # rank silhouette score in a descending manner
    arrange(desc(sil_pair)) %>% 
    
    # extract all the new markers from each contrast
    mutate(markers_new = map(markers, ~ .x %>% pull(symbol) %>% unique()))
}

# pairwise selection 1 marker at a time
single_marker_pw_select <- 
  
  function(.contrast, .level, .target_size=NULL, .discard_num=NULL, .method="PCA") {
    
    # initialize variables
    
    contrast_copy <- .contrast %>% 
      mutate(ancestor = !!as.symbol(pre(.level)))
    
    # initialise a signature tibble to store signature markers for each cell type in each iteration
    signature <- tibble(
      ancestor = contrast_copy %>% 
        pull(ancestor)
    ) %>% 
      mutate(markers_cumu = map(ancestor, ~ vector()))
    
    condition <- 
      purrr::when(
        !is.null(.target_size) ~ any(map_int(signature$markers_cumu, ~ length(.x)) < .target_size),
        !is.null(.discard_num) ~ any(j < .discard_num),
        ~ stop("please supply either target signature size or the numer of genes to be discarded")
      )
    
    # initialise an output tibble containing all results of interest
    contrast_summary_tb <- vector()
    
    # set the base markers
    contrast_pair_tb0 <- 
      
      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(.level, 1) %>% 
      
      mutate(ancestor = !!as.symbol(pre(.level))) %>% 
      
      nest(data = - ancestor) %>% 
      
      # enables markers from each contrast pair to be entitled to every other contrast so that 
      # all the base genes rather than one base gene per contrast can be incorporated in each iteration
      mutate(data = map(data, ~ .x %>% 
                          nest(markers = - contrast_pretty) %>% 
                          expand(contrast_pretty, markers) %>% 
                          unnest(markers)
      )) %>% 
      # base markers
      mutate(base_markers = map(data, ~ .x %>% pull(symbol) %>% unique() ))
    
    # remove base markers from contrast_copy input before further selection
    contrast_copy <- contrast_copy %>%
      left_join(contrast_pair_tb0, by= "ancestor", suffix=c("", ".y")) %>% 
      select(-ends_with(".y")) %>% 
      mutate(markers = map2(markers, base_markers, ~.x %>% 
                              filter(!symbol %in% .y) ))
    
    sil_score <- tibble(ancestor = contrast_copy %>% 
                          pull(ancestor)) %>% 
      mutate(sil_pre = map_dbl(ancestor, ~ 0))
    
    # counter for number of genes discarded
    j <- map_int(contrast_pair_tb0$base_markers, ~ length(.x))
    
    while (condition & 
           (!identical(map_int(contrast_copy$markers, ~ nrow(.x)), integer(0))) ) {
      
      contrast_pair_tb <- 
        
        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        contrast_copy %>% 
        
        # select top 1 markers from each contrast, output is an unnested tibble
        sig_select(.level, 1) %>% 
        
        mutate(ancestor = !!as.symbol(pre(.level))) %>% 
        
        # combine markers at the base level before iteration
        bind_rows(contrast_pair_tb0 %>% 
                    select(- base_markers) %>% 
                    unnest(data) ) %>% 
        
        nest(data = - ancestor) %>% 
        
        mutate(data = map(data, ~ .x %>% sil_for_each_contrast(.level, .method) )) %>% 
        
        mutate(is_sil_greater = 
                 map2_lgl(data, ancestor, 
                          ~ if (.x$sil_pair[1]> with(sil_score, sil_pre[ancestor==.y])) 
                          {TRUE} else {FALSE} )) %>% 
        
        mutate(sil_current = map_dbl(data, ~ .x[[1, "sil_pair"]] )) %>% 
        
        mutate(markers_to_filter = map2(data, is_sil_greater, ~
                                          if(.y == TRUE) {
                                            .x %>% 
                                              # select the first marker (with highest silhouette score)
                                              slice(1) %>% 
                                              pull(markers_new) %>% 
                                              unlist()
                                          } else {
                                            .x %>% 
                                              # select all the unique markers from all contrasts for removal
                                              pull(markers_new) %>% 
                                              unlist() %>% 
                                              unique()
                                          } )) %>% 
        # if we want the marker, keep contrast information, else, we don't care if the contrast is correct
        mutate(contrast = map_chr(data, ~ .x$contrast_pretty[[1]])) %>% 
        
        mutate(markers_cumu = map(ancestor, ~ vector())) %>% 
        
        select(-data)
      
      sil_score <- sil_score %>% 
        mutate(sil_pre = map2_dbl(sil_pre, ancestor, ~ 
                                    if (with(contrast_pair_tb, is_sil_greater[ancestor == .y]) )
                                    {with(contrast_pair_tb, sil_current[ancestor == .y])} else{.x} ))
      
      # append the base + 1 markers that result in highest silhouette score
      signature <- signature %>% 
        mutate(markers_cumu = map2(markers_cumu, ancestor, ~ 
                                     if (with(contrast_pair_tb, is_sil_greater[ancestor == .y]) ) 
                                     {.x %>% 
                                         append(
                                           with(contrast_pair_tb, markers_to_filter[ancestor == .y][[1]])
                                         ) %>% 
                                         unique()} else {.x} ))
      
      contrast_summary_tb <- contrast_summary_tb %>% 
        append(
          contrast_pair_tb %>% 
            mutate(markers_cumu = map2(markers_cumu, ancestor,
                                       ~ .x %>% 
                                         append(
                                           with(signature, markers_cumu[ancestor==.y][[1]])
                                         ))) %>% 
            nest(data = - is_sil_greater) %>% 
            with(data[is_sil_greater==TRUE])
        )
      
      contrast_copy <- contrast_copy %>% 
        mutate(markers = map2(markers, ancestor,
                              ~ if(with(contrast_pair_tb, is_sil_greater[ancestor==.y])) {
                                .x %>% 
                                  filter(!symbol %in% with(signature, markers_cumu[ancestor==.y][[1]]))
                              } else {
                                .x %>% 
                                  filter(!symbol %in% with(contrast_pair_tb, markers_to_filter[ancestor==.y][[1]]))
                              } 
        ))
      
      # number of genes discarded for each node
      j <- j + 
        (map_int(contrast_pair_tb$markers_to_filter, ~ length(.x)) -
           map_int(contrast_pair_tb0$base_markers, ~ length(.x)) ) * 
        (!contrast_pair_tb$is_sil_greater)
      
      cat("genes discarded for each node: ", j, "\n")
      cat("genes selected for each node: ", map_int(signature$markers_cumu, ~ length(.x)),  "\n")
    }
    
    return(contrast_summary_tb %>% 
             bind_rows() %>% 
             nest(sig_df = - ancestor))
  }

# calculates silhouette score for each set of signature (cumulative markers) at a signature size
sil_score_for_markers <-function(.sig_df, .contrast, .level, .method) {
  
  .contrast %>%
    
    mutate(ancestor = !!as.symbol(pre(.level))) %>% 
    
    left_join(.sig_df, by="ancestor") %>% 
    
    unnest(sig_df) %>% 
    
    # filter markers that are in the signature
    mutate(markers = map2(markers, markers_cumu, 
                          ~.x %>% 
                            filter(symbol %in% .y))) %>% 
    
    # format statistics from pairwise contrast
    mutate(markers  = map(markers, 
                          ~ .x %>% 
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
                            mutate(stat_df = map(stat_df, ~.x %>% 
                                                   pivot_wider(names_from = stats, 
                                                               values_from = .value))) %>% 
                            
                            unnest(stat_df) )) %>% 
    
    # Add original data data to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>% 
    
    # select only columns needed
    select(-data) %>% 
    
    nest(sil_df = c(!!as.symbol(pre(.level)), markers)) %>% 
    
    mutate(sil_df = map(sil_df, ~ .x %>% unnest(markers))) %>% 
    
    # Calculate silhouette score for PCA plot resulted from the markers selected
    mutate(sil_df = map(sil_df, ~ sil_func(.x, .level, .method))) %>% 
    
    mutate(sil = map_dbl(sil_df, ~ .x$sil)) %>% 
    mutate(real_size = map_int(sil_df, ~ .x$real_size)) %>% 
    
    nest(data = - ancestor)
}

# =======================================================================
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

# signature size optimation functions:
ratio <- function(.plot_data) {
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    filter(real_size > 10) %>%
    
    # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
    mutate(rs_rescaled = rescale(real_size, c(0, 1))) %>% 
    
    # calculate the difference between rescaled size of the max sil point with that of all other points
    # the bigger different the better (even for negative numbers)
    mutate(diff_size = rs_rescaled[which.max(sil)] - rs_rescaled) %>%
    
    # calculate the difference between silhouette score of the max sil point and that of all other points
    # the smaller the better
    mutate(diff_sil = max(sil) - sil) %>%
    
    # calculate the ratio between diff_size/diff_sil, the bigger the better
    mutate(ratio = ifelse(diff_sil==0, diff_sil, diff_size / diff_sil))
}

optim_size <- function(.plot_data) {
  
  # if the highest silhouette score is the first point then select this point (no left points)
  if (which.max(.plot_data$sil) == 1) {
    op_size <- c(sig_size = .plot_data$sig_size[1], 
                 real_size = .plot_data$real_size[1])
    
    # if the highest silhouette score is the last point then select the best left point (no right points)
  } else if (which.max(.plot_data$sil) == nrow(.plot_data)) {
    
    # find the index of the point that gives the best ratio to the left of max silhouette score point
    lop_index <- which.max(.plot_data$ratio[.plot_data$ratio > 0])
    
    op_size <- c(sig_size = .plot_data$sig_size[lop_index],
                 real_size = .plot_data$real_size[lop_index])
    
  } else {
    
    # find the index that gives the optimal point to the left of max silhouette score point
    lop_index <- which.max(.plot_data$ratio[.plot_data$ratio > 0])
    
    # find the index that gives the optimal point to the right of max silhouette score point
    rop_index <- which.max(.plot_data$sil) +
      which.max(.plot_data$ratio[.plot_data$ratio < 0])
    
    # choose the sizes of the optimal left point if its silhouette score is bigger than the optimal right point
    if (.plot_data$sil[lop_index] >= .plot_data$sil[rop_index]) {
      op_size <- c(sig_size = .plot_data$sig_size[lop_index],
                   real_size = .plot_data$real_size[lop_index])
    } else {
      # choose the max silhouette point if silhouette score of the optim right point is bigger than that of the left
      op_size <- c(sig_size = .plot_data$sig_size[which.max(.plot_data$sil)],
                   real_size = .plot_data$real_size[which.max(.plot_data$sil)])
    }
  }
  return(op_size)
}

penalized_sil <- function(.plot_data) {
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    filter(real_size > 10) %>%
    
    # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
    mutate(rs_rescaled = rescale(real_size, c(0, 1))) %>% 
    
    mutate(penalized_sil = sil - 0.1 * rs_rescaled)
}

optim_size2 <- function(.plot_data) {
  index <- which.max(.plot_data$penalized_sil)
  op_size <- c(sig_size = .plot_data$sig_size[index], 
               real_size = .plot_data$real_size[index])
  return(op_size)
}
```

# 1 Pre-analysis =============================

## 1.1 Preprocess data
```{r preprocess, message=F, warning = FALSE}
# 1 Setup data frame & preprocessing

tt_L1 <- counts %>%
  # create an ancestor node for cell types on level_1
  mutate(level_0 = "cell") %>%
  preprocess(LEVEL)

saveRDS(tt_L1, "tt_L1.rds")
```

```{r cell_types, message=F, warning=F}
# View cell types on ancestor level_0
tt_L1 %>% 
  unnest(data) %>% 
  select(!!as.symbol(pre(LEVEL))) %>% 
  distinct()

```

## 1.2 Cluster without marker selection for contrast

```{r naive, results = FALSE, warning = FALSE, message=F}
# 2 No selection of markers
tt_naive_L1 <-
  tt_L1 %>%
  
  # Scale and reduce dimensions
  mutate(data = map(
    data,
    ~ .x %>%
      reduce_dimensions(method="PCA") %>%
      reduce_dimensions(method="MDS") %>%
      reduce_dimensions(method="tSNE")
  ))

```

```{r PCA_naive_cell, message=F}
PCA_naive_cell <- tt_naive_L1 %>% 
  pluck("data", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_cell
```

# 2 Hierarchy + Pairwise Analysis =======================

## 2.0 Generate pairwise contrast

```{r pairwise_contrast, message=F, results = FALSE, warning = FALSE}
contrast_PW_L1 <- tt_L1 %>% 
  contrast_PW(LEVEL)

```

## 2.1 Ellipse Area Analysis============

```{r ellipse, message=F, results = FALSE, warning = FALSE,}
sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
ellip_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(ellip = map(sig_size, ~ sig_select(contrast_PW_L1, LEVEL, .x))) %>%
  mutate(ellip = map(ellip, ~ ellip_func(.x, LEVEL, METHOD) ))

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level

ellip_data <- ellip_tb %>% 
  ellip_scale(LEVEL)

```

### 2.1.1 Signature size selection & PCA plot for cell

Signature size selection plot for cell:
  
  ```{r cell_elli, message=F}
cell_elli <- ellip_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

cell_elli
```

The elbow point indicates optimal sig_size is $7$. PCA at sig_size = $7$: 
  
  ```{r PCA1_cell_elli, message=F, warning = FALSE}
sig_size <- 7

PCA1_cell_elli <- ellip_tb %>% 
  pluck("ellip", 7) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA1_cell_elli
```

## 2.2 Silhouette Analysis===============

```{r silhouette, message=F, warning = FALSE}
sig_size <- 60

sil_tb_pw <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ sig_select(contrast_PW_L1, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL, METHOD)))

pw_naive_sil_L1 <- sil_tb_pw %>% 
  unnest(sil_df) %>%
  mutate(ancestor = !!as.symbol(pre(LEVEL))) %>% 
  mutate(cumulative_signature = map(
    rdim,
    ~ .x %>% 
      pluck("data_symbol", 1) %>% 
      pull(symbol)
  )) %>% 
  mutate(level = LEVEL) %>% 
  select(level, ancestor, sig_size, real_size, silhouette = sil, cumulative_signature) %>% 
  nest(plot_data = - c(level, ancestor))

saveRDS(pw_naive_sil_L1, "pw_naive_sil_L1.rds")


sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(plot_data, ~ .x %>% ratio()) ) %>% 
  mutate(optim_size = map(plot_data, ~ .x %>% optim_size()))
```
sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(plot_data, ~ .x %>% penalized_sil()) ) %>% 
  mutate(optim_size = map(plot_data, ~ .x %>% optim_size2()))

sil_data %>% mutate(optim_size = map_df(optim_size, ~.x %>% unlist()))

### 2.2.1 Signature size selection & PCA plot for cell

Signature size selection plot for cell:
  
  ```{r cell_sil, warning=FALSE, message=F}
cell_sil <- sil_data %>%
  # select ancestor cell type
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point() +
  ggtitle("level_1 cell PW")

cell_sil

ggsave("cell_sil_PW.png", cell_sil)
```

```{r max_sig, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  which.max()
```

sil_data %>% 
  pluck("plot_data", 1) %>%
  slice(44)

The peak is reached when sig_size = $44$. PCA at sig_size = $44$:
  
  ```{r PCA5_tCD4_memory_sil, message=F, warning = FALSE}
# sig_size = 44; real_size = 

PCA1_cell_sil <- sil_data %>% 
  # select ancestor cell type
  pluck("plot_data", 1) %>% 
  # select data at optimal size
  pluck("rdim", 44) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("level_1 cell sig_size = 44")

PCA1_cell_sil

ggsave("level1_cell_size_44.png", PCA1_cell_sil)

# best sig_size = 39, real_size = 342 here

sil_data %>% 
  # select ancestor cell type
  pluck("plot_data", 1) %>% 
  # select data at optimal size
  pluck("rdim", 39) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("level_1 cell sig_size = 39")

```

# 3 Hierarchy + mean contrast analysis ============

## 3.0 Generate mean contrast

```{r mean_contrast, message=F, warning=F}
contrast_MC_L1 <- tt_L1 %>% 
  contrast_MC(LEVEL)
```

## 3.1 Ellipse Area Analysis==============

```{r ellipse_3.1, warning = FALSE, message=F}
sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
ellip_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(ellip = map(sig_size, ~ sig_select(contrast_MC_L1, LEVEL, .x))) %>%
  mutate(ellip = map(ellip, ~ ellip_func(.x, LEVEL, METHOD)))

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
ellip_data <- ellip_tb %>% ellip_scale(LEVEL)

```

### 3.1.1 Signature size selection & PCA plot for cell

Signature size selection plot for cell:
  
  ```{r cell_elli_3.1.1, message=F}
cell_elli <- ellip_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

cell_elli
```

The elbow point indicates optimal sig_size is $7$. PCA at sig_size = $7$: 
  
  ```{r PCA5_tCD4_memory_3.1.1, warning = FALSE, message=F}
sig_size <- 7

PCA1_cell_elli <- ellip_tb %>% 
  pluck("ellip", 7) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA1_cell_elli
```

## 3.2 Silhouette Analysis==================

```{r silhouette_3.2, warning = FALSE, message=F}
sig_size <- 60

sil_tb_mc <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ sig_select(contrast_MC_L1, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL, METHOD)))

mc_naive_sil_L1 <- sil_tb_mc %>% 
  unnest(sil_df) %>%
  mutate(ancestor = !!as.symbol(pre(LEVEL))) %>% 
  mutate(cumulative_signature = map(
    rdim,
    ~ .x %>% 
      pluck("data_symbol", 1) %>% 
      pull(symbol)
  )) %>% 
  mutate(level = LEVEL) %>% 
  select(level, ancestor, sig_size, real_size, silhouette = sil, cumulative_signature) %>% 
  nest(plot_data = - c(level, ancestor))

saveRDS(mc_naive_sil_L1, "mc_naive_sil_L1.rds")

sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(plot_data, ~ .x %>% ratio()) ) %>% 
  mutate(optim_size = map(plot_data, ~ .x %>% optim_size()))

```
sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(plot_data, ~ .x %>% penalized_sil()) ) %>% 
  mutate(optim_size = map(plot_data, ~ .x %>% optim_size2()))

sil_data %>% mutate(optim_size = map_df(optim_size, ~.x %>% unlist()))

### 3.2.1 Signature size selection & PCA plot for cell

Signature size selection plot for cell:
  
  ```{r cell_sil_3.2.1, warning = FALSE, message=F}
cell_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point() +
  ggtitle("cell_sil_MC")

cell_sil

ggsave("cell_sil_MC.png", cell_sil)
```

```{r max_sig, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  which.max()
```
sil_data %>% 
  pluck("plot_data", 1) %>%
  pull(sil) %>% 
  max()

The peak is reached when sig_size = $60$. PCA at sig_size = $60$:
  
  ```{r PCA1_cell_sil_3.2.1, warning = FALSE, message=F}
# sig_size = 60 real_size = 240

PCA1_cell_sil <- sil_data %>% 
  pluck("plot_data", 1) %>% 
  pluck("rdim", 60) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("level1_cell_MC_sig_size = 60")

PCA1_cell_sil

# here optimal sig_size = 49, real_size = 196

PCA1_cell_sil <- sil_data %>% 
  pluck("plot_data", 1) %>% 
  pluck("rdim", 49) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("level1_cell_MC_sig_size = 49")

ggsave("level1_cell_MC_sig_size_49.png", PCA1_cell_sil)

```
# New pairwise signature selection method 

## selecting a single marker at a time with all base markers in each iteration

# data output from new method (base marker = 1 from each contrast, no reduced_dimensions) ============

# pw_markers_L1 gives the output without PC values from reduced_dimensions cuz the file would be too big otherwise
pw_markers_L1 <- contrast_PW_L1 %>% 
  single_marker_pw_select(LEVEL, .discard_num = 1000, METHOD)

saveRDS(pw_markers_L1, "pw_markers_L1.rds")

mc_markers_L1 <- contrast_MC_L1 %>% 
  single_marker_pw_select(LEVEL, .discard_num = 1000, METHOD)

saveRDS(mc_markers_L1, "mc_markers_L1.rds")

# pw_markers_L1 trend plots & PCA

# trend plot cell
pw_markers_L1 %>% 
  pluck("signature_data", 1) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(10)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("cell, 186 marker selected")

ggsave("cell.png", x)


# PCA plot cell
pw_markers_L1 %>% 
  slice(1) %>% 
  unnest(signature_data) %>% 
  tail(1) %>%
  mutate(silhouette_data = map2(
    cumulative_signature, ancestor,
    ~ silhouette_for_markers(.x, .y, contrast_PW_L1, LEVEL, METHOD)
  )) %>% 
  pluck("silhouette_data", 1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA cell, 186 markers selected, silhouette = 0.874")

# data output from new method (base marker = 2 from each contrast, no reduced_dimensions) ===============

# pw_markers_L1.2 gives the output without PC values from reduced_dimensions cuz it'd be too big otherwise
pw_markers_L1.2 <- contrast_PW_L1 %>% 
  single_marker_pw_select(LEVEL, .discard_num = 1000, METHOD)


saveRDS(pw_markers_L1.2, "pw_markers_L1.2.rds")

# pw_markers_L1.2 trend plots & PCA

# trend plot cell
pw_markers_L1.2 %>% 
  pluck("signature_data", 1) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("cell, 209 marker selected")


# PCA plot cell
pw_markers_L1.2 %>% 
  unnest("signature_data") %>% 
  tail(1) %>%
  mutate(silhouette_data = map2(
    cumulative_signature, ancestor,
    ~ silhouette_for_markers(.x, .y, contrast_PW_L1, LEVEL, METHOD)
  )) %>% 
  pluck("silhouette_data", 1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA cell, 209 markers selected, silhouette = 0.870")

