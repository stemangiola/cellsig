devtools::install_github("stemangiola/nanny@convert-to-S3", force = TRUE)
devtools::install_github("stemangiola/tidybulk@dev", force = TRUE)
#============================================================================================================

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
library(limma)

load("dev/orig_counts.rda")

LEVEL = 3

# Setup data frame
tt <- 
  
  # Load dataset
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

# Inspect cell types on level 2
tt %>% filter(level==LEVEL) %>% unnest(data) %>% select(cell_type) %>% distinct()

tt %>% filter(level==LEVEL) %>% unnest(data)

# No markers
tt_naive <-  
  tt %>%
  
  # Scale and reduce dimensions
  mutate(data = map(
    data,
    ~ .x %>%
      scale_abundance() %>%
      reduce_dimensions(method="PCA") %>% 
      reduce_dimensions(method="MDS") %>%
      reduce_dimensions(method="tSNE")
  )) %>%
  
  # Cluster
  mutate(data = map(data, ~ cluster_elements(.x, method="SNN")))


# For debugging you have 2 options
## call debugonce(get_constrasts_from_df), run the function that calls this function, which is constrast(tt), 
## when get_constrasts_from_df is called, you can just run any line to see the output
## insert "browser()" inside the function, load the function, run the function



# Functions
## 1
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



## 2
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

contrast31 <- c(
  "cell_typeb_cell - (cell_typegranulocyte + cell_typemast_cell + cell_typemono_derived + cell_typenatural_killer + cell_typet_cell)/5",
  "cell_typegranulocyte - (cell_typeb_cell + cell_typemast_cell + cell_typemono_derived + cell_typenatural_killer + cell_typet_cell)/5",
  "cell_typemast_cell - (cell_typeb_cell + cell_typegranulocyte + cell_typemono_derived + cell_typenatural_killer + cell_typet_cell)/5",
  "cell_typemono_derived - (cell_typeb_cell + cell_typegranulocyte + cell_typemast_cell + cell_typenatural_killer + cell_typet_cell)/5",
  "cell_typenatural_killer - (cell_typeb_cell + cell_typegranulocyte + cell_typemast_cell + cell_typemono_derived + cell_typet_cell)/5",
  "cell_typet_cell - (cell_typeb_cell + cell_typegranulocyte + cell_typemast_cell + cell_typemono_derived + cell_typenatural_killer)/5"
)

contrast32 <-  c(
  "cell_typeb_cell - (cell_typegranulocyte + cell_typemast_cell + cell_typemono_derived + cell_typenatural_killer + cell_typet_cell)",
  "cell_typegranulocyte - (cell_typeb_cell + cell_typemast_cell + cell_typemono_derived + cell_typenatural_killer + cell_typet_cell)",
  "cell_typemast_cell - (cell_typeb_cell + cell_typegranulocyte + cell_typemono_derived + cell_typenatural_killer + cell_typet_cell)",
  "cell_typemono_derived - (cell_typeb_cell + cell_typegranulocyte + cell_typemast_cell + cell_typenatural_killer + cell_typet_cell)",
  "cell_typenatural_killer - (cell_typeb_cell + cell_typegranulocyte + cell_typemast_cell + cell_typemono_derived + cell_typet_cell)",
  "cell_typet_cell - (cell_typeb_cell + cell_typegranulocyte + cell_typemast_cell + cell_typemono_derived + cell_typenatural_killer)"
)

contrast <- function(tt, LEVEL){
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
    
    # Select markers from each cell type
    # mutate(markers = map(markers, ~ select_markers_for_each_cell_type(.x, sig_size))) %>% 
    
    # Add original data data to markers
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
    
    # obtain covariance matrix for each cell type
    mutate(cov = map(PC, ~ cov(.x))) %>% # pluck("cov", 1)
    
    # calculate the eigenvalues for the covariance matrix of each cell type
    mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% # pluck("eigval", 1)
    
    # transformation
    mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>% # unnest(area)
    
    # below is the actual area 
    mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% # unnest(area)
    
    # log10 transformation to make values bigger
    # mutate(area = map_dbl(area, ~ 1/ log10(.x))) %>% 
    
    # summation of areas across all cell types
    select(area) %>% 
    sum()
}


# Ellipse area==================================================================================================

sig_size = 60

# for single sig_size calculation
tt %>% 
  contrast(., LEVEL) %>% 
  ellip_area()

# for serial sig_size calculation

area_sum <- 1: sig_size

for (i in 1: sig_size) {
  sig_size = i
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

# when using get_contrasts_from_df(), sig_size =  gives the minimum total_area: ;
# when using contrast1, sig_size = 27  gives the minimum total_area: 8.713;
# when using contrast2, sig_size =  gives the minimum total_area: 

sig_size
sig_size 

# Silhouette analysis=============================================================================

# for single sig_size calculation

# calculate the dissimilarity matrix
distance <- PCA_level3 %>%
  select(contains("PC")) %>% 
  dist()

# calculate silhouette score for each data point
sil <- silhouette(as.numeric(PCA_level3$cell_type), distance)
summary(sil)

# plot all silhouette scores and store the data in sil_data
sil_data <- fviz_silhouette(sil)

# find the average silhouette width (or score) for all data points
sil_score <- sil_data$data$sil_width %>% mean()


# identify samples that have negative silhouette values
neg_sil_index <- which(sil[, "sil_width"] < 0)
neg_sil_index
sil[neg_sil_index, , drop = FALSE]

# for serial sig_size calculation

sig_size = 60

scores <- 1:sig_size

for (i in 1: sig_size){
  sig_size = i
  all_contrasts <- contrast(tt, LEVEL)
  
  PCA_level3 <- all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p)
  
  distance <- PCA_level3 %>% 
    select(contains("PC")) %>% 
    dist()
  
  sil <- silhouette(as.numeric(PCA_level3$cell_type), distance)
  sil_data <- fviz_silhouette(sil)
  scores[i] <- sil_data$data$sil_width %>% mean()
}

# convert scores from ggplot2 object into a vector
sil_scores <- 1: sig_size

for (i in 1: sig_size){
  sil_scores[i] <- scores[[i]]
}

tb2 <- tibble(size = 1 :sig_size, sil_scores = sil_scores )
ggplot(tb2, aes(size, sil_scores)) +
  geom_line() +
  geom_point()

max(sil_scores)
which.max(sil_scores)

# when using get_contrasts_from_df(), sig_size =  gives maximum silhouette score: ; 
# when using contrast1, sig_size =  gives maximum silhouette score: ; 
# when using contrast2, sig_size =  gives maximum silhouette score:

sig_size
sig_size



# Plots==============================================================================================
## Plot Markers
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


## Plot PCA

all_contrasts <- contrast(tt, LEVEL)


PCA_level3 <- all_contrasts %>%
  distinct(sample, symbol, count_scaled, cell_type) %>%
  reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p)

PCA_level3

PCA3 <- PCA_level3 %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3

PCA3 %>% ggplotly(tooltip = c("label", "cell_type"))





# 7 NEW =============

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


```{r data}
# Load data
tt_all <- readRDS("intermediate_data/tt_all.rds")

tt_L3 <- tt_all %>% 
  pluck("tt", 3)
```

```{r literal}
# select level of interest
LEVEL = "level_3"

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
    
    # Add original data data to the markers selected
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
  
  while (any(map_int(signature$markers_cumu, ~ length(.x)) < .target_size) & 
         any(j < .discard_num) & 
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
    
    # option3: log transformation of ratio to apply weight
    # mutate(ratio = ifelse(diff_sil==0, diff_sil, 0.5*log(diff_size) - log(diff_sil) ))
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


# 1 Pre-analysis ==================

## 1.1 Preprocess data
```{r preprocess, message=F, results = FALSE, warning = FALSE}
# 1 Setup data frame & preprocessing

tt_L3 <- preprocess(counts, LEVEL)

saveRDS(tt_L3, "tt_L3.rds")
```

```{r cell_types, message=F, warning=F}
# View cell types on level_1
tt_L3 %>% 
  unnest(data) %>% 
  select(!!as.symbol(pre(LEVEL))) %>% 
  distinct()
```

## 1.2 Cluster without marker selection for contrast

```{r naive, results = FALSE, warning = FALSE, message=F}
# 2 No selection of markers
tt_naive_L3 <-
  tt_L3 %>%
  
  # Scale and reduce dimensions
  mutate(data = map(
    data,
    ~ .x %>%
      reduce_dimensions(method="PCA")
    # reduce_dimensions(method="MDS") %>%
    # reduce_dimensions(method="tSNE")
  ))

```

```{r PCA_naive_mono, warning=F, message=F}
PCA_naive_mono <- tt_naive_L3 %>% 
  pluck("data", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_mono
```

```{r PCA_naive_t, warning=F, message=F}
PCA_naive_t <- tt_naive_L3 %>% 
  pluck("data", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_t
```

```{r PCA_naive_granulocyte, warning=F, message=F}
PCA_naive_granulocyte <- tt_naive_L3 %>% 
  pluck("data", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_granulocyte
```

```{r PCA_naive_b, warning=F, message=F}
PCA_naive_b <- tt_naive_L3 %>% 
  pluck("data", 4) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_b
```

```{r PCA_naive_NK, warning=F, message=F}
PCA_naive_NK <- tt_naive_L3 %>% 
  pluck("data", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_NK
```

# 2 Hierarchy + Pairwise Analysis ============================

## 2.0 Generate pairwise contrast

```{r pairwise_contrast, message=F, warning = FALSE}
contrast_PW_L3 <- tt_L3 %>% 
  contrast_PW(LEVEL)

```

## 2.1 Ellipse Area Analysis ============================

```{r ellipse, message=F, results = FALSE, warning = FALSE,}
sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
ellip_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(ellip = map(sig_size, ~ sig_select(contrast_PW_L3, LEVEL, .x))) %>%
  mutate(ellip = map(ellip, ~ ellip_func(.x, LEVEL, METHOD) ))

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
ellip_data <- ellip_tb %>% ellip_scale(LEVEL)

```

### 2.1.1 Signature size selection & PCA plot for monocyte_derived cell

Signature size selection plot for mono_derived cell:
  
  ```{r mono_elli, warning=F, message=F}
mono_elli <- ellip_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

mono_elli
```

The elbow point indicates optimal sig_size is $5$. PCA at sig_size = $5$: 
  
  ```{r PCA3_mono_elli, message=F, message=F, warning = FALSE}
sig_size <- 5

PCA3_mono_elli <- ellip_tb %>% 
  pluck("ellip", 5) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_mono_elli
```

### 2.1.2 Signature size selection & PCA plot for t cell

Signature size selection plot for t cell:
  
  ```{r t_elli, warning=F, message=F}
t_elli <- ellip_data %>% 
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

t_elli
```

The elbow point indicates optimal sig_size is $4$. PCA at sig_size = $4$: 
  
  ```{r PCA3_t_elli, message=F, warning = FALSE}
sig_size <- 4

PCA3_t_elli <- ellip_tb %>% 
  pluck("ellip", 4) %>% 
  pluck("rdim", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_t_elli
```

### 2.1.3 Signature size selection & PCA plot for granulocyte

Signature size selection plot for granulocyte:
  
  ```{r granulocyte_elli, warning=F, message=F}
granulocyte_elli <- ellip_data %>% 
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

granulocyte_elli
```

The elbow point indicates optimal sig_size is $2$. PCA at sig_size = $2$: 
  
  ```{r PCA3_granulocyte, message=F, warning = FALSE}
sig_size <- 2

PCA3_granulocyte <- ellip_tb %>% 
  pluck("ellip", 2) %>% 
  pluck("rdim", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_granulocyte
```

### 2.1.4 Signature size selection & PCA plot for b cells

Signature size selection plot for b cell:
  
  ```{r b_elli, warning=F, message=F}
b_elli <- ellip_data %>% 
  pluck("plot_data", 4) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

b_elli
```

The elbow point indicates optimal sig_size is $4$. PCA at sig_size = $4$: 
  
  ```{r PCA3_b, message=F, warning = FALSE}
sig_size <- 4

PCA3_b <- ellip_tb %>% 
  pluck("ellip", 4) %>% 
  pluck("rdim", 4) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_b
```

### 2.1.5 Signature size selection & PCA plot for NK cell

Signature size selection plot for NK:
  
  ```{r NK_elli, warning=F, message=F}
NK_elli <- ellip_data %>% 
  pluck("plot_data", 5) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

NK_elli
```

The elbow point indicates optimal sig_size is $4$. PCA at sig_size = $4$: 
  
  ```{r PCA3_NK, message=F, warning = FALSE}
sig_size <- 4

PCA3_NK <- ellip_tb %>% 
  plcuk("ellip", 4) %>% 
  pluck("rdim", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_NK
```

## 2.2 Silhouette Analysis =======================================

```{r silhouette, message=F, results = FALSE, warning = FALSE}
sig_size <- 60

sil_tb_pw <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ sig_select(contrast_PW_L3, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL, METHOD)))

pw_naive_sil_L3 <- sil_tb_pw %>% 
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

saveRDS(pw_naive_sil_L3, "pw_naive_sil_L3.rds")

sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(plot_data, ~ .x %>% ratio())) %>% 
  mutate(optim_size = map(plot_data, ~ .x %>% optim_size()))

```
sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(plot_data, ~ .x %>% penalized_sil())) %>% 
  mutate(optim_size = map(plot_data, ~ .x %>% optim_size2()))

sil_data %>% mutate(optim_size = map_df(optim_size, ~.x %>% unlist()))

### 2.2.1 Signature size selection & PCA plot for monocyte_derived cell

Signature size selection plot for mono_derived cell:
  
  ```{r mono_sil, warning=F, message=F}
mono_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

mono_sil
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

The peak is reached when sig_size = $60$. PCA at sig_size = $60$:
  
  ```{r PCA3_mono_sil, message=F, warning = FALSE}
# sig_size = 60, real_size =

PCA3_mono_sil <- sil_data %>% 
  pluck("plot_data", 1) %>% 
  pluck("rdim", 60) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_mono_sil

# sig_size = 57, real_size = 254
PCA3_mono_sil <- sil_data %>% 
  pluck("plot_data", 1) %>% 
  pluck("rdim", 57) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()
```

### 2.2.2 Signature size selection & PCA plot for t cell

Signature size selection plot for t cell:
  
  ```{r t_sil, warning=F, message=F}
t_sil <- sil_data %>%
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

t_sil
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

The peak is reached when sig_size = $5$. PCA at sig_size = $5$:
  
  ```{r PCA3_t_sil, message=F, warning = FALSE}
# sig_size = 5, real_size = 26

PCA3_t_sil <- sil_data %>% 
  pluck("plot_data", 2) %>% 
  pluck("rdim", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_t_sil
```

### 2.2.3 Signature size selection & PCA plot for granulocyte

Signature size selection plot for granulocyte:
  
  ```{r granulocyte_sil, warning=F, message=F}
granulocyte_sil <- sil_data %>%
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

granulocyte_sil
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

The peak is reached when sig_size = $1$. PCA at sig_size = $1$:
  
  ```{r PCA3_granulocyte_sil, message=F, warning = FALSE}
# sig_size = 1, real_size =

PCA3_granulocyte_sil <- sil_data %>% 
  pluck("plot_data", 3) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_granulocyte_sil

# sig_size = 51, real_size = 102

PCA3_granulocyte_sil <- sil_data %>% 
  pluck("plot_data", 3) %>% 
  pluck("rdim", 51) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()
```

### 2.2.4 Signature size selection & PCA plot for b cell

Signature size selection plot for b cell:
  
  ```{r b_sil, warning=F, message=F}
b_sil <- sil_data %>%
  pluck("plot_data", 4) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

b_sil
```

```{r max_sig4, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 4) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 4) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $4$. PCA at sig_size = $4$:
  
  ```{r PCA3_b_sil, message=F, warning = FALSE}
# sig_size = 4, real_size = 

PCA3_b_sil <- sil_data %>% 
  pluck("plot_data", 4) %>% 
  pluck("rdim", 4) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_b_sil

# sig_size = 7, real_size = 14
PCA3_b_sil <- sil_data %>% 
  pluck("plot_data", 4) %>% 
  pluck("rdim", 7) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()
```

### 2.2.5 Signature size selection & PCA plot for NK

Signature size selection plot for NK:
  
  ```{r NK_sil, warning=F, message=F}
NK_sil <- sil_data %>%
  pluck("plot_data", 5) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

NK_sil
```

```{r max_sig5, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 5) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 5) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $15$. PCA at sig_size = $15$:
  
  ```{r PCA3_NK_sil, message=F, warning = FALSE}
# sig_size = 15, real_size =30

PCA3_NK_sil <- sil_tb %>% 
  pluck("sil_df", 15) %>% 
  pluck("rdim", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_NK_sil
```

# 3 Hierarchy + mean contrast analysis ========================

## 3.0 Generate mean contrast

```{r mean_contrast, message=F, warning=F}
contrast_MC_L3 <- tt_L3 %>% 
  contrast_MC(LEVEL)
```

## 3.1 Ellipse Area Analysis ===================================

```{r ellipse_3.1, warning = FALSE, message=F}
sig_size <- 20

# create a tibble that stores the confidence ellipse area output for each signature size
ellip_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(ellip = map(sig_size, ~ sig_select(contrast_MC_L3, LEVEL, .x))) %>%
  mutate(ellip = map(ellip, ~ ellip_func(.x, LEVEL, METHOD)))

# rescale areas and plot total areas vs the total number of markers selected from cell types in a level
ellip_data <- ellip_tb %>% ellip_scale(LEVEL)
```

### 3.1.1 Signature size selection & PCA plot for monocyte derived cells

Signature size selection plot for mono_derived cell:
  
  ```{r mono_elli_3.1.1, warning=F, message=F}
mono_elli <- ellip_data %>% 
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # scale_x_continuous(sec.axis = sec_axis(as.factor())) +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

mono_elli
```

The elbow point indicates optimal sig_size is $7$. PCA at sig_size = $7$: 
  
  ```{r PCA3_mono_elli_3.1.1, warning = FALSE, message=F}
sig_size <- 7

PCA3_mono_elli <- ellip_tb %>% 
  pluck("ellip", 7) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_mono_elli
```

### 3.1.2 Signature size selection & PCA plot for t cell

Signature size selection plot for t cell:
  
  ```{r t_elli_3.1.2, warning = FALSE, message=F}
t_elli <- area_data %>% 
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

t_elli
```

The elbow point indicates optimal sig_size is $7$. PCA at sig_size = $7$: 
  
  ```{r PCA3_t_elli_3.1.2, warning = FALSE, message=F}
sig_size <- 7

PCA3_t_elli <- ellip_tb %>% 
  pluck("ellip", 7) %>% 
  pluck("rdim", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_t_elli
```

### 3.1.3 Signature size selection & PCA plot for granulocytes

Signature size selection plot for granulocytes:
  
  ```{r granulocyte_elli_3.1.3, warning = FALSE, message=F}
granulocyte_elli <- ellip_data %>% 
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

granulocyte_elli
```

The elbow point indicates optimal sig_size is $2$. PCA at sig_size = $2$: 
  
  ```{r PCA3_granulocyte_3.1.3, warning = FALSE, message=F}
sig_size <- 2

PCA3_granulocyte__elli <- ellip_tb %>% 
  pluck("ellip", 2) %>% 
  pluck("rdim", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_granulocyte__elli
```

### 3.1.4 Signature size selection & PCA plot for b cell

Signature size selection plot for b cells:
  
  ```{r b_elli_3.1.4, warning = FALSE, message=F}
b_elli <- ellip_data %>% 
  pluck("plot_data", 4) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

b_elli
```

The elbow point indicates optimal sig_size is $4$. PCA at sig_size = $4$: 
  
  ```{r PCA3_b_3.1.4, warning = FALSE, message=F}
sig_size <- 4

PCA3_b__elli <- ellip_tb %>% 
  pluck("ellip", 4) %>% 
  pluck("rdim", 4) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_b__elli
```

### 3.1.5 Signature size selection & PCA plot for NK cells

Signature size selection plot for NK:
  
  ```{r NK_elli_3.1.5, warning = FALSE, message=F}
NK_elli <- ellip_data %>% 
  pluck("plot_data", 5) %>% 
  ggplot(aes(real_size, area_value, colour=area_type)) + 
  geom_line() +
  geom_point() +
  # facet_grid(rows = vars(area_type), scales = "free_y")
  facet_wrap(~ area_type, scales = "free_y")

NK_elli
```

The elbow point indicates optimal sig_size is $4$. PCA at sig_size = $4$: 
  
  ```{r PCA3_NK_3.1.5, warning = FALSE, message=F}
sig_size <- 4

PCA3_NK_elli <- ellip_tb %>% 
  pluck("ellip", 4) %>% 
  pluck("rdim", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_NK_elli
```

## 3.2 Silhouette Analysis ======================================

```{r silhouette_3.2, warning = FALSE, message=F}
sig_size <- 60

sil_tb_mc <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ sig_select(contrast_MC_L3, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL, METHOD)))

mc_naive_sil_L3 <- sil_tb_mc %>% 
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

saveRDS(mc_naive_sil_L3, "mc_naive_sil_L3.rds")


sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(plot_data, ~ .x %>% ratio())) %>% 
  mutate(optim_size = map(plot_data, ~ .x %>% optim_size()))
```
sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL))) %>% 
  mutate(plot_data = map(plot_data, ~ .x %>% penalized_sil())) %>% 
  mutate(optim_size = map(plot_data, ~ .x %>% optim_size2()))

sil_data %>% mutate(optim_size = map_df(optim_size, ~.x %>% unlist()))
### 3.2.1 Signature size selection & PCA plot for monocyte derived cells

Signature size selection plot for mono_derived cell:
  
  ```{r mono_sil_3.2.1, warning = FALSE, message=F}
mono_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

mono_sil
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

The peak is reached when sig_size = $29$. PCA at sig_size = $29$, \\
but silhouette score has only minimal improvement than that when sig_size = $17$:
  
  ```{r PCA3_mono_sil_3.2.1, warning = FALSE, message=F}
sig_size <- 29

PCA3_mono_sil <- sil_tb %>%
  pluck("sil_df", 29) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_mono_sil
```

### 3.2.2 Signature size selection & PCA plot for t cells

Signature size selection plot for t cell:
  
  ```{r t_sil_3.2.2, warning = FALSE, message=F}
t_sil <- sil_data %>%
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

t_sil
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

The peak is reached when sig_size = $26$. PCA at sig_size = $26$:
  
  ```{r PCA3_t_sil_3.2.2, warning = FALSE, message=F}
sig_size <- 26

PCA3_t_sil <- sil_tb %>% 
  pluck("sil_df", 26) %>% 
  pluck("rdim", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_t_sil
```

### 3.2.3 Signature size selection & PCA plot for granulocytes

Signature size selection plot for granulocytes:
  
  ```{r granulocyte_sil_3.2.3, warning = FALSE, message=F}
granulocyte_sil <- sil_data %>%
  pluck("plot_data", 3) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

granulocyte_sil
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

The peak is reached when sig_size = $1$. PCA at sig_size = $1$:
  
  ```{r PCA3_granulocyte_sil_3.2.3, warning = FALSE, message=F}
sig_size <- 1

PCA3_granulocyte_sil <- sil_tb %>% 
  pluck("sil_df", 1) %>% 
  pluck("rdim", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_granulocyte_sil
```

### 3.2.4 Signature size selection & PCA plot for b cells

Signature size selection plot for b cells:
  
  ```{r b_sil_3.2.4, warning = FALSE, message=F}
b_sil <- sil_data %>%
  pluck("plot_data", 4) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

b_sil
```

```{r max_sig_3.2.4, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 4) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 4) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $4$. PCA at sig_size = $4$:
  
  ```{r PCA3_b_sil_3.2.4, warning = FALSE, message=F}
sig_size <- 4

PCA3_b_sil <- sil_tb %>% 
  pluck("sil_df", 4) %>% 
  pluck("rdim", 4) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_b_sil
```

### 3.2.5 Signature size selection & PCA plot for granulocytes

Signature size selection plot for NK cells:
  
  ```{r NK_sil_3.2.5, warning = FALSE, message=F}
NK_sil <- sil_data %>%
  pluck("plot_data", 5) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point()

NK_sil
```

```{r max_sig_3.2.5, warning=F, message=F}
sil_data %>% 
  pluck("plot_data", 5) %>%
  pull(sil) %>% 
  which.max()
```

```{r}
sil_data %>% 
  pluck("plot_data", 5) %>%
  pull(sil) %>% 
  max()
```

The peak is reached when sig_size = $15$. PCA at sig_size = $15$:
  
  ```{r PCA3_NK_sil_3.2.5, warning = FALSE, message=F}
sig_size <- 15

PCA3_NK_sil <- sil_tb %>% 
  pluck("sil_df", 15) %>% 
  pluck("rdim", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_NK_sil
```

# New pairwise selection method

## selecting a single marker at a time with all base markers in each iteration =====================

# pw_markers_L3 gives the full output including PC values from reduced_dimensions but this file is too big
pw_markers_L3 <- contrast_PW_L3 %>% 
  single_marker_pw_select(LEVEL, .discard_num = 1000, METHOD)

saveRDS(pw_markers_L3, "pw_markers_L3.rds")

mc_markers_L3 <- contrast_MC_L3 %>% 
  single_marker_pw_select(LEVEL, .discard_num = 1000, METHOD)

saveRDS(mc_markers_L3, "mc_markers_L3.rds")

# trend plots & PCA ==================
# trend plot mono derived
pw_markers_L3 %>% 
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
  ggtitle("mono_derived, 171 marker selected")


# PCA plot mono derived
pw_markers_L3 %>% 
  pluck("signature_data", 1) %>% 
  tail(1) %>%
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA mono_derived, 171 markers selected, silhouette = 0.855")


# trend plot t_cell
pw_markers_L3 %>% 
  pluck("signature_data", 2) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1) %>% 
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("t_cell, 38 marker selected")

## PCA plot t_cell
pw_markers_L3 %>% 
  pluck("signature_data", 2) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA t_cell, 38 markers selected, sil_score=0.787")

# trend plot granulocyte
pw_markers_L3 %>% 
  pluck("signature_data", 3) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1) %>%
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("granulocyte, 3 marker selected")


## PCA plot granulocyte
pw_markers_L3 %>% 
  pluck("signature_data", 3) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA granulocyte, 3 markers selected, sil_score=0.808")

# trend plot b_cell
pw_markers_L3 %>% 
  pluck("signature_data", 4) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1) %>%
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("b_cell, 43 marker selected")

## PCA plot b_cell
pw_markers_L3 %>% 
  pluck("signature_data", 4) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA b cell, 43 markers selected, sil_score=0.855")

# trend plot NK
pw_markers_L3 %>% 
  pluck("signature_data", 5) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1) %>%
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("nk_cell, 61 marker selected")

## PCA plot NK
pw_markers_L3 %>% 
  pluck("signature_data", 5) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA nk cell, 61 markers selected, sil_score=0.92")

# pw_markers_L3.2 trend plots & PCA ==================

# trend plot mono derived
pw_markers_L3.2 %>% 
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
  ggtitle("mono_derived, 167 marker selected")


# PCA plot mono derived
pw_markers_L3.2 %>% 
  pluck("signature_data", 1) %>% 
  tail(1) %>%
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA mono_derived, 167 markers selected, silhouette = 0.859")


# trend plot t_cell
pw_markers_L3.2 %>% 
  pluck("signature_data", 2) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("t_cell, 97 marker selected")

## PCA plot t_cell
pw_markers_L3.2 %>% 
  pluck("signature_data", 2) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA t_cell, 97 markers selected, sil_score=0.71")

# trend plot granulocyte
pw_markers_L3.2 %>% 
  pluck("signature_data", 3) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("granulocyte, 11 marker selected")


## PCA plot granulocyte
pw_markers_L3.2 %>% 
  pluck("signature_data", 3) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA granulocyte, 11 markers selected, sil_score=0.689")

# trend plot b_cell
pw_markers_L3.2 %>% 
  pluck("signature_data", 4) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("b_cell, 46 marker selected")

## PCA plot b_cell
pw_markers_L3.2 %>% 
  pluck("signature_data", 4) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA b cell, 46 markers selected, sil_score=0.836")

# trend plot NK
pw_markers_L3.2 %>% 
  pluck("signature_data", 5) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("nk_cell, 61 marker selected")

## PCA plot NK
pw_markers_L3.2 %>% 
  pluck("signature_data", 5) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA nk cell, 61 markers selected, sil_score=0.92")

