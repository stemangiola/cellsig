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


# 7 NEW ====================================================================================================
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
LEVEL = "level_3"
```


```{r functions, message=F}
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

```


# 1 Pre-analysis

## 1.1 Preprocess data
```{r preprocess, message=F, results = FALSE, warning = FALSE}
# 1 Setup data frame & preprocessing

tt <- preprocess(counts, LEVEL)

```

```{r cell_types, message=F, warning=F}

# View cell types on level_1
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

```{r PCA_naive_tCD4_memory, message=F}
PCA_naive_mono <- tt_naive %>% 
  pluck("data", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_mono
```

```{r PCA_naive_t, message=F}
PCA_naive_t <- tt_naive %>% 
  pluck("data", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_t
```

```{r PCA_granuo, message=F}
PCA_naive_granulocyte <- tt_naive %>% 
  pluck("data", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_granulocyte
```
PCA_naive_b <- tt_naive %>% 
  pluck("data", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_b

PCA_naive_NK <- tt_naive %>% 
  pluck("data", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA_naive_NK


# 2 Hierarchy + Pairwise Analysis

## 2.1 Ellipse Area Analysis

```{r ellipse, message=F, results = FALSE, warning = FALSE,}
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

### 2.1.1 Signature size selection & PCA plot for monocyte_derived cell

Signature size selection plot for t_CD4_memory cell:
  
  ```{r mono_elli, warning=F, message=F}
mono_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_mono_elli <- PCA_level3 %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_mono_elli
```

### 2.1.2 Signature size selection & PCA plot for t cell

Signature size selection plot for t cell:
  
  ```{r t_elli, warning=F, message=F}
t_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_t_elli <- PCA_level3 %>% 
  pluck("pca", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_t_elli

```


### 2.1.3 Signature size selection & PCA plot for granulocyte

Signature size selection plot for granulocyte:
  
  ```{r granulocyte_elli, warning=F message=F}
granulocyte_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_granulocyte <- PCA_level3 %>% 
  pluck("pca", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_granulocyte
```

### 2.1.4 Signature size selection & PCA plot for b cells

Signature size selection plot for b cell:
  
  ```{r b_elli, warning=F, message=F}
b_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_b <- PCA_level3 %>% 
  pluck("pca", 4) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_b
```

### 2.1.5 Signature size selection & PCA plot for NK cell

Signature size selection plot for NK:
  
  ```{r NK_elli, warning=F message=F}
NK_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_NK <- PCA_level3 %>% 
  pluck("pca", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_NK
```


## 2.2 Silhouette Analysis

```{r silhouette, message=F, results = FALSE, warning = FALSE}
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


### 2.2.1 Signature size selection & PCA plot for monocyte_derived cell

Signature size selection plot for mono_derived cell:
  
  ```{r mono_sil, warning=F, message=F}
mono_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

mono_sil
```


The peak is reached when sig_size = $11$. PCA at sig_size = $11$:
  
  ```{r PCA3_mono_sil, message=F, warning = FALSE}
sig_size <- 11

PCA_level3 <- tt %>% 
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

PCA3_mono_sil <- PCA_level3 %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_mono_sil

```

### 2.2.2 Signature size selection & PCA plot for t cell

Signature size selection plot for t cell:
  
  ```{r t_sil, warning=F, message=F}
t_sil <- sil_data %>%
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

t_sil
```

The peak is reached when sig_size = $5$. PCA at sig_size = $5$:
  
  ```{r PCA3_granulocyte_sil, message=F, warning = FALSE}
sig_size <- 5

PCA_level3 <- tt %>% 
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

PCA3_t_sil <- PCA_level3 %>% 
  pluck("pca", 2) %>% 
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
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

granulocyte_sil
```

The peak is reached when sig_size = $1$. PCA at sig_size = $1$:
  
  ```{r PCA3_granulocyte_sil, message=F, warning = FALSE}
sig_size <- 1

PCA_level3 <- tt %>%
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

PCA3_granulocyte_sil <- PCA_level3 %>% 
  pluck("pca", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_granulocyte_sil
```

### 2.2.4 Signature size selection & PCA plot for b cell

Signature size selection plot for b cell:
  
  ```{r b_sil, warning=F, message=F}
b_sil <- sil_data %>%
  pluck("plot_data", 4) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

b_sil
```

The peak is reached when sig_size = $4$. PCA at sig_size = $4$:
  
  ```{r PCA3_b_sil, message=F, warning = FALSE}
sig_size <- 4

PCA_level3 <- tt %>%
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

PCA3_b_sil <- PCA_level3 %>% 
  pluck("pca", 4) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_b_sil
```

### 2.2.5 Signature size selection & PCA plot for NK

Signature size selection plot for NK:
  
  ```{r NK_sil, warning=F, message=F}
NK_sil <- sil_data %>%
  pluck("plot_data", 5) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

NK_sil
```

The peak is reached when sig_size = $15$. PCA at sig_size = $15$:
  
  ```{r PCA3_granulocyte_sil, message=F, warning = FALSE}
sig_size <- 15

PCA_level3 <- tt %>%
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

PCA3_NK_sil <- PCA_level3 %>% 
  pluck("pca", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_NK_sil
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

### 3.1.1 Signature size selection & PCA plot for monocyte derived cells

Signature size selection plot for mono_derived cell:
  
  ```{r mono_elli_3.1.1, warning=F, message=F}
mono_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_mono_elli <- PCA_level3 %>% 
  pluck("pca", 1) %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_t_elli <- PCA_level3 %>% 
  pluck("pca", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_t_elli

```


### 3.1.3 Signature size selection & PCA plot for granulocytes

Signature size selection plot for granulocytes:
  
  ```{r granulocyte_elli_3.1.3, warning = FALSE, message=F}
granulocyte_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_granulocyte__elli <- PCA_level3 %>% 
  pluck("pca", 3) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_granulocyte__elli
```

### 3.1.4 Signature size selection & PCA plot for b cell

Signature size selection plot for b cells:
  
  ```{r b_elli_3.1.4, warning = FALSE, message=F}
b_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_b__elli <- PCA_level3 %>% 
  pluck("pca", 4) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_b__elli
```

### 3.1.5 Signature size selection & PCA plot for NK cells

Signature size selection plot for NK:
  
  ```{r NK_elli_3.1.5, warning = FALSE, message=F}
NK_elli <- area_data %>% 
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

PCA_level3 <- tt %>% 
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

PCA3_NK_elli <- PCA_level3 %>% 
  pluck("pca", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_NK_elli
```


## 3.2 Silhouette Analysis

```{r silhouette_3.2, warning = FALSE, message=F}
sig_size <- 30

sil_tb <-
  tibble(sig_size = 1:sig_size) %>%
  # slice(1) %>%
  mutate(sil_df = map(sig_size, ~ contrast(tt, LEVEL, .x))) %>%
  mutate(sil_df = map(sil_df, ~.x %>% sil_func(LEVEL)))


sil_data <- sil_tb %>%
  unnest(sil_df) %>%
  nest(plot_data = - !!as.symbol(pre(LEVEL)))

```


### 3.2.1 Signature size selection & PCA plot for monocyte derived cells

Signature size selection plot for mono_derived cell:
  
  ```{r mono_sil_3.2.1, warning = FALSE, message=F}
mono_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

mono_sil
```


The peak is reached when sig_size = $29$. PCA at sig_size = $29$, \\
but silhouette score has only minimal improvement than that when sig_size = $17$:
  
  ```{r PCA3_mono_sil_3.2.1, warning = FALSE, message=F}
sig_size <- 17

PCA_level3 <- tt %>% 
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

PCA3_mono_sil <- PCA_level3 %>% 
  pluck("pca", 1) %>% 
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
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

t_sil
```

The peak is reached when sig_size = $26$. PCA at sig_size = $26$:

  ```{r PCA3_t_sil_3.2.2, warning = FALSE, message=F}
sig_size <- 26

PCA_level3 <- tt %>% 
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

PCA3_t_sil <- PCA_level3 %>% 
  pluck("pca", 2) %>% 
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
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

granulocyte_sil
```

The peak is reached when sig_size = $3$. PCA at sig_size = $1$:
  
  ```{r PCA3_granulocyte_sil_3.2.3, warning = FALSE, message=F}
sig_size <- 1

PCA_level3 <- tt %>%
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

PCA3_granulocyte_sil <- PCA_level3 %>% 
  pluck("pca", 3) %>% 
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
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

b_sil
```

The peak is reached when sig_size = $4$. PCA at sig_size = $4$:
  
  ```{r PCA3_b_sil_3.2.4, warning = FALSE, message=F}
sig_size <- 4

PCA_level3 <- tt %>%
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

PCA3_b_sil <- PCA_level3 %>% 
  pluck("pca", 4) %>% 
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
  ggplot(aes(real_size, sil_score)) +
  geom_line() +
  geom_point()

NK_sil
```

The peak is reached when sig_size = $15$. PCA at sig_size = $15$:
  
  ```{r PCA3_NK_sil_3.2.5, warning = FALSE, message=F}
sig_size <- 15

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

PCA3_NK_sil <- PCA_level3 %>% 
  pluck("pca", 5) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA3_NK_sil
```
