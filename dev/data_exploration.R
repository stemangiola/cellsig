devtools::install_github("stemangiola/nanny@convert-to-S3", force = TRUE)
devtools::install_github("stemangiola/tidybulk@dev", force = TRUE)

library(tidyverse)
library(plotly)
library(nanny)
library(ggrepel)
library(GGally)
library(tidyHeatmap)
library(furrr)
plan(multiprocess, workers=5)

# To be loaded after all libraries
library(tidybulk)
library(cluster)
library(proxy)
library(factoextra)
load("data/counts.rda")




# Setup data frame
tt <- 
  
  # Load dataset
  counts %>%
  tidybulk(sample, symbol, count) %>%

  # Group by level because otherwise samples are duplicated
  nest(data = -level) %>%
  
  # Redefine factors inside each level
  mutate(data = future_map(data, ~ droplevels(.x))) %>%
  
  # Fill missing data. There are many genes that
  # are not shared by the majority of samples
  mutate(data = future_map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
  
  # Scale for future PCA plotting
  mutate(data = future_map(data, ~ .x %>% identify_abundant() %>% scale_abundance()))

# No markers
tt_naive <-  
  tt %>%

  # Scale and reduce dimensions
  mutate(data = future_map(
    data,
    ~ .x %>%
      scale_abundance() %>%
      reduce_dimensions(method="MDS") %>%
      reduce_dimensions(method="PCA") %>%
      reduce_dimensions(method="tSNE")
  )) %>%

  # Cluster
  mutate(data = future_map(data,~ cluster_elements(.x, method="SNN")))


# Functions
## 1
get_constrasts_from_df = function(.data){
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
select_markers_for_each_contrast = function(.data, sigsize){
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
                           slice(1:sigsize)     
    )) %>%
    unnest(stat_df)
}


## 3
contrast <- function(tt){
  tt %>%
    
    # Investigate one level
    filter(level==1) %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = get_constrasts_from_df(.x),
                                    action="only"
      ) 
      
    )) %>%
    
    # Select rank from each contrast
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, sigsize))) %>%
    
    # Add marker info to original data
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
}

## 4 calculate the area of confidence ellipses and the sum of their areas
ellip_area <- function(all_contrasts){
  all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    nanny::reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p) %>% 
    select(-sample) %>% 
    nest(PC = -cell_type) %>% # pluck("PC", 1) %>% select(contains('PC')) %>% cov()
    mutate(cov = map(PC, ~ cov(.x))) %>% # pluck("cov", 1)
    mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% # pluck("eigval", 1)
    mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>% # unnest(area)
    # one line below is the actual area, followed by log10 transformation to make values bigger
    mutate(area = map(area, ~ prod(.x)*pi)) %>% # unnest(area)
    mutate(area = map_dbl(area, ~ 1/-log10(.x))) %>% 
    pluck("area") %>% 
    sum()
}


# Ellipse area==================================================================================================
sigsize = 17
all_contrasts <- contrast(tt)


PCA_level1 <- all_contrasts %>%
  distinct(sample, symbol, count_scaled, cell_type) %>%
  nanny::reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p)

PCA_level1

n = 20
area_sum <- 1:n

for (i in 1:n) {
  sigsize = i
 area_sum[i] <- tt %>% 
   contrast() %>% 
   ellip_area()
}

tb <- tibble(size = 1:n, total_area = area_sum)
ggplot(tb, aes(size, total_area)) +
  geom_line() +
  geom_point()

which.min(area_sum)  # 17 genes gives the minimum area_sum 

# Silhouette analysis=============================================================================

distance <- PCA_level1 %>%
  select(contains("PC")) %>% 
  dist()


sil <- silhouette(as.numeric(PCA_level1$cell_type), distance)
summary(sil)

sil_data <- fviz_silhouette(sil)
sil_score <- sil_data$data$sil_width %>% mean()

sil_scores <- 1:n

# identify samples that have negative silhouette values
neg_sil_index <- which(sil[, "sil_width"] < 0)
neg_sil_index
sil[neg_sil_index, , drop = FALSE]

for (i in 1:n){
  sigsize = i
  all_contrasts <- contrast(tt)
  PCA_level1 <- all_contrasts %>% 
    distinct(sample, symbol, count_scaled, cell_type) %>%
    nanny::reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p)
  distance <- PCA_level1 %>% 
    select(contains("PC")) %>% 
    dist()
  sil <- silhouette(as.numeric(PCA_level1$cell_type), distance)
  sil_data <- fviz_silhouette(sil)
  sil_scores[i] <- sil_data$data$sil_width %>% mean()
}

scores <- 1:n
for (i in 1:n){
  scores[i] <- sil_scores[[i]]
}
tb2 <- tibble(size = 1:n, scores = scores )
ggplot(tb2, aes(size, scores)) +
  geom_line() +
  geom_point()

max(scores)
which.max(scores) # 9 genes leads to maximum silhouette score

sigsize = 9
all_contrasts <- contrast(tt)


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
PCA1 <- PCA_level1 %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type, label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw()

PCA1

PCA1 %>% ggplotly(tooltip = c("label", "cell_type"))


# Incorrect====================================================================================================
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

  mutate(a_dis = map_dbl(data, ~ sqrt(.x$PC1 - .x$x_center)^2 + (.x$PC2 - .x$y_center)^2)) %>% 