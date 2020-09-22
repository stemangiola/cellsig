# devtools::install_github("stemangiola/nanny@convert-to-S3")
# devtools::install_github("stemangiola/tidybulk@dev")

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

load("data/counts.rda")




# Setup data frame
tt <- 
  
  # Load dataset
  cellsig::counts %>%
  tidybulk(sample, symbol, count) %>%

  # Group by level because otherwise samples are duplicated
  nest(data = -level) %>%
  
  # Redefine factors inside each level
  mutate(data = future_map(data, ~ droplevels(.x))) %>%
  
  # Fill missing data. There are many genes that
  # are not shared by the majority of samples
  mutate(data = future_map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
  
  # Scale for future PCA plotting
  mutate(data = future_map(
    data, ~ .x %>% 
       identify_abundant(factor_of_interest = cell_type) %>% 
       scale_abundance()
  ))

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

select_markers_for_each_contrast = function(.data){
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
                           slice(1:10)        
    )) %>%
    unnest(stat_df)
}

# Level 1
all_contrasts <-
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
  mutate(markers                = map( markers               , ~ .x %>% select_markers_for_each_contrast)) %>%

  # Add marker info to original data
  mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
  select(markers) %>%
  unnest(markers) %>%
  
  # make contrasts pretty
  mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", "")) 
  

# Plot Markers
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

# Plot PCA
all_contrasts %>%
  distinct(sample, symbol, count_scaled,cell_type) %>%
  nanny::reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="get", transform = log1p) %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + 
  geom_point() +
  theme_bw()


