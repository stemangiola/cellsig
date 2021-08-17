# How to call this file
# From bash
# Rscript database_plus_tree_as_input_benchmark_plot_as_output.R
#
# Output -> a PDF file
# Input nothing
#
# What does it do
# Produce a pdf running the whole pipeline, including all methods, and testing with sihuette and deconvolution


# Silhouette

silhouette_score <- function(.reduced_dimensions, .distance, .level){
  
  .reduced_dimensions %>% 
    
    pull(!!as.symbol(.level)) %>% 
    
    as.factor() %>% 
    
    as.numeric() %>% 
    
    silhouette(.distance) %>% 
    
    summary()
  
}

silhouette_function <- function(.selected, .reduction_method){
  
  .selected %>% 
    
    # reduce dimensions
    mutate(reduced_dimensions = map2(
      markers, level, 
      ~ dimension_reduction(.x, .y, .reduction_method)
    )) %>% 
    
    # calculate distance matrix using PC1 & PC2
    mutate(distance = map(
      reduced_dimensions,
      ~ distance_matrix(.x, .reduction_method)
    )) %>% 
    
    # calculate silhouette score
    mutate(silhouette = pmap(
      list(reduced_dimensions, distance, level),
      ~ silhouette_score(..1, ..2, ..3)
    )) %>% 
    
    # remove unnecessary columns
    select(-c(markers, distance))
  
}

yy <- full_df %>% 
  nest(signature = -method) %>% 
  mutate(signature = map(signature, ~.x %>% pull(signature) %>% unlist() %>% unique())) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score for all signatures combined in each method
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

cibersortx <- readRDS("dev/topInf_scaleFALSE/cibersortx.new.rds")
cibersort_signature <- cibersortx$signature[[1]]

cibersortx <- tibble(method = "cibersortx") %>% 
  mutate(signature = list(cibersort_signature)) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

yy <- yy %>% bind_rows(cibersortx)

yy %>% 
  mutate(cluster.silhouette = map(silhouette, ~ .x$clus.avg.widths)) %>% 
  mutate(avg.silhouette = map_dbl(silhouette, ~ .x$avg.width)) %>% 
  select(-c(reduced_dimensions, silhouette)) %>% 
  unnest(cluster.silhouette) %>% 
  ggplot(aes(x=reorder(method, avg.silhouette), y=cluster.silhouette, colour=method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  theme(axis.text.x = element_blank())



all_methods_silhouette <- full_df %>%
  mutate(data = map(data, ~ if ("cumulative_signature" %in% names(.x))
  {rename(.x, signature=cumulative_signature)}else(.x))) %>% 
  nest(signature = -method) %>% 
  mutate(signature = map(signature, ~ .x %>% do_optimisation("penalised", 0.4))) %>% 
  mutate(signature = map(signature, ~.x %>% pull(signature) %>% unlist() %>% unique())) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score for all signatures combined in each method
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

cibersortx <- readRDS("topInf_scaleFALSE/cibersortx.new.rds")

all_methods_comparison <- all_methods_silhouette %>% 
  bind_rows(cibersortx) %>% 
  arrange(desc(silhouette))

all_methods_comparison

all_methods_comparison %>% 
  ggplot(aes(reorder(method, silhouette), silhouette, fill = method)) +
  geom_col() +
  geom_text(aes(label = round(silhouette, 3)), vjust = 1.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("All methods comparison using silhouette score, penalty_rate=1.2")
