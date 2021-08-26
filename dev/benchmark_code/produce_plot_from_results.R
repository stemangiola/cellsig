argvs

input_directory = arg[1]

dir(input_director) %>%
  map(~ .x %>% readRDS(.x)) %>%
  map(.... to produce the plot)


cibersortx <- 
read_delim("dev/jian_R_files/cibersortx/CIBERSORTx_Job21_phenoclass_1.CIBERSORTx_Job21_reference_1.bm.K999.txt", 
           "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  pull(NAME) %>% 
  list() %>% 
  tibble(stream = "cibersortx", signature = .)



evaluation <- function(.signature, .reduction_method, .preprocessed_non_hierarchy){
  
  # modified silhouette_function and silhouette_score for evaluation
  
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
  
  .preprocessed_non_hierarchy %>% 
    unnest(tt) %>% 
    unnest(data) %>% 
    filter(symbol %in% .signature) %>% 
    nest(markers = -c(level, ancestor)) %>% 
    # calculate silhouette score for all signatures combined in each method
    silhouette_function(.reduction_method = .reduction_method) %>% 
    select(reduced_dimensions, silhouette)
}

x <- 
  signature_hierarchial_methods %>%
  mutate(stream = "mean_contrast_edgR_logFC_naive_penalty") %>%
  select(benchmark, stream) %>%
  mutate(signature = map(benchmark, ~ .x$signature %>% unlist() %>% unique())) %>%
  bind_rows(cibersortx) %>% 
  mutate(silhouette = map(signature, ~ evaluation(.x, "PCA", counts_non_hierarchy))) %>% 
 
  unnest(silhouette) %>% 
  select(-reduced_dimensions) %>% 
  mutate(cluster_silhouette = map(silhouette, ~ .x$clus.avg.widths)) %>% 
  mutate(avg_silhouette = map_dbl(silhouette, ~ .x$avg.width)) %>% 
  unnest(cluster_silhouette) %>% 
  
  ggplot(aes(x=reorder(stream, avg_silhouette), y=cluster_silhouette, colour=stream)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  theme(axis.text.x = element_blank())

  




