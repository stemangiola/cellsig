# dendritic have no genes that have FDR < 0.05, min(FDR)= 0.1538621
debugonce(select_markers_for_each_contrast)

contrast_PW_L4 %>% 
  pluck("markers", 4) %>% 
  select_markers_for_each_contrast(1)

y <- contrast_PW_L4 %>% 
  pluck("markers", 4) %>% 
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
  mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value)))

# no genes with FDR < 0.05, minimum is 0.154
y %>% 
  pluck("stat_df", 1) %>% 
  filter(FDR < 0.05)

y %>% 
  pluck("stat_df", 1) %>% 
  filter(logFC > 2)

y %>% 
  pluck("stat_df", 1) %>% 
  select(FDR) %>% 
  min()

y %>% 
  pluck("stat_df", 1) %>% 
  select(FDR) %>% 
  max()
