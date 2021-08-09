# This file show tSNE error at the first and second node on level_3

# Generate input data
sil_df3 <- contrast_all %>% 
  slice(3) %>% 
  mutate(markers = map2(contrast_MC, level, ~ sig_select(.x, .y, 3))) %>% 
  select(markers) %>% 
  unnest(markers)

saveRDS(sil_df3, "sil_df3.rds")

sil_df3 <- readRDS("intermediate_data/sil_df3.rds")

# Error in Rtsne.default(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.01084476806697,  : 
# Remove duplicates before running TSNE.

## Error 1 at mono_derived
sil_df3 %>% 
  nest(rdim = -level_2) %>% 
  pluck("rdim", 1) %>% 
  distinct(sample, symbol, count_scaled, level_3) %>% 
  reduce_dimensions(sample, symbol, count_scaled,
                    method = "tSNE",
                    action = "add",
                    transform = log1p)
# Error 2 at t_cell
sil_df3 %>% 
  nest(rdim = -level_2) %>% 
  pluck("rdim", 2) %>% 
  distinct(sample, symbol, count_scaled, level_3) %>% 
  reduce_dimensions(sample, symbol, count_scaled,
                    method = "tSNE",
                    action = "add",
                    transform = log1p)

sil_df3 %>% 
  nest(rdim = -level_2) %>% 
  pluck("rdim", 3) %>% 
  distinct(sample, symbol, count_scaled, level_3) %>% 
  reduce_dimensions(sample, symbol, count_scaled,
                    method = "tSNE",
                    action = "add",
                    transform = log1p)
# Contrast below shows that there aren't any duplicates
sil_df3 %>% 
  nest(rdim = -level_2) %>% 
  pluck("rdim", 1)

sil_df3 %>% 
  nest(rdim = -level_2) %>% 
  pluck("rdim", 1) %>% 
  distinct(sample, symbol, count_scaled, level_3)


