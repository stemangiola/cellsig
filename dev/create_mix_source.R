
load("dev/counts_imputed.rda")


mix_base = 
  counts_imputed %>%
  
  # keep genes present in all
  filter(count_scaled %>% is.na %>% `!`) %>%
  nest(data = -c(cell_type, symbol)) %>%
  add_count(symbol) %>%
  filter(n == max(n)) %>%
  unnest(data) %>% 
  select(-n)

S = 1

proportions = 
  gtools::rdirichlet(S, rep(1, length(unique(mix_base$cell_type)))) %>%
  as.data.frame() %>%
  setNames(unique(mix_base$cell_type))

mix = generate_mixture_from_proportion_matrix(mix_base, proportions) 


tidybulk::deconvolve_cellularity(
  mix,
  replicate, symbol, count_mix,
  #reference = ...,
  method = "llsr", 
  prefix = "llrs_", 
  action = "get"
  ) %>%
  pivot_longer(cols=starts_with("llrs_"), names_prefix ="llrs_", names_to="cell_type", values_to="estimated_proportion") %>%
  left_join(
    mix %>% 
      unnest(data_samples) %>% 
      distinct(replicate, cell_type, proportion)
  )
