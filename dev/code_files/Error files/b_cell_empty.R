# Error file, why b cell has become empty

b_cell_empty <- 
  
  # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
  contrast_copy %>% 
  
  slice(4) %>% 
  
  # select top 1 markers from each contrast, output is an unnested tibble
  sig_select(LEVEL, 1)

saveRDS(b_cell_empty, "b_cell_empty.rds")
