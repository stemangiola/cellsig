# create new bayes data frame based on the new tree


readRDS("dev/counts_bayes_imputed_old_tree.rds") %>% 
  select(-contains("level")) %>% 
  rename(feature = .feature, sample = .sample) %>% 
  left_join(new_tree %>% ToDataFrameTypeColFull(fill=NA), by = "cell_type") %>% 
  saveRDS("dev/intermediate_data/counts_bayes_imputed.rds", compress = "xz")



