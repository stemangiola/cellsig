source("dev/jian_R_files/cellsig_for_kamran.R")


load("dev/counts.rda") 

counts = counts %>% filter(cell_type %in% c("epithelial","nk_resting", "nk_primed", "eosinophil", "neutrophil")) %>%
  select(-c('level_1', 'level_2', 'level_3', 'level_4', 'level_5'))


tree = yaml.load_file("dev/tree_test.yaml") %>% 
  as.Node

sig_de_rank = counts %>% 
  main(.sample = sample, .symbol = symbol, .count = count, .cell_type = cell_type,
       .is_hierarchy=TRUE, 
       .tree = tree,
       .contrast_method = pairwise_contrast, .ranking_method = rank_edgR_robust_likelihood_ratio, .rank_stat="PValue",
       .selection_method = "silhouette", .reduction_method = "PCA", .dims=4, .discard_number = 1000,
       .optimisation_method = "penalty",
       .is_complete = TRUE) %>%
  unnest(children) %>% unnest(enriched) %>% unnest(signature) %>% as_tibble()



