source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")

main <- function(.transcript, .is_hierarchy=TRUE, .level=NULL, 
                 .contrast_method, .ranking_method, .rank_stat=NULL, .bayes=NULL, 
                 .selection_method, .kmax=60, .discard_number=2000, .reduction_method = "PCA",
                 .optimisation_method, .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
                 .is_complete = FALSE) {
  
  
  if ( (!.is_hierarchy)|(.selection_method == "naive")) {
    
    .transcript %>%
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      scale_input_counts(.is_hierarchy=.is_hierarchy, .level=.level) %>%
      
      # Input: data.frame columns_1 <int> | ...
      do_ranking(.ranking_method=.ranking_method, .contrast_method=.contrast_method, .rank_stat=.rank_stat, .bayes=.bayes) %>%  
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      do_selection(.selection_method=.selection_method, .reduction_method=.reduction_method, .discard_number=.discard_number, .kmax=.kmax) %>% 
      
      # do_optimisation(.optimisation_method=.optimisation_method) %>% 
      
      format_output(.is_complete = .is_complete)
    
  } else {
    
    .transcript %>%
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      scale_input_counts(.is_hierarchy=.is_hierarchy, .level=.level) %>%
      
      do_ranking(.ranking_method=.ranking_method, .contrast_method=.contrast_method, .rank_stat=.rank_stat, .bayes=.bayes) %>% 
      
      mutate(level.copy = level) %>% 
      nest(data = -level.copy) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          
          do_selection(.selection_method=.selection_method, .reduction_method=.reduction_method, .discard_number=.discard_number, .kmax=.kmax) %>%
          
          # do_optimisation(.optimisation_method = .optimisation_method) %>%
          
          format_output(.is_complete = .is_complete)
        
      )) %>% 
      
      unnest(data) %>% 
      select(-level.copy)
    
  }
}

# BEGINNING
# mean_contrast, silhouette, hierarchy, curvature 

readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/counts_stefano_tree.rds") %>% 
  main(.contrast_method=mean_contrast, .ranking_method=rank_edgR_robust_likelihood_ratio, 
       .rank_stat = "logFC", .selection_method = "silhouette", .optimisation_method = NULL,
       .is_complete = TRUE) %>% 
  saveRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/signature_stefano_unoptimised.rds", 
          compress = "xz")
  


 # counts %>% 
 #  
 #  # tree_and_signatures_to_database(tree = tree, 
 #  #                                 .sample = sample, 
 #  #                                 .cell_type = cell_type,
 #  #                                 .symbol = symbol,
 #  #                                 .count = count) %>% 
 # 
 #  scale_input_counts(.is_hierarchy=TRUE) %>% 


# # debug
# counts %>% 
#   
#   distinct(sample, cell_type) %>% 
#   nest(s = -sample) %>%
#   mutate(n = map_int(s,~ nrow(.x))) %>%
#   filter(n > 1)  %>%
#   unnest(s)
# arrange(desc(n))
# 
# all(is.na(counts$count))
# 
# # find which function produce duplicates
# counts_first_db_raw %>% 
#   tidybulk(sample, symbol, count) %>%
#   
#   # aggregate duplicate sample/gene pairs in the data
#   # aggregate_duplicates(sample, symbol, count) %>%
#   
#   # rectangularise data
#   nest(data = -c(symbol, cell_type)) %>%
#   add_count(symbol) %>%
#   filter(n == max(n)) %>%
#   unnest(data) %>% 
#   
#   # Imputation of missing data
#   impute_missing_abundance(~ cell_type) %>%
#     
#   group_by(sample, cell_type, symbol) %>% 
#   count() %>% 
#   arrange(desc(n))
#   
#   # scale counts
#   identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
#   scale_abundance()
