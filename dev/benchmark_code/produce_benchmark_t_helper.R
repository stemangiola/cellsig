source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")

rank_bayes <- function(.hierarchical_counts, .sample, .symbol, .cell_type, 
                       .contrast_method, .bayes, .tree, 
                       .rank_stat=NULL){
  
  # Args:
  # .bayes is the bayes data that have been imputed and obtained gene imputation ratio by the code:
  # counts_bayes_imputed %>%
  #   select(-level) %>%
  #   scale_input_counts(.is_hierarchy = TRUE)
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)
  
  .bayes %>%
    
    # force the column names of bayes data to be consistent with input expression data
    dplyr::rename(!!.symbol := feature, !!.sample := sample, !!.cell_type := cell_type) %>%
    do_hierarchy(.is_hierarchy = all(.hierarchical_counts$level != "root"), .tree = .tree,
                 .sample=!!.sample, .symbol=!!.symbol, .cell_type= !!.cell_type) %>%
    
    unnest(tt) %>%
    
    mutate(contrast = map2(data, level, ~ .contrast_method(.x, .y))) %>%
    
    unnest(contrast) %>%
    mutate(contrast = str_remove_all(contrast, glue("{level}"))) %>%
    
    mutate(lower_quantile = map2(
      data, contrast,
      ~ .x %>%
        # filter for target cell type in the contrast
        filter(!!.cell_type == str_extract(.y, ".*(?=\\s\\-)")) %>%
        # filter out genes with imputation ratio greater than 0.2 (only used for user pipeline not benchmark)
        # filter(ratio_imputed_samples < 0.2) %>%
        select(!!.symbol, lower_quantile='25%')
    )) %>%
    
    mutate(mean_upper_quantile = map2(
      data, contrast,
      ~ {
        # obtain background cell type(s)
        background <- (.y) %>%
          str_extract("(?<=\\-\\s).*") %>%
          str_split("\\+") %>%
          unlist() %>%
          str_remove_all("(?<=\\/).*|\\W")
        
        (.x) %>%
          # calculate the mean 75% quantile of each gene over all background cell types
          filter(!!.cell_type %in% background) %>%
          group_by(!!.symbol) %>%
          summarise(!!.symbol, mean_upper_quantile = mean(`75%`)) %>%
          distinct() %>%
          ungroup()
      }
    )) %>%
    
    mutate(stat_df = map2(
      lower_quantile, mean_upper_quantile,
      ~ left_join(.x, .y, by= quo_name(.symbol))
    )) %>%
    select(-c(lower_quantile, mean_upper_quantile)) %>%
    
    mutate(stat_df = map(
      stat_df,
      ~ .x %>%
        mutate(difference = lower_quantile - mean_upper_quantile) %>%
        arrange(desc(difference))
    )) %>%
    
    select(-data) %>%
    
    nest(markers = -c(level, ancestor)) %>%
    
    right_join(.hierarchical_counts %>% unnest(tt), by = c("level", "ancestor")) %>%
    
    unnest(markers)
  
}

main <- function(.input, .sample, .symbol, .count=NULL, .cell_type,
                 .is_hierarchy=TRUE, .level=NULL, 
                 .tree, .node=NULL,
                 .contrast_method, .ranking_method, .rank_stat=NULL, .bayes=NULL, 
                 .selection_method, .kmax=60, .discard_number=2000, .reduction_method = "PCA", .dims=2,
                 .optimisation_method, .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
                 .is_complete = TRUE) {
  
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  .cell_type = enquo(.cell_type)
  
  # subtree = tree_subset(.tree=.tree, .node=.node)
  
  if ( (!.is_hierarchy)|(.selection_method == "naive")) {
    
    .input %>%
      
      # adapt_tree(.tree = .tree, .node = .node) %>%
      # 
      # tree_and_signatures_to_database(tree=subtree, ., .sample=!!.sample, .cell_type=!!.cell_type,
      #                                .symbol=!!.symbol, .count=!!.count) %>%
      # 
      # do_scaling(.sample = !!.sample, .symbol= !!.symbol , .count= !!.count, .cell_type=!!.cell_type) %>% 
      #   
      # do_imputation(.sample = !!.sample, .symbol=!!.symbol, .cell_type=!!.cell_type) %>% 
      
      
    do_hierarchy(.sample=!!.sample,
                 .symbol=!!.symbol,
                 .cell_type = !!.cell_type,
                 .tree = .tree,
                 .is_hierarchy=.is_hierarchy,
                 .level=.level) %>%
    
    # Input: data.frame columns_1 <int> | ...
    do_ranking(.sample=!!.sample, 
               .symbol=!!.symbol,
               .cell_type = !!.cell_type,
               .ranking_method=.ranking_method, 
               .contrast_method=.contrast_method, 
               .rank_stat=.rank_stat, 
               .bayes=.bayes,
               .tree = .tree) %>%  
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      do_selection(.sample=!!.sample, 
                   .symbol=!!.symbol,
                   .selection_method=.selection_method, 
                   .reduction_method=.reduction_method, 
                   .discard_number=.discard_number, 
                   .kmax=.kmax,
                   .dims=.dims) %>% 
      
      do_optimisation(.optimisation_method=.optimisation_method, .symbol=!!.symbol) %>% 
      
      format_output(.is_complete = .is_complete)
    
  } else {
    
    .input %>%
      
      # adapt_tree(.tree = .tree, .node = .node) %>% 
      # 
      # tree_and_signatures_to_database(tree=subtree, ., .sample=!!.sample, .cell_type=!!.cell_type, 
      #                                .symbol=!!.symbol, .count=!!.count) %>% 
      # 
      # do_scaling(.sample = !!.sample, .symbol= !!.symbol , .count= !!.count, .cell_type=!!.cell_type) %>% 
      #   
      # do_imputation(.sample = !!.sample, .symbol=!!.symbol, .cell_type=!!.cell_type) %>% 
      
      
    do_hierarchy(.sample=!!.sample,
                 .symbol=!!.symbol,
                 .cell_type = !!.cell_type,
                 .tree = .tree,
                 .is_hierarchy=.is_hierarchy,
                 .level=.level) %>%
    
    do_ranking(.sample=!!.sample, 
               .symbol=!!.symbol,
               .ranking_method=.ranking_method, 
               .contrast_method=.contrast_method, 
               .rank_stat=.rank_stat, 
               .bayes=.bayes,
               .tree = .tree) %>% 
      
      mutate(level.copy = level) %>% 
      nest(data = -level.copy) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          
          do_selection(.sample=!!.sample, 
                       .symbol=!!.symbol,
                       .selection_method=.selection_method, 
                       .reduction_method=.reduction_method, 
                       .discard_number=.discard_number, 
                       .kmax=.kmax,
                       .dims=.dims) %>%
          
          do_optimisation(.optimisation_method = .optimisation_method, .symbol=!!.symbol) %>%
          
          format_output(.is_complete = .is_complete)
        
      )) %>% 
      
      unnest(data) %>% 
      select(-level.copy)
    
  }
}


args <- commandArgs(trailingOnly = TRUE)

is_hierarchy = args[1]
contrast = args[2]
ranking_method = args[3] 
ranking_stat = args[4]
bayes = args[5] %>% as.integer()
selection = args[6] 
optimisation = args[7] 
output_file = args[8]

.is_hierarchy = is_hierarchy %>% when(
  (.)=="hierarchical" ~ TRUE,
  (.)=="non_hierarchical" ~ FALSE
)

.contrast_method = contrast %>% when(
  (.)=="mean_contrast" ~ mean_contrast,
  (.)=="pairwise_contrast" ~ pairwise_contrast
)

.ranking_method = ranking_method %>% when(
  (.)=="edgR" ~ rank_edgR_quasi_likelihood,
  (.)=="bayes" ~ rank_bayes
)

.rank_stat = ranking_stat %>% when(
  (.)=="_" ~ NULL,
  ~ (.)
)

.bayes  =
  if(bayes==1L){
    readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/benchmark_results_t_helper/counts_bayes_imputed_t_helper_tree.rds")
  } else {NULL}

readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/benchmark_results_t_helper/counts_imputed_t_helper_tree.rds") %>%
  dplyr::rename(symbol = feature) %>% 
  main(.sample = sample, 
       .symbol = symbol,
       .count = NULL,
       .tree = NULL,
       .cell_type = cell_type,
       .is_hierarchy=.is_hierarchy,
       .contrast_method=.contrast_method, 
       .ranking_method=.ranking_method, 
       .rank_stat=.rank_stat, 
       .bayes=.bayes, 
       .selection_method=selection, .kmax=60, .discard_number=2000, .reduction_method = "PCA", .dims=2,
       .optimisation_method=optimisation, .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
       .is_complete = TRUE) %>% 
  saveRDS(file = output_file, compress = "xz")
