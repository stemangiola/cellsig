# Modularisation ===============================================

# input tibble_of transcript abundance
# output tibble(node = c("root", "immuno", "CD4", ...), signature_genes)

# do_selection_naive_pairwise_hierarchical = function() 
# do_selection_naive_mean_contrast_hierarchical = function()
# do_selection_pairwise_silouette_hierarchical = function()
# do_selection_naive_pairwise_non_hierarchical = function() 
# do_selection_naive_mean_contrast_non_hierarchical = function()
# do_selection_pairwise_silouette_non_hierarchical = function()
# do_selection_cibersort = function()

# do_selection = function(input_tibble_of_transcript_abundance, method_ranking = function(), method_optimisation = function(),  method_selection = function(...)){
#
# do_ranking -> method_ranking(input_tibble_of_transcript_abundance)
# feature selection -> method_selection(input_tibble_of_transcript_abundance)
# feature optimisation -> method_optimisation(...)
# return feature union
# } 

# EXAMPLE
# do_selection(input_tibble_of_transcript_abundance, do_selection_naive_pairwise_hierarchical) %>%
#   distinct(gene) %>%
#   left_join(input_tibble_of_transcript_abundance) %>%
#   recuce_dimensions() %>%
#   calculate_siluette()

counts %>% 
  scale_input_counts(.level = "level_4") %>% 
  generate_contrast(.contrast_method = "pairwise") %>% 
  do_ranking(.ranking_method = "logFC") %>% 
  do_selection("naive", .kmax = 60, METHOD) %>% 
  do_optimisation("penalised") %>% 
  format_output()

mean_contrast.naive.non_hierarchy <- contrast_MC_NH %>% 
  do_ranking(.ranking_method = "logFC") %>%
  do_selection("naive", .kmax = 60, METHOD) %>% 
  do_optimisation("penalised") %>% 
  format_output(.is_complete = TRUE)

pairwise.naive.non_hierarchy <- contrast_PW_NH %>% 
  do_ranking(.ranking_method = "logFC") %>%
  do_selection("naive", .kmax = 60, METHOD) %>% 
  do_optimisation("penalised") %>% 
  format_output(.is_complete = TRUE)

mean_contrast.naive.hierarchy <- contrast_MC_H %>% 
  do_ranking(.ranking_method = "logFC") %>%
  do_selection("naive", .kmax = 60, METHOD) %>% 
  do_optimisation("penalised") %>% 
  format_output(.is_complete = TRUE)

pairwise.naive.hierarchy <- contrast_PW_H %>% 
  do_ranking(.ranking_method = "logFC") %>%
  do_selection("naive", .kmax = 60, METHOD) %>% 
  do_optimisation("penalised") %>% 
  format_output(.is_complete = TRUE)

mean_contrast.silhouette.non_hierarchy.before_optimisation <- contrast_MC_NH %>% 
  do_ranking(.ranking_method = "logFC") %>%
  do_selection("silhouette", .discard_number = 10000, .reduction_method = METHOD) %>% 
  do_optimisation("penalised") %>% 
  format_output(.is_complete = TRUE)

# currently running (before do_optimisation() step)
# pairwise is computationally unfeasible for silhouette selection
pairwise.silhouette.non_hierarchy <- contrast_PW_NH %>% 
  do_ranking(.ranking_method = "logFC") %>%
  do_selection("silhouette", .discard_number = 10000, .reduction_method = METHOD) %>% 
  do_optimisation("penalised") %>% 
  format_output(.is_complete = TRUE)


ranked_H <- contrast_MC_H %>% 
  do_ranking(.ranking_method = "logFC")

mean_contrast.silhouette.hierarchy <- contrast_MC_H %>% 
  do_ranking(.ranking_method = "logFC") %>%
  mutate(level.copy = level) %>% 
  nest(data = -level.copy) %>% 
  mutate(data = map(
    data,
    ~ .x %>% 
      do_selection("silhouette", .discard_number = 2000, .reduction_method = METHOD) %>%
      do_optimisation("penalised") %>% 
      format_output(.is_complete = TRUE)
  )) %>% 
  unnest(data) %>% 
  select(-level.copy)

# unoptimised data
saveRDS(mean_contrast.silhouette.hierarchy.unOP, 
        "mean_contrast.silhouette.hierarchy.unOP.rds",
        compress = "xz")

pairwise.silhouette.hierarchy <- contrast_PW_H %>% 
  do_ranking(.ranking_method = "logFC") %>%
  mutate(level.copy = level) %>% 
  nest(data = -level.copy) %>% 
  mutate(data = map(
    data,
    ~ .x %>% 
      do_selection("silhouette", .discard_number = 2000, .reduction_method = METHOD) %>%
      do_optimisation("penalised") %>% 
      format_output(.is_complete = TRUE)
  )) %>% 
  unnest(data) %>% 
  select(-level.copy)

# unoptimised data
saveRDS(pairwise.silhouette.hierarchy.unOP, 
        "pairwise.silhouette.hierarchy.unOP.rds",
        compress = "xz")

saveRDS(mean_contrast.naive.non_hierarchy, "mean_contrast.naive.non_hierarchy.rds", compress = "xz")
saveRDS(pairwise.naive.non_hierarchy, "pairwise.naive.non_hierarchy.rds", compress = "xz")
saveRDS(mean_contrast.naive.hierarchy, "mean_contrast.naive.hierarchy.rds", compress = "xz")
saveRDS(pairwise.naive.hierarchy, "pairwise.naive.hierarchy.rds", compress = "xz")
saveRDS(mean_contrast.silhouette.hierarchy, "mean_contrast.silhouette.hierarchy.rds", compress = "xz")
saveRDS(pairwise.silhouette.hierarchy, "pairwise.silhouette.hierarchy.rds", compress = "xz")
saveRDS(mean_contrast.silhouette.non_hierarchy, "mean_contrast.silhouette.non_hierarchy.rds", compress = "xz")
saveRDS()

# Main ============================

main <- function(.tree, .transcript, .sample, .cell_type, .symbol, .count,
                 .is_hierarchy=TRUE, .level=NULL, 
                 .contrast_method, .ranking_method, .selection_method,
                 .kmax = NULL, .discard_number = NULL, .reduction_method = "PCA",
                 .optimisation_method="penalised_silhouette", .penalty_rate = 0.2,
                 .is_complete = FALSE) {
  
  # mao the given developmental tree to a data frame
  counts <- .transcript %>% 
    tree_and_signatures_to_database(tree = tree, 
                                    .sample = sample, 
                                    .cell_type = cell_type,
                                    .symbol = symbol,
                                    .count = count)
  
  
  if ( (!.is_hierarchy)|(.selection_method == "naive")) {
    
    .counts %>% 
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      scale_input_counts(.is_hierarchy, .level) %>% 
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      generate_contrast(.contrast_method) %>% 
      
      # Input: data.frame columns_1 <int> | ...
      do_ranking(.ranking_method) %>% 
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      do_selection(.selection_method, .kmax, .reduction_method) %>% 
      
      do_optimisation(.optimisation_method, .penalty_rate) %>% 
      
      format_output(.is_complete)
    
  } else {
    
    counts %>% 
      
      scale_input_counts(.is_hierarchy, .level) %>% 
      
      generate_contrast(.contrast_method) %>% 
      
      do_ranking(.ranking_method) %>% 
      
      mutate(level.copy = level) %>% 
      nest(data = -level.copy) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          
          do_selection(.selection_method, .discard_number, .reduction_method) %>%
          
          do_optimisation(.optimisation_method, .penalty_rate) %>%
          
          format_output(.is_complete)
        
      )) %>% 
      
      unnest(data) %>% 
      select(-level.copy)
    
  }
}

# Map hierarchy to a dataframe =============================
tree_and_signatures_to_database = function(tree, signatures, .sample, .cell_type, .symbol, .count){
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  
  signatures %>%
    
    # Add tree info
    left_join(
      tree %>%
        data.tree::Clone() %>%
        ToDataFrameTypeColFull(fill=NA) ,
      by = quo_name(.cell_type)
    ) %>%
    filter(level_1 %>% is.na %>% `!`) %>%
    
    # Reduce size
    mutate_if(is.character, as.factor) %>% 
    droplevels %>% 
    mutate(!!.count := !!.count %>% as.integer) %>%
    
    # Filter only symbol existing
    filter(!!.symbol %>% is.na %>% `!`) %>%
    
    # Aggregate
    aggregate_duplicates(!!.sample, !!.symbol, !!.count) %>% 
    select(-one_of("merged_transcripts"))
}

#Test
counts2 <- transcriptome %>% 
  tree_and_signatures_to_database(tree = tree, 
                                  .sample = sample, 
                                  .cell_type = cell_type,
                                  .symbol = symbol,
                                  .count = count)

# Input scale abundance ================

scale_input_counts <- function(.transcript, .is_hierarchy = TRUE, .level = NULL){
  
  if (.is_hierarchy) { # scaling under each ancestor node at each level
    
    if (is.null(.level) == TRUE) { # scale counts for all levels present
      
      tt_hierarchy <- 
        
        tibble(level = names(.transcript) %>% str_subset("level")) %>% 
        
        mutate(tt = map(level, ~ .transcript %>% 
                          
                          mutate(level_0 = "root") %>% 
                          
                          preprocess(.x))) %>% 
        
        mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename("ancestor" = pre(.y))))
      
    } else { # scale counts for the level specified by .level
      
      tt_hierarchy <- 
        
        tibble(level = .level) %>% 
        
        mutate(tt = map(level, ~ .transcript %>% 
                          
                          mutate(level_0 = "root") %>% 
                          
                          preprocess(.x))) %>% 
        
        mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename("ancestor" = pre(.y))))
    }
    
    return(tt_hierarchy)
    
  } else { # non-hierarchical: scaling all cell types under the root node
    
    .level <- "root"
    
    tt_non_hierarchy <- 
      
      tibble(level = .level) %>% 
      
      mutate(tt = map(level, ~ .transcript %>%
                        
                        # create a root column for pre(.level)
                        mutate(!!as.symbol(.level) := cell_type) %>% 
                        mutate(!!as.symbol(pre(.level)) := .level) %>% 
                        
                        preprocess(.x))) %>% 
      
      mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename("ancestor" = pre(.y))))
    
    return(tt_non_hierarchy)
  }
}

preprocess <- function(.transcript, .level) {
  
  # load data
  .transcript %>%
    
    tidybulk(sample, symbol, count) %>%
    
    # aggregate duplicate sample/gene pairs in the data
    # aggregate_duplicates(sample, symbol, count) %>% 
    
    # rectangularise data
    nest(data = -c(symbol, cell_type)) %>%
    add_count(symbol) %>%
    filter(n == max(n)) %>%
    unnest(data) %>% 
    
    # Imputation of missing data
    impute_missing_abundance(~ cell_type) %>%
    
    # scale counts
    identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
    scale_abundance() %>% 
    
    # filter for cells at the level of interest
    filter(is.na(!!as.symbol(.level))==FALSE) %>%
    
    # nest by ancestor
    nest(data = - !!as.symbol(pre(.level)))
  
  # # scale count for further analysis
  # mutate(data=map(data, ~ .x %>%
  #                   identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
  #                   scale_abundance()
  # ))
}

pre <- function(.level) {
  
  if (.level == "root" ) {
    
    return("level_0")
    
  } else {
    
    ancestor <- 
      
      .level %>% 
      
      # str_split returns a list
      str_split("_") %>% 
      
      {as.numeric(.[[1]][2])-1} %>% 
      
      paste("level", ., sep = "_")
    
    return(ancestor)
  }
}

# Test
tt_non_hierarchy <- scale_input_counts(counts, .is_hierarchy = FALSE)
tt_L4_new <- scale_input_counts(counts, .is_hierarchy = TRUE, LEVEL)
tt_hierarchy <- scale_input_counts(counts, .is_hierarchy = TRUE)

LEVEL <- "level_4"

saveRDS(tt_non_hierarchy, "tt_non_hierarchy.rds", compress = "xz")
saveRDS(tt_L4, "tt_L4.rds", compress = "xz")
saveRDS(tt_hierarchy, "tt_hierarchy.rds", compress = "xz")
tt_L4

# summary output of preprocessed counts ===============
counts_preprocessed <- counts_preprocessed %>% 
  unnest(data) %>% 
  select(-c(level.copy, markers))

counts_summary <- counts_preprocessed %>% 
  summarise_counts()

counts_summary %>% 
  nest(data = -cell_type)

summarise_counts <- function(.preprocessed){
  .preprocessed %>% 
    unnest(data) %>% 
    dplyr::rename("gene" = "symbol") %>% 
    group_by(cell_type, gene) %>% 
    summarise(cell_type, gene, median_scaled_abundance = median(count_scaled) %>% ceiling()) %>% 
    distinct(cell_type, gene, median_scaled_abundance) %>% 
    ungroup()
}
  

# Generate contrast for ranking ===================

hypothesis_test_edgeR

hypothesis_test_bayes

pairwise_contrast = function(.data, .level){
  
  .data %>% 
    distinct(!!as.symbol(.level)) %>% 
    mutate(!!as.symbol(.level) := paste0(.level, !!as.symbol(.level))) %>% 
    
    # Permute
    mutate(cell_type2 = !!as.symbol(.level)) %>% 
    tidyr::expand(!!as.symbol(.level), cell_type2) %>% 
    filter(!!as.symbol(.level) != cell_type2) %>% 
    
    # Create contrasts
    mutate(contrast = sprintf("%s - %s", !!as.symbol(.level), cell_type2)) %>%
    pull(contrast)
  
}

mean_contrast <- function(.data, .level){
  
  cell_types <- .data %>% 
    distinct(!!as.symbol(.level)) %>%
    pull() %>% 
    paste0(.level, .)
  
  .data %>% 
    distinct(!!as.symbol(.level)) %>% 
    mutate(!!as.symbol(.level) := paste0(.level, !!as.symbol(.level))) %>% 
    mutate(background = map(
      !!as.symbol(.level), ~ cell_types[cell_types != .x])) %>% 
    mutate(contrast = map2_chr(
      !!as.symbol(.level), background,
      ~ sprintf("%s - (%s)/%s", .x, paste(.y, collapse="+"), length(.y))
    )) %>% 
    pull(contrast)
}

# OBSOLETE, combined with do_ranking now
generate_contrast <- function(.preprocessed, .contrast_method){
  
  .preprocessed %>%
    unnest(tt) %>% 
    
    # Differential transcription
    mutate(markers = map2(
      data, level,
      ~ .x %>% 
        test_differential_abundance(
          as.formula(sprintf("~ 0 + %s", .y)),
          .contrasts = .contrast_method(.x, .y),
          action="only") 
    ))
}

# Test
tt_L4 %>% 
  select(-level) %>% 
  unnest(tt) %>% 
  pluck("data", 1) %>% 
  pairwise_contrast(LEVEL)

contrast_MC_L4 <- tt_L4 %>% 
  generate_contrast(.contrast_method = mean_contrast)

contrast_MC_NH <- tt_non_hierarchy %>% 
  generate_contrast(.contrast_method = "mean_contrast")

contrast_PW_NH <- tt_non_hierarchy %>% 
  generate_contrast(.contrast_method = "pairwise")

contrast_PW_H <- tt_hierarchy %>% 
  generate_contrast(.contrast_method = "pair")

contrast_MC_H <- tt_hierarchy %>% 
  generate_contrast(.contrast_method = "mean_contrast")

saveRDS(contrast_PW_NH, "contrast_PW_NH.rds", compress = "xz")
saveRDS(contrast_MC_NH, "contrast_MC_NH.rds", compress = "xz")
saveRDS(contrast_MC_H, "contrast_MC_H.rds", compress = "xz")
saveRDS(contrast_PW_H, "contrast_PW_H.rds", compress = "xz")

# Rank ==========================================

do_ranking <- function(.preprocessed, .ranking_method, .contrast_method, .rank_stat=NULL, 
                       .bayes=NULL){
  
  # rank_stat takes either "Pvalue" or "logFC" 
  .preprocessed %>%
    
    # .bayes = .cellsig_theoretical_transcript_abundace_distribution
    .ranking_method(.contrast_method, .rank_stat, .bayes) %>% 
    
    # filter out potential nodes in which no genes are considered significant by rank_by_logFC
    filter(map_int(markers, ~ .x %>% unnest(stat_df) %>% nrow()) != 0)
    
}

rank_edgR_quasi_likelihood <- function(.preprocessed, .contrast_method, .rank_stat, .bayes=NULL){
  
  .preprocessed %>%
    unnest(tt) %>% 
    
    # Differential transcription: generate contrast
    mutate(markers = map2(
      data, level,
      ~ .x %>% 
        test_differential_abundance(
          as.formula(sprintf("~ 0 + %s", .y)),
          .contrasts = .contrast_method(.x, .y),
          method = "edgeR_quasi_likelihood",
          action="only") 
    )) %>% 
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ rank_by_stat(.x, .rank_stat) )) %>% 
    
    # remove prefixes from contrast expressions
    mutate(markers = map2(
      markers, level,
      ~ .x %>% 
        mutate(contrast = map2_chr(contrast, .y, ~ str_replace_all(.x, .y, "")))
    ))
}

rank_edgR_robust_likelihood_ratio <- function(.preprocessed, .contrast_method, .rank_stat="PValue",
                                              .bayes=NULL){
  
  .preprocessed %>%
    unnest(tt) %>% 
    
    # Differential transcription: generate contrast
    mutate(markers = map2(
      data, level,
      ~ .x %>% 
        test_differential_abundance(
          as.formula(sprintf("~ 0 + %s", .y)),
          .contrasts = .contrast_method(.x, .y),
          method = "edger_robust_likelihood_ratio",
          test_above_log2_fold_change = 1,
          action="only") 
    )) %>% 
    
    # Select markers from each contrast by rank of Pvalue
    mutate(markers = map(markers, ~ rank_by_stat(.x, "PValue") )) %>% 
    
    # remove prefixes from contrast expressions
    mutate(markers = map2(
      markers, level,
      ~ .x %>% 
        mutate(contrast = map2_chr(contrast, .y, ~ str_replace_all(.x, .y, "")))
    ))
}

rank_by_stat <-  function(.markers, .rank_stat){
  
  # .rank_stat = enquo(.rank_stat)
  
  .markers %>%
    
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
    mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>%
    
    # Filter out insignificant genes and rank the significant ones
    
    # THIS WILL HAVE TO CHANGE
    mutate(stat_df = map(
      stat_df, 
      ~ .x %>%
        filter(FDR < 0.05 & logFC > 2) %>%
        filter(logCPM > mean(logCPM)) )) %>% 
    
    mutate(stat_df = map(
      stat_df,
      ~ if(.rank_stat == "logFC"){
        .x %>% dplyr::arrange(desc(logFC))
      }else{
        .x %>% dplyr::arrange(PValue)
        }
    ))
  
}


rank_bayes <- function(.preprocessed, .contrast_method, .rank_stat=NULL, .bayes){
  
  .preprocessed %>% 
    
    unnest(tt) %>% 
    
    mutate(descendants = map2(data, level, ~ .contrast_method(.x, .y))) %>% 
    
    unnest(descendants) %>% 
    mutate(descendants = str_remove_all(descendants, "level_\\d")) %>% 
    separate(descendants, into = c("target", "background"), sep = " - ") %>% 
    
    mutate(background = str_split(background, "\\+")) %>% 
    mutate(background = map(background, ~ .x %>% str_remove_all("\\W|\\d$"))) %>% 
    
    mutate(lower_quantile = map(
      target,
      ~ .bayes %>% 
        filter(cell_type == .x) %>% 
        select(symbol, lower_quantile='25%') %>% 
        arrange(symbol)
    )) %>% 
    
    mutate(mean_upper_quantile = map(
      background,
      ~ .bayes %>% 
        # calculate the mean 75% quantile of each gene over all background cell types
        filter(cell_type %in% .x) %>% 
        group_by(symbol) %>% 
        summarise(symbol, mean_upper_quantile = mean(`75%`)) %>% 
        distinct() %>% 
        ungroup()
    )) %>% 
    
    mutate(stat_df = map2(
      lower_quantile, mean_upper_quantile,
      ~ inner_join(.x, .y, by= "symbol")
    )) %>% 
    select(-c(lower_quantile, mean_upper_quantile)) %>% 
    
    mutate(stat_df = map(
      stat_df,
      ~ .x %>% 
        mutate(difference = lower_quantile - mean_upper_quantile) %>% 
        arrange(desc(difference))
    )) %>% 
    
    nest(markers = -c(level, ancestor, data))
  
}

# OBSOLETE
do_ranking <- function(.contrast, .ranking_method){
  
  if (.ranking_method == "logFC") {
    
    .contrast %>% 
      
      # Select markers from each contrast by rank of stats
      mutate(markers = map(markers, ~ rank_by_logFC(.x) )) %>% 
      
      # filter out potential nodes in which no genes are considered significant by rank_by_logFC
      filter(map_int(markers, ~ .x %>% unnest(stat_df) %>% nrow()) != 0) %>% 
      
      # remove prefixes from contrast expressions
      mutate(markers = map2(
        markers, level,
        ~ .x %>% 
          mutate(contrast = map2_chr(contrast, .y, ~ str_replace_all(.x, .y, "")))
      ))
    
  }
  
}

# Test

ranked_PW_L4 <- tt_L4 %>% 
  do_ranking(.ranking_method = rank_edgR_quasi_likelihood, 
             .contrast_method = pairwise_contrast,
             .rank_stat = "logFC")

# Selection =======================================

do_selection <- 
  function(.ranked, .selection_method, .kmax=NULL, .discard_number=NULL, .reduction_method="PCA") {
    
    # .k is the number of genes selected from each cell_type contrast
    
    if (.selection_method == "naive") {
      
      .ranked %>% 
        
        do_naive_selection(.kmax, .reduction_method)
      
    } else {
      
      .ranked %>% 
        
        single_marker_pw_selection_using_silhouette(.discard_number, .reduction_method)
      
    }
    
  }

## Naive selection =============================

# calculate silhouette score for a series of sig_sizes

do_naive_selection <- function(.ranked, .kmax, .reduction_method) {
  
  # Args:
  # .ranked: output from do_ranking
  # .kmax: maximum number of markers selected from each cell type contrast
  # .reduction_method: method used to reduce dimensions such as "PCA", "tSNE", "MSA"
  
  tibble(number_of_markers_from_each_contrast = 1: .kmax) %>% 
    
    # select signature and calculate silhouette score 
    mutate(data = map(
      number_of_markers_from_each_contrast,
      ~ naive_selection(.ranked, .x) %>% 
        silhouette_function(.reduction_method)
    )) %>% 
    
    # nest by ancestor nodes/cell types
    unnest(data) %>%
    nest(data = - c(level, ancestor))
  
}

naive_selection <- function(.ranked, .k) {
  
  # Args:
  # .ranked: output from do_ranking()
  # .k: the number of genes selected from each cell_type contrast
  
  .ranked %>% 
    
    # selection markers from each contrast
    mutate(markers = map(
      markers,
      ~ .x %>% 
        mutate(stat_df = map(stat_df, ~ .x %>% dplyr::slice(1: .k))) %>% 
        unnest(stat_df)
    )) %>% 
    
    # Add original data info to the markers selected, 
    # use inner_join to ensure symbols are present in both markers and data
    mutate(markers = map2(markers, data, ~ inner_join(.x, .y, by="symbol"))) %>%
    
    # remove unnecessary column
    select(-data) %>% 
    
    # collect from which contrasts signature genes are extracted
    # mutate(contrast = map(markers, ~ .x %>% distinct(contrast, symbol))) %>% 
    
    # collect signature genes selected
    mutate(signature = map(markers, ~ .x$symbol %>% unique())) %>% 
    
    # number of the signature genes
    mutate(real_size = map_int(signature, ~ length(.x)))
  
}

# calculate silhouette score

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
    mutate(silhouette = pmap_dbl(
      list(reduced_dimensions, distance, level),
      ~ silhouette_score(..1, ..2, ..3)
    )) %>% 
    
    # remove unnecessary columns
    select(-c(markers, distance))
  
}

dimension_reduction <- function(.markers, .level, .reduction_method) {
  
  .markers %>% 
    
    distinct(sample, symbol, count_scaled, !!as.symbol(.level)) %>% 
    
    reduce_dimensions(sample, symbol, count_scaled, 
                      action = "get",
                      method = .reduction_method,
                      # .dims = 2,
                      transform = log1p,
                      top = Inf,
                      scale = FALSE,
                      check_duplicates = FALSE)
}

distance_matrix <- function(.reduced_dimensions, .reduction_method){
  
  .reduced_dimensions %>% 
    
    select(contains(str_sub(.reduction_method, end = -2L))) %>% 
    
    factoextra::get_dist(method = "euclidean")
}

silhouette_score <- function(.reduced_dimensions, .distance, .level){
  
  .reduced_dimensions %>% 
    
    pull(!!as.symbol(.level)) %>% 
    
    as.factor() %>% 
    
    as.numeric() %>% 
    
    silhouette(.distance) %>% 
    
    summary() %>% 
    
    .$avg.width
  
}


# Test

naive_PW_L4 <- ranked_PW_L4 %>% 
  naive_selection(1)

xx <- naive_PW_L4 %>% 
  silhouette_function(METHOD)

xx <- ranked_PW_L4 %>% 
  do_naive_selection(5, METHOD)

naive_MC_L4 <- bayes_ranked_L4 %>% 
  naive_selection(5)

bayes_naive_selected <- naive_MC_L4 %>% 
  silhouette_function(METHOD)

## Silhouette selection =========================================

single_marker_pw_selection_using_silhouette <- 
  function(.ranked, .discard_number=NULL, .reduction_method="PCA") {
    
    # initialize variables
    
    # ranked_copy is created as a pool of markers for selection, 
    # which continuously decrease with each iterative selection,
    # input .rank is used for calculating silhouette score for the selected markers
    ranked_copy <- .ranked 
    
    # initialise a signature tibble to store signature markers for each cell type in each iteration
    signature <- .ranked %>% 
      select(level, ancestor) %>% 
      mutate(signature = map(ancestor, ~ vector())) %>% 
      mutate(last_silhouette = 0)
    
    # initialise an output tibble containing all results of interest
    summary_tb <- tibble(
      level = character(),
      ancestor = character(),
      new_challengers = list(),
      winner = list(),
      winning_contrast = list(),
      signature = list(),
      # reduced_dimensions = list(),
      silhouette = double()
    )
    
    # set the base markers
    contrast_pair_tb0 <- 
      
      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      .ranked %>%
      
      # select top 1 markers from each contrast
      naive_selection(1) %>%
      
      dplyr::rename("new_challengers" = "signature") %>% 
      
      mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
      
      mutate(winning_contrast = map(
        markers,
        ~ .x %>% pull(contrast) %>% unique()
          # distinct(contrast, symbol) %>% 
          # mutate(contrast = contrast %>% str_extract(".*(?=\\s\\-)")) %>% 
          # mutate(contrast_symbol = map2_chr(contrast, symbol, ~ paste(.x, .y, sep = "."))) %>% 
          # pull(contrast_symbol)
      )) %>% 
      
      mutate(signature = winner) %>% 
      
      silhouette_function(.reduction_method) %>% 
      
      select(-c(reduced_dimensions, real_size))
    
    
    signature <- signature %>%
      
      # append cumulative markers
      mutate(signature = map2(
        signature, ancestor,
        ~ .x %>% 
          append(with(contrast_pair_tb0, signature[ancestor==.y][[1]]))
      )) %>% 
      
      # append silhouette scores for these markers
      mutate(last_silhouette = map_dbl(
        ancestor,
        ~ with(contrast_pair_tb0, silhouette[ancestor==.x])
      ))
    
    summary_tb <- summary_tb %>% 
      bind_rows(contrast_pair_tb0)
    
    # remove base markers from contrast_copy input before further selection
    ranked_copy <- ranked_copy %>%
      mutate(markers = map2(
        markers, ancestor, 
        ~ .x %>%
          unnest(stat_df) %>% 
          filter(!symbol %in% with(signature, signature[ancestor==.y][[1]])) %>% 
          nest(stat_df = - contrast)
      ))
    
    # counter for number of genes discarded
    j <- map_int(signature$signature, ~ length(.x))
    
    # count the number of iterations
    i <- 0
    while (any(j < .discard_number) &
           # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
           all(map_int(ranked_copy$markers, 
                       # hence the boundary should be the number of satisfactory genes selected
                       ~ .x %>% unnest(stat_df) %>% nrow()) > 0)) {
      
      contrast_pair_tb <- 
        
        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        ranked_copy %>% 
        
        # select top 1 markers from each contrast, ignore the signature output
        naive_selection(1) %>% 
        
        # pick the one new challenger from each contrast
        mutate(markers = map(
          markers,
          ~ .x %>% 
            nest(new_challenger = - contrast) %>% 
            mutate(new_challenger = map_chr(new_challenger, ~.x %>% distinct(symbol) %>% pull()))
        )) %>% 
        unnest(markers) %>% 
        select(-c(signature, real_size)) %>% 
        
        # append the new challenger from each contrast to the base markers for that ancestor node
        mutate(challengers_for_silhouette = map2(
          new_challenger, ancestor, 
          ~ with(signature, signature[ancestor==.y][[1]]) %>% 
            append(.x)
        )) %>% 
        
        # calculate silhouette score for the challengers from each contrast
        mutate(silhouette = map2_dbl(
          challengers_for_silhouette, ancestor, 
          ~ silhouette_for_markers(.ranked, .x, .y, .reduction_method) %>% 
            pull(silhouette)
        )) %>% 
        
        # arrange silhouette score in a descending manner within each ancestor node
        group_by(ancestor) %>% 
        arrange(desc(silhouette), .by_group = TRUE) %>% 
        ungroup() %>% 
        
        # check if the silhouette score for the challengers is greater than previous silhouette score
        mutate(is_greater = map2_lgl(
          silhouette, ancestor, 
          ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
        )) %>% 
        
        # nest under ancestor node to select markers that is TRUE for is_greater
        nest(data = - c(level, ancestor)) %>% 
        
        # record new_challengers
        mutate(new_challengers = map(data, ~ .x %>% pull(new_challenger))) %>% 
        
        # check if the biggest silhouette score is greater than previous score, if true we have a winner, else no winner
        mutate(winner = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ]$new_challenger
        } else {NA}
        )) %>% 
        
        # record which contrast the winner comes from
        mutate(winning_contrast = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ]$contrast # %>% str_extract(".*(?=\\s\\-)")
        } else {NA}
        )) %>% 
        
        # cummulative signature: winner + previously selected
        mutate(signature = pmap(
          list(data, winner, ancestor),
          ~ if(!is.na(..2)) {
            with(..1[1, ], challengers_for_silhouette[[1]])
          } else {
            with(signature, signature[ancestor==..3][[1]])
          }
        )) %>% 
        
        # silhouette score
        mutate(silhouette = map_dbl(data, ~ .x[[1, "silhouette"]]))
      
      
      # append the base + 1 markers that result in highest silhouette score
      signature <- signature %>% 
        
        mutate(signature = map(
          ancestor,
          ~ with(contrast_pair_tb, signature[ancestor==.x][[1]])
        )) %>% 
        
        mutate(last_silhouette = map2_dbl(
          ancestor, last_silhouette,
          ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x]))) {
            with(contrast_pair_tb, silhouette[ancestor==.x])
          } else {.y}
        )) 
      
      # append the winning signatures into the output summary table
      summary_tb <- summary_tb %>% 
        bind_rows(
          contrast_pair_tb %>% 
            filter(!is.na(winner)) %>% 
            # mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>% 
            select(-data)
        )
      
      # remove the signatures and unsuccessful genes from the selection list(ranked_copy)
      ranked_copy <- ranked_copy %>% 
        mutate(markers = map2(
          markers, ancestor, 
          ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y]))){
            .x %>% 
              unnest(stat_df) %>% 
              filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]])) %>% 
              nest(stat_df = -contrast)
          } else {
            .x %>% 
              unnest(stat_df) %>% 
              filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]])) %>% 
              nest(stat_df = -contrast)
          }
        ))
      
      # number of genes discarded for each node
      j <- j + 
        
        # unsuccessful candidates
        map_int(contrast_pair_tb$new_challengers, ~length(.x)) *
        is.na(contrast_pair_tb$winner) +
        
        # winning candidates
        map_int(contrast_pair_tb$winner, ~length(.x)) *
        !is.na(contrast_pair_tb$winner)
      
      cat("genes discarded for each node: ", j, "\n")
      cat("genes selected for each node: ", map_int(signature$signature, ~ length(.x)),  "\n")
      
      i <- i + 1
      cat("iteration: ", i, "\n")
      
    }
    
    # format output for optimisation
    output <- summary_tb %>% 
      mutate(real_size = map_int(signature, ~ length(.x))) %>% 
      nest(data = - c(level, ancestor))
    
    return(output)
  }

silhouette_for_markers <-function(.ranked, .signature, .ancestor, .reduction_method) {
  
  .ranked %>%
    
    filter(ancestor == .ancestor) %>% 
    
    select(-markers) %>%
    
    # filter markers that are in the signature
    mutate(data = map(data, ~.x %>% 
                        filter(symbol %in% .signature))) %>% 
    
    # format input
    dplyr::rename("markers" = "data") %>% 
    
    silhouette_function(.reduction_method)
  
}

## CIBERSORTx selection ==========================================
# use "cibersort_signature" from cibersortx.R

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

saveRDS(cibersortx, "cibersortx.rds", compress = "xz")

# Optimisation =========================================

do_optimisation <- function(.selected, 
                            .optimisation_method, 
                            .penalty_rate=0.2,
                            .kernel = "normal", 
                            .bandwidth = 0.05, 
                            .gridsize = 100){
  
  if(.optimisation_method == "penalty") {
    
    .selected %>% 
      
      mutate(optimal_size = map_int(data, ~ penalised_silhouette(.x, .penalty_rate))) %>% 
      
      unnest(data) %>% 
      
      filter(real_size == optimal_size)
    
  } else if (.optimisation_method == "curvature") {
    
    .selected %>% 
      
      curvature_of_kernel_smoothed_trend(.kernel, .bandwidth, .gridsize) %>% 
      
      unnest(data) %>% 
      
      filter(real_size <= optimal_size) %>% 
      
      select(-c(size.rescaled, smoothed)) %>% 
      
      nest(data = -c(level, ancestor, signature)) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          unnest(winner, winning_contrast) %>% 
          nest(enriched = -winning_contrast) %>% 
          mutate(enriched = map(enriched, ~ .x %>% pull(winner) %>% unique()))
      ))
  }
  
}

curvature <- function(.drv1, .drv2){
  abs(.drv2) / (1 + .drv1^2)^(3/2)
}

curvature_of_kernel_smoothed_trend <- function(.plot_data, 
                                               .kernel = "normal", 
                                               .bandwidth = 0.05, 
                                               .gridsize = 100){
  .plot_data %>% 
    
    mutate(data = map(
      data,
      ~ .x %>% 
        mutate(size.rescaled = rescale(real_size))
    )) %>% 
    
    mutate(smoothed.estimate = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette, 
                drv = 0L, degree=2, kernel = .kernel, 
                bandwidth = .bandwidth, gridsize = .gridsize) %>% 
        as_tibble() %>% 
        `colnames<-`(c("grid", "estimate"))
    )) %>% 
    
    mutate(first.derivative = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette, 
                drv = 1L, degree=2, kernel = .kernel, 
                bandwidth = .bandwidth, gridsize = .gridsize) %>% 
        as_tibble() %>% 
        `colnames<-`(c("grid", "deriv1"))
    )) %>% 
    
    mutate(second.derivative = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette, 
                drv = 2L, degree=2, kernel = .kernel, 
                bandwidth = .bandwidth, gridsize = .gridsize) %>% 
        as_tibble() %>% 
        `colnames<-`(c("grid", "deriv2"))
    )) %>% 
    
    mutate(smoothed = pmap(
      list(smoothed.estimate, first.derivative, second.derivative),
      ~..1 %>% 
        left_join(..2, by = "grid") %>% 
        left_join(..3, by = "grid")
    )) %>% 
    
    mutate(smoothed = map(
      smoothed,
      ~ .x %>% 
        mutate(curvature = map2_dbl(
          deriv1, deriv2,
          ~ curvature(.x, .y)
        ))
    )) %>% 
    
    select(-c(smoothed.estimate, first.derivative, second.derivative)) %>% 
    
    # optimal_size
    mutate(optimal_size = map_dbl(
      smoothed,
      ~ with(.x, grid[which(peaks(curvature))[which.max(curvature[peaks(curvature)])]]) 
    )) %>% 
    
    mutate(optimal_size = map2_dbl(
      data, optimal_size,
      ~ with(.x, size.rescaled[which.min(abs(size.rescaled - .y))])
    )) %>% 
    
    mutate(optimal_size = map2_int(
      data, optimal_size,
      ~ .x %>% 
        with(real_size[size.rescaled == .y])
    )) %>% 
    
    mutate(optimal_size = ifelse(optimal_size<10, 10, optimal_size))
  
}

penalised_silhouette <- function(.plot_data, .penalty_rate=0.2) {
  
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    # filter(real_size > 10) %>%
    
    # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
    mutate(size_rescaled = rescale(real_size)) %>% 
    
    mutate(penalised_silhouette = silhouette - .penalty_rate * size_rescaled) %>% 
    
    filter(penalised_silhouette == max(penalised_silhouette)) %>% 
    
    pull(real_size)
}

# Test

mean_contrast.silhouette.hierarchy.unOP %>% 
  do_optimisation("penalised", .penalty_rate = 0.61)

# Format output ===================================

format_output <- function(.optimised, .is_complete=FALSE){
  
  if (!.is_complete) {
    
    .optimised %>% 
      
      select(node = ancestor, signature)
    
  } else {
    
    .optimised
    
  }
  
  
}

# Test

x %>% format_output()


# Compare all methods ============================
# use the complete output from format_output as input

# import summary data from all methods
naive <- list.files("dev/topInf_scaleFALSE/", pattern = ".*naive\\..*\\..*")
silhouette <- list.files("dev/topInf_scaleFALSE/", pattern = ".*silhouette\\..*\\..*")

naive_df <- map_dfr(naive, ~ readRDS(paste0("dev/topInf_scaleFALSE/", .x))) %>% 
  select(level, ancestor, real_size, signature, silhouette) %>% 
  mutate(method = rep(str_replace_all(naive, '\\.rds', ''), c(14, 1, 14, 1)))

o <- rep(str_replace_all(silhouette, '\\.rds', ''), c(14, 1, 14))
silhouette_df <- map_dfr(silhouette, ~ readRDS(paste0("dev/topInf_scaleFALSE/", .x))) %>% 
  select(level, ancestor, real_size, signature=cumulative_signature, silhouette) %>% 
  mutate(method = o)
rm(o)

full_df <- silhouette_df %>% 
  bind_rows(naive_df)

all_methods_silhouette <- full_df %>% 
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


# summary table comparing all methods using silhouette score

all_methods_comparison <- all_methods_silhouette %>% 
  bind_rows(cibersortx) %>% 
  mutate(method = str_remove_all(method, "\\.new"))
  arrange(desc(silhouette))

all_methods_comparison

saveRDS(all_methods_comparison, "all_methods_comparison.new.rds", compress = 'xz')

# summary bar plot comparing all methods using silhouette score

all_methods_comparison %>% 
  
  ggplot(aes(reorder(method, silhouette), silhouette, fill = method)) +
  geom_col() +
  geom_text(aes(label = round(silhouette, 3)), vjust = 1.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("All methods comparison using silhouette score")

# PCA plot 
all_methods_comparison %>% 
  pluck("reduced_dimensions", 8) %>% 
  ggplot(aes(PC1, PC2, color = root), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("cibersortx")
ggtitle("mean_contrast.silhouette.non_hierarchy")


