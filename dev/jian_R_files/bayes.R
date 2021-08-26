
# final function (input .bayes is the dataframe processed as what's below the function)
rank_by_bayes <- function(.preprocessed, .contrast_method, .bayes){
  
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

do_ranking <- function(.preprocessed, .contrast_method, .ranking_method, .bayes){
  
  if (.ranking_method == "logFC") {
    
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
      )) %>% 
      
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
    
  }else{
    
    .preprocessed %>% 
      
      rank_by_bayes(.contrast_method, .bayes)
    
  }
  
}

bayes <- readRDS("dev/cellsig_theoretical_transcript_abundance_distribution.rds") %>% 
  
  # take the first row of the duplicated data so that each gene in a cell type has only one set of quantiles
  nest(data = - c(symbol, cell_type)) %>% 
  mutate(data = map(data, ~.x[1, ])) %>% 
  
  # ensure all cell types have the same set of genes
  add_count(symbol) %>%
  filter(n == max(n)) %>%

  unnest(data)
  

bayes_ranked_L4 <- tt_L4 %>% 
  rank_by_bayes(mean_contrast, bayes)


bayes_tibble_L4 <- tt_L4 %>% 
  
  unnest(tt) %>% 
  
  mutate(descendants = map2(data, level, ~ mean_contrast(.x, .y))) %>% 
  
  unnest(descendants) %>% 
  mutate(descendants = str_remove_all(descendants, "level_\\d")) %>% 
  separate(descendants, into = c("target", "background"), sep = " - ") %>% 
  
  mutate(background = str_split(background, "\\+")) %>% 
  mutate(background = map(background, ~ .x %>% str_remove_all("\\W|\\d$")))

xx <- "( t_CD4_memory )/1"
xx %>% str_remove_all("\\W|\\d$")

x_L4 <- bayes_tibble_L4 %>% 
  
  mutate(lower_quantile = map(
    target,
    ~ bayes %>% 
      filter(cell_type == .x) %>% 
      select(symbol, lower_quantile='25%') %>% 
      arrange(symbol)
  )) %>% 
  
  mutate(mean_upper_quantile = map(
    background,
    ~ bayes %>% 
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
  )) 

x_L4 <- x_L4 %>%  
  nest(markers = -c(level, ancestor, data))


y <- y %>% 
  mutate(quantiles = map2(
    lower_quantile, mean_upper_quantile,
    ~ inner_join(.x, .y, by= "symbol")
  )) %>% 
  select(-c(lower_quantile, mean_upper_quantile))
  
bayes %>% 
  filter(cell_type %in% c("epithelial", "endothelial", "fibroblast")) %>% 
  # group_by(symbol) %>% 
  # summarise(symbol, mean_upper_quantile = mean(`75%`))
  pivot_wider(id_cols = c(symbol, level), names_from = cell_type, values_from = `75%`)
  
# problem: there are duplicate genes for some cell types
bayes %>% 
  filter(cell_type %in% c("epithelial", "endothelial", "fibroblast")) %>% 
  nest(data = - c(symbol, cell_type)) %>% 
  mutate(data = map(data, ~.x[1, ])) %>% 
  unnest(data) %>% 
  group_by(symbol) %>% 
  summarise(symbol, mean_upper_quantile = mean(`75%`))
  
  pivot_wider(id_cols = c(symbol, level), names_from = cell_type, values_from = `75%`)
  


  filter(cell_type == "fibroblast" & symbol == "ABCB6")

bayes %>% 
  distinct(symbol) %>% 
  filter(cell_type == "fibroblast" & symbol == "ADCY5")

# problem: some genes are present in one cell type but not in other cell types of the same level
bayes %>% 
  # rectangularise data
  nest(data = -c(symbol, cell_type)) %>%
  add_count(symbol) %>%
  filter(n == max(n)) %>%
  unnest(data)
