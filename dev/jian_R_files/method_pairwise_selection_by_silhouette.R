# new method

sig_select <- function(.contrast, .level, .sig_size) {
  .contrast %>% 
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, .sig_size))) %>%
    
    # Add original data info to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(level, ancestor, markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, .level, "") %>% str_replace(.level, ""))
}

# debug workplace
single_marker_pw_select <- function(.contrast, .level, .discard_num, .method) {
  
  .discard_num = 200
  
  # initialize variables
  
  contrast_copy <- contrast_PW_L3 %>% 
    mutate(ancestor = !!as.symbol(pre(LEVEL)))
  
  # initialise a signature tibble to store signature markers for each cell type in each iteration
  signature <- tibble(
    ancestor = contrast_copy %>% 
      pull(ancestor)
  ) %>% 
    mutate(cumulative_signature = map(ancestor, ~ vector())) %>% 
    mutate(last_silhouette = 0)
  
  # initialise an output tibble containing all results of interest
  summary_tb <- tibble(
    ancestor = character(),
    new_challengers = list(),
    winner = list(),
    cumulative_signature = list(),
    reduced_dimensions = list(),
    winning_silhouette = double()
  )
  
  # set the base markers
  contrast_pair_tb0 <-
    
    # contrast_copy contains all the statistics of all cell_type contrasts for each gene
    contrast_PW_L4 %>%
    
    # select top 1 markers from each contrast, output is an unnested tibble
    sig_select(LEVEL, 1) %>%
    
    # mutate(ancestor = !!as.symbol(pre(LEVEL))) %>%
    
    nest(data = - c(level, ancestor)) %>%
    
    mutate(data = map(data, ~ .x %>% 
                        nest(new_challenger = - contrast_pretty) %>% 
                        
                        mutate(new_challenger = map(
                          new_challenger, ~ .x %>% pull(symbol) %>% unique()))
                      
                      # mutate(challengers_for_silhouette = list(new_challenger))
    )) %>% 
    
    mutate(new_challengers = map(data, ~.x %>% 
                                   pull(new_challenger) %>% 
                                   unlist())) %>% 
    
    mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
    
    mutate(cumulative_signature = winner) %>% 
    
    select(-data) %>% 
    
    mutate(silhouette = map2(
      new_challengers, ancestor,
      ~ silhouette_for_markers(.x, .y, contrast_PW_L3, LEVEL, METHOD) %>% 
        select(reduced_dimensions, winning_silhouette = silhouette)
    )) %>% 
    
    unnest(silhouette)
  
  # unnest(data) %>% 
  # 
  # mutate(silhouette = map2(
  #   challengers_for_silhouette, ancestor,
  #   ~ silhouette_for_markers(.x, .y, contrast_PW_L3, LEVEL, METHOD) %>% 
  #     select(-ancestor)
  # )) %>% 
  # 
  # unnest(silhouette) %>% 
  # 
  # mutate(is_greater = map2_lgl(
  #   silhouette, ancestor,
  #   ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
  # )) %>% 
  # 
  # nest(data = - ancestor) %>% 
  # 
  # mutate(new_challengers = map(data, ~.x %>% pull(new_challenger))) %>% 
  # 
  # mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
  # 
  # mutate(cumulative_signature = winner)
  
  
  signature <- signature %>%
    mutate(cumulative_signature = map2(
      cumulative_signature, ancestor,
      ~ .x %>% 
        append(with(contrast_pair_tb0, cumulative_signature[ancestor==.y][[1]]))
    )) %>% 
    mutate(last_silhouette = map_dbl(
      ancestor,
      ~ with(contrast_pair_tb0, winning_silhouette[ancestor==.x])
    ))
  
  # mutate(last_silhouette = map_dbl(
  #   ancestor,
  #   ~ with(with(contrast_pair_tb0, data[ancestor==.x])[[1]], unique(silhouette))
  # ))
  
  summary_tb <- summary_tb %>% 
    bind_rows(contrast_pair_tb0)
  
  # remove base markers from contrast_copy input before further selection
  contrast_copy <- contrast_copy %>%
    mutate(markers = map2(
      markers, ancestor, 
      ~ .x %>%
        filter(!symbol %in% with(signature, cumulative_signature[ancestor==.y][[1]]))
    ))
  
  # counter for number of genes discarded
  j <- map_int(signature$cumulative_signature, ~ length(.x))
  
  while (any(j < .discard_num) & 
         all(map_int(contrast_copy$markers, 
                     ~ .x %>% select_markers_for_each_contrast(1) %>% nrow()) > 0)) {
    
    contrast_pair_tb <- 
      
      # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(LEVEL, 1) %>% 
      
      mutate(ancestor = !!as.symbol(pre(LEVEL))) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(data = map(data, ~ .x %>% 
                          nest(new_challenger = - contrast_pretty) %>% 
                          mutate(new_challenger = map_chr(
                            new_challenger, 
                            ~.x %>% pull(symbol) %>% unique()
                          ))
      )) %>% 
      
      unnest(data) %>% 
      
      mutate(challengers_for_silhouette = map2(
        new_challenger, ancestor, 
        ~ with(signature, cumulative_signature[ancestor==.y][[1]]) %>% 
          append(.x)
      )) %>% 
      
      mutate(silhouette = map2(
        challengers_for_silhouette, ancestor, 
        ~ silhouette_for_markers(.x, .y, contrast_PW_L3, LEVEL, METHOD) %>% 
          select(-ancestor)
      )) %>% 
      
      unnest(silhouette) %>% 
      
      group_by(ancestor) %>% 
      
      arrange(desc(silhouette), .by_group = TRUE) %>% 
      
      ungroup() %>% 
      
      mutate(is_greater = map2_lgl(
        silhouette, ancestor, 
        ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
      )) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(new_challengers = map(
        data,
        ~ .x %>% 
          pull(new_challenger)
      )) %>% 
      
      mutate(winner = map(
        data,
        ~ if(.x[1, ]$is_greater){
          .x[1, ]$new_challenger
        } else {NA}
      )) %>% 
      
      mutate(cumulative_signature = pmap(
        list(data, winner, ancestor),
        ~ if(!is.na(..2)) {
          with(..1[1, ], challengers_for_silhouette[[1]])
        } else {
          with(signature, cumulative_signature[ancestor==..3][[1]])
        }
      ))
    
    # append the base + 1 markers that result in highest silhouette score
    signature <- signature %>% 
      
      mutate(cumulative_signature = map(
        ancestor,
        ~ with(contrast_pair_tb, cumulative_signature[ancestor==.x][[1]])
      )) %>% 
      
      mutate(last_silhouette = map2_dbl(
        ancestor, last_silhouette,
        ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x][[1]]))) {
          with(contrast_pair_tb, data[ancestor==.x][[1]][[1, "silhouette"]])
        } else {.y}
      ))
    
    summary_tb <- summary_tb %>% 
      bind_rows(
        contrast_pair_tb %>% 
          filter(!is.na(winner)) %>% 
          mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>% 
          mutate(winning_silhouette = map_dbl(data, ~ .x$silhouette[1])) %>% 
          select(-data)
      )
    
    contrast_copy <- contrast_copy %>% 
      mutate(markers = map2(
        markers, ancestor, 
        ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y][[1]]))){
          .x %>% 
            filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))
        } else {
          .x %>% 
            filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]]))
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
    cat("genes selected for each node: ", map_int(signature$cumulative_signature, ~ length(.x)),  "\n")
  }
  
  return(summary_tb %>% 
           nest(signature_data = - ancestor))
}



# below are real functions ===================================

# Hierarchical

# pairwise selection 1 marker at a time
single_marker_pw_select <- function(.contrast, .level, .discard_num, .method) {
  
  # initialize variables
  
  contrast_copy <- .contrast %>% 
    mutate(ancestor = !!as.symbol(pre(.level)))
  
  # initialise a signature tibble to store signature markers for each cell type in each iteration
  signature <- tibble(
    ancestor = contrast_copy %>% 
      pull(ancestor)
  ) %>% 
    mutate(cumulative_signature = map(ancestor, ~ vector())) %>% 
    mutate(last_silhouette = 0)
  
  # initialise an output tibble containing all results of interest
  summary_tb <- tibble(
    ancestor = character(),
    new_challengers = list(),
    winner = list(),
    cumulative_signature = list(),
    # reduced_dimensions = list(),
    winning_silhouette = double()
  )
  
  # set the base markers
  contrast_pair_tb0 <-
    
    # contrast_copy contains all the statistics of all cell_type contrasts for each gene
    contrast_copy %>%
    
    # select top 1 markers from each contrast, output is an unnested tibble
    sig_select(.level, 1) %>%
    
    mutate(ancestor = !!as.symbol(pre(.level))) %>%
    
    nest(data = - ancestor) %>%
    
    mutate(data = map(data, ~ .x %>% 
                        nest(new_challenger = - contrast_pretty) %>% 
                        
                        mutate(new_challenger = map(
                          new_challenger, ~ .x %>% pull(symbol) %>% unique()))
    )) %>% 
    
    mutate(new_challengers = map(data, ~.x %>% pull(new_challenger) %>% unlist())) %>% 
    
    mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
    
    mutate(cumulative_signature = winner) %>% 
    
    select(-data) %>% 
    
    mutate(silhouette = map2(
      new_challengers, ancestor,
      ~ silhouette_for_markers(.x, .y, .contrast, .level, .method) %>% 
        select(
          # reduced_dimensions, 
          winning_silhouette = silhouette)
    )) %>% 
    
    unnest(silhouette)
  
  
  signature <- signature %>%
    mutate(cumulative_signature = map2(
      cumulative_signature, ancestor,
      ~ .x %>% 
        append(with(contrast_pair_tb0, cumulative_signature[ancestor==.y][[1]]))
    )) %>% 
    mutate(last_silhouette = map_dbl(
      ancestor,
      ~ with(contrast_pair_tb0, winning_silhouette[ancestor==.x])
    ))
  
  summary_tb <- summary_tb %>% 
    bind_rows(contrast_pair_tb0)
  
  # remove base markers from contrast_copy input before further selection
  contrast_copy <- contrast_copy %>%
    mutate(markers = map2(
      markers, ancestor, 
      ~ .x %>%
        filter(!symbol %in% with(signature, cumulative_signature[ancestor==.y][[1]]))
    ))
  
  # counter for number of genes discarded
  j <- map_int(signature$cumulative_signature, ~ length(.x))
  
  # count the number of iterations
  i <- 0
  while (any(j < .discard_num) &
         # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
         all(map_int(contrast_copy$markers, 
                     # hence the boundary should be the number of satisfactory genes selected
                     ~ .x %>% select_markers_for_each_contrast(1) %>% nrow()) > 0)) {
    
    contrast_pair_tb <- 
      
      # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(.level, 1) %>% 
      
      mutate(ancestor = !!as.symbol(pre(.level))) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(data = map(data, ~ .x %>% 
                          nest(new_challenger = - contrast_pretty) %>% 
                          mutate(new_challenger = map_chr(
                            new_challenger, 
                            ~.x %>% pull(symbol) %>% unique()
                          ))
      )) %>% 
      
      unnest(data) %>% 
      
      mutate(challengers_for_silhouette = map2(
        new_challenger, ancestor, 
        ~ with(signature, cumulative_signature[ancestor==.y][[1]]) %>% 
          append(.x)
      )) %>% 
      
      mutate(silhouette = map2(
        challengers_for_silhouette, ancestor, 
        ~ silhouette_for_markers(.x, .y, .contrast, .level, .method) %>% 
          select(-ancestor)
      )) %>% 
      
      unnest(silhouette) %>% 
      
      group_by(ancestor) %>% 
      
      arrange(desc(silhouette), .by_group = TRUE) %>% 
      
      ungroup() %>% 
      
      mutate(is_greater = map2_lgl(
        silhouette, ancestor, 
        ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
      )) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(new_challengers = map(
        data,
        ~ .x %>% 
          pull(new_challenger)
      )) %>% 
      
      mutate(winner = map(
        data,
        ~ if(.x[1, ]$is_greater){
          .x[1, ]$new_challenger
        } else {NA}
      )) %>% 
      
      mutate(cumulative_signature = pmap(
        list(data, winner, ancestor),
        ~ if(!is.na(..2)) {
          with(..1[1, ], challengers_for_silhouette[[1]])
        } else {
          with(signature, cumulative_signature[ancestor==..3][[1]])
        }
      ))
    
    
    # append the base + 1 markers that result in highest silhouette score
    signature <- signature %>% 
      
      mutate(cumulative_signature = map(
        ancestor,
        ~ with(contrast_pair_tb, cumulative_signature[ancestor==.x][[1]])
      )) %>% 
      
      mutate(last_silhouette = map2_dbl(
        ancestor, last_silhouette,
        ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x][[1]]))) {
          with(contrast_pair_tb, data[ancestor==.x][[1]][[1, "silhouette"]])
        } else {.y}
      ))
    
    summary_tb <- summary_tb %>% 
      bind_rows(
        contrast_pair_tb %>% 
          filter(!is.na(winner)) %>% 
          # mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>% 
          mutate(winning_silhouette = map_dbl(data, ~ .x$silhouette[1])) %>% 
          select(-data)
      )
    
    contrast_copy <- contrast_copy %>% 
      mutate(markers = map2(
        markers, ancestor, 
        ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y][[1]]))){
          .x %>% 
            filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))
        } else {
          .x %>% 
            filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]]))
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
    cat("genes selected for each node: ", map_int(signature$cumulative_signature, ~ length(.x)),  "\n")
    
    i <- i + 1
    cat("iteration: ", i, "\n")
    
  }
  
  output <- summary_tb %>% 
    nest(signature_data = - ancestor)
  
  return(output)
}

# calculates silhouette score for each set of signature (cumulative markers) at a signature size
silhouette_for_markers <-function(.signature, .ancestor, .contrast, .level, .method) {
  
  .contrast %>%
    
    mutate(ancestor = !!as.symbol(pre(.level))) %>% 
    
    filter(ancestor == .ancestor) %>% 
    
    # filter markers that are in the signature
    mutate(markers = map2(markers, ancestor, ~.x %>% 
                            filter(symbol %in% .signature))) %>% 
    
    # format statistics from pairwise contrast
    mutate(markers  = map(markers, 
                          ~ .x %>% 
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
                            mutate(stat_df = map(stat_df, ~.x %>% 
                                                   pivot_wider(names_from = stats, 
                                                               values_from = .value))) %>% 
                            
                            unnest(stat_df) )) %>% 
    
    # Add original data data to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>% 
    
    # select only columns needed
    select(-data) %>% 
    
    nest(silhouette_data = c(!!as.symbol(pre(.level)), markers)) %>% 
    
    mutate(silhouette_data = map(silhouette_data, ~ .x %>% unnest(markers))) %>% 
    
    # Calculate silhouette score for PCA plot resulted from the markers selected
    mutate(silhouette_data = map(silhouette_data, ~ sil_func(.x, .level, .method))) %>% 
    
    mutate(reduced_dimensions = map(silhouette_data, ~ .x$rdim[[1]])) %>% 
    mutate(silhouette = map_dbl(silhouette_data, ~ .x$sil)) %>% 
    select(-silhouette_data)
  
}

# Non-hierarchical

single_marker_pw_select0 <- function(.contrast, .discard_num, .method) {
  
  # initialize variables
  
  contrast_copy <- .contrast %>% 
    mutate(ancestor = level_0)
  
  # initialise a signature tibble to store signature markers for each cell type in each iteration
  signature <- tibble(
    ancestor = contrast_copy %>% 
      pull(ancestor)
  ) %>% 
    mutate(cumulative_signature = map(ancestor, ~ vector())) %>% 
    mutate(last_silhouette = 0)
  
  # initialise an output tibble containing all results of interest
  summary_tb <- tibble(
    ancestor = character(),
    new_challengers = list(),
    winner = list(),
    cumulative_signature = list(),
    # reduced_dimensions = list(),
    winning_silhouette = double()
  )
  
  # set the base markers
  contrast_pair_tb0 <-
    
    # contrast_copy contains all the statistics of all cell_type contrasts for each gene
    contrast_copy %>%
    
    # select top 1 markers from each contrast, output is an unnested tibble
    sig_select0(1) %>%
    
    mutate(ancestor = level_0) %>%
    
    nest(data = - ancestor) %>%
    
    mutate(data = map(data, ~ .x %>% 
                        nest(new_challenger = - contrast_pretty) %>% 
                        
                        mutate(new_challenger = map(
                          new_challenger, ~ .x %>% pull(symbol) %>% unique()))
    )) %>% 
    
    mutate(new_challengers = map(data, ~.x %>% pull(new_challenger) %>% unlist())) %>% 
    
    mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
    
    mutate(cumulative_signature = winner) %>% 
    
    select(-data) %>% 
    
    mutate(silhouette = map2(
      new_challengers, ancestor,
      ~ silhouette_for_markers0(.x, .y, .contrast, .method) %>% 
        select(
          # reduced_dimensions, 
          winning_silhouette = silhouette)
    )) %>% 
    
    unnest(silhouette)
  
  
  signature <- signature %>%
    mutate(cumulative_signature = map2(
      cumulative_signature, ancestor,
      ~ .x %>% 
        append(with(contrast_pair_tb0, cumulative_signature[ancestor==.y][[1]]))
    )) %>% 
    mutate(last_silhouette = map_dbl(
      ancestor,
      ~ with(contrast_pair_tb0, winning_silhouette[ancestor==.x])
    ))
  
  summary_tb <- summary_tb %>% 
    bind_rows(contrast_pair_tb0)
  
  # remove base markers from contrast_copy input before further selection
  contrast_copy <- contrast_copy %>%
    mutate(markers = map2(
      markers, ancestor, 
      ~ .x %>%
        filter(!symbol %in% with(signature, cumulative_signature[ancestor==.y][[1]]))
    ))
  
  # counter for number of genes discarded
  j <- map_int(signature$cumulative_signature, ~ length(.x))
  
  # count the number of iterations
  i <- 0
  while (any(j < .discard_num) &
         # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
         all(map_int(contrast_copy$markers, 
                     # hence the boundary should be the number of satisfactory genes selected
                     ~ .x %>% select_markers_for_each_contrast(1) %>% nrow()) > 0)) {
    
    contrast_pair_tb <- 
      
      # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select0(1) %>% 
      
      mutate(ancestor = level_0) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(data = map(data, ~ .x %>% 
                          nest(new_challenger = - contrast_pretty) %>% 
                          mutate(new_challenger = map_chr(
                            new_challenger, 
                            ~.x %>% pull(symbol) %>% unique()
                          ))
      )) %>% 
      
      unnest(data) %>% 
      
      mutate(challengers_for_silhouette = map2(
        new_challenger, ancestor, 
        ~ with(signature, cumulative_signature[ancestor==.y][[1]]) %>% 
          append(.x)
      )) %>% 
      
      mutate(silhouette = map2(
        challengers_for_silhouette, ancestor, 
        ~ silhouette_for_markers0(.x, .y, .contrast, .method) %>% 
          select(-ancestor)
      )) %>% 
      
      unnest(silhouette) %>% 
      
      group_by(ancestor) %>% 
      
      arrange(desc(silhouette), .by_group = TRUE) %>% 
      
      ungroup() %>% 
      
      mutate(is_greater = map2_lgl(
        silhouette, ancestor, 
        ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
      )) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(new_challengers = map(
        data,
        ~ .x %>% 
          pull(new_challenger)
      )) %>% 
      
      mutate(winner = map(
        data,
        ~ if(.x[1, ]$is_greater){
          .x[1, ]$new_challenger
        } else {NA}
      )) %>% 
      
      mutate(cumulative_signature = pmap(
        list(data, winner, ancestor),
        ~ if(!is.na(..2)) {
          with(..1[1, ], challengers_for_silhouette[[1]])
        } else {
          with(signature, cumulative_signature[ancestor==..3][[1]])
        }
      ))
    
    
    # append the base + 1 markers that result in highest silhouette score
    signature <- signature %>% 
      
      mutate(cumulative_signature = map(
        ancestor,
        ~ with(contrast_pair_tb, cumulative_signature[ancestor==.x][[1]])
      )) %>% 
      
      mutate(last_silhouette = map2_dbl(
        ancestor, last_silhouette,
        ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x][[1]]))) {
          with(contrast_pair_tb, data[ancestor==.x][[1]][[1, "silhouette"]])
        } else {.y}
      ))
    
    summary_tb <- summary_tb %>% 
      bind_rows(
        contrast_pair_tb %>% 
          filter(!is.na(winner)) %>% 
          # mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>% 
          mutate(winning_silhouette = map_dbl(data, ~ .x$silhouette[1])) %>% 
          select(-data)
      )
    
    contrast_copy <- contrast_copy %>% 
      mutate(markers = map2(
        markers, ancestor, 
        ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y][[1]]))){
          .x %>% 
            filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))
        } else {
          .x %>% 
            filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]]))
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
    cat("genes selected for each node: ", map_int(signature$cumulative_signature, ~ length(.x)),  "\n")
    
    i <- i + 1
    cat("iteration: ", i, "\n")
    
  }
  
  output <- summary_tb %>% 
    nest(signature_data = - ancestor)
  
  return(output)
}

# calculates silhouette score for each set of signature (cumulative markers) at a signature size
silhouette_for_markers0 <-function(.signature, .ancestor, .contrast, .method) {
  
  .contrast %>%
    
    mutate(ancestor = level_0) %>% 
    
    filter(ancestor == .ancestor) %>% 
    
    # filter markers that are in the signature
    mutate(markers = map2(markers, ancestor, ~.x %>% 
                            filter(symbol %in% .signature))) %>% 
    
    # format statistics from pairwise contrast
    mutate(markers  = map(markers, 
                          ~ .x %>% 
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
                            mutate(stat_df = map(stat_df, ~.x %>% 
                                                   pivot_wider(names_from = stats, 
                                                               values_from = .value))) %>% 
                            
                            unnest(stat_df) )) %>% 
    
    # Add original data data to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>% 
    
    # select only columns needed
    select(-data) %>% 
    
    nest(silhouette_data = c(level_0, markers)) %>% 
    
    mutate(silhouette_data = map(silhouette_data, ~ .x %>% unnest(markers))) %>% 
    
    # Calculate silhouette score for PCA plot resulted from the markers selected
    mutate(silhouette_data = map(silhouette_data, ~ sil_func0(.x, .method))) %>% 
    
    mutate(reduced_dimensions = map(silhouette_data, ~ .x$rdim[[1]])) %>% 
    mutate(silhouette = map_dbl(silhouette_data, ~ .x$sil)) %>% 
    select(-silhouette_data)
  
}

# Testing =====

# pw_markers_L3 gives the full output including PC values from reduced_dimensions but this file is too big
pw_markers_L2 <- contrast_PW_L2 %>% 
  single_marker_pw_select(LEVEL, .discard_num = 1000, METHOD)

# pw_markers_L3_sub removes the reduced_dimensions column from the output. 
# Use silhouette_for_markers function to calculate PC values from any specific signature for PCA plot.
pw_markers_L1_sub <- pw_markers_L1 %>% 
  mutate(signature_data = map(
    signature_data,
    ~ .x %>% select(-reduced_dimensions)
  ))

saveRDS(pw_markers_L2, "pw_markers_L2.rds")

# pw_markers_L1 trend plots & PCA ==================

# trend plot t_CD4
pw_markers_L4 %>% 
  pluck("signature_data", 1) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(10)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("t_CD4, 3 marker selected")

ggsave("t_CD4.png", x)

# PCA plot t_CD4
x <- pw_markers_L4 %>% 
  slice(1) %>% 
  unnest(signature_data) %>% 
  tail(1) %>%
  mutate(silhouette_data = map2(
    cumulative_signature, ancestor,
    ~ silhouette_for_markers(.x, .y, contrast_PW_L4, LEVEL, METHOD)
  )) %>% 
  pluck("silhouette_data", 1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA t_CD4, 3 markers selected, silhouette = 0.952")


# trend plot t_cell
pw_markers_L3.2 %>% 
  pluck("signature_data", 2) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("t_cell, 97 marker selected")

## PCA plot t_cell
pw_markers_L3.2 %>% 
  pluck("signature_data", 2) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA t_cell, 97 markers selected, sil_score=0.71")

# trend plot granulocyte
pw_markers_L3.2 %>% 
  pluck("signature_data", 3) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("granulocyte, 11 marker selected")


## PCA plot granulocyte
pw_markers_L3.2 %>% 
  pluck("signature_data", 3) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA granulocyte, 11 markers selected, sil_score=0.689")

# trend plot b_cell
pw_markers_L3.2 %>% 
  pluck("signature_data", 4) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("b_cell, 46 marker selected")

## PCA plot b_cell
pw_markers_L3.2 %>% 
  pluck("signature_data", 4) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA b cell, 46 markers selected, sil_score=0.836")

# trend plot NK
pw_markers_L3.2 %>% 
  pluck("signature_data", 5) %>% 
  mutate(signature_size = map_int(
    cumulative_signature,
    ~ length(.x)
  )) %>% 
  # tail(1)
  ggplot(aes(signature_size, winning_silhouette)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("nk_cell, 61 marker selected")

## PCA plot NK
pw_markers_L3.2 %>% 
  pluck("signature_data", 5) %>% 
  tail(1) %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("PCA nk cell, 61 markers selected, sil_score=0.92")

