get_alpha = function(slope, which_changing, cell_types){
  
  # Get the alpha matrix
  
  intercept = rep(0, length(cell_types))
  slope_arr = rep(0, length(cell_types))
  
  slope_arr[which_changing] = slope
  matrix(intercept %>%	c(slope_arr), ncol = 2)
  
}

get_survival_X = function(S, PFI_all_cancers){
 
    PFI_all_cancers %>%
    filter(PFI.2 == 1 & !is.na(PFI.time.2) & PFI.time.2 > 0) %>%
    select(real_days = PFI.time.2 ) %>%
    mutate(real_days = real_days %>% scale(center = F) %>% as.numeric) %>%
    sample_n(S) %>%
    mutate(sample = sprintf("S%s", 1:n())) %>%
    mutate(alive = sample(0:1, n(), replace = T)) %>%
    mutate(days = ifelse(alive==1, real_days/2, real_days) ) %>%
    mutate(intercept = 1)
}

generate_mixture_from_model = function(.data, X_df, alpha, foreign_prop = 0, foreign_sample = NULL) {
  add_attr = function(var, attribute, name) {
    attr(var, name) <- attribute
    var
  }
  
  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }
  
  softmax <- function (x) {
    exp(x - logsumexp(x))
  }
  
  # Regress on the log days
  X = X_df %>% mutate(real_days = log(real_days)) %>% select(intercept, real_days) %>% nanny::as_matrix()
  
  samples_per_run =
    map_dfr(
      1:nrow(X), ~ 
        .data %>%
        distinct(cell_type, sample) %>%
        group_by(cell_type) %>%
        sample_n(1) %>%
        ungroup() %>%
        mutate(run = .x)
    )
  
  ct_names = .data %>% distinct(cell_type) %>% pull(1)
  
  alpha_df = alpha %>% as.data.frame %>% setNames(sprintf("alpha_%s", 1:2)) %>% mutate(cell_type  = ct_names)
  
  ct_changing = alpha_df %>% filter(alpha_2 != 0) %>% pull(cell_type)
  
  cell_type_proportions =
    # Choose samples
    samples_per_run %>%
    
    # Choose proportions
    left_join(
      # Decide theoretical, noise-less proportions for each sample
      X %*% t(alpha) %>%
        apply(1, softmax) %>%
        t %>%
        `*` (40) %>%
        as.data.frame() %>%
        as_tibble() %>%
        setNames(ct_names) %>%
        mutate(run = 1:n()) %>%
        gather(cell_type, alpha, -run)
    ) %>%
    
    # Add X
    left_join(X_df %>% select(-sample) %>% mutate(run = 1:n())) %>%
    
    # Add alpha
    left_join(alpha_df) %>%
    
    group_by(run) %>%
    mutate(p = gtools::rdirichlet(1, alpha) %>% as.numeric()) %>%
    ungroup()
  
  # Add fold increase decrease
  fold_change = 
    matrix(c(rep(1, 2), c(0, max(X,2))), ncol = 2)  %*% t(alpha) %>%
    apply(1, softmax) %>%
    t %>%
    `*` (40) %>%
    apply(1, softmax) %>%
    .[ct_names == ct_changing,] %>%
    {	max(.) / min(.)	} %>%
    { slope = alpha[,2][ alpha[,2]!=0]; ifelse(slope<0, -(.), (.)) }
  
  # Add counts
  dirichlet_source =
    cell_type_proportions %>%
    left_join(.data, by = c("cell_type", "sample"))
  
  # Make mix
  dirichlet_source %>%
    mutate(c = `count_scaled` * p) %>%
    group_by(run, symbol) %>%
    summarise(`count_mix` = c %>% sum) %>%
    ungroup %>%
    
    left_join(dirichlet_source %>% nanny::subset(run) ) %>%
    
    mutate(fold_change = fold_change) %>%
    
    # Add neuron
    when(
      !is.null(foreign_sample) ~  left_join(., foreign_sample, by = "symbol") %>%
        mutate(prop_neural = foreign_prop) %>%
        mutate(`count_mix` = (`count_mix` * (1-prop_neural))+(count_neuro*prop_neural)),
      ~ (.)
    ) %>%
    
    # Add proportions
    add_attr(cell_type_proportions, "proportions") 
  
  
}

#' @export
generate_mixture_from_proportion_matrix = function(.data, proportions, foreign_prop = 0, foreign_sample = NULL) {

    proportions %>% 
    as_tibble(rownames = "replicate") %>% 
    gather(cell_type, proportion, -replicate) %>%
  
    mutate(sample = map_chr(
      cell_type, 
      ~ !!.data %>% 
        filter(cell_type == .x) %>% 
        distinct(sample) %>% 
        pull(sample) %>% 
        sample(size = 1)
    )) %>%
    
    left_join(.data, by = c("cell_type", "sample")) %>%
    
    nest(data_samples = -c(replicate, symbol)) %>%
    mutate(count_mix = map_dbl(  data_samples,  ~ sum( .x$count_scaled * .x$proportion ) )) %>%
    
    # Add neuron
    when(
      !is.null(foreign_sample) ~  left_join(., foreign_sample, by = "symbol") %>%
        mutate(prop_neural = foreign_prop) %>%
        mutate(`count_mix` = (`count_mix` * (1-prop_neural))+(count_neuro*prop_neural)),
      ~ (.)
    ) 
  
  
}


#' @export
ref_plus_model_to_mixture = function(ref, S, slope, which_changing, PFI_all_cancers, foreign_prop = 0, foreign_sample = NULL){
  
  cell_types =  ref %>% pull(cell_type) %>% unique
  
  
  alpha = get_alpha(slope, which_changing, cell_types)
  X_df = get_survival_X(S, PFI_all_cancers)
  
  ref %>%
    generate_mixture_from_model(X_df, alpha, foreign_prop = foreign_prop, foreign_sample = foreign_sample) %>%
    mutate(`count_mix` = as.integer(`count_mix`), run = as.character(run)) %>%
    rename(sample = run) 
}


