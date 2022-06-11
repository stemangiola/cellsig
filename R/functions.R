#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom tibble enframe
generate_quantities_standalone = function(fit, G){
  
  
  rstan::gqs(
    stanmodels$generated_quantities,
    #rstan::stan_model("inst/stan/generated_quantities.stan"),
    draws =  as.matrix(fit),
    data = list(G)
  ) %>%
    
    rstan::extract("counts") %$% counts %>%
    as.data.frame() %>%
    setNames(1:G) %>%
    as_tibble(rownames = "draw") %>%
    gather(G, generated_quantity, -draw) %>%
    mutate(G = as.integer(G)) %>%
    nest(data = -G) %>%
    mutate(quantiles = map(
      data, 
      ~ quantile(
        .x$generated_quantity, 
        probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
      ) %>% 
        enframe() %>% 
        spread(name, value)
    )) %>%
    unnest(quantiles) %>%
    select(-data)
  
}

#' ref_intercept_only
#'
#' @description This function calls the stan model.
#'
#'
#' @importFrom tibble tibble
#'
#' @import dplyr
#' @import tidyr 
#' @import purrr
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes median_qi
#' @import tidybayes
#'
#' @import stringr
#'
#' @import data.tree
#'
#' @param mix A matrix
#' @param my_design A matrix
#' @param cov_to_test A character string
#' @param fully_bayesian A boolean
#' @param is_mix_microarray A boolean
#' @param verbose A boolean
#' @param save_report A boolean
#' @param custom_ref A matrix
#' @param multithread A boolean
#' @param do_debug A boolean
#' @param cell_type_root A character string
#' @param choose_internal_ref A design matrix
#' @param omit_regression A boolean
#' @param save_fit A boolean
#' @param seed An integer
#'
#' @return An ARMET object
#'
#' @export
#'
ref_intercept_only = function(reference,
                              sample_abundance_multiplier,
                              cores = 8,
                              approximate_posterior = F
) {

  exposure_rate_col = enquo(sample_abundance_multiplier)

  # Non centred
  lambda_mu_prior = c(8, 2)
  lambda_sigma_prior =  c(log(3.3) , 1)
  lambda_skew_prior =  c(-2.7, 2)
  sigma_intercept_prior = c(1.9 , 0.5)
  
  res1 = 
    
    run_model_ref(
      reference,
      !!exposure_rate_col,
      cores,
      approximate_posterior
      #iterations = iterations,
      #sampling_iterations = sampling_iterations
    )

  G = res1[[1]] %>% distinct(cell_type, .feature, G) %>% nrow
  
  generated_quantities = 
    
    # Run model
    generate_quantities_standalone( res1[[2]], G  ) 
  
  res1[[1]] %>% 
    nanny::subset(c(.feature, cell_type)) %>%
    select(-starts_with("level_")) %>%

    # Attach lambda sigma
    left_join(  generated_quantities,    by = c("G")  ) %>%
    select(-GM) 

}

#' Add attribute to abject
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}


#' @importFrom rstan vb
#' @importFrom rstan sampling
#'
#' @export
#'
run_model_ref = function(
                         reference_filtered,
                         exposure_rate_col,
                         shards,
                         approximate_posterior,
                         iterations = 250,
                         sampling_iterations = 100) {

  
  
  lambda_mu_prior = c(8, 2)
  lambda_sigma_prior =  c(log(3.3) , 1)
  lambda_skew_prior =  c(-2.7, 2)
  sigma_intercept_prior = c(1.9 , 0.5)
  
  exposure_rate_col = enquo(exposure_rate_col)
  
  df = ref_format(reference_filtered ) 

  G = df %>% distinct(G) %>% nrow()
  GM = df %>% distinct(.feature) %>% nrow()
  S = df %>% distinct(.sample) %>% nrow()
  CL = df %>% nrow

  counts_linear = df %>%  pull(count)
  G_to_counts_linear = df %>% pull(G)
  exposure_rate = df %>% pull(!!exposure_rate_col)
  

  # library(rstan)
  # fileConn<-file("~/.R/Makevars")
  # writeLines(c( "CXX14FLAGS += -O2","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
  # close(fileConn)
  # ARMET_tc_model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_ref.stan", auto_write = F)

  if(shards > 1) Sys.setenv("STAN_NUM_THREADS" = shards)
  
  list(df,
       approximate_posterior %>%
         when(
           # Variational
           (.) == TRUE ~  vb_iterative(
             cellsig:::stanmodels$ARMET_ref,
             #rstan::stan_model("inst/stan/ARMET_ref.stan"),
             output_samples = 5000,
             iter = 50000,
             tol_rel_obj = 0.01,
             algorithm = "meanfield"
             #,
             #save_warmup = FALSE
           ),  
             
          # HMC
           ~ sampling(
             cellsig:::stanmodels$ARMET_ref,
           #rstan::stan_model("inst/stan/ARMET_ref.stan"),
           chains = 3,
           cores = 3,
           iter = iterations,
           warmup = iterations - sampling_iterations,
           save_warmup = FALSE
         ) %>%
         {
           (.)  %>% rstan::summary() %$% summary %>% as_tibble(rownames = "par") %>% arrange(Rhat %>% desc) %>% print
           (.)
         }
       ))

}

#' @importFrom rstan vb
#' @importFrom rstan sampling
#'
#' @export
#'
get_mean_overdispersion_association = function(
  reference_filtered,
  exposure_rate_col,
  shards,
  approximate_posterior,
  iterations = 250,
  sampling_iterations = 100) {
  
  
  exposure_rate_col = enquo(exposure_rate_col)
  
  # Make cells all the same
  reference_filtered = 
    reference_filtered %>%
    
    # Subset samples
    nest(data = -.sample) %>% 
    sample_n(min(30, n())) %>% 
    unnest(data) %>% 
    
    # Subset genes
    nest(data = -.feature) %>% 
    sample_n(min(10000, n())) %>% 
    unnest(data) 
  
  
  df = ref_format(reference_filtered ) 
  
  G = df %>% distinct(G) %>% nrow()
  GM = df %>% distinct(.feature) %>% nrow()
  S = df %>% distinct(.sample) %>% nrow()
  CL = df %>% nrow
  
  counts_linear = df %>%  pull(count)
  G_to_counts_linear = df %>% pull(G)
  exposure_rate = df %>% pull(!!exposure_rate_col)
  
  
  # library(rstan)
  # fileConn<-file("~/.R/Makevars")
  # writeLines(c( "CXX14FLAGS += -O2","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
  # close(fileConn)
  # ARMET_tc_model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_ref.stan", auto_write = F)
  
  if(shards > 1) Sys.setenv("STAN_NUM_THREADS" = shards)
  
  list(df,
       approximate_posterior %>%
         when(
           # Variational
           (.) == TRUE ~  vb_iterative(
             stanmodels$ARMET_ref,
             #rstan::stan_model("inst/stan/ARMET_ref.stan"),
             output_samples = 5000,
             iter = 50000,
             tol_rel_obj = 0.01,
             algorithm = "meanfield"
             #,
             #save_warmup = FALSE
           ),  
           
           # HMC
           ~ sampling(
             stanmodels$ARMET_ref,
             #rstan::stan_model("inst/stan/ARMET_ref.stan"),
             chains = 3,
             cores = 3,
             iter = iterations,
             warmup = iterations - sampling_iterations,
             save_warmup = FALSE
           ) %>%
             {
               (.)  %>% rstan::summary() %$% summary %>% as_tibble(rownames = "par") %>% arrange(Rhat %>% desc) %>% print
               (.)
             }
         ))
  
}


fit_mixed_effect = function(df,mod_estimate = stanmodels$mixed_effect, mod_rng = stanmodels$mixed_effect_generate,  assoc_intercept, assoc_slope, assoc_sd_sd, assoc_sd_shape, lambda_mu, lambda_sigma,lambda_skew,
                            iterations = 250,
                            sampling_iterations = 100, vb = FALSE){
  
  
  
  data = list(
    N = nrow(df),
    Y = df$count_scaled,
    grouping_gene_idx_N  = df %>% mutate(database_for_cell_type = factor(database_for_cell_type), .feature = factor(.feature)) %>% select(database_for_cell_type, .feature),
    
    G = length(unique(df$.feature)),
    D = length(unique(df$database_for_cell_type)),
    grouping_gene_idx_D =
      df %>%
      mutate(database_for_cell_type = factor(database_for_cell_type), .feature = factor(.feature)) %>%
      select(database_for_cell_type, .feature) %>%
      distinct(),
    grainsize=1
  )
  
  
  
  if(vb) fit = (
    vb_iterative(
      mod_estimate,
      iter = 10000,
      tol_rel_obj = 0.01,
      data = data
    )
    
  )
  else fit = (sampling(mod_estimate, data = data, cores = 4, iter = 600, warmup = 300))
  
  rng =  rstan::gqs(
    mod_rng,
    draws =  as.matrix(fit),
    data = data
  )
  
  return(list(fit = fit, rng = rng))
  
}



#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @export
#'
ref_format = function(ref) {
  # Get reference based on mix genes
  ref %>% 

    # Add marker symbol indexes
    nest(data = -c(.feature, cell_type)) %>%
    rowid_to_column(var = "G") %>%
    unnest(data) %>%

    # Add sample indeces
    nest(data = -.sample) %>%
    rowid_to_column(var = "S") %>%
    unnest(data) %>%
    
    # Add .sample indeces
    nest(data = -.feature) %>%
    rowid_to_column(var = "GM") %>%
    unnest(data) 

}

#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @export
infer_sequencing_depth_bias = function(counts, shards = 10, hk600){
  
  model_input = 
    counts %>%
    
    filter(level_1 %>% is.na %>% `!`) %>%
    mutate(cell_type = level_1) %>%
    
    # Get only house keeping genes
    filter(symbol %in% hk600) %>%
    
    # Add sample indexes
    nest(data = -sample) %>%
    rowid_to_column(var = "S") %>%
    unnest(data) %>%
    
    # Add feature indexes
    nest(data = -symbol) %>%
    rowid_to_column(var = "GM") %>%
    unnest(data) 
  
  # Input data
  GM = model_input %>% distinct(symbol) %>% nrow()
  S = model_input %>% distinct(sample) %>% nrow()
  CL = model_input %>% nrow
  
  counts_linear = model_input %>% pull(count)
  GM_to_counts_linear = model_input %>% pull(GM)
  S_linear = model_input %>% pull(S)
  
  # Non centered
  lambda_mu_prior = c(8, 2)
  lambda_sigma_prior =  c(log(3.3) , 1)
  lambda_skew_prior =  c(-2.7, 2)
  sigma_intercept_prior = c(1.9 , 0.5)
  
  Sys.setenv("STAN_NUM_THREADS" = shards)
  
  
  fit = 
    vb_iterative(
      #stanmodels$ARMET_ref,
      rstan::stan_model("inst/stan/infer_exposure.stan"),
      output_samples = 500,
      iter = 50000,
      tol_rel_obj = 0.01,
      algorithm = "meanfield", data = list(
        shards = 10,
        
        GM = model_input %>% distinct(symbol) %>% nrow(),
        S = model_input %>% distinct(sample) %>% nrow(),
        CL = model_input %>% nrow,
        
        counts_linear = model_input %>% pull(count),
        GM_to_counts_linear = model_input %>% pull(GM),
        S_linear = model_input %>% pull(S),
        
        # Non centered
        lambda_mu_prior = c(8, 2),
        lambda_sigma_prior =  c(log(3.3) , 1),
        lambda_skew_prior =  c(-2.7, 2),
        sigma_intercept_prior = c(1.9 , 0.5)
        )
    )
  

  fit %>%
    summary_to_tibble("exposure_rate", "S") %>%
    filter(.variable != "exposure_rate_minus_1") %>%
    left_join(
      model_input %>% distinct(sample, S),
      by="S"
    ) %>% 
    select(sample, exposure_rate = mean) %>%
    right_join(counts, by="sample")
}

#' @export
impute_abundance_using_levels = function(.data, .abundance){
  
# Input dataset
#   Rows: 28,455,434
#   Columns: 11
#   $ sample        <chr> "ENCFF060YNO", "ENCFF060YNO", "ENCFF060YNO", "ENCFF060YNO", "ENCFF060YNO", "ENCFF06…
#   $ exposure_rate <dbl> -0.5029584, -0.5029584, -0.5029584, -0.5029584, -0.5029584, -0.5029584, -0.5029584,…
#   $ cell_type     <chr> "dendritic_myeloid", "dendritic_myeloid", "dendritic_myeloid", "dendritic_myeloid",…
#   $ count         <int> 20, 26, 0, 27635, 1, 0, 0, 0, 85, 0, 0, 78, 65, 0, 0, 0, 0, 0, 0, 0, 0, 398, 146, 1…
#   $ symbol        <chr> "A1BG", "A1BG-AS1", "A1CF", "A2M", "A2M-AS1", "A2ML1", "A2MP1", "A3GALT2", "A4GALT"…
#   $ database      <chr> "ENCODE", "ENCODE", "ENCODE", "ENCODE", "ENCODE", "ENCODE", "ENCODE", "ENCODE", "EN…
#   $ level_1       <chr> "immune_cell", "immune_cell", "immune_cell", "immune_cell", "immune_cell", "immune_…
#   $ level_2       <chr> "mono_derived", "mono_derived", "mono_derived", "mono_derived", "mono_derived", "mo…
#   $ level_3       <chr> "dendritic_myeloid", "dendritic_myeloid", "dendritic_myeloid", "dendritic_myeloid",…
#   $ level_4       <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
#   $ level_5       <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…

  
  # counts_imputed = 
  #   counts %>% 
  #   mutate(count_scaled = count / exp(exposure_rate)) %>% 
  #   impute_abundance_using_levels(count_scaled)
  # .abundance = enquo(.abundance)
  
  .abundance = enquo(.abundance)
  
  .data %>%
  #mutate(count_scaled = count * exp(exposure_rate)) %>%
  complete(nesting(sample, exposure_rate, cell_type, database, level_1 , level_2 , level_3 ,  level_4 ,level_5), symbol) %>%
  
  # Complete last level
  nest(data = -c(symbol, cell_type)) %>%
  mutate(count_median = map_dbl(data, ~  .x %>% pull(!!.abundance) %>%  median(na.rm = T) )) %>%
  mutate(count_median = case_when(!is.na(cell_type) ~ count_median)) %>%
  unnest(data) %>%
  mutate(!!.abundance := case_when(is.na(!!.abundance) ~ count_median, TRUE ~ !!.abundance)) %>%
  select(-count_median) %>%
  
  # Complete level 5
  nest(data = -c(symbol, level_5 )) %>%
  mutate(count_median = map_dbl(data, ~  .x %>% pull(!!.abundance) %>%  median(na.rm = T) )) %>%
  mutate(count_median = case_when(!is.na(level_5 ) ~ count_median)) %>%
  unnest(data) %>%
  mutate(!!.abundance := case_when(is.na(!!.abundance) ~ count_median, TRUE ~ !!.abundance)) %>%
  select(-count_median) %>%
  
  # Complete level 4
  nest(data = -c(symbol, level_4 )) %>%
  mutate(count_median = map_dbl(data, ~  .x %>% pull(!!.abundance) %>%  median(na.rm = T) )) %>%
  mutate(count_median = case_when(!is.na(level_4 ) ~ count_median)) %>%
  unnest(data) %>%
  mutate(!!.abundance := case_when(is.na(!!.abundance) ~ count_median, TRUE ~ !!.abundance)) %>%
  select(-count_median) %>%
  
  # Complete level 3
  nest(data = -c(symbol, level_3 )) %>%
  mutate(count_median = map_dbl(data, ~  .x %>% pull(!!.abundance) %>%  median(na.rm = T) )) %>%
  mutate(count_median = case_when(!is.na(level_3 ) ~ count_median)) %>%
  unnest(data) %>%
  mutate(!!.abundance := case_when(is.na(!!.abundance) ~ count_median, TRUE ~ !!.abundance)) %>%
  select(-count_median) %>%
  
  # Complete level 2
  nest(data = -c(symbol, level_2 )) %>%
  mutate(count_median = map_dbl(data, ~  .x %>% pull(!!.abundance) %>%  median(na.rm = T) )) %>%
  mutate(count_median = case_when(!is.na(level_2 ) ~ count_median)) %>%
  unnest(data) %>%
  mutate(!!.abundance := case_when(is.na(!!.abundance) ~ count_median, TRUE ~ !!.abundance)) %>%
  select(-count_median) %>%
  
  # Complete level 1
  nest(data = -c(symbol, level_1 )) %>%
  mutate(count_median = map_dbl(data, ~  .x %>% pull(!!.abundance) %>%  median(na.rm = T) )) %>%
  mutate(count_median = case_when(!is.na(level_1 ) ~ count_median)) %>%
  unnest(data) %>%
  mutate(!!.abundance := case_when(is.na(!!.abundance) ~ count_median, TRUE ~ !!.abundance)) %>%
  select(-count_median)
}

#' @export
ToDataFrameTypeColFull = function(tree, fill = T, ...) {
  t = tree %>% data.tree::Clone()
  
  tree_df = 
    1:(t %$% Get("level") %>% max) %>%
    map_dfr(
      ~ data.tree::Clone(t) %>%
        {
          data.tree::Prune(., function(x)
            x$level <= .x +1)
          .
        } %>%
        data.tree::ToDataFrameTypeCol() %>%
        as_tibble
      
    ) %>%
    distinct() 
  
  tree_df_filled = 
    tree_df %>%
    
    purrr::when(
      1 & ("level_2" %in% colnames(.)) ~ mutate(., level_2 = ifelse(level_2 %>% is.na, level_1, level_2)),
      TRUE ~ (.)
    ) %>%
    purrr::when(
      1 & ("level_3" %in% colnames(.)) ~ mutate(., level_3 = ifelse(level_3 %>% is.na, level_2, level_3)),
      TRUE ~ (.)
    ) %>%
    purrr::when(
      1 & ("level_4" %in% colnames(.)) ~ mutate(., level_4 = ifelse(level_4 %>% is.na, level_3, level_4)),
      TRUE ~ (.)
    ) %>%
    purrr::when(
      1 & ("level_5" %in% colnames(.)) ~ mutate(., level_5 = ifelse(level_5 %>% is.na, level_4, level_5)),
      TRUE ~ (.)
    ) %>%
    purrr::when(
      1 & ("level_6" %in% colnames(.)) ~ mutate(., level_6 = ifelse(level_6 %>% is.na, level_5, level_6)),
      TRUE ~ (.)
    ) %>%
    dplyr::select(..., everything())
  
  tree_df %>%
    select(-1) %>%
    setNames(tree_df %>% colnames %>% .[-ncol(tree_df)]) %>%
    mutate(cell_type = tree_df_filled %>% pull(ncol(tree_df)))
  
}

#' @export
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

#' @details This function takes a data frame as input, with a column containing cell type and outputs a one-level tree
#' @export
from_dataframe_to_one_level_tree = function(counts, cell_type_column){
  # load count data
  
  cell_type_column = enquo(cell_type_column)
  
  pseudo_tree_from_count_data = 
    counts %>%
    select(!!cell_type_column) %>%
    rename(cell_type = !!cell_type_column) %>% 
    unique() 
  
  pseudo_tree_from_count_data$pathString = 
    paste("Tissue",
          pseudo_tree_from_count_data$cell_type,
          sep = "/"
    )
  
  as.Node(pseudo_tree_from_count_data)
}