

#' sccomp_glm main
#'
#' @description The function for linear modelling takes as input a table of cell counts with three columns containing a cell-group identifier, sample identifier, integer count and the covariates (continuous or discrete). The user can define a linear model with an input R formula, where the first covariate is the factor of interest. Alternatively, sccomp accepts single-cell data containers (Seurat, SingleCellExperiment44, cell metadata or group-size). In this case, sccomp derives the count data from cell metadata.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#' @importFrom cmdstanr cmdstan_model
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column (optional depending on the input class) | covariate columns.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each covariates.
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .count A column name as symbol. The cell_group abundance (read count). Used only for data frame count output. The variable in this column should be of class integer.
#'
#' @param prior_mean_variable_association A list of the form list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(5,8)). Where for intercept and slope parameters, we specify mean and standard deviation, while for standard deviation, we specify shape and rate. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param bimodal_mean_variability_association A boolean. Whether to model the mean-variability as bimodal, as often needed in the case of single-cell RNA sequencing data, and not usually for CyTOF and microbiome data. The plot summary_plot()$credible_intervals_2D can be used to assess whether the bimodality should be modelled.
#'
#' @param percent_false_positive A real between 0 and 100 non included. This used to identify outliers with a specific false positive rate.
#' @param exclude_priors A boolean. Whether to run a prior-free model, for benchmarking purposes.
#' @param use_data A booelan. Whether to sun the model data free. This can be used for prior predictive check.
#' @param max_sampling_iterations An integer. This limit the maximum number of iterations in case a large dataset is used, for limiting the computation time.
#' @param pass_fit A boolean. Whether to pass the Stan fit as attribute in the output. Because the Stan fit can be very large, setting this to FALSE can be used to lower the memory imprint to save the output.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param verbose A boolean. Prints progression.
#' @param noise_model A character string. The two noise models available are multi_beta_binomial (default) and dirichlet_multinomial.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#'
#' @return A nested tibble `tbl`, with the following columns
#' \itemize{
#'   \item cell_group - column including the cell groups being tested
#'   \item parameter - The parameter being estimated, from the design matrix dscribed with the input formula_composition and formula_variability
#'
#'   \item c_lower - lower (2.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_effect - mean of the posterior distribution for a composition (c) parameter.
#'   \item c_upper - upper (97.5%) quantile of the posterior distribution fo a composition (c)  parameter.
#'   \item c_pH0 - Probability of the null hypothesis (no difference) for  a composition (c). This is not a p-value.
#'   \item c_FDR - False-discovery rate of the null hypothesis (no difference) for  a composition (c).
#'
#'   \item v_lower - (optional, present if variability is modelled dependent on covariates) lower (2.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_effect - (optional, present if variability is modelled dependent on covariates) mean of the posterior distribution for a variability (v) parameter
#'   \item v_upper - (optional, present if variability is modelled dependent on covariates) upper (97.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_pH0 - (optional, present if variability is modelled dependent on covariates) Probability of the null hypothesis (no difference) for a variability (v). This is not a p-value.
#'   \item v_FDR - (optional, present if variability is modelled dependent on covariates) False-discovery rate of the null hypothesis (no difference), for a variability (v).
#' }
#'
#' @examples
#'
#' data("counts_obj")
#'
#' estimate =
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type,
#'    ~1,
#'    sample,
#'    cell_group,
#'    count,
#'     approximate_posterior_inference = "all",
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#' @export
#'
#'
cellsig_multilevel_varing_intercept <- function(.data,
                                                .sample,
                                                .feature,
                                                .abundance,
                                                .cell_group,
                                                .scaling_multiplier,
                                                .multilevel_grouping,
                                                
                                                # Other parameters
                                                cores = detectCores(),
                                                priors = list(
                                                  assoc_intercept_mean = 1,
                                                  assoc_slope_mean = -0.55,
                                                  assoc_sd_sd_mean = 1.22,
                                                  assoc_sd_shape_mean = 1.14,
                                                  lambda_mu_mean = 9.8,
                                                  lambda_sigma_mean = 2,
                                                  lambda_skew_mean = -5
                                                ),
                                                iterations_warmup = 250,
                                                iterations_sampling = 300,
                                                pass_fit = FALSE,
                                                use_cmdstanr = FALSE) {
  UseMethod("cellsig_multilevel_varing_intercept", .data)
}

#' @importFrom rstan rstan_options
#' @importFrom forcats fct_relevel
#' @export
cellsig_multilevel_varing_intercept.data.frame = function(
                            .data,
                             .sample,
                             .feature,
                             .abundance,
                            .cell_group,
                             .scaling_multiplier,
                             .multilevel_grouping,
                             
                             # Other parameters
                            cores = detectCores(),
                            priors = list(
                              assoc_intercept_mean = 1,
                              assoc_slope_mean = -0.55,
                              assoc_sd_sd_mean = 1.22,
                              assoc_sd_shape_mean = 1.14,
                              lambda_mu_mean = 9.8,
                              lambda_sigma_mean = 2,
                              lambda_skew_mean = -5
                            ),
                             iterations_warmup = 250,
                             iterations_sampling = 300,
                            pass_fit = FALSE,
                            use_cmdstanr = FALSE) {
  
  # Use quotation for column names
  .sample = enquo(.sample)
  .feature = enquo(.feature)
  .abundance = enquo(.abundance)
  .cell_group = enquo(.cell_group)
  .scaling_multiplier = enquo(.scaling_multiplier)
  .multilevel_grouping = enquo(.multilevel_grouping)
  
  # Chains and cores
  chains = 3
  cores = max(cores, chains)
  
  # Iterations per chain
  iterations_sampling_per_chain = iterations_sampling / chains
  
  # Add utility columns
  .data = 
    .data %>% 
    
    # Unique database cell type
    unite("database_for_cell_type_feature", c(!!.multilevel_grouping, !!.cell_group, !!.feature), remove = FALSE ) %>% 
    unite("feature_cell_type", c(!!.feature, !!.cell_group), remove = FALSE) %>% 
    
    # !!! Needed for new model formatting - DO NOT CHANGE OTHERWISE THE MODEL WILL PRODUCE MEANINGLESS RESULTS
    arrange(feature_cell_type, database_for_cell_type_feature, !!.sample) %>% 
    mutate(feature_cell_type = factor(feature_cell_type)) %>% 
    mutate(database_for_cell_type_feature = fct_relevel(database_for_cell_type_feature, unique(as.character(database_for_cell_type_feature)))) %>% 
    mutate(!!.sample := as.factor(!!.sample)) %>% 
    
    
    # Count scaled
    mutate(count_scaled = !!.abundance * !!.scaling_multiplier) %>% 
    
    # Exposure rate
    mutate(exposure_rate = -log(!!.scaling_multiplier) ) 
    
  # Build model input
  model_data = list(
    N = nrow(.data),
    Y = .data %>% pull(!!.abundance),
    exposure_rate =
      .data %>%
      distinct(!!.sample, exposure_rate) %>%
      arrange(!!.sample) %>%
      pull(exposure_rate),
    
    grouping_gene_idx_N  =
      .data %>%
      select(database_for_cell_type_feature, feature_cell_type,!!.sample),
    
    G = .data %>% pull(feature_cell_type) %>% unique() %>%  length(),
    D = length(unique(.data$database_for_cell_type_feature)),
    S = length(unique(.data %>% pull(!!.sample))),
    
    grouping_gene_idx_D =
      .data %>%
      select(database_for_cell_type_feature, feature_cell_type) %>% 
      distinct(),
    
    grainsize = 1
  ) %>%
    
    # Add priors
    c(priors) %>%
    
    # Add offset
    c(
      list(
        # Offsets
        gene_mean_offset =
          .$grouping_gene_idx_N %>%
          bind_cols(Y = .data %>% pull(count_scaled)) %>%
          with_groups(feature_cell_type, ~ summarise(.x, m = mean(log1p(
            Y
          )))) %>%
          pull(m)
      )
    )
  
  init_fun <- function(...) list(
    #gene_mean = data$grouping_gene_idx_N %>% bind_cols(Y=.data$count_scaled) %>% with_groups(.feature_cell_type, ~ summarise(.x, m = mean(log1p(Y))))%>% pull(m) 
    gene_sd_alpha = 3,
    gene_sd_beta = 3
  )
  
  init_fun_vb = function(...) list(
    gene_mean = fit_vb %>% summary("gene_mean") %$% summary %>% .[,1],
    gene_sd = fit_vb %>% summary("gene_sd") %$% summary %>% .[,1],
    gene_sd = fit_vb %>% summary("shape") %$% summary %>% .[,1]
  )
  
  vb = FALSE
  
  # RSTAN
  if(!use_cmdstanr) {
    
    # Set global options
    rstan_options(threads_per_chain = ceiling(cores / chains))
    
    # Sample
    if(vb) fit = (
      vb_iterative(
        stanmodels$mixed_effect,
        iter = 10000,
        tol_rel_obj = 0.01,
        data = model_data,
        init = init_fun
      )
      
    )
    else fit = sampling(
      stanmodels$mixed_effect,
      data = model_data ,
      cores = chains,
      chains = chains,
      iter = iterations_sampling_per_chain + iterations_warmup,
      warmup = iterations_warmup,
      init = init_fun
    )
    
    rng =  rstan::gqs(
      stanmodels$mixed_effect_generate,
      draws =  as.matrix(fit),
      data = model_data
    )
    
    rng_summary = 
      rng %>% 
      rstan::summary("Y_gen", c(0.1, 0.5, 0.9)) %$%
      summary %>% 
      as_tibble() %>% 
      rowid_to_column(var = ".feature_idx")
  }
 else {
   
   # Lad model code
   if(file.exists("mixed_effect_cmdstanr.rds"))
     mod = readRDS("mixed_effect_cmdstanr.rds")
   else {
     write_file(mixed_effect_cmdstanr, "mixed_effect_cmdstanr.stan")
     mod = cmdstan_model( "mixed_effect_cmdstanr.stan", cpp_options = list(stan_threads = TRUE) )
     mod  %>% saveRDS("mixed_effect_cmdstanr.rds")
   }
   
   fit = 
     mod$sample(
     data = model_data ,
     chains = chains,
     parallel_chains = chains,
     threads_per_chain = ceiling(cores / chains),
     iter_warmup = iterations_warmup,
     iter_sampling = iterations_sampling_per_chain,
     #refresh = ifelse(verbose, 1000, 0),
     save_warmup = FALSE,
     init = list(list(gene_sd_alpha = 3, gene_sd_beta = 3), list(gene_sd_alpha = 3, gene_sd_beta = 3), list(gene_sd_alpha = 3, gene_sd_beta = 3)),
     output_dir = "."
   ) %>%
     suppressWarnings()
   
   # Lad model code
   if(file.exists("mixed_effect_generate_cmdstanr.rds"))
     mod_generate = readRDS("mixed_effect_generate_cmdstanr.rds")
   else {
     write_file(mixed_effect_generate_cmdstanr, "mixed_effect_generate_cmdstanr.stan")
     mod_generate = cmdstan_model( "mixed_effect_generate_cmdstanr.stan" )
     mod_generate %>% saveRDS("mixed_effect_generate_cmdstanr.rds")
   }
   
   rng = mod_generate$generate_quantities(
     fit,
     data = model_data,
     parallel_chains = chains
   )
   
   rng_summary = 
     rng$summary("Y_gen", ~quantile(.x, probs = c(0.1, 0.5, 0.9))) %>% 
     rowid_to_column(var = ".feature_idx")
   
 }
  

  .data %>% 
    #mutate(.feature_idx = as.integer(.feature_cell_type)) %>% 
    mutate(.feature_idx = as.integer(feature_cell_type)) %>% 
    distinct(.feature_idx, !!.feature, !!.cell_group) %>% 
    left_join(rng_summary,  by = ".feature_idx"  ) %>% 
    
    left_join(
      rng$summary("Y_gen_log") %>% 
        .[,c("mean", "sd")] %>% 
        setNames(c("log_mean", "log_sd")) %>% 
        rowid_to_column(var = ".feature_idx")  ,
      by = ".feature_idx"
    ) %>% 
    
    # Add attributes
    add_attr(
      .data %>%
        mutate(feature_cell_type_idx = as.integer(feature_cell_type), database_for_cell_type_feature_idx = as.integer(database_for_cell_type_feature))%>% 
        select(!!.feature, !!.cell_group, !!.multilevel_grouping, feature_cell_type_idx, database_for_cell_type_feature_idx) %>% 
        distinct() ,
      "indeces"
    ) %>% 
    
    # Add attributes
    when(
      pass_fit ~ add_attr(., fit, "fit") %>% add_attr(rng, "rng"),
      ~ (.)
    )
  
}