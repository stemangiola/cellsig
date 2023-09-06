# ~/third_party_sofware/cctools-7.2.0-x86_64-centos7/bin/makeflow_monitor makefile.makeflow.makeflowlog
# ~/third_party_sofware/cctools-7.2.0-x86_64-centos7/bin/makeflow -T slurm -J 200  dev/TCGA_makeflow_pipeline/makefile_ARMET_TCGA.makeflow


library(rstan)
library(tidyverse)
#library(cmdstanr)

#library(cellsig)

# mod <- stan_model(stanc_ret = stanc_builder('~/PostDoc/cellsig/dev/mixed_effect.stan'),auto_write = T,verbose = T)
# gen_mod = stan_model(stanc_ret = stanc_builder("~/PostDoc/cellsig/dev/mixed_effect_generate.stan"),auto_write = T,verbose = T)
# 
# mod %>% saveRDS("~/PostDoc/cellsig/dev/benchmark_code/mod.rds")
# gen_mod %>% saveRDS("~/PostDoc/cellsig/dev/benchmark_code/gen_mod.rds")


#mod_cmdstanr = cmdstan_model( "~/PostDoc/cellsig/dev/mixed_effect.stan",  )

rstan_options(threads_per_chain = 5)

fit_mixed_effect = function(df,mod_estimate = stanmodels$mixed_effect,
                            mod_rng = stanmodels$mixed_effect_generate,
                            feature,
                            counts,
                            assoc_intercept,
                            assoc_slope,
                            assoc_sd_sd,
                            assoc_sd_shape,
                            lambda_mu,
                            lambda_sigma,
                            lambda_skew,
                            iterations = 250,
                            sampling_iterations = 100,
                            chains = 4,
                            vb = FALSE){
  feature = enquo(feature)
  counts = enquo(counts)
  
  data = list(
    N = nrow(df),
    Y = df %>% pull(!!counts),
    exposure_rate = df %>% distinct(.sample, exposure_rate) %>% arrange(.sample) %>% pull(exposure_rate),
    grouping_gene_idx_N  =
      df %>% 
      mutate(database_for_cell_type = factor(database_for_cell_type), !!feature := factor(!!feature)) %>%
      select(database_for_cell_type, !!feature, .sample),
    
    G = df %>% pull(!!feature) %>% unique() %>%  length(),
    D = length(unique(df$database_for_cell_type)),
    S = length(unique(df$.sample)),
    
    grouping_gene_idx_D =
      df %>%
      mutate(database_for_cell_type = factor(database_for_cell_type), !!feature := factor(!!feature)) %>%
      select(database_for_cell_type, !!feature) %>%
      distinct(),
    
    # grouping_sample_idx_N = df %>% pull(.sample),
    assoc_intercept = assoc_intercept,
    assoc_slope = assoc_slope,
    assoc_sd_sd= assoc_sd_sd,
    assoc_sd_shape= assoc_sd_shape,
    lambda_mu= lambda_mu,
    lambda_sigma = lambda_sigma,
    lambda_skew = lambda_skew,
    
    grainsize=1
  )
  
  init_fun <- function(...) list(
    #gene_mean = data$grouping_gene_idx_N %>% bind_cols(Y=df$count_scaled) %>% with_groups(.feature_cell_type, ~ summarise(.x, m = mean(log1p(Y))))%>% pull(m) 
    gene_sd_alpha = 3,
    gene_sd_beta = 3
  )
  
  init_fun_vb = function(...) list(
    gene_mean = fit_vb %>% summary("gene_mean") %$% summary %>% .[,1],
    gene_sd = fit_vb %>% summary("gene_sd") %$% summary %>% .[,1],
    gene_sd = fit_vb %>% summary("shape") %$% summary %>% .[,1]
  )
  
  if(vb) fit = (
    vb_iterative(
      mod_estimate,
      iter = 10000,
      tol_rel_obj = 0.01,
      data = data,
      init = init_fun
    )
    
  )
  else fit = sampling(
    mod_estimate,
    data = data %>%
      c(list(gene_mean_offset = data$grouping_gene_idx_N %>% bind_cols(Y=df$count_scaled) %>% with_groups(.feature_cell_type, ~ summarise(.x, m = mean(log1p(Y))))%>% pull(m)  )),
    cores = chains,
    chains = chains,
    iter = iterations,
    warmup = iterations - sampling_iterations,
    init = init_fun
  )
  
  rng =  rstan::gqs(
    mod_rng,
    draws =  as.matrix(fit),
    data = data
  )
  
  return(list(fit = fit, rng = rng))
  
}

# fit_mixed_effect = function(df,mod_estimate = stanmodels$mixed_effect, mod_rng = stanmodels$mixed_effect_generate,  assoc_intercept, assoc_slope, assoc_sd_sd, assoc_sd_shape, lambda_mu, lambda_sigma,lambda_skew,
#                             iterations = 250,
#                             sampling_iterations = 100, chains = 4, vb = FALSE){
#   
#   
#   
#   data = list(
#     N = nrow(df),
#     Y = df$count,
#     exposure_rate = df$exposure_rate,
#     grouping_gene_idx_N  = df %>% mutate(database_for_cell_type = factor(database_for_cell_type), .feature = factor(.feature)) %>% select(database_for_cell_type, .feature),
#     
#     G = length(unique(df$.feature)),
#     D = length(unique(df$database_for_cell_type)),
#     grouping_gene_idx_D =
#       df %>%
#       mutate(database_for_cell_type = factor(database_for_cell_type), .feature = factor(.feature)) %>%
#       select(database_for_cell_type, .feature) %>%
#       distinct(),
#     
#     assoc_intercept = assoc_intercept, 
#     assoc_slope = assoc_slope, 
#     assoc_sd_sd= assoc_sd_sd, 
#     assoc_sd_shape= assoc_sd_shape, 
#     lambda_mu= lambda_mu, 
#     lambda_sigma = lambda_sigma,
#     lambda_skew = lambda_skew,
#     
#     grainsize=1
#   )
#   
#   
#   
#   if(vb) fit = (
#     vb_iterative(
#       mod_estimate,
#       iter = 10000,
#       tol_rel_obj = 0.01,
#       data = data
#     )
#     
#   )
#   else fit = (sampling(mod_estimate, data = data, cores = chains, chains = chains, iter = iterations, warmup = iterations - sampling_iterations))
#   
#   rng =  rstan::gqs(
#     mod_rng,
#     draws =  as.matrix(fit),
#     data = data
#   )
#   
#   return(list(fit = fit, rng = rng))
#   
# }

# Read arguments
args = commandArgs(trailingOnly=TRUE)
file_in = args[1]
file_out = args[2]
cores = as.integer(args[3])


# input_df =
#   readRDS(file_in) %>% 
#   unite("database_for_cell_type", c(database, cell_type, .feature), remove = FALSE ) %>%
#   unite(".feature", c(.feature, cell_type), remove = FALSE ) %>%
# 
#   mutate(database_for_cell_type = factor(database_for_cell_type), .feature = factor(.feature)) %>%
#   mutate(count = as.integer(count), count_scaled = as.integer(count_scaled)) %>%
#   mutate(exposure_rate = -log(multiplier)) %>%
#   mutate(exposure_multiplier = exp(exposure_rate))

input_df =
  readRDS(file_in) %>%
  
  unite("database_for_cell_type", c(database, cell_type, .feature), remove = FALSE ) %>%
  mutate(
    database_for_cell_type = factor(database_for_cell_type),
    .feature = factor(.feature),
    .sample = factor(.sample)
  ) %>%
  
  unite(".feature_cell_type", c(.feature, cell_type), remove = FALSE ) %>%
  mutate(.feature_cell_type = factor(.feature_cell_type)) %>%
  
  mutate(exposure_rate = -log(multiplier)) %>%
  mutate(exposure_multiplier = exp(exposure_rate)) %>%
  
  # Exposure omit
  #mutate(exposure_rate = 0) %>%
  mutate(count = as.integer(count), count_scaled = as.integer(count_scaled))

output_list =
  input_df %>%
  
  fit_mixed_effect(
    mod_estimate  = readRDS("~/PostDoc/cellsig/dev/mod.rds"),
    mod_rng  = readRDS("~/PostDoc/cellsig/dev/gen_mod.rds"),
    feature = .feature_cell_type,
    counts = count,
    assoc_intercept = 4.2,
    assoc_slope = -0.55,
    assoc_sd_sd = 1.22,
    assoc_sd_shape = 1.14,
    lambda_mu = 9.8,
    lambda_sigma = 1.69,
    lambda_skew = -4.58,
    chains = 3,
    vb=FALSE,
    iterations = 350,
    sampling_iterations = 200
  )

# output_list =
#   fit_mixed_effect(
#     df = input_df,
#     mod_estimate  = readRDS("~/PostDoc/cellsig/dev/mod.rds"), #stan_model("~/PostDoc/cellsig/dev/mixed_effect.stan"), # readRDS("~/PostDoc/cellsig/dev/mod.rds"),
#     mod_rng  = readRDS("~/PostDoc/cellsig/dev/gen_mod.rds"),
#     assoc_intercept = 4.2, 
#     assoc_slope = -0.55, 
#     assoc_sd_sd = 1.22, 
#     assoc_sd_shape = 1.14, 
#     lambda_mu = 9.8, 
#     lambda_sigma = 1.69,
#     lambda_skew = -4.58,
#     vb=FALSE,
#     iterations = 1000,
#     sampling_iterations = 200,
#     chains = 3
#   )

output_df = 
  input_df %>% 
  #mutate(.feature_idx = as.integer(.feature_cell_type)) %>% 
  mutate(.feature_idx = as.integer(.feature)) %>% 
  distinct(.feature_idx, .feature, cell_type, level) %>% 
  left_join(
    output_list$rng %>% 
      rstan::summary("Y_gen", c(0.1, 0.5, 0.9)) %$%
      summary %>% 
      as_tibble() %>% 
      rowid_to_column(var = ".feature_idx")
  ) 


output_list %>% 
  c(list(output_df = output_df)) %>%
  saveRDS(file_out)


