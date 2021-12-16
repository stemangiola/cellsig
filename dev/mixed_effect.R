library(rstan)
library(brms)
library(tidyverse)

counts = readRDS("~/PostDoc/cellsig/dev/counts.rds")

df = 
  counts %>% 
  unite("database_for_cell_type", c(database, level_4), remove = FALSE ) %>% 
  filter(.feature=="A2ML1") %>% 
  filter(level_4=="t_CD4_memory") %>%
  mutate(count_scaled = as.integer(count_scaled)) 

# counts %>% filter(.feature=="LRRC36") %>% filter(level_4=="t_CD4_memory") %>% ggplot(aes(level_4,count_scaled, fill=interaction(database, cell_type))) + geom_boxplot()
# 
# fit = df %>%  
#   brm( count_scaled ~ (1 | database), data = .,  family = negbinomial())

data = list(
  N = nrow(df),
  Y = df$count_scaled,
  num_dataset = length(unique(df$database)),
  num_coefficient = 1,
  grouping_idx = factor(df$database) %>% as.integer
)

#mod = stan_model("~/PostDoc/cellsig/dev/mixed_effect.stan")
mod <- stan_model(stanc_ret = stanc_builder('~/PostDoc/cellsig/dev/mixed_effect.stan'),auto_write = F,verbose = T)


my_fit = sampling(mod, data = data, cores = 4)


df = 
  counts %>% 
  unite("database_for_cell_type", c(database, level_4, .feature), remove = FALSE ) %>% 
  nest(data = -.feature) %>% 
  sample_n(100) %>% 
  unnest(data) %>% 
  filter(level_4=="t_CD4_memory") %>%
  mutate(count_scaled = as.integer(count_scaled)) 


data = list(
  N = nrow(df),
  Y = df$count_scaled,
  grouping_gene_idx_N  = df %>% mutate(database_for_cell_type = factor(database_for_cell_type), level_4 = factor(level_4)) %>% select(database_for_cell_type, level_4),

  G = length(unique(df$.feature)),
  D = length(unique(df$database_for_cell_type)),  
  grouping_gene_idx_D =  
    df %>% 
    mutate(database_for_cell_type = factor(database_for_cell_type), level_4 = factor(level_4)) %>% 
    select(database_for_cell_type, level_4) %>% 
    distinct()
)

mod = stan_model("dev/mixed_effect.stan")

my_fit = sampling(mod, data = data, cores = 4, iter = 600, warmup = 300)


# Other analysis

.libPaths("R/library_for_cell_sig_with_new_rstan/rstudio")

x = 
  counts %>% 
  with_groups(c(cell_type, .feature, database), ~ summarise(.x, count = median(count_scaled))) %>% 
  nest(data = -.feature) %>% 
  sample_n(1000) %>% 
  unnest(data) %>% 
  mutate(exposure_rate = 0) %>% 
  unite(".sample", c(database, cell_type), remove = FALSE) %>% 
  mutate(count = as.integer(count))

# Non centred
fit = x %>% 
  run_model_ref(exposure_rate,  shards = 1,  approximate_posterior = FALSE)

fit %>% saveRDS("dev/run_for_mean_overdispersion_association.rds")


# non aggregated
# mean     se_mean          sd       2.5%        25%        50%        75%      97.5%     n_eff     Rhat
# sigma_intercept  4.2116786 0.015305152 0.075132774  4.0675466  4.1587226  4.2135158  4.2662779  4.3447396  24.09814 1.134510
# sigma_slope     -0.5485913 0.001637523 0.009478366 -0.5655248 -0.5551634 -0.5490307 -0.5423468 -0.5292065  33.50367 1.109947
# sigma_sigma      1.1406230 0.003475751 0.023706466  1.0986771  1.1238373  1.1407187  1.1585306  1.1848475  46.51965 1.026858
# lambda_mu        9.8004146 0.004860596 0.044021864  9.7247181  9.7712344  9.7987880  9.8304879  9.8894398  82.02718 1.038594
# lambda_sigma     1.6889501 0.001995169 0.014935331  1.6614980  1.6791389  1.6879354  1.6974082  1.7214793  56.03644 1.074574
# lambda_skew     -4.5779882 0.017743853 0.184550666 -4.9373235 -4.6977757 -4.5689116 -4.4447598 -4.2416715 108.17711 1.022169
# 
# aggregated
# mean      se_mean          sd       2.5%        25%        50%        75%     97.5%     n_eff      Rhat
# sigma_intercept  4.6142189 0.0030368850 0.039431625  4.5458458  4.5866034  4.6129284  4.6439842  4.685032 168.59033 0.9953181
# sigma_slope     -0.6118821 0.0003802697 0.005012327 -0.6222671 -0.6149592 -0.6117695 -0.6085909 -0.603095 173.73824 0.9924648
# sigma_sigma      1.2266322 0.0017470746 0.013392337  1.1994317  1.2179779  1.2263259  1.2364188  1.252214  58.76109 1.0532557
# lambda_mu       10.1356283 0.0017001326 0.025571430 10.0838454 10.1180580 10.1362845 10.1540440 10.182078 226.22700 1.0006050
# lambda_sigma     2.0733526 0.0007319888 0.008497813  2.0552918  2.0678701  2.0733129  2.0797962  2.087900 134.77380 1.0009592
# lambda_skew     -8.2240221 0.0151203963 0.218904867 -8.6750264 -8.3628650 -8.2267583 -8.0828402 -7.813414 209.59672 0.9948057

df = 
  counts %>% 
  unite("database_for_cell_type", c(database, level_4, .feature), remove = FALSE ) %>% 
  nest(data = -.feature) %>% 
  sample_n(100) %>% 
  unnest(data) %>% 
  filter(level_4=="t_CD4_memory") %>%
  mutate(count_scaled = as.integer(count_scaled)) 

# mod = stan_model("dev/mixed_effect.stan")

fit_mix =
  fit_mixed_effect(
    df,
    mod,
    assoc_intercept = 4.2, 
    assoc_slope = -0.55, 
    assoc_sd_sd = 1.22, 
    assoc_sd_shape = 1.14, 
    lambda_mu = 9.8, 
    lambda_sigma = 1.69,
    lambda_skew = -4.58
  )

fit_mix %>% spread_draws(gene_mean[G], gene_sd[G]) %>% mutate(mu = rnorm(n(), gene_mean, gene_sd)) %>% arrange(desc(mu))

plot(summary(fit_mix, "gene_mean")$summary[,1], summary(fit_mix, "gene_sd")$summary[,1])

gen_mod = stan_model("~/PostDoc/cellsig/dev/mixed_effect_generate.stan")


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
    distinct()
)

rng =  rstan::gqs(
  gen_mod,
  draws =  as.matrix(fit_mix),
  data = data
)

summary(rng, "Y_gen", c(0.1, 0.9))$summary %>% 
  as_tibble() %>% 
  rowid_to_column(var = ".feature_idx") %>% 
  left_join(
    df %>% 
      mutate(database_for_cell_type = factor(database_for_cell_type), .feature = factor(.feature)) %>% 
      select(database_for_cell_type, cell_type, .feature, count_scaled) %>% 
      mutate(.feature_idx = as.integer(.feature))
  ) %>% 
  
  ggplot(aes(.feature, count_scaled + 1, color = database_for_cell_type)) + 
  geom_violin()  +
  geom_errorbar(aes(ymin = `10%`+1, ymax=`90%`+1), color="black",  alpha = 0.4, size = 0.5 ) +
  scale_y_log10() +
  guides(color="none") +
  theme_bw()

rng %>% tidybayes::gather_draws(Y_gen[G])


df %>% 
  with_groups(c(cell_type, .feature, database), ~ summarise(.x, count_scaled = median(count_scaled))) %>% 
  with_groups(c(cell_type, .feature), ~ summarise(.x, count_scaled = median(count_scaled))) %>% 
  pull(count_scaled) %>% 
  log1p %>% 
  hist

df %>% 
  filter(.feature=="AP3B1") %>% 
  filter(level_4=="t_CD4_memory") %>% 
  ggplot(aes(level_4,count_scaled, fill=interaction(database, cell_type))) + 
  geom_boxplot()



# Parallel
.libPaths("~/R/library_for_cell_sig_with_new_rstan/rstudio")

# install.packages("RcppParallel", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# Sys.setenv(TBB_LIB = RcppParallel::tbbLibraryPath())
# Sys.setenv(TBB_INC  = RcppParallel::CxxFlags())


unloadNamespace( asNamespace("rstan"))
unloadNamespace( asNamespace("StanHeaders"))
unloadNamespace( asNamespace("RcppParallel"))

library(rstan)
rstan_options(threads_per_chain = 4)

d <- read.csv("~/PostDoc/cellsig/dev/RedcardData.csv", stringsAsFactors = FALSE)
d2 <- d[!is.na(d$rater1),]
redcard_data <- list(n_redcards = d2$redCards, n_games = d2$games, rating = d2$rater1)
redcard_data$N <- nrow(d2)
redcard_data$grainsize <- 1

gmod <- stan_model(stanc_ret = stanc_builder('~/PostDoc/cellsig/dev/logistic1.stan'),auto_write = T,verbose = T)

# Sys.setenv(STAN_NUM_THREADS=3)
f = sampling(gmod, redcard_data,
       chains = 4,
       cores = 4,
       refresh = 1000)
