library(rstan)

df = 
  counts %>% 
  filter(.feature=="A2ML1") %>% 
  filter(level_4=="t_CD4_memory") %>%
  mutate(count_scaled = as.integer(count_scaled))

fit = df %>%  
  brm( count_scaled ~ (1 | database), data = .,  family = negbinomial())

data = list(
  N = nrow(df),
  Y = df$count_scaled,
  num_dataset = length(unique(df$database)),
  num_coefficient = 1,
  grouping_idx = factor(df$database) %>% as.integer
)

mod = stan_model("dev/mixed_effect.stan")

my_fit = sampling(mod, data = data, cores = 4)
