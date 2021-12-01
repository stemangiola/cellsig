library(rstan)

df = 
  counts %>% 
  unite("database_for_cell_type", c(database, level_4), remove = FALSE ) %>% 
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



data = list(
  N = nrow(df),
  Y = df$count_scaled,
  grouping_gene_idx_N  = df %>% mutate(database_for_cell_type = factor(database_for_cell_type), level_4 = factor(level_4)) %>% select(database_for_cell_type, level_4),

  G = length(unique(df$level_4)),
  D = length(unique(df$database_for_cell_type)),  
  grouping_gene_idx_D =  
    df %>% 
    mutate(database_for_cell_type = factor(database_for_cell_type), level_4 = factor(level_4)) %>% 
    select(database_for_cell_type, level_4) %>% 
    distinct()
)

mod = stan_model("dev/mixed_effect.stan")

my_fit = sampling(mod, data = data, cores = 4)
