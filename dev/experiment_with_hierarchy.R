library(rstan)
library(tidyverse)
library(furrr)
plan(multiprocess, workers=10)

mod = "
data{
int N;
int K;
vector[N] y;
int idx[N];
}
parameters{
  real mu_mu;
  real<lower=0> mu_sigma;
  
  real<lower=0> sigma_mu;
  real<lower=0> sigma_sigma;
  
  vector<offset=mu_mu, multiplier=mu_sigma>[K] mu;
  vector<lower=0>[K] sigma;
  

}
model{
y ~ normal(mu[idx], sigma[idx]);

mu ~ normal(mu_mu, mu_sigma);
sigma~ normal(sigma_mu, sigma_sigma);

sigma_sigma ~ student_t(3, 0, 5);
mu_sigma ~ student_t(3, 0, 5);
}
"
my_mod = stan_model(model_code = mod)


.data = 
  tibble(idx = c(1,2, 3,4), .mean = c(1, 2, 1.5, 1.8), .sigma = c(1, 1, 1,1)) %>%
  mutate(y = map2(.mean, .sigma, ~ rnorm(100, .x, .y))) %>%
  unnest(y) 

.data = 
  tibble(idx = c(1,2), .mean = c(1, 2), .sigma = c(1, 1)) %>%
  mutate(y = map2(.mean, .sigma, ~ rnorm(100, .x, .y))) %>%
  unnest(y) 


fit = sampling(my_mod, data = list(N=nrow(.data), y = .data$y, idx = .data$idx, K=(.data$idx %>% unique %>% length)))


lv4 = 
  counts_fit %>%
  filter(level==4) %>%
  unnest(lv_data) %>%
  nanny::subset(c(cell_type, symbol))

fit_df = 
  lv4 %>%
  left_join( cellsig:::get_ancestor_child(cellsig::tree), by = "cell_type" ) %>%
  group_by(ancestor, symbol) %>%
  summarise(
    sigma_mu = sd(lambda_log ),
    sigma_sigma = sd(sigma_inv_log)
  )  %>%
  ungroup() %>%
  pivot_longer(names_to = ".which", values_to = ".sd", cols = c(sigma_mu, sigma_sigma)) %>%
  nest(data = - c(ancestor, .which)) %>%
  mutate(fit = future_map(data, ~ gamlss::fitDist(.x$.sd, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)))

fit_df %>%
  filter(.which == "sigma_mu") %>%
  pull(fit)

fit_df %>%
  filter(.which == "sigma_sigma") %>%
  pull(fit)

fit_df %>%
  unnest(data) %>%
  ggplot(aes(.sd, color = ancestor)) + geom_density() + facet_wrap(~.which)

fit_df %>%
  dplyr::select(-fit) %>%
  unnest(data) %>%
  pivot_wider(names_from = .which, values_from = c(.sd)) %>%
  filter(ancestor=="t_CD8") %>%
  tidygate::gate(.element = symbol, .dim1 = sigma_mu , .dim2 = sigma_sigma)


fit = 
  lv4 %>%
  filter(grepl("nk", cell_type)) %>%
  group_by(symbol) %>%
  summarise(
    sigma_mu = sd(lambda_log ),
    sigma_sigma = sd(sigma_inv_log)
  ) %>%
  pull(sigma_mu) %>%
  

summary(fit)
