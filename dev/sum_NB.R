n = 1000
prop = c(0.6, 0.3, 0.1)
mu = c(20, 200, 1000)
sigma = exp(c(1.3420415 + log(mu) * -0.3386389))


mix_df = 
  tibble(
    prop = prop,
    mu = mu
  ) %>%
  mutate(
    population = 1:n(),
    sigma = sigma
  ) %>%
  
  
  mutate(samples = map(population, ~ 1:10000)) %>%
  unnest(samples) %>%
  mutate(draws = pmap(list(mu, sigma, prop), ~ rnbinom(n * ..3, mu = ..1, size = ..2))) %>%
  unnest(draws) %>%
  nest(draws = -samples) %>%
  mutate(mix = map_dbl(draws, ~ sum(.x)))
  
  
  mutate(draws = map2(
    mu, sigma, 
    ~ rnbinom(10000, mu = .x, size = .y) %>% 
      as_tibble() %>%
      rownames_to_column(var = "draw")
  )) %>%
  unnest(draws) %>%
  nest(data = -draw) %>%
  mutate(mix = map_dbl(data, ~ sum(.x$value * .x$prop))) %>%
  
  # Print parameters
  {
    fitdistrplus::fitdist(data = as.integer((.)$mix), distr = "nbinom", method = "mle") %>%
      print();
    (.)
  }


mix_df %>%
  ggplot(aes(mix)) + 
  geom_density() + 
  geom_density(aes(value, color=population, data = mix_df %>% unnest(data))) +
  scale_x_log10()


 gamlss::fitDist(as.integer(mix_df$mix), k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)


 
sum_NB = function(lambda_mat, sigma_mat, prop_mat){
   
  lambda_sum = prop_mat %*% lambda_mat; 
   sigma_sum =                        
   (lambda_sum^2) /
     (
       (prop_mat) %*%
         (
           (lambda_mat^2) /
             sigma_mat
         )
     ) ;
   
  c(lambda_sum, sigma_sum)
  }
  
sum_NB(matrix(mu), matrix(sigma), matrix(prop, nrow = 1))


neg_binomial_sum_moments = function( mus, phis) {
  mu_approx = sum(mus)
  
  list(
    mu_approx = mu_approx,
    phi_approx = (mu_approx^2) / sum( (mus^2) / phis )
  )
 
  

}

neg_binomial_sum_moments(mu, sigma)
           