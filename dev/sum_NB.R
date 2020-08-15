library(tidyverse)

# Number of potential cells This number is relative and does not really matter
n = 1000

# Proportion of the 3 cell type in evey sample/person
prop = c(0.9, 0.05, 0.05)

# Parameters
mu = c(20, 200, 1000)
sigma = c(5, 25, 150)

# Reparametrisation
rgamma_rep = function(n, mu, phi){
  
  phi = phi * mu
  
  alpha <- mu * mu / phi; 
  beta <- mu / phi;
  
  rgamma(n, shape = alpha, rate = beta)
}


mix_df = 
  
  # Set up dataframe
  tibble(  prop = prop,  mu = mu  ) %>%
  mutate(  population = 1:n(),  sigma = sigma ) %>%
  mutate(n = 1000) %>%
  
  # Create 1000 samples/people
  mutate(sample = map(population, ~ 1:1000)) %>%
  unnest(sample) %>%
  
  # Draw the abundance expected value for each sample/person for each cell type 
  mutate(gamma_draw = rgamma_rep(n(), mu, sigma)) %>%
  
  # Generate trancripts for 1000 cells from each cell type as reference, and sum because it is what we observe through sequencing
  mutate(sum_cell_pure = map2_int(gamma_draw, n, ~ rpois(.y, .x) %>% sum)) %>%
  
  # Generate transcripts for 1000 total cells, divided among 3 cell types according to proportions, and sum the total gene abundance
  mutate(sum_cell_in_mix = map2_int(gamma_draw, n*prop, ~ rpois(.y, .x) %>% sum)) %>%
  nest(sample_df = -sample) %>%
  
  # Sum all outputs for each sample/person
  mutate(mix_gamma_draw = map_dbl(sample_df, ~ sum(.x$gamma_draw))) %>%
  mutate(mix = map_dbl(sample_df, ~ sum(.x$sum_cell_in_mix))) %>%

  # Print estimated parameters of the mix
  {
    fitdistrplus::fitdist(data = as.integer((.)$mix), distr = "nbinom", method = "mle") %>%
      print();
    (.)
  } %>%
  
  # Plot
  {
    ((.) %>% unnest( sample_df) %>%
       ggplot(aes(mix)) + 
       geom_density() + 
       geom_density(aes(sum_cell_pure, color=as.factor(population))) +
       scale_x_log10() +
       ggtitle("the pure components and resulting mix in black")) %>%
      plot()
    
    ((.) %>% unnest( sample_df) %>%
      ggplot(aes(mix)) + 
      geom_histogram() + 
       geom_histogram(aes(sum_cell_in_mix, color=as.factor(population)), alpha=0.2) +
      #geom_density() +
      scale_x_log10() +
       ggtitle("sum of the actual mixes")) %>%
      plot()
    (.)
    
  }





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
           

LL = function(a1, a2, a3, b1, b2, b3){
  -sum(
    log(
    .Call("_coga_dcoga_approx", PACKAGE = "coga", x, c(a1, a2, a3), c(b1, b2, b3)/100)
    )
  )
}

x <- mix_df$mix_gamma_draw

mle(
  minuslogl = LL, 
  start = list(a1 = 4, a2 = 8, a3 =6, b1 = 20, b2 =4, b3 = 0.6),
  lower = c( 0, 0, 0, 0, 0, 0),
  upper = c(Inf, Inf, Inf, Inf, Inf, Inf),
  method = "L-BFGS-B"
)
