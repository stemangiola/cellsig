library(job)


  mixed_effect_cmdstanr =  read_file("inst/stan/mixed_effect.stan")
  
  save(mixed_effect_cmdstanr, file="data/mixed_effect_cmdstanr.rda", compress = "xz")

  mixed_effect_generate_cmdstanr =  read_file("inst/stan/mixed_effect_generate.stan")
  
  save(mixed_effect_generate_cmdstanr, file="data/mixed_effect_generate_cmdstanr.rda", compress = "xz")
