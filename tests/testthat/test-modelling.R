test_that("model runs", {

  library(cellsig)

  # debugonce(cellsig_multilevel_varing_intercept)
  
 x =
   test_random_intercept |>
    
    #mutate(multiplier = exp(exposure_rate)) |> 
    cellsig_multilevel_varing_intercept(
      .sample, 
      .feature,
      count, 
      cell_type,
      multiplier, 
      database, 
      pass_fit = TRUE
<<<<<<< HEAD
      #, 
      #use_cmdstanr = TRUE
=======
>>>>>>> b9c107af17b9250647f20b346108dc0faf7f4272
    )
})
