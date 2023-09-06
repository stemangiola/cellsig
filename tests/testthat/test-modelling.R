test_that("model runs", {


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
      pass_fit = TRUE, 
      cores = 1
    )
})
