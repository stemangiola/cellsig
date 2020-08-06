test_that("model runs", {

  cellsig::counts %>%
    filter(level ==1 & cell_type == "epithelial") %>%
    
    # subset
    nest(data = -c(symbol, house_keeping)) %>%
    nest(hk = -house_keeping) %>%
    arrange(house_keeping ) %>%
    mutate(how_many = c(5, 50)) %>%
    mutate(hk = map2(hk, how_many, ~ .x %>% sample_n(.y))) %>%
    unnest(hk) %>%
    unnest(data) %>%
    
    ref_intercept_only(1, cores = 1, approximate_posterior = T)
})
