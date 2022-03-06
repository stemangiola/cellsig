test_that("model runs", {

  readRDS("dev/modeling_results/level_1_cell_type_endothelial_partition_1_input.rds") %>% 
    
    nest(data = -.feature) %>%
    sample_n(50) %>%
    unnest(data) %>%
    
    mutate(multiplier = exp(exposure_rate)) %>% 
    cellsig_multilevel_varing_intercept(.sample, .feature, count, cell_type, multiplier, database)
})
