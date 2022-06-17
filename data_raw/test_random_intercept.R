test_random_intercept = 
  readRDS("dev/modeling_results_TEMP/level_1_cell_type_endothelial_partition_1_input.rds") |> 
  
  nest(data = -.feature) |>
  sample_n(50) |>
  unnest(data) %>% 
  select(.sample, .feature, count, cell_type, multiplier, database)

save(test_random_intercept, file="data/test_random_intercept.rda", compress = "xz")
