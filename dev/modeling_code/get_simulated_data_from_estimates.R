library(tidyverse)
library(glue)
library(cmdstanr)
library(posterior)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
directory_path = args[1]

# batch = 3
# directory_path = glue("dev/benchmark_database_crossvalidation_batch_{batch}/")
# output_path = glue("dev/benchmark_database_crossvalidation_batch_{batch}/simulated_data_crossvalidation_batch_{batch}.rds")

output_path = glue("{directory_path}/simulated_data.rds")

sample_cell_type_count = 
  directory_path |> 
  dir( pattern = "input.rds", full.names = TRUE) |> 
  map_df(~ readRDS(.x) |> distinct(.sample, cell_type)) |> 
  distinct(.sample, cell_type) |> 
  count(cell_type)

directory_path |> 
  dir( pattern = "result.rds", full.names = TRUE) |> 
  map_df(~{
    result = readRDS(.x)
    
    samples_to_sample =
      sample_cell_type_count |> 
      filter(cell_type == (pull(result, cell_type) |> unique() |> as.character())) |> 
      pull(n)
    
    attr(result, "indeces") |> 
      distinct(.feature, cell_type, feature_cell_type_idx) |> 
      left_join(sample_cell_type_count, by = "cell_type") |> 
      left_join(
        
        attr(result, "rng")$draws("Y_gen") |> 
          as_draws_df() |> 
          sample_n(samples_to_sample, replace=TRUE) |> 
          rowid_to_column(var = ".sample") |> 
          tidyr::pivot_longer(-.sample) |> 
          mutate(value = as.integer(value)) |> 
          tidyr::extract(name, "feature_cell_type_idx", "Y_gen\\[([0-9]+)\\]", convert = TRUE),
        by = "feature_cell_type_idx"
        
      )
    
  }) |> 
  select(-feature_cell_type_idx) |> 
  mutate_if(is.character, as.factor) |> 
  saveRDS(output_path, compress = "xz")



