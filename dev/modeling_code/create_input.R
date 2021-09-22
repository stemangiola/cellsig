library(tidyverse)
library(magrittr)
library(tidybulk)
library(tidySummarizedExperiment)

# library(cellsig)
# library(future)
# library("future.batchtools")
# library(furrr)

# slurm <- future::tweak(batchtools_slurm,
#                        template = sprintf("%s/dev/modeling_results/slurm_batchtools.tmpl", local_dir),
#                        resources=list(
#                          cores = cores,
#                          memory_mb = 5000,
#                          time = "48:00:00"
#                        )
# )
# plan(slurm)

# create_partitions = function(.data, .level, .partitions = 30){
#   .data %>%
#     
#     nest(data = -c(cell_type, .feature)) %>%
#     mutate(partition = sample(1:.partitions, size = n(), replace = T)) %>%
#     unnest(data) %>%
#     nest(data = -partition)
# }

# Save files
create_partition_files = function(.data, .level, .partitions = 30){
  .data %>%
    mutate(level = .level) %>%
    nest(data = -c(cell_type, .feature)) %>%
    
    # For testing
    #sample_frac(0.01) %>%
    mutate(partition = sample(1:.partitions, size = n(), replace = T)) %>%
    unnest(data) %>%
    nest(data = -partition) %>%
    mutate(saved = map2_lgl(
      data, partition,
      ~ {
        .x %>% saveRDS(sprintf("%s/dev/modeling_results/level_%s_patition_%s_input.rds", local_dir, .level, .y))
        TRUE
      }
    ))
}

# Read counts
readRDS("dev/counts.rds") %>%
  select(-cell_type) %>%
  pivot_longer(
    contains("level_"), names_prefix="level_", 
    names_to = "level", values_to="cell_type",
    names_transform=list(level=as.integer)
  ) %>%
  filter(cell_type %>% is.na %>% `!`) %>%
  mutate(count = as.integer(count)) %>%
  
  # Create partition files
  nest(data = -level) %>%
  mutate(partitions = map2(
    data, level,
    ~ create_partition_files(.x, .y, 15)
  ))
