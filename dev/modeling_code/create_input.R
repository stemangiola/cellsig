library(tidyverse)
library(magrittr)
library(tidybulk)
library(tidySummarizedExperiment)
library(glue)

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

local_dir = "~/PostDoc/cellsig"

# Save files
create_partition_files = function(.data, .level, .partitions = 30){
  .data %>%
    nest(data = -c(.feature)) %>%
    
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
  nest(data = -c(level, cell_type, .feature)) %>%
  nest(data = -c(level, cell_type)) %>% 
  
  mutate(data = map(
    data, 
    ~ mutate(.x, partition = sample(1:100, size = n(), replace = T))
  )) %>%
  unnest(data) %>% 
  unnest(data) %>% 
  
  # Save
  nest(data = -c(level, cell_type, partition)) %>%
  mutate(saved = pmap_lgl(
    list(data, level, cell_type, partition),
    ~ {
      ..1 %>% 
        mutate(
          level = ..2,
          cell_type = ..3,
          partition= ..4
        ) %>% 
      saveRDS(glue("{local_dir}/dev/modeling_results/level_{..2}_cell_type_{..3}_partition_{..4}_input.rds") )
      TRUE
    }
  ))
