
# srun --job-name "InteractiveJob" --cpus-per-task 8 --mem-per-cpu 55000 --time 48:00:00 --pty bash
# Rscript dev/modeling_code/create_input.R dev/benchmark_database_crossvalidation_batch_2/training_data_2_parsed.rds dev/benchmark_database_crossvalidation_batch_2
# ~/third_party_sofware/cctools-7.2.0-x86_64-centos7/bin/makeflow -T slurm -J 200  dev/benchmark_database_crossvalidation/run_model.makeflow

library(tidyverse)
library(magrittr)
library(tidybulk)
library(tidySummarizedExperiment)
library(glue)
library(cellsig)
library(yaml)
library(data.tree)

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

# Read arguments
args = commandArgs(trailingOnly=TRUE)
dataset_in_path = args[1]
directory_out = args[2]
tree_in_path = args[3]

local_dir = "~/PostDoc/cellsig"

dir.create(file.path(local_dir, directory_out), showWarnings = FALSE)

# Load data
dataset_in = readRDS(dataset_in_path) 

# If no tree provided, just une one level
tree_in = 
  tree_in_path %>% 
  when(
    is.na(.) ~ from_dataframe_to_one_level_tree(dataset_in, cell_type),
    ~ read_yaml(.) %>% as.Node
  )

# PARSE
dataset_in %>%
  
  # Imputation
  as_SummarizedExperiment(sample, symbol, count) %>% 
  impute_missing_abundance(~ cell_type, force_scaling = TRUE) %>% 
  as_tibble() %>% 
  filter(count %>% is.na %>% `!`) %>% 
  
  # Parsing
  tree_and_signatures_to_database(
    tree_in,
    .,
    .sample,
    cell_type,
    .feature,
    count
  ) %>%
  identify_abundant(.sample, .feature, count) %>%
  scale_abundance(.sample, .feature, count) %>%
  dplyr::select(-count_scaled) %>%
  filter(!.imputed) %>% 
  select(-.imputed) %>% 

  # CREATE INPUTS
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
  
  mutate(number_of_partitions = if_else(cell_type=="immune_cell", 100, 40)) %>%
  mutate(data = map2(
    data, number_of_partitions,
    ~ mutate(.x, partition = sample(1:.y, size = n(), replace = T))
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
      droplevels() %>% 
      saveRDS(glue("{local_dir}/{directory_out}/level_{..2}_cell_type_{..3}_partition_{..4}_input.rds"), compress=FALSE )
      TRUE
    }
  ))



# CREATE MAKEFILE

cores = 15
tab = "\t"

sprintf("CATEGORY=create_input\nMEMORY=20024\nCORES=%s", cores) %>%
  
  c(
    dir(glue("{local_dir}/{directory_out}/"), pattern = "input.rds") %>%
      
      enframe(value = "file") %>%
      mutate(cores = !!cores) %>%
      mutate(command = map2_chr(
        file, cores,
        ~{
          my_basename = basename(.x) %>%  sub("^([^.]*).*", "\\1", .)
          my_basename = glue("{my_basename}_result.rds")
          glue("{directory_out}/{my_basename}: {directory_out}/{.x}\n{tab}Rscript dev/modeling_code/run_model.R {directory_out}/{.x} {directory_out}/{my_basename} {.y}" )
          }
        )
      ) %>%
      pull(command) %>%
      unlist()
  ) %>%
  write_lines(glue("{local_dir}/{directory_out}/run_model.makeflow")) 

