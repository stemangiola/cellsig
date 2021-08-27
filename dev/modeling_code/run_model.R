library(unixtools)
dir.create(t <- paste(sprintf("~/.Rtemp/%s", basename(tempdir())), Sys.getpid(), sep='-'), FALSE, TRUE, "0700")
set.tempdir(t)

library(tidyverse)
library(magrittr)
# library(cellsig)
# library(future)
# library("future.batchtools")
# library(furrr)
library(tidybulk)
library(tidySummarizedExperiment)

local_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/cellsig"
cores = 8

# slurm <- future::tweak(batchtools_slurm,
#                        template = sprintf("%s/dev/modeling_results/slurm_batchtools.tmpl", local_dir),
#                        resources=list(
#                          cores = cores,
#                          memory_mb = 5000,
#                          time = "48:00:00"
#                        )
# )
# plan(slurm)



# Save files
create_partition_files = function(.data, .level, .partitions = 30){
  .data %>%
    
    nest(data = -c(cell_type, .feature)) %>%
    
    # For testing
    #sample_frac(0.01) %>%
    mutate(partition = sample(1:.partitions, size = n(), replace = T)) %>%
    unnest(data) %>%
    nest(data = -partition) %>%
    mutate(saved = map2_lgl(
      data, partition,
      ~ {
        .x %>% saveRDS(sprintf("%s/dev/modeling_results/level_%s_patition_%s.rds", local_dir, .level, .y))
        TRUE
      }
    ))
}

create_partitions = function(.data, .level, .partitions = 30){
  .data %>%
    
    nest(data = -c(cell_type, .feature)) %>%
    mutate(partition = sample(1:.partitions, size = n(), replace = T)) %>%
    unnest(data) %>%
    nest(data = -partition)
}

# Read counts
counts =
  readRDS("dev/counts.rds") %>%
  select(-cell_type) %>%
  pivot_longer(
    contains("level_"), names_prefix="level_", 
    names_to = "level", values_to="cell_type",
    names_transform=list(level=as.integer)
  ) %>%
  filter(cell_type %>% is.na %>% `!`)



#sys("rm modeling_results/*rds")

# Create partition files
counts %>%
  nest(data = -level) %>%
  mutate(partitions = map2(
    data, level,
    ~ create_partition_files(.x, .y, 15)
  ))

# tibble(level=1:5) %>%
#   mutate(partitions = map(
#     level,
#     ~ counts %>% 
#       filter(level==.x) %>% 
#       create_partition_files(.x, 15)
#   ))

# Create input
sprintf("CATEGORY=create_input\nMEMORY=20024\nCORES=%s\nWALL_TIME=14000", cores) %>%
  
  c(
  dir(sprintf("%s/dev/modeling_results/", local_dir), pattern = ".rds") %>%
    grep("result", ., invert = T, value = T) %>%
    enframe(value = "file") %>%
    mutate(cores = !!cores) %>%
    mutate(command = map2_chr(
      file, cores,
      ~sprintf(
          "dev/modeling_results/%s: dev/modeling_results/%s\n\tRscript dev/modeling_results/core_run_model.R dev/modeling_results/%s dev/modeling_results/%s %s",
          sprintf("%s_result.rds", basename(.x) %>%  sub("^([^.]*).*", "\\1", .)),
          .x,
          .x,  
          sprintf("%s_result.rds", basename(.x) %>%  sub("^([^.]*).*", "\\1", .)) ,
          .y
      ))
    ) %>%
    pull(command) %>%
    unlist()
  ) %>%
  write_lines("dev/modeling_code/run_model.makeflow") 
#%>%
#  saveRDS(sprintf("%s/dev/temp.rds", local_dir))


# dir(sprintf("%s/dev/modeling_results/", local_dir), pattern = ".rds", full.names = T) %>%
#   grep("result", ., value = T) %>%
#   map_dfr(~ readRDS(.x)) %>%
#   saveRDS("dev/cellsig_theoretical_transcript_abundance_distribution.rds", compress = "xz")
