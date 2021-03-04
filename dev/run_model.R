library(tidyverse)
library(magrittr)
library(cellsig)
library(future)
library("future.batchtools")
library(furrr)
library(tidyr)

local_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/cellsig"

slurm <- future::tweak(batchtools_slurm,
                       template = sprintf("%s/dev/modeling_files/slurm_batchtools.tmpl", local_dir),
                       resources=list(
                         cores = 1,
                         memory_mb = 20000
                       )
)

plan(slurm)

# counts =
#   readRDS(file="~/PhD/deconvolution/ARMET/dev/counts_infer_NB.rds") %>%
#   #left_join( (.) %>% distinct(`symbol original`) %>% mutate(run = sample(1:5, n(), replace = T))) %>%
#   #filter(`Cell type category` == "house_keeping" | run == my_run) %>%
#   mutate(symbol = `symbol original`) %>%
#
#   # Replace the category house_keeping
#   mutate(`house keeping` = `Cell type category` == "house_keeping") %>%
#   rename(temp = `Cell type category`) %>%
#   left_join(
#     (.) %>%
#       filter(temp != "house_keeping") %>%
#       distinct(sample, temp, level) %>%
#       rename(`Cell type category` = temp)
#   ) %>%
#   select(-temp) %>%
#   filter(`Cell type category` %>% is.null %>% `!`) %>%
#
#   # Adapt it as input
#   select(sample, symbol, count, `Cell type category`, level, `count scaled`, `house keeping`)


# #
load("dev/counts.rda")

# Save files
create_partition_files = function(.data, .level, .partitions = 30){
  .data %>%

    # IMPORTANT! To change for every level
    mutate(level=.level, cell_type = !!as.symbol(sprintf("level_%s", .level))) %>%

    # Process
    filter(cell_type %>% is.na %>% `!`) %>%
    nest(data = -c(cell_type, symbol)) %>%
    #sample_frac(0.01) %>%
    mutate(partition = sample(1:.partitions, size = n(), replace = T)) %>%
    unnest(data) %>%
    nest(data = -partition) %>%
    mutate(saved = map2_lgl(
      data, partition,
      ~ {
        .x %>% saveRDS(sprintf("dev/modeling_files/level_%s_patition_%s.rds", .level, .y))
        TRUE
      }
    ))
}


create_partitions = function(.data, .level, .partitions = 30){
  .data %>%
    
    # IMPORTANT! To change for every level
    mutate(cell_type = !!as.symbol(sprintf("level_%s", .level))) %>%
    
    # Process
    filter(cell_type %>% is.na %>% `!`) %>%
    nest(data = -c(cell_type, symbol)) %>%
    mutate(partition = sample(1:.partitions, size = n(), replace = T)) %>%
    unnest(data) %>%
    nest(data = -partition)
}


tibble(level=1:5) %>%
  mutate(partitions = map(level, ~ create_partitions(counts, .x, 3000))) %>%
  unnest(partitions) %>%

  # FOR TESTING
  slice(1:5) %>%
  mutate(inference = future_imap(
    data,
    # ~ 5
    ~ .x %>%
      ref_intercept_only(
        exposure_rate,
        cores = 1,
        approximate_posterior = T
    ),
    .options = furrr_options(packages = c("tidyverse", "magrittr", "cellsig"))
  )) %>%
  saveRDS(sprintf("%s/dev/temp.rds", local_dir))

# tibble(level=1) %>%
#   mutate(partitions = map(level, ~ create_partitions(counts, .x, 3000))) %>%
#   unnest(partitions) %>%
#   
#   # FOR TESTING
#   slice(1) %>%
#   pull(data) %>%
#   .[[1]] %>%
#   ref_intercept_only(
#     exposure_rate,
#     cores = 1,
#     approximate_posterior = TRUE
#   ) 

