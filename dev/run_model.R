library(unixtools)
dir.create(t <- paste(sprintf("~/.Rtemp/%s", basename(tempdir())), Sys.getpid(), sep='-'), FALSE, TRUE, "0700")
set.tempdir(t)

library(tidyverse)
library(magrittr)
library(cellsig)
library(future)
library("future.batchtools")
library(furrr)

cores = 8

local_dir = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/cellsig"

slurm <- future::tweak(batchtools_slurm,
                       template = sprintf("%s/dev/modeling_files/slurm_batchtools.tmpl", local_dir),
                       resources=list(
                         cores = cores,
                         memory_mb = 5000,
                         time = "48:00:00"
                       )
)

plan(slurm)



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
        .x %>% saveRDS(sprintf("%s/dev/modeling_files/level_%s_patition_%s.rds", local_dir, .level, .y))
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

# Create files
# load("dev/counts.rda")
# sys("rm modeling_files/*rds")
# tibble(level=1:5) %>%
#   mutate(partitions = map(level, ~ create_partition_files(counts, .x, 30))) 

dir(sprintf("%s/dev/modeling_files/", local_dir), pattern = ".rds", full.names = T) %>%
  enframe(value = "file") %>%
  mutate(cores = !!cores) %>%
  mutate(inference = future_map2(
    file, cores,
    ~ 11,
    # ~ {
    #   
    #   # library(unixtools)
    #   # dir.create(t <- paste(sprintf("~/.Rtemp/%s", basename(tempdir())), Sys.getpid(), sep='-'), FALSE, TRUE, "0700")
    #   # set.tempdir(t)
    #   
    #   library(tidyverse)
    #   library(magrittr)
    #   library(cellsig)
    #   
    #   .x %>%
    #   readRDS() %>%
    #   ref_intercept_only(
    #     exposure_rate,
    #     cores = .y,
    #     approximate_posterior = T
    #   )
    #   }
    #,
    #.options = furrr_options(packages = c("tidyverse", "magrittr", "cellsig"))
  )) 
#%>%
#  saveRDS(sprintf("%s/dev/temp.rds", local_dir))


# On the fly - PROBABLY TOO MUCH
# counts_partitions =
#   tibble(level=1:5) %>%
#   mutate(partitions = map(level, ~ create_partition_files(counts, .x, 30))) %>%
#   unnest(partitions)
# 
# counts_partitions %>%
# 
#   mutate(inference = future_map(
#     data,
#     ~ .x %>%
#       ref_intercept_only(
#         exposure_rate,
#         cores = 1,
#         approximate_posterior = T
#     ),
#     .options = furrr_options(packages = c("tidyverse", "magrittr", "cellsig"))
#   )) %>%
#   saveRDS(sprintf("%s/dev/temp.rds", local_dir))


