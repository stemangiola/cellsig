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
options("tidybulk_do_validate"= FALSE) 



#sys("rm modeling_files/*rds")
tibble(level=1:5) %>%
  mutate(partitions = map(
    level,
    ~ counts %>% 
      filter(level==.x) %>% 
      create_partition_files(.x, 15)
  ))

# Create input
sprintf("CATEGORY=create_input\nMEMORY=20024\nCORES=%s\nWALL_TIME=14000", 12) %>%
  
  c(
  dir(sprintf("%s/dev/modeling_files/", local_dir), pattern = ".rds") %>%
    grep("result", ., invert = T, value = T) %>%
    enframe(value = "file") %>%
    mutate(cores = !!cores) %>%
    mutate(command = map2_chr(
      file, cores,
      ~sprintf(
          "dev/modeling_files/%s: dev/modeling_files/%s\n\tRscript dev/modeling_files/core_run_model.R dev/modeling_files/%s dev/modeling_files/%s %s",
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
  write_lines("dev/modeling_files/run_model.makeflow") 
#%>%
#  saveRDS(sprintf("%s/dev/temp.rds", local_dir))


dir(sprintf("%s/dev/modeling_files/", local_dir), pattern = ".rds", full.names = T) %>%
  grep("result", ., value = T) %>%
  map_dfr(~ readRDS(.x)) %>%
  saveRDS("dev/cellsig_theoretical_transcript_abundance_distribution.rds", compress = "xz")
