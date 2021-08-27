# ~/third_party_sofware/cctools-7.2.0-x86_64-centos7/bin/makeflow_monitor makefile.makeflow.makeflowlog
# ~/third_party_sofware/cctools-7.2.0-x86_64-centos7/bin/makeflow -T slurm -J 200  dev/TCGA_makeflow_pipeline/makefile_ARMET_TCGA.makeflow

library(tidyverse)
library(magrittr)
library(cellsig)


# Read arguments
args = commandArgs(trailingOnly=TRUE)
file_in = args[1]
file_out = args[2]
cores = as.integer(args[3])

readRDS(file_in) %>%

  # TEMPORARY
  mutate(count = as.integer(count)) %>%
  
ref_intercept_only(
  exposure_rate,
  cores = cores,
  approximate_posterior = T
) %>%
  saveRDS(file_out)
