# ~/third_party_sofware/cctools-7.2.0-x86_64-centos7/bin/makeflow_monitor makefile.makeflow.makeflowlog
# ~/third_party_sofware/cctools-7.2.0-x86_64-centos7/bin/makeflow -T slurm -J 200  dev/TCGA_makeflow_pipeline/makefile_ARMET_TCGA.makeflow


library(rstan)
library(tidyverse)
library(cellsig)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
file_in = args[1]
file_out = args[2]
cores = as.integer(args[3])

# file_in = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/cellsig/dev/modeling_results/level_2_cell_type_mono_derived_partition_18_input.rds"
readRDS(file_in) %>%
  
  cellsig_multilevel_varing_intercept(
    .sample,
    .feature,
    count, 
    cell_type,
    multiplier, 
    database, 
    cores = 15, 
    pass_fit = TRUE
  ) %>%
  saveRDS(file_out)




