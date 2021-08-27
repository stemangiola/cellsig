
library(tidyverse)
library(magrittr)
library(cellsig)


# Read arguments
args = commandArgs(trailingOnly=TRUE)
file_in = args[1]
file_out = args[2]
cores = as.integer(args[3])

readRDS(file_in) %>%
ref_intercept_only(
  exposure_rate,
  cores = cores,
  approximate_posterior = T
) %>%
  saveRDS(file_out)
