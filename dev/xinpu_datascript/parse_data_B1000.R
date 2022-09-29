library(Seurat)
library(dplyr)
library(tidyseurat)
library(purrr)

# Load the raw data
cellsig_B1000 <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/cellsigDB/cellsig_B1000.rds")
