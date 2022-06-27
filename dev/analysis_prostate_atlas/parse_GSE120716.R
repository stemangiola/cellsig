library(Seurat)
library(tidyverse)
library(tidyseurat)
library(SeuratDisk)
library(scater)

library(zellkonverter)
GSE120716 = readH5AD(  file = "PPCG_deconvolution_signatures_RAW_DATA/GSE120716/henry18_0.processed.h5ad")

names(assays(GSE120716)) = "counts"
GSE120716 %>% 
  logNormCounts() %>% 
  as.Seurat() %>% 
  saveRDS("PPCG_deconvolution_signatures_RAW_DATA/GSE120716/GSE120716_parsed.rds")



