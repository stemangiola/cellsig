library(Seurat)
library(tidyverse)
library(tidyseurat)
library(SeuratDisk)

Convert("./dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/henry18_0.processed.h5ad", ".h5seurat")
seuratObject <- LoadH5Seurat("./dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/henry18_0.processed.h5Seurat")

henry = Seurat::Read10X_h5("./dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/henry18_0.processed.h5ad")

Tuong_et_al_2021_300921 <- readRDS("~/master_project/cellsig/dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/Tuong_et_al_2021/Tuong_et_al_2021_300921.RDS")
