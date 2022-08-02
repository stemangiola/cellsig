library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyseurat)
library(purrr)
library(stringr)

tuong<-readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/Tuong_et_al_2021/Tuong_et_al_2021_300921.RDS")
tuong.seurat<- as.Seurat(tuong)
assay<-tuong.seurat%>%mutate(dataset='EGAS00001005787')
#cannot directly use rename function here
#use immune cell type as cell_type; the data also contains non-immune cell types.
assay$cell_type<-assay@meta.data[["celltype.immune"]] 
assay@meta.data[["celltype.immune"]]<-NULL

saveRDS(assay,file='dev/xinpu_datascript/parsed_data/EGAS00001005787_final.rds')
#assay<-readRDS(file='dev/xinpu_datascript/parsed_data/EGAS00001005787_final.rds')

#list unique cell type
assay%>%distinct(dataset,sample,cell_type)
# cell number in each sample
assay%>%group_by(dataset,sample)%>%summarise(count=n())
