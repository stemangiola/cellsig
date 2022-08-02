library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyseurat)
library(purrr)
library(readxl)

# load meta data which contains cell type
meta_data<- read_excel("dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/GSE176031/41467_2021_27322_MOESM2_ESM.xlsx",
                       sheet=4)
# load data: (note) cannot use raw data here since UMI cannot match cell name
GSE176031_seurat <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/GSE176031/GSE176031_seurat.rds")

# merge meta data and seurat object
assay<- 
  GSE176031_seurat %>% 
  left_join(meta_data, by = c(".cell" = "Cell_name")) %>%
  filter(!is.na(SingleR_annotation))

assay<-assay%>%
  #the cell type column
  rename(cell_type=SingleR_annotation)%>%
  mutate(dataset='GSE176031')

assay$sample<-assay$orig.ident
saveRDS(assay,file='dev/xinpu_datascript/parsed_data/GSE176031_final.rds')

#assay<-readRDS(file='dev/xinpu_datascript/parsed_data/GSE176031_final.rds')

#list unique cell type
assay%>%distinct(dataset,sample,cell_type)
# cell number in each sample
assay%>%group_by(dataset,sample)%>%summarise(count=n())
