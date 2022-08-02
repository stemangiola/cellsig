library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyseurat)
library(purrr)
library(readxl)



meta_data <-readxl::read_excel("dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/GSE137829/42003_2020_1476_MOESM4_ESM.xlsx",sheet='metadata of 6 CRPC')

# ALTERNATIVE WAY USING TIDY PARADIGM

assay =
  dir("dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/GSE137829", pattern = "txt.gz", full.names = TRUE) %>%
  map(~{
    counts<-read.table(.x,sep = "\t",header=TRUE)
    rownames(counts)<-make.names(counts[,2],unique=TRUE)
    #remove Gene_id, remove cell marker names
    counts[,3:length(counts)] %>%
    # Create seurat container
    CreateSeuratObject(.)
  }) 

for (i in 1:length(assay)){
  assay[[i]]<-AddMetaData(
    object=assay[[i]],
    metadata=paste('patient #',i,sep = ""),
    col.name = 'patient'
  )
  assay[[i]]<-assay[[i]]%>%
    left_join(meta_data,by=c('.cell'='...1','patient'='orig.ident'))%>%
    filter(!is.na(CellType))
}

assay<-assay%>%purrr::reduce(tidyseurat::bind_rows)

assay<-assay%>%mutate(dataset='GSE137829',sample=case_when(
  patient == "patient #1"~ "GSM4089151",
  patient == "patient #2" ~ "GSM4089152",
  patient == "patient #3" ~ "GSM4089153",
  patient == "patient #4" ~ "GSM4089154",
  patient == "patient #5" ~ "GSM4711414",
  patient == "patient #6" ~ "GSM4711415"))%>%
  rename(cell_type=CellType)


saveRDS(assay,file='dev/xinpu_datascript/parsed_data/GSE137829_final.rds')
#assay<-readRDS('dev/xinpu_datascript/parsed_data/GSE137829_final.rds')


#list unique cell type
assay%>%distinct(dataset,sample,cell_type)
# cell number in each sample
assay%>%group_by(dataset,sample)%>%summarise(count=n())







              
