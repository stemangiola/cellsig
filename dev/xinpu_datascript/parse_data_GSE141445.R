library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyseurat)
library(purrr)
library(stringr)

raw_data<- read.table("./dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/GSE141445/GSM4203181_data.matrix.txt"
                      ,sep="\t",header=TRUE) 
meta_data<-read_csv("dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/GSE141445/celltypes_NCB_PCa13.csv")
#process meta_data's UMI make it the same as our raw_data
umi_new<-as.vector(str_replace(meta_data$...1,'-','.'))
meta_data['...1']<-umi_new
assay<- CreateSeuratObject(raw_data) 
#integrate meta_data
assay<-assay%>%
  left_join(meta_data, by = c(".cell" = "...1"))%>%
  filter(!is.na(type))%>%
  mutate(orig.ident.x='GSM4203181',dataset='GSE141445')

assay<-assay%>%rename(sample=orig.ident.x,cell_type=type)
saveRDS(assay,file='dev/xinpu_datascript/parsed_data/GSE141445_final.rds')
#assay<-readRDS(file='dev/xinpu_datascript/parsed_data/GSE141445_final.rds')

#list unique cell type
assay%>%distinct(dataset,sample,cell_type)
# cell number in each sample
assay%>%group_by(dataset,sample)%>%summarise(count=n())

