library(Seurat)
library(dplyr)
library(tidyseurat)
library(purrr)
library(job)

#parsed seurat data
dir="dev/xinpu_datascript/parsed_data/"
files<-list.files(path=dir,pattern='_final.rds')
folder<-paste0(dir,files)
assay_list<-lapply(folder,readRDS)
prostate_samples<-assay_list%>%purrr::reduce(tidyseurat::bind_rows)
#put into background
job::job({saveRDS(prostate_samples,file='dev/xinpu_datascript/parsed_data/5datasets.rds')})

prostate_samples<-readRDS(file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/5datasets.rds')
#get cell number per sample
prostate_samples%>%count(dataset,sample)%>% print(n = Inf)
#get cell type per sample
prostate_samples%>%distinct(dataset,cell_type)%>% print(n = Inf)

# Get the cell numbers per sample
sample.combined<-readRDS(file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')
sample.combined%>%count(sample,dataset)%>%print(n=Inf)

sample.combined%>%count(dataset,cell_type)%>%print(n=Inf)

old_celltype<-sample.combined%>%distinct(cell_type)%>%print(n=Inf)                        
#write.csv(old_celltype,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/rename_celltype.csv')                                                

# Check sample which contains cancer cells
sample.combined%>%filter(cell_type=='cancer_associated_fibroblast')%>%select(dataset,sample,cell_type)%>%count(dataset,sample)

