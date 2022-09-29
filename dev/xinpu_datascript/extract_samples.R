library(tidyverse)
library(Seurat)
library(tidyseurat)
library(dplyr)
library(ggplot2)

sample.combined<-readRDS('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')
samples<- SplitObject(sample.combined,split.by='sample')

set.seed(111)
sample500<-lapply(samples,function(t){
  if(length(Idents(t))>500){
    cell.name<-sample(x=names(t@active.ident),size=500,replace=F)
    t<-subset(t,cells=cell.name)
  }else{
    t<-t
  }
})

sample.small500<-sample500%>%purrr::reduce(tidyseurat::bind_rows)
#sample.small500<-readRDS('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/sample_500.rds')
sample.small500[['integrated']]@var.features<-sample.combined[['integrated']]@var.features
saveRDS(sample.small500,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/sample_500.rds')
#job::job({saveRDS(sample.small500,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/sample_500.rds')})

