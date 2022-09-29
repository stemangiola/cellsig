# We check if there are samples within a dataset has a particular poor quality.
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(dplyr)
library(ggplot2)

#EGAS00001005115
#EGAS00001005115_final <-readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/EGAS00001005115_final.rds")
EGAS00001005115_final<-NormalizeData(EGAS00001005115_final)
EGAS00001005115_final<-ScaleData(EGAS00001005115_final,verbose=FALSE)
EGAS00001005115_final<- FindVariableFeatures(EGAS00001005115_final, selection.method = "vst", nfeatures = 2000)
EGAS00001005115_final<-RunPCA(EGAS00001005115_final,npc=30,verbose=FALSE)
EGAS00001005115_final<-RunUMAP(EGAS00001005115_final,reduction='pca',dims=1:30)
EGAS00001005115_final<-FindNeighbors(EGAS00001005115_final,reduction='pca',dims=1:30)
EGAS00001005115_final<-FindClusters(EGAS00001005115_final,resolution=0.6)
p1<-DimPlot(EGAS00001005115_final,reduction='umap',split.by = 'sample')
p6<-DimPlot(EGAS00001005115_final,reduction='umap',group.by = 'cell_type',label=TRUE)
p11<-DimPlot(EGAS00001005115_final,reduction='umap',group.by = 'sample',shuffle = TRUE)


#EGAS00001005787
#EGAS00001005787_final <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/EGAS00001005787_final.rds")
EGAS00001005787_final <-NormalizeData(EGAS00001005787_final)
EGAS00001005787_final<-ScaleData(EGAS00001005787_final,verbose=FALSE)
EGAS00001005787_final<- FindVariableFeatures(EGAS00001005787_final, selection.method = "vst", nfeatures = 2000)
EGAS00001005787_final<-RunPCA(EGAS00001005787_final,npc=30,verbose=FALSE)
EGAS00001005787_final<-RunUMAP(EGAS00001005787_final,reduction='pca',dims=1:30)
EGAS00001005787_final<-FindNeighbors(EGAS00001005787_final,reduction='pca',dims=1:30)
EGAS00001005787_final<-FindClusters(EGAS00001005787_final,resolution=0.6)
p2<-DimPlot(EGAS00001005787_final,reduction='umap',group.by = 'sample',shuffle = TRUE)
p7<-DimPlot(EGAS00001005787_final,reduction='umap',group.by = 'cell_type',label=TRUE)

#GSE137829
#GSE137829_final <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE137829_final.rds")
GSE137829_final <-NormalizeData(GSE137829_final)
GSE137829_final<-ScaleData(GSE137829_final,verbose=FALSE)
GSE137829_final<- FindVariableFeatures(GSE137829_final, selection.method = "vst", nfeatures = 2000)
GSE137829_final<-RunPCA(GSE137829_final,npc=30,verbose=FALSE)
GSE137829_final<-RunUMAP(GSE137829_final,reduction='pca',dims=1:30)
GSE137829_final<-FindNeighbors(GSE137829_final,reduction='pca',dims=1:30)
GSE137829_final<-FindClusters(GSE137829_final,resolution=0.6)
p3<-DimPlot(GSE137829_final,reduction='umap',split.by = 'sample')

#GSE141445
#GSE141445_final <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE141445_final.rds")
GSE141445_final <-NormalizeData(GSE141445_final)
GSE141445_final<-ScaleData(GSE141445_final,verbose=FALSE)
GSE141445_final<- FindVariableFeatures(GSE141445_final, selection.method = "vst", nfeatures = 2000)
GSE141445_final<-RunPCA(GSE141445_final,npc=30,verbose=FALSE)
GSE141445_final<-RunUMAP(GSE141445_final,reduction='pca',dims=1:30)
GSE141445_final<-FindNeighbors(GSE141445_final,reduction='pca',dims=1:30)
GSE141445_final<-FindClusters(GSE141445_final,resolution=0.6)
p4<-DimPlot(GSE141445_final,reduction='umap',split.by = 'sample')
p8<-DimPlot(GSE141445_final,reduction='umap',group.by = 'sample',shuffle = TRUE)
p9<-DimPlot(GSE141445_final,reduction='umap',group.by = 'cell_type')

#GSE176031
#GSE176031_final <- readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE176031_final.rds")
GSE176031_final<-NormalizeData(GSE176031_final)
GSE176031_final<-ScaleData(GSE176031_final,verbose=FALSE)
GSE176031_final<- FindVariableFeatures(GSE176031_final, selection.method = "vst", nfeatures = 2000)
GSE176031_final<-RunPCA(GSE176031_final,npc=30,verbose=FALSE)
GSE176031_final<-RunUMAP(GSE176031_final,reduction='pca',dims=1:30)
GSE176031_final<-FindNeighbors(GSE176031_final,reduction='pca',dims=1:30)
GSE176031_final<-FindClusters(GSE176031_final,resolution=0.6)
p5<-DimPlot(GSE176031_final,reduction='umap',split.by = 'sample')
p10<-DimPlot(GSE176031_final,reduction='umap',group.by = 'sample',shuffle = TRUE)

# check mast cells
sample.combined<-sample.combined%>%mutate(is_mast_cell = grepl('mast', cell_type, ignore.case=TRUE ))
DimPlot(sample.combined,reduction='umap',group.by = 'is_mast_cell',shuffle = TRUE)

# Check normal sample and tumor sample
# EGAS00001005787
EGAS00001005787_normal<-EGAS00001005787_final%>%filter(group=='normal') %>% distinct(sample)
EGAS00001005787_tumor<-EGAS00001005787_final%>%filter(group=='tumor') %>% distinct(sample)

# EGAS00001005115
EGAS00001005115_tumor<-EGAS00001005115_final%>%filter(disease__ontology_label=='prostate cancer')%>%distinct(sample)

# GSE137829
GSE137829_tumor<-GSE137829_final%>%distinct(sample)

#GSE141445
GSE141445_tumor<-GSE141445_final%>%distinct(sample)

#GSE176031
GSE176031_normal<-GSE176031_final%>%filter(grepl('N',sample))%>%distinct(sample)
GSE176031_tumor<-GSE176031_final%>%filter(grepl('T',sample))%>%distinct(sample)

# bind the tumor sample from different datasets together
tumor_set<-bind_rows(EGAS00001005787_tumor,GSE141445_tumor,GSE176031_tumor,EGAS00001005115_tumor,GSE137829_tumor)%>%print(n=Inf)
#saveRDS(tumor_set,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/tumor_set')
# bind the normal sample from different datasets together
normal_set<-bind_rows(EGAS00001005787_normal,GSE176031_normal)
#saveRDS(normal_set,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/normal_set')


