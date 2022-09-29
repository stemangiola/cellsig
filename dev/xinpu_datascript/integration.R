library(tidyverse)
library(Seurat)
library(tidyseurat)
library(dplyr)
library(ggplot2)

# Get the rds files of the samples
readfile<-function(x){
  dir='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/results/preprocessing_results/preprocessing_output/'
  readRDS(paste0(dir,x))
}

dir='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/results/preprocessing_results/preprocessing_output/'
name.list<- list.files(path=dir, pattern='.rds')
sample.list<-lapply(name.list,readfile)

# Integration for SCT assay
# Use the default feature number 2000
features<- SelectIntegrationFeatures(object.list = sample.list) 
sample.list.int<-PrepSCTIntegration(object.list = sample.list,anchor.features = features)

# anchore with reference
anchor.ref<-FindIntegrationAnchors(object.list=sample.list.int,reference=6,normalization.method = 'SCT',anchor.features = features)
# The integration result
sample.combined<-IntegrateData(anchorset=anchor.ref,normalization.method='SCT',features.to.integrate = )
#saveRDS(sample.combined,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')
sample.combined<-readRDS(file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')

# Integration analysis
sample.combined<-ScaleData(sample.combined,verbose=FALSE)
sample.combined<-RunPCA(sample.combined,npc=50,verbose=FALSE)
sample.combined<-RunUMAP(sample.combined,reduction='pca',dims=1:50)
sample.combined<-FindNeighbors(sample.combined,reduction='pca',dims=1:50)
sample.combined<-FindClusters(sample.combined,resolution=2)
#saveRDS(sample.combined,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')

#sample.combined<-readRDS('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')
#visualization
p1<-DimPlot(sample.combined,reduction='umap',label=TRUE,repel=TRUE)
print(p1)
#group by dataset
p2<-DimPlot(sample.combined,reduction='umap',group.by = 'dataset')
p3<-DimPlot(sample.combined,reduction='umap',split.by = 'dataset')
#group by sample
p4<-DimPlot(sample.combined,reduction='umap',group.by='sample')+ggplot2::theme(legend.position = 'bottom')
p4
p5<-DimPlot(sample.combined,reduction='umap',group.by='cell_type')+ggplot2::theme(legend.position = 'bottom')
# saveRDS(sample.combined,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')


#cell type with more than 30 cells
sample.combined%>%count(seurat_clusters,cell_type)
sample.combined%>%count(seurat_clusters,cell_type)%>%filter(n>30)%>%with_groups(seurat_clusters,~ .x %>%summarise(cell_types = paste(cell_type,collapse=',')))%>%print(n=Inf)
