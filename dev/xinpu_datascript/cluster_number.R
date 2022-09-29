source('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/third_party_software/sc-SHC/significance_analysis.R')
library(Seurat)
library(tidyseurat)

# Get the original cluster label
#sample.combined<-readRDS(file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')

# try to sample cells from the seurat object and reduce the size
set.seed(111)
cell.name<-sample(x=names(sample.combined@active.ident),size=10000,replace=F)
sample.combined.small<-subset(sample.combined,cells=cell.name)

# load assay data
data.combined<-GetAssayData(sample.combined.small@assays[['SCT']],slot='data')
#data.combined<-GetAssayData(sample.combined)
varfeature<-sample.combined.small@assays[["integrated"]]@var.features
data.combined<-data.combined[rownames(data.combined) %in% varfeature,]
data.combined<-as.matrix(data.combined)

#########################
# Get 500 cells from each sample
#sample.combined.small<-readRDS('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/sample_500.rds')
data.combined<- GetAssayData(sample.combined.small@assays[["SCT"]],slot='data')
varfeature<-sample.combined.small@assays[["integrated"]]@var.features
data.combined<-data.combined[rownames(data.combined) %in% varfeature,]
data.combined<-as.matrix(data.combined)

# Get new clusters
new_cluster<-testClusters(data.combined,as.character(Idents(sample.combined.small)),num_features=2000)
table(new_cluster,Idents(sample.combined.small))

#################################
# Try with one sample's SCT assay
Idents(sample.combined) <- "sample"
single.sample<-subset(sample.combined,idents='AUG_PB1')
single.data<-GetAssayData(single.sample@assays[['SCT']],slot='counts')
single.var<-single.sample@assays[["integrated"]]@var.features
single.data<-single.data[rownames(single.data) %in% single.var,]
single.data<-as.matrix(single.data)
Idents(single.sample)<-"integrated_snn_res.0.6"
new<-testClusters(single.data,as.character(Idents(single.sample)),num_features=2000)
table(new,Idents(single.sample))
