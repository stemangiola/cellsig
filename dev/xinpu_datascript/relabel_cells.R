library(tidyverse)
library(Seurat)
library(tidyseurat)
library(dplyr)
library(ggplot2)
library(xlsx)
source('xinpu_datascript/functions.R')

sample.combined<-readRDS('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')
# Get the cell types from each cluster
cluster_cell<-sample.combined%>%group_by(seurat_clusters)%>%count(seurat_clusters,cell_type)%>%
  filter(n>30)%>%arrange(desc(n),.by_group = TRUE)%>%print(n=Inf)
write_csv(cluster_cell,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/cell_type.csv')

#remove doublet from the datasets
sample.combined<-sample.combined%>%filter(cell_type!='Doublets')

# format the cell_type names
cell_type_old<-sample.combined@meta.data[["cell_type"]]
cell_type_new<-sample.combined@meta.data[["cell_type"]]
new<- read_csv('/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/rename_celltype.csv',col_names =TRUE,col_select=c('cell_type','rename'))
old_name<-new$cell_type
new_name<-new$rename
for (i in 1:length(old_name)){
  cell_type_new<-cell_type_new%>%replace(cell_type_new==old_name[i],new_name[i])
}
sample.combined@meta.data[["cell_type_old"]]<-cell_type_old
sample.combined@meta.data[["cell_type"]]<-cell_type_new
#saveRDS(sample.combined,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')

# Visualize umap of sample.combined to filter cells for further clustering
DimPlot(sample.combined,reduction='umap',label=TRUE,repel=TRUE)

#sample.combined<-readRDS('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')
DefaultAssay(sample.combined) <- 'SCT'
FeaturePlot(sample.combined, features = c("CD14", "FCGR3A", "CD79A", "CD3G", "EPCAM", 'LUM',"VIM", "PLVAP",
                                          "CD68"))
FeaturePlot(sample.combined, features = c('LUM',"VIM",'DCN','THY1','COL1A1','COL14A1','FN1','CAV1'))

DefaultAssay(sample.combined) <- 'integrated'
# Get cell marker in clusters
cluster.markers<-FindAllMarkers(sample.combined,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
var.heatmap<-DoHeatmap(sample.combined, features = top10$gene) + NoLegend()
saveRDS(cluster.markers,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/cluster_marker.rds')


#Sanity check
# Get percentage of mitochondria
mt<-sample.combined@meta.data[["percent.mt"]]
FeaturePlot(sample.combined,reduction='umap',features='percent.mt',order = TRUE)
# Color by total RNA
FeaturePlot(sample.combined,reduction='umap',features='nCount_RNA',order = TRUE)
# Color by sample
DimPlot(sample.combined,reduction='umap',group.by = 'sample')+ggplot2::theme(legend.position = 'bottom')
# 10X vs SMART-Seq (GSE176031)
DimPlot(sample.combined,reduction='umap',group.by = 'dataset',cols=c('GSE176031'='red'))+ggplot2::ggtitle('SMART-Seq')
# cancer vs normal改
sample<-sample.combined@meta.data[["sample"]]
cancer<-c()
for (i in 1:length(sample)){
  cell_type_new<-cell_type_new%>%replace(cell_type_new==old_name[i],new_name[i])
}
DimPlot(sample.combined,reduction='umap',group.by = 'sample',cols=c(tumor_set='red'))+ggplot2::ggtitle('cancer sample')
DimPlot(sample.combined,reduction='umap',group.by = 'sample',cols=c(normal_set='grey'))+ggplot2::ggtitle('cancer sample')

#markers<-ComputeMarkers
# #filter epithelial and non-epithelial cells
# immune_cluster<-c('0','10','12','13','14','17','18','20','21','23','26','27','28','31','32','34')
# epithelial<-sample.combined%>%filter(!as.character(sample.combined$seurat_clusters)%in%immune_cluster)
# other_cell<-sample.combined%>%filter(as.character(sample.combined$seurat_clusters)%in%immune_cluster)
# # delete the old integrated assay
# DefaultAssay(epithelial)<-'RNA'
# DefaultAssay(other_cell)<-'RNA'
# epithelial[['integrated']] <- NULL
# other_cell[['integrated']] <- NULL
# 
# # split epithelial and non-epithelial by sample and rerun SCT to find the highly variable genes
# epi.list<-SplitObject(epithelial,split.by='sample')
# epi.list<-lapply(epi.list,function(x){
#   x<-SCTransform(x)
# })
# 
# other.list<-SplitObject(other_cell,split.by='sample')
# other.list<-lapply(other.list,function(x){
#   x<-SCTransform(x)
# })
# 
# # re-cluster the cell by epithelial and non-epithelial type separately 
# # epithelial,luminal and cancer
# epi.features<- SelectIntegrationFeatures(object.list = epi.list) 
# epi.list.int<-PrepSCTIntegration(object.list = epi.list,anchor.features = epi.features)
# # anchore with reference; reference is 'JD1800156SL' largest cell number
# epithelial%>% count(sample,sort=TRUE)%>%print(n=Inf)
# which(names(epi.list)=='JD1800156SL')
# epi.anchor.ref<-FindIntegrationAnchors(object.list=epi.list.int,reference=11,normalization.method = 'SCT',anchor.features = epi.features)
# # The integration result
# epithelial<-IntegrateData(anchorset=epi.anchor.ref,normalization.method='SCT')
# 
# epithelial<-ScaleData(epithelial,verbose=FALSE)
# epithelial<-RunPCA(epithelial,npc=50,verbose=FALSE)
# epithelial<-RunUMAP(epithelial,reduction='pca',dims=1:30,spread  = 0.5,min.dist  = 0.01, n.neighbors = 10)
# epithelial<-FindNeighbors(epithelial,reduction='pca',dims=1:30)
# epithelial<-FindClusters(epithelial,resolution=2)
# p1<-DimPlot(epithelial,reduction='umap',label=TRUE,repel=TRUE)
# #saveRDS(epithelial,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/epithelial.rds')
# 
# # print the cell names in the epithelial clusters
# epi<-epithelial%>%group_by(seurat_clusters)%>%count(seurat_clusters,cell_type)%>%
#   filter(n>30)%>%arrange(desc(n),.by_group = TRUE)%>%print(n=Inf)
# #write_csv(epi,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/epithelial.csv')
# 
# epithelial%>%count(seurat_clusters,cell_type)%>%filter(n>30)%>%with_groups(seurat_clusters,~ .x %>%summarise(cell_types = paste(cell_type,collapse=',')))%>%print(n=Inf)
# # non-epithelial
# # remove sample with less than 100 cells for integration
# other_cell%>%count(sample,sort=TRUE)%>%print(n=Inf)
# other.list$GSM4089152<-NULL
# other.list$PR5254_T<-NULL
# other.list$`SLX-15736SIGAE8HTHM2BBXXs_4`<-NULL
# other.features<- SelectIntegrationFeatures(object.list = other.list) 
# other.list.int<-PrepSCTIntegration(object.list = other.list,anchor.features = other.features)
# # anchore with reference; reference is GSM4711414, largest set
# other.anchor.ref<-FindIntegrationAnchors(object.list=other.list.int,reference=6,normalization.method = 'SCT',anchor.features = other.features)
# # The integration result
# other_cell<-IntegrateData(anchorset=other.anchor.ref,normalization.method='SCT')
# 
# other_cell<-ScaleData(other_cell,verbose=FALSE)
# other_cell<-RunPCA(other_cell,npc=50,verbose=FALSE)
# other_cell<-RunUMAP(other_cell,reduction='pca',dims=1:50)
# other_cell<-FindNeighbors(other_cell,reduction='pca',dims=1:50)
# other_cell<-FindClusters(other_cell,resolution=2)
# p2<-DimPlot(other_cell,reduction='umap',label=TRUE,repel=TRUE)
# #saveRDS(other_cell,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/other_cell.rds')
# 
# other<-other_cell%>%group_by(seurat_clusters)%>%count(seurat_clusters,cell_type)%>%
#   filter(n>30)%>%arrange(desc(n),.by_group = TRUE)%>%print(n=Inf)
# #write_csv(other,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/other_type.csv')
# 
# other_cell%>%count(seurat_clusters,cell_type)%>%filter(n>30)%>%with_groups(seurat_clusters,~ .x %>%summarise(cell_types = paste(cell_type,collapse=',')))%>%print(n=Inf)
