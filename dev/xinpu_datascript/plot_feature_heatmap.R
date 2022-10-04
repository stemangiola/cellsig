library(tidyverse)
library(Seurat)
library(tidyseurat)
library(dplyr)
library(ggplot2)
library(xlsx)

sample.combined<-readRDS('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')
# Plot certain marker gene in umap
DefaultAssay(sample.combined) <- 'SCT'
#monocyte, monocyte(NK), b cell, T cell
FeaturePlot(sample.combined, features = c("CD14", "FCGR3A", "CD79A", "CD3G", "EPCAM", 'LUM',"VIM", "PLVAP",
                                          "CD68"))
# fibroblast choose LUMN
FeaturePlot(sample.combined, features = c('LUM',"VIM",'DCN','THY1','COL1A1','COL14A1','FN1','CAV1'))
# Basal
FeaturePlot(sample.combined, features ='KRT17' )
DefaultAssay(sample.combined) <- 'integrated'

# Get cell marker in clusters
cluster_T<-sample.combined%>%filter(seurat_clusters==c(0,12,23,20,12))
cluster_epithelial<-sample.combined%>%filter(seurat_clusters==c(1,2,3,4,5,6,7,8,9,11,15,16,19,22,24,25,28,29,30,35,36,37,38))
cluster_endothelial<-sample.combined%>%filter(seurat_clusters==c(17,26))
cluster_fibroblast<-sample.combined%>%filter(seurat_clusters==c(14,18,21))
# B cell, monocyte, mast cell
cluster_other<-sample.combined%>%filter(seurat_clusters==c(10,27,31,32,34))

# Get the markers
cluster.markers_T<-FindAllMarkers(cluster_T,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers_epi<-FindAllMarkers(cluster_epithelial,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers_endo<-FindAllMarkers(cluster_endothelial,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers_fibro<-FindAllMarkers(cluster_fibroblast,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers_other<-FindAllMarkers(cluster_other,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Find top10 gene
top10_T<-cluster.markers_T %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

top10_epi<-cluster.markers_epi %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

top10_endo<-cluster.markers_endo %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

top10_fibro<-cluster.markers_fibro %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

top10_other<-cluster.markers_other %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


#Produce heatmap and save as pdf
var.heatmap_T<-DoHeatmap(cluster_T, features = top10_T$gene,size=2) + NoLegend()+
  theme(text=element_text(size =5))
ggsave('/home/users/allstaff/chen.x/master_project/cellsig/feature_heatmap_T.pdf',height=9, plot=var.heatmap_T)

var.heatmap_epi<-DoHeatmap(cluster_epithelial, features = top10_epi$gene,size=2) + NoLegend()+
  theme(text=element_text(size =5))
ggsave('/home/users/allstaff/chen.x/master_project/cellsig/feature_heatmap_epi.pdf',height=9, plot=var.heatmap_epi)

var.heatmap_endo<-DoHeatmap(cluster_endothelial, features = top10_endo$gene,size=2) + NoLegend()+
  theme(text=element_text(size =5))
ggsave('/home/users/allstaff/chen.x/master_project/cellsig/feature_heatmap_endo.pdf',height=9, plot=var.heatmap_endo)

var.heatmap_fibro<-DoHeatmap(cluster_fibroblast, features = top10_fibro$gene,size=2) + NoLegend()+
  theme(text=element_text(size =5))
ggsave('/home/users/allstaff/chen.x/master_project/cellsig/feature_heatmap_fibro.pdf',height=9, plot=var.heatmap_fibro)

var.heatmap_other<-DoHeatmap(cluster_other, features = top10_other$gene,size=2) + NoLegend()+
  theme(text=element_text(size =5))
ggsave('/home/users/allstaff/chen.x/master_project/cellsig/feature_heatmap_other.pdf',height=9, plot=var.heatmap_other)



cluster.markers<-FindAllMarkers(sample.combined,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#cluster.markers<-readRDS(file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/cluster_marker.rds')
#saveRDS(cluster.markers,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/cluster_marker.rds')
cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
#Produce heatmap and save as pdf
var.heatmap<-DoHeatmap(sample.combined, features = top10$gene,size=2) + NoLegend()
ggsave('feature_heatmap.pdf',plot=var.heatmap)
# using ggsave to save the heatmap as a pdf

