library(tidyverse)
library(Seurat)
library(tidyseurat)
library(dplyr)
library(ggplot2)
library(patchwork)

sample.combined<-readRDS('/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/integrated_sample.rds')
# Test different resolution of clustering
sample_01<-sample.combined<-FindClusters(sample.combined,resolution=0.1)
sample_02<-sample.combined<-FindClusters(sample.combined,resolution=0.2)
sample_03<-sample.combined<-FindClusters(sample.combined,resolution=0.3)
sample_04<-sample.combined<-FindClusters(sample.combined,resolution=0.4)
sample_05<-sample.combined<-FindClusters(sample.combined,resolution=0.5)
sample_06<-sample.combined<-FindClusters(sample.combined,resolution=0.6)
sample_07<-sample.combined<-FindClusters(sample.combined,resolution=0.7)
sample_08<-sample.combined<-FindClusters(sample.combined,resolution=0.8)
sample_09<-sample.combined<-FindClusters(sample.combined,resolution=0.9)
sample_10<-sample.combined<-FindClusters(sample.combined,resolution=1)

# We can have label,but repel(repel label) need to be False. 
# Or patchwork cannot arrange figure size.
p1<-DimPlot(sample_01,reduction='umap',label=TRUE)+ggtitle('resolution 0.1')
p2<-DimPlot(sample_02,reduction='umap',label=TRUE)+ggtitle('resolution 0.2')
p3<-DimPlot(sample_03,reduction='umap',label=TRUE)+ggtitle('resolution 0.3')
p4<-DimPlot(sample_04,reduction='umap',label=TRUE)+ggtitle('resolution 0.4')
p5<-DimPlot(sample_05,reduction='umap',label=TRUE)+ggtitle('resolution 0.5')
p6<-DimPlot(sample_06,reduction='umap',label=TRUE)+ggtitle('resolution 0.6')
p7<-DimPlot(sample_07,reduction='umap',label=TRUE)+ggtitle('resolution 0.7')
p8<-DimPlot(sample_08,reduction='umap',label=TRUE)+ggtitle('resolution 0.8')
p9<-DimPlot(sample_09,reduction='umap',label=TRUE)+ggtitle('resolution 0.9')
p10<-DimPlot(sample_10,reduction='umap',label=TRUE)+ggtitle('resolution 1')

(p1+p2+p3)/(p4+p5+p6)
(p7+p8)/(p9+p10)

wrap_plots(p1,p2,p3,p4)
