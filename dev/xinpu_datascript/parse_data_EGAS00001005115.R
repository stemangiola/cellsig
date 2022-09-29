library(Seurat)
library(dplyr)
library(tidyseurat)
library(tidyr)
library(purrr)
library(readxl)
library(stringr)

counts<-Matrix::readMM('dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/EGAS00001005115/EGAS00001005115/Wu_etal_2021_allcells_raw_counts.mtx')
#get row names and column names of the matrix
row<- read.table("dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/EGAS00001005115/EGAS00001005115/Wu_etal_2021_allcells_genes.tsv", header=F, sep="\t") 
coln<-read.table("dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/EGAS00001005115/EGAS00001005115/Wu_etal_2021_allcells_barcodes.tsv", header=F, sep="\t")

#create Seurat object
assay<-matrix(as.numeric(as.matrix(counts)),nrow=nrow(counts),
       dimnames=list(t(row['V1']),t(coln['V1'])))%>%
  CreateSeuratObject(.,project='prostate')
assay$umi_code<-colnames(assay)

meta_data<-read.delim('dev/xinpu_datascript/PPCG_deconvolution_signatures_RAW_DATA/EGAS00001005115/EGAS00001005115/Wu_etal_2021_metadata.txt',sep='\t',header=TRUE)%>%
  .[2:nrow(.),]
assay<-assay%>%left_join(meta_data, by = c("umi_code" = "NAME",'orig.ident'='donor_id'))%>%
  filter(!is.na(CellType))%>%
  mutate(dataset='EGAS00001005115')
assay<-assay%>%rename(cell_type=CellType)
assay<-assay%>%rename(sample=biosample_id)

#remove breast cancer and keep prostate cancer only
assay<-EGAS00001005115_final%>%filter(disease__ontology_label=='prostate cancer')
saveRDS(assay,file='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/parsed_data/EGAS00001005115_final.rds')

#assay<-readRDS(file='dev/xinpu_datascript/parsed_data/EGAS00001005115_final.rds')

#list unique cell type
assay%>%distinct(dataset,sample,cell_type)
# cell number in each sample
assay%>%group_by(dataset,sample)%>%summarise(count=n())



