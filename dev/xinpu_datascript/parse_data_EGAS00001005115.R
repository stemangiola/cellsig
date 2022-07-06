library(Seurat)
library(dplyr)
library(tidyseurat)
library(tidyr)
library(purrr)
library(readxl)
library(stringr)

counts<-Matrix::readMM('EGAS00001005115_RAW/EGAS00001005115/Wu_etal_2021_allcells_raw_counts.mtx')
#get row names and column names of the matrix
row<- read.table("EGAS00001005115_RAW/EGAS00001005115/Wu_etal_2021_allcells_genes.tsv", header=F, sep="\t") 
coln<-read.table("EGAS00001005115_RAW/EGAS00001005115/Wu_etal_2021_allcells_barcodes.tsv", header=F, sep="\t")

#create Seurat object
assay<-matrix(as.numeric(as.matrix(counts)),nrow=nrow(counts),
       dimnames=list(t(row['V1']),t(coln['V1'])))%>%
  CreateSeuratObject(.,project='prostate')
assay$umi_code<-colnames(assay)

meta_data<-read.delim('EGAS00001005115_RAW/EGAS00001005115/Wu_etal_2021_metadata.txt',sep='\t',header=TRUE)%>%
  .[2:nrow(.),]
assay<-assay%>%left_join(meta_data, by = c("umi_code" = "NAME",'orig.ident'='donor_id'))%>%
  filter(!is.na(CellType))%>%
  mutate(database='EGAS00001005115')
assay<-assay%>%rename(cell_type=CellType)
saveRDS(assay,file='parsed_data/EGAS00001005115.rds')






