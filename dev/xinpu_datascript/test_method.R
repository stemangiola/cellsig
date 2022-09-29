library(BiocFileCache)
library(Seurat)
set.seed(73122)

# Load in 10X data
bfc <- BiocFileCache(ask=FALSE)
path <- bfcrpath(bfc,'https://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/293t_filtered_gene_bc_matrices.tar.gz')
tmp <- tempfile()
untar(path,exdir=tmp)
data1 <- Read10X(data.dir = file.path(tmp,'filtered_matrices_mex','hg19')) 

data1 <- as.matrix(data1[,1:200])
top50 <- order(rowSums(data1),decreasing=T)[1:50]
bottom50 <- order(rowSums(data1),decreasing=F)[1:50]
data1[c(top50,bottom50),1:100] <- data1[sample(c(top50,bottom50),replace=F),1:100]


data.seurat <- CreateSeuratObject(data1)
data.seurat <- NormalizeData(data.seurat)
data.seurat <- FindVariableFeatures(data.seurat,nfeatures=100)
data.seurat <- ScaleData(data.seurat)
data.seurat <- RunPCA(data.seurat,npcs=6)
data.seurat <- RunUMAP(data.seurat,features=VariableFeatures(data.seurat))
data.seurat <- FindNeighbors(data.seurat,dims=1:6)
data.seurat <- FindClusters(data.seurat,resolution=1)

new_seurat <- testClusters(data1,as.character(Idents(data.seurat)),num_PCs=6,
                           num_features=100,parallel=F)
table(new_seurat,Idents(data.seurat))


