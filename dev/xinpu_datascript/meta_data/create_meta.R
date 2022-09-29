library(Seurat)
library(dplyr)
library(tidyseurat)

#EGAS00001005115_meta
EGAS00001005115_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/EGAS00001005115_final.rds")
EGAS00001005115_meta<-tibble(.cell=names(Idents(EGAS00001005115_final)))
saveRDS(EGAS00001005115_meta,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/meta_data/EGAS00001005115_meta.rds')

#EGAS00001005787_meta
EGAS00001005787_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/EGAS00001005787_final.rds")
EGAS00001005787_meta<-tibble(.cell=names(Idents(EGAS00001005787_final)))
saveRDS(EGAS00001005787_meta,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/meta_data/EGAS00001005787_meta.rds')

#GSE137829_meta
GSE137829_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE137829_final.rds")
GSE137829_meta<-tibble(.cell=names(Idents(GSE137829_final)))
saveRDS(GSE137829_meta,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/meta_data/GSE137829_meta.rds')

#GSE141445_meta
GSE141445_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE141445_final.rds")
GSE141445_meta<-tibble(.cell=names(Idents(GSE141445_final)))
saveRDS(GSE141445_meta,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/meta_data/GSE141445_meta.rds')

#GSE176031_meta
GSE176031_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE176031_final.rds")
GSE176031_meta<-tibble(.cell=names(Idents(GSE176031_final)))
saveRDS(GSE176031_meta,file='/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/meta_data/GSE176031_meta.rds')





