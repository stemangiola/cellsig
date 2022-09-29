library(Seurat)
library(tidyseurat)

fun.save<-function(x){
  #/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript
  file.name<-paste0('/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/samples/',x$sample[[1]])
  file.name<-paste0(file.name,'.rds')
  #print(file.name)
  saveRDS(x,file=file.name)
}

#EGAS00001005787_final
EGAS00001005787_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/EGAS00001005787_final.rds")
EGAS00001005787_final <-RenameAssays(object = EGAS00001005787_final, originalexp = 'RNA')
EGAS00001005787.list<-SplitObject(EGAS00001005787_final,split.by='sample')
lapply(EGAS00001005787.list,fun.save)

#EGAS00001005115_final
EGAS00001005115_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/EGAS00001005115_final.rds")
EGAS00001005115.list<-SplitObject(EGAS00001005115_final,split.by='sample')
lapply(EGAS00001005115.list,fun.save)

#GSE137829
GSE137829_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE137829_final.rds")
GSE137829.list<-SplitObject(GSE137829_final,split.by='sample')
lapply(GSE137829.list,fun.save)

#GSE141445
GSE141445_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE141445_final.rds")
GSE141445.list<-SplitObject(GSE141445_final,split.by='sample')
lapply(GSE141445.list,fun.save)

#GSE176031
GSE176031_final <- readRDS("/stornext/Home/data/allstaff/c/chen.x/master_project/cellsig/dev/xinpu_datascript/parsed_data/GSE176031_final.rds")
GSE176031.list<-SplitObject(GSE176031_final,split.by='sample')
lapply(GSE176031.list,fun.save)





