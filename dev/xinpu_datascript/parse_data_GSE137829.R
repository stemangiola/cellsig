library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyseurat)
library(purrr)
library(readxl)


# GOALS OF PARSING
# Create an rds file for each dataset (total fo 14 rds files).
# Each rds file should include the following columns
# cell_type
# sample
# database

# FEEDBACK
# USE here package to refer to paths, rather than absolute paths
# always give meaningfull names for: variables, file, anmes, directory names. Dont use achronims or abbreviations. Imagine e collegue of yours takes your projects and should understand what is what without asking you any question.
# Do not use capitals or spaces in your variable, file, directory names, instead of using CellType, we use cell_type

# assays<-dir('/Users/euphemia/Desktop/xinpu/project/papenfuss/scrna_data/GSE137829_RAW/')
# wd<-'/Users/euphemia/Desktop/xinpu/project/papenfuss/scrna_data/GSE137829_RAW/'
# dirs<-paste0(wd,assays)
# samples_name=c('GSM4089151','GSM4089152','GSM4089153','GSM4089154','GSM4711414',
#                'GSM4711415')
# 

# 
# scRNAlist<-list()
# for (i in 1:length(dirs)){
#   raw_data<-read.csv(dirs[i],sep = "\t")
#   counts<-raw_data[,3:length(raw_data)] #remove Gene_id, remove cell marker names
#   rownames(counts)<-make.names(raw_data[,2],unique=TRUE) #set genes names to be the row names
#   scRNAlist[[i]]<-CreateSeuratObject(counts,project=samples_name[i]) #no quality control
#   
#   # Merge the metadata
#   scRNAlist[[i]] %>% 
#     left_join(X42003_2020_1476_MOESM4_ESM, by = c(".cell" = "...1", "sample = orig.ident")) %>%
#     filter(!is.na(CellType))
#   #scRNAlist[[i]]<- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])  #add prefix to barcodes to avoid duplication of name
# }
# 
# # merge different dataset together
# scrna<-merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
# 
# MY_MERGED_DATA = scRNAlist[[i]] %>% 
#   left_join(X42003_2020_1476_MOESM4_ESM, by = c(".cell" = "...1")) %>%
#   filter(!is.na(CellType))
# 
# MY_MERGED_DATA %>%
#   rename(sample = orig.ident.y, cell_type = CellType) %>%
#   mutate(database = "XXXX")



meta_data <-readxl::read_excel("GSE137829_RAW/42003_2020_1476_MOESM4_ESM.xlsx",sheet='metadata of 6 CRPC')

# ALTERNATIVE WAY USING TIDY PARADIGM
assay =
  dir("GSE137829_RAW/", pattern = "txt.gz", full.names = TRUE) %>%
  map(~{
    counts<-read.csv(.x,sep = "\t") 
    rownames(counts)<-make.names(counts[,2],unique=TRUE)
    #remove Gene_id, remove cell marker names
    counts<-counts[,3:length(counts)] %>%
    # Create seurat container
    CreateSeuratObject(.) %>% 

    # Merge the metadata
    left_join(meta_data, by = c(".cell" = "...1")) %>%
    filter(!is.na(CellType))
  }) %>%
  # Reduce applies one function to all elements of a list, returned by map.
  # bind_rows is the function that merges two of more dataset together.
  purrr::reduce(tidyseurat::bind_rows)

assay<-assay%>%mutate(database='GSE137829',sample=case_when(
  orig.ident.y == "patient #1"~ "GSM4089151",
  orig.ident.y == "patient #2" ~ "GSM4089152",
  orig.ident.y == "patient #3" ~ "GSM4089153",
  orig.ident.y == "patient #4" ~ "GSM4089154",
  orig.ident.y == "patient #5" ~ "GSM4711414",
  orig.ident.y == "patient #6" ~ "GSM4711415"))

saveRDS(scrna,file='/da1/home/penggongxin/cxp/project/parsed_data/GSE137829.rds')

              
