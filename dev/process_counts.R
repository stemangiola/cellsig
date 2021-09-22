library(tidyverse)
library(cellsig)
library(tidybulk)
library(tidySummarizedExperiment)

counts_first_db_raw = readRDS("dev/counts_first_db_raw.rds")
load("dev/counts_second_db_raw.rda")
counts_second_db_raw = new_data %>% select(sample, symbol, count, dataset, cell_type)

data("tree")

options("tidybulk_do_validate"= FALSE) 

# Join datasets
counts_first_db_raw %>%
  select(-cell_type_original) %>%
  bind_rows(counts_second_db_raw %>% rename(database = dataset)) %>% 
  select(sample, cell_type, symbol, count) %>%
  
  # Parse into hierarchical dataset
  tree_and_signatures_to_database(tree, ., sample, cell_type, symbol, count)  %>%
  
  # Remove redundant samples
  remove_redundancy(sample, symbol, count, correlation_threshold = 0.999, top = 500, method = "correlation") %>%
  droplevels() %>% 
  
  # Eliminate suspicious samples
  filter(!grepl("GSM3722278|GSM3722276|GSM3722277", sample)) %>%
  
  # eliminate genes that are not in all cell types level 1
  nest(data = -c(level_1, symbol)) %>%
  add_count( symbol) %>%
  filter(n==4) %>%
  select(-n) %>%
  unnest(data) %>%
  
  # Convert to SE
  as_SummarizedExperiment(sample, symbol, count)  %>%
  
  # Scale with first degree imputation. 
  # This because there are no common genes to all samples
  impute_missing_abundance(~ cell_type, suffix="") %>%
  identify_abundant() %>%
  scale_abundance() %>%
  filter(!.imputed) %>%
  select(-count_imputed, -.imputed )  %>%
  
  # Just needed for the old version
  select(-one_of("exposure_rate")) %>%
  
  # Calculate exposure for Bayes model
  mutate(exposure_rate = -log(multiplier)) %>%
  
  # Save
  saveRDS("dev/counts.rds", compress = "xz")



# symbol   level_1         n
# <chr>    <chr>       <int>
#   1 AATK-AS1 immune_cell     1
# 2 ACN9     immune_cell     1
# 3 ACPT     immune_cell     1
# 4 ACRC     immune_cell     1
# 5 ADCK3    immune_cell     1
# 6 ADCK4    immune_cell     1
# 7 ADRBK1   immune_cell     1
# 8 ADRBK2   immune_cell     1
# 9 AGPAT6   immune_cell     1
# 10 AGPAT9   immune_cell     1



# filter(cell_type %>% is.na %>% `!`) %>%

# # Eliminate genes that are not in all cell types
# inner_join( (.) %>% distinct(symbol, cell_type) %>% count(symbol) %>% filter(n == max(n)) ) %>%
# 
# # Setup house keeping genes
# mutate(cell_type = ifelse(symbol %in% (read_csv("dev/database/ARMET/dev/hk_600.txt", col_names = FALSE) %>% pull(1)), "house_keeping", cell_type)) %>%
# 
# Select just needed columns
#select(sample, `Cell type`, cell_type_formatted, cell_type, level, count, symbol, ensembl_gene_id, `Data base`,  `Sample name`) %>%

# # Median redundant
# group_by(level, `Data base`) %>%
# multidplyr::partition(cluster) %>%
# do(
#   ttBulk::aggregate_duplicates((.), .sample = sample, .transcript = symbol, .abundance = count)
# ) %>%
# collect() %>%
# ungroup %>%

# do_parallel_start(n_cores, "symbol") %>%
# do({
#   `%>%` = magrittr::`%>%`
#   library(tidyverse)
#   library(magrittr)
#
#   (.) %>%
#     group_by(sample, symbol, `Cell type`, cell_type, `Cell type formatted`,  `Data base`) %>%
#     summarise(`count` = `count` %>% median(na.rm = T)) %>%
#     ungroup()
# }) %>%
# do_parallel_end() %>%

# # Normalise
# group_by(level) %>%
#   multidplyr::partition(cluster) %>%
#   do(
#     ttBulk::scale_abundance((.), sample, symbol, `count`)
#   ) %>%
#   collect() %>%
#   ungroup %>%
#   
#   mutate(`count scaled` = `count scaled` %>% as.integer) %>%
#   mutate(`count scaled log` = `count scaled` %>% `+` (1) %>% log) %>%
#   
#   # mutate symbol
#   mutate(`symbol original` = symbol) %>%
#   unite(symbol, c("Cell type category", "symbol"), remove = F) %>%
#   
#   # Remove redundant samples
#   group_by(level) %>%
#   do({
#     
#     nr =
#       (.) %>%
#       filter(cell_type != "house_keeping") %>%
#       group_by( cell_type) %>%
#       do({
#         
#         threshold = (.) %>% distinct(cell_type) %>% left_join( ct_to_correlation_threshold ) %>% pull(threshold)
#         
#         
#         # Remove redundant samples
#         (.) %>%
#           anti_join(
#             (.) %>%
#               distinct(symbol, sample, `count`) %>%
#               ttBulk::filter_variable(sample, symbol, count, top = 300) %>%
#               spread(sample, `count`) %>%
#               drop_na %>%
#               gather(sample, `count`, -symbol) %>%
#               rename(rc = `count`) %>%
#               mutate_if(is.factor, as.character) %>%
#               widyr::pairwise_cor(sample, symbol, rc, sort=T, diag = FALSE, upper = F) %>%
#               filter(correlation > threshold) %>%
#               distinct(item1) %>%
#               rename(sample = item1)
#           ) %>%
#           
#           # Sample homogeneous population
#           mutate( `threshold contribution` = (.) %>%  distinct(sample, `Cell type formatted`) %>% count(`Cell type formatted`) %>% pull(n) %>% quantile(0.8) %>% floor) %>%
#           group_by(`Cell type formatted`) %>%
#           inner_join( (.) %>% distinct(sample, `threshold contribution`) %>%  filter(row_number(sample) <= `threshold contribution`)) %>%
#           ungroup() %>%
#           select(- `threshold contribution`)
#         
#       }) %>%
#       ungroup
#     
#     (.) %>%
#       filter(cell_type == "house_keeping") %>%
#       inner_join( nr %>% distinct(sample)) %>%
#       bind_rows( nr )
#     
#   } )%>%
#   ungroup %>%
#   
#   mutate(symbol = symbol %>% as.factor) %>%
#   mutate(sample = sample %>% as.factor) %>%
#   saveRDS(file="dev/database/ARMET/dev/counts_infer_NB.rds")


# # Plots and Study
# library(ttBulk)
# (
#   readRDS(file="counts_infer_NB.rds") %>%
#     dplyr::distinct(sample, symbol, `count`, `Cell type formatted`, `Data base`) %>%
#     ttBulk::ttBulk(sample, symbol, `count`) %>%
#     distinct() %>%
#     aggregate_duplicates() %>%
#     ttBulk::scale_abundance() %>%
#     distinct(`Data base`, sample, symbol, `Cell type formatted`, `count scaled`)  %>%
#     reduce_dimensions(sample, symbol, `count scaled`, method = "tSNE", .dims=2) %>%
#     select(contains("tSNE"), `Data base`, `Cell type formatted`) %>%
#     distinct %>%
#     ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=`Cell type formatted`)) + geom_point(size=2) +
#     theme_bw() +
#     theme(
#       panel.border = element_blank(),
#       axis.line = element_line(),
#       panel.grid.major = element_line(size = 0.2),
#       panel.grid.minor = element_line(size = 0.1),
#       text = element_text(size=12),
#       legend.position="bottom",
#       aspect.ratio=1,
#       #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       strip.background = element_blank(),
#       axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
#       axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
#     )
# ) %>% plotly::ggplotly()
#
# (
#   ARMET::ARMET_ref %>%
#     filter(level==3) %>%
#     distinct(`Data base`, sample, symbol, cell_type, `read count scaled bayes`)  %>%
#      inner_join(rr$signatures[[3]] %>% distinct(symbol)) %>%
#     reduce_dimensions(sample, symbol, `read count scaled bayes`, method = "tSNE", .dims=2) %>%
#     select(contains("tSNE"), `Data base`, cell_type) %>%
#     distinct %>%
#     ggplot(aes(x = `tSNE 1`, y = `tSNE 2`, color=cell_type, shape=`Data base`)) + geom_point(size=2) +
#     theme_bw() +
#     theme(
#       panel.border = element_blank(),
#       axis.line = element_line(),
#       panel.grid.major = element_line(size = 0.2),
#       panel.grid.minor = element_line(size = 0.1),
#       text = element_text(size=12),
#       legend.position="bottom",
#       aspect.ratio=1,
#       #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       strip.background = element_blank(),
#       axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
#       axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
#     )
# ) %>% plotly::ggplotly()


#
# (
#   ARMET::ARMET_ref %>%
#     distinct(`Data base`, sample, symbol, `Cell type formatted`, `read count`)  %>%
#     ttBulk(sample, symbol, `read count`) %>%
#     aggregate_duplicates() %>%
#     scale_abundance() %>%
#     #inner_join(rr$signatures[[3]] %>% distinct(symbol)) %>%
#     reduce_dimensions(sample, symbol, `read count scaled`, method = "tSNE", .dims=2) %>%
#     select(contains("tSNE"), `Data base`, `Cell type formatted`) %>%
#     distinct %>%
#     ggplot(aes(x = `tSNE 1`, y = `tSNE 2`, color=`Cell type formatted`, shape=`Data base`)) + geom_point(size=2) +
#     theme_bw() +
#     theme(
#       panel.border = element_blank(),
#       axis.line = element_line(),
#       panel.grid.major = element_line(size = 0.2),
#       panel.grid.minor = element_line(size = 0.1),
#       text = element_text(size=12),
#       legend.position="bottom",
#       aspect.ratio=1,
#       #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       strip.background = element_blank(),
#       axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
#       axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
#     )
# ) %>% plotly::ggplotly()
#
# counts_annot =
#   counts %>%
#   filter(`read count` %>% is.na %>% `!`) %>%
#   filter(symbol %>% is.na %>% `!`) %>%
#   left_join(
#     ARMET::tree %>%
#       ARMET::ToDataFrameTypeColFull(fill = F) %>%
#       rename(`Cell type formatted` = level_5)
#   ) %>%
#   ttBulk(sample, symbol, `read count`) %>%
#   aggregate_duplicates() %>%
#   scale_abundance() %>%
#   reduce_dimensions(method = "MDS" ) %>%
#   reduce_dimensions(method = "tSNE")



# raw_df = readRDS("../../ARMET_dev/dev/counts_infer_NB.rds")
# 
# sample_ct = 
#   raw_df %>%
#   distinct(sample, cell_type, level) %>%
#   filter(cell_type != "house_keeping")
# 
# raw_df_ct = 
#   raw_df %>%
#   mutate(house_keeping = cell_type == "house_keeping" ) %>%
#   select(-cell_type) %>%
#   left_join(sample_ct) %>%
#   mutate(symbol = gsub("house_keeping_", "", symbol))
# 
# counts = 
#   raw_df_ct %>%
#   select(
#     database = `Data base`,
#     sample, 
#     level, 
#     cell_type = cell_type, 
#     cell_type_formatted =  `Cell type formatted`,
#     cell_tyoe_original = `Cell type`,
#     symbol,
#     count,
#     house_keeping
#   ) %>%
#   mutate_if(is.character, as.factor)
# 
# 
# save(counts, file="data/counts.rda", compress = "xz")



