library(tidyverse)
library(cellsig)
library(tidybulk)

counts_first_db_raw = readRDS("dev/counts_first_db_raw.rds")
load("dev/counts_second_db_raw.rda")
counts_second_db_raw = new_data %>% select(sample, symbol, count, dataset, cell_type)

#' @export
ToDataFrameTypeColFull = function(tree, fill = T, ...) {
  t = tree %>% data.tree::Clone()
  
  tree_df = 
    1:(t %$% Get("level") %>% max) %>%
    map_dfr(
      ~ data.tree::Clone(t) %>%
        {
          data.tree::Prune(., function(x)
            x$level <= .x +1)
          .
        } %>%
        data.tree::ToDataFrameTypeCol() %>%
        as_tibble
      
    ) %>%
    distinct() 
  
  tree_df_filled = 
    tree_df %>%
    
    purrr::when(
      1 & ("level_2" %in% colnames(.)) ~ mutate(., level_2 = ifelse(level_2 %>% is.na, level_1, level_2)),
      TRUE ~ (.)
    ) %>%
    purrr::when(
      1 & ("level_3" %in% colnames(.)) ~ mutate(., level_3 = ifelse(level_3 %>% is.na, level_2, level_3)),
      TRUE ~ (.)
    ) %>%
    purrr::when(
      1 & ("level_4" %in% colnames(.)) ~ mutate(., level_4 = ifelse(level_4 %>% is.na, level_3, level_4)),
      TRUE ~ (.)
    ) %>%
    purrr::when(
      1 & ("level_5" %in% colnames(.)) ~ mutate(., level_5 = ifelse(level_5 %>% is.na, level_4, level_5)),
      TRUE ~ (.)
    ) %>%
    purrr::when(
      1 & ("level_6" %in% colnames(.)) ~ mutate(., level_6 = ifelse(level_6 %>% is.na, level_5, level_6)),
      TRUE ~ (.)
    ) %>%
    dplyr::select(..., everything())
  
  tree_df %>%
    select(-1) %>%
    setNames(tree_df %>% colnames %>% .[-ncol(tree_df)]) %>%
    mutate(cell_type = tree_df_filled %>% pull(ncol(tree_df)))
  
}


counts = 
  counts_first_db_raw %>%
  select(-cell_type_original) %>%
  bind_rows(counts_second_db_raw %>% rename(database = dataset)) %>%
  
  # Add tree info
  left_join(
    tree %>%
      data.tree::Clone() %>%
      ToDataFrameTypeColFull(fill=NA) 
  ) %>%
  filter(level_1 %>% is.na %>% `!`) %>%
  
  # Reduce size
  mutate_if(is.character, as.factor) %>% 
  droplevels %>% 
  mutate(count = count %>% as.integer) %>%
  
  # Filter only symbol existing
  filter(symbol %>% is.na %>% `!`) %>%

  # Aggregate
  aggregate_duplicates(sample, symbol, count) %>%
  
  select(-ensembl_gene_id, -`merged transcripts`)  %>%

  # Infer exposure rate  
  infer_sequencing_depth_bias()



save(counts, file="dev/counts.rda", compress = "xz")

counts_imputed =
  counts %>%
  mutate(count_scaled = count / exp(exposure_rate)) %>%
  impute_abundance_using_levels(count_scaled)

save(counts_imputed, file="dev/counts_imputed.rda", compress = "xz")

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


# Plots and Study
library(ttBulk)
(
  readRDS(file="counts_infer_NB.rds") %>%
    dplyr::distinct(sample, symbol, `count`, `Cell type formatted`, `Data base`) %>%
    ttBulk::ttBulk(sample, symbol, `count`) %>%
    distinct() %>%
    aggregate_duplicates() %>%
    ttBulk::scale_abundance() %>%
    distinct(`Data base`, sample, symbol, `Cell type formatted`, `count scaled`)  %>%
    reduce_dimensions(sample, symbol, `count scaled`, method = "tSNE", .dims=2) %>%
    select(contains("tSNE"), `Data base`, `Cell type formatted`) %>%
    distinct %>%
    ggplot(aes(x = `tSNE1`, y = `tSNE2`, color=`Cell type formatted`)) + geom_point(size=2) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(),
      panel.grid.major = element_line(size = 0.2),
      panel.grid.minor = element_line(size = 0.1),
      text = element_text(size=12),
      legend.position="bottom",
      aspect.ratio=1,
      #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_blank(),
      axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
      axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
    )
) %>% plotly::ggplotly()
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



