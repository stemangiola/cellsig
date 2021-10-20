# devtools::install_github("stemangiola/nanny@convert-to-S3", force = TRUE)
# devtools::install_github("stemangiola/tidybulk@dev", force = TRUE)
# devtools::install_github("stemangiola/tidybulk@for-jian", force = TRUE)
# devtools::install_github("stemangiola/cellsig@dev", force = TRUE)

library(yaml)
library(tidytext)
library(data.tree)
library(tidytree)
library(ape)
library(glue)
library(rlang)
library(factoextra)
library(stringr)
library(scales)
library(KernSmooth)
library(splus2R)
library(data.tree)
library(cluster)
library(tidyverse)
library(tidybulk)
library(cellsig)
library(patchwork)
library(tidySummarizedExperiment)


# # OLD Functions for data of old format===============================================================================
# ## 1 preprocess data
# 
# 
# preprocess <- function(.data) {
#   # Load dataset
#   .data %>%
#     
#     # Remove entries with NA in gene symbol
#     filter(symbol %>% is.na %>% `!`) %>% 
#     
#     tidybulk(sample, symbol, count) %>%
#     
#     # Group by level because otherwise samples are duplicated
#     nest(data = -level) %>%
#     # filter(level ==3) %>%
#     
#     # Redefine factors inside each level
#     mutate(data = map(data, ~ droplevels(.x))) %>%
#     
#     # Remove redundancy
#     mutate(data = map(data, ~aggregate_duplicates(.x))) %>% 
#     
#     # Fill missing data. There are many genes that are not shared by the majority of samples
#     mutate(data = map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
#     
#     
#     # Scale for future PCA plotting
#     mutate(data = map(
#       data, ~ .x %>% 
#         identify_abundant(factor_of_interest = cell_type) %>%
#         scale_abundance() 
#     ))
#   
# }
# 
# 
# 
# ## 2 contrast functions
# 
# ### 2.1 pairwise comparisons
# get_contrasts_from_df = function(.data){
#   .data %>% 
#     
#     distinct(cell_type) %>% 
#     
#     # Permute
#     mutate(cell_type2 = cell_type) %>% 
#     expand(cell_type, cell_type2) %>% 
#     filter(cell_type != cell_type2) %>% 
#     
#     # Create contrasts
#     mutate(contrast = sprintf("cell_type%s - cell_type%s", cell_type, cell_type2)) %>%
#     pull(contrast)
#   
# }
# 
# ### 2.2 create a contrast vector for limma::makeContrasts() or tidybulk::test_differential abundance()
# 
# make_contrasts <- function(tt, LEVEL){
#   
#   prefix = "cell_type"
#   
#   # find all cell types
#   cell_types <- tt %>% 
#     filter(level == LEVEL) %>% 
#     unnest(data) %>% 
#     distinct(cell_type) %>% 
#     pull() 
#   
#   # format cell_types with prefix
#   cell_types <- paste(prefix, cell_types, sep="")
#   
#   # initialise a vector called contrasts
#   contrasts <- 1: length(cell_types)
#   
#   # create all contrasts and store them in contrasts
#   for(i in 1: length(cell_types) ){
#     background = paste(cell_types[-i], collapse = "+")
#     divisor = length(cell_types[-i])
#     contrasts[i] <- sprintf("%s-(%s)/%s", cell_types[i], background, divisor)
#   }
#   
#   return(contrasts)
# }
# 
# 
# ### 2.3 contrast with no average background
# 
# make_contrasts2 <- function(tt, LEVEL){
#   
#   prefix = "cell_type"
#   
#   # find all cell types
#   cell_types <- tt %>% 
#     filter(level == LEVEL) %>% 
#     unnest(data) %>% 
#     distinct(cell_type) %>% 
#     pull() 
#   
#   # format cell_types with prefix
#   cell_types <- paste(prefix, cell_types, sep="")
#   
#   # initialise a vector called contrasts
#   contrasts <- 1: length(cell_types)
#   
#   # create all contrasts and store them in contrasts
#   for(i in 1: length(cell_types) ){
#     background = paste(cell_types[-i], collapse = "+")
#     contrasts[i] <- sprintf("%s-(%s)", cell_types[i], background)
#   }
#   
#   return(contrasts)
# }
# 
# ## 3 marker selection
# 
# select_markers_for_each_contrast = function(.data, sig_size){
#   .data %>%
#     
#     # Group by contrast. Comparisons both ways.
#     pivot_longer(
#       cols = contains("___"),
#       names_to = c("stats", "contrast"), 
#       values_to = ".value", 
#       names_sep="___"
#     ) %>% 
#     
#     # Markers selection
#     nest(stat_df = -contrast) %>%
#     
#     # Reshape inside each contrast
#     mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>%
#     
#     # Rank
#     mutate(stat_df = map(stat_df, ~.x %>%
#                            filter(FDR < 0.05 & logFC > 2) %>%
#                            filter(logCPM > mean(logCPM)) %>%
#                            arrange(logFC %>% desc()) %>%
#                            slice(1:sig_size)
#                          
#     )) %>%
#     
#     unnest(stat_df)
# }
# 
# 
# ## 4
# 
# contrast <- function(tt, LEVEL, sig_size){
#   tt %>%
#     
#     # Investigate one level
#     filter(level == LEVEL) %>%
#     
#     # Differential transcription
#     mutate(markers = map(
#       data,
#       ~ test_differential_abundance(.x,
#                                     ~ 0 + cell_type, 
#                                     .contrasts = get_contrasts_from_df(.x),
#                                     action="only") 
#     )) %>%
#     
#     # Select rank from each contrast
#     mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, sig_size))) %>%
#     
#     # Add original data info to markers
#     mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
#     select(markers) %>%
#     unnest(markers) %>%
#     
#     # make contrasts pretty
#     mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
# }
# 
# 
# ## 5 calculate the area of confidence ellipses and the sum of their areas
# PCA_level5
# 
# ellip_area <- function(all_contrasts){
#   # reduce dimension
#   PCA <- all_contrasts %>% 
#     distinct(sample, symbol, count_scaled, cell_type) %>%
#     reduce_dimensions(sample, symbol, count_scaled,  method = "PCA", action="add", transform = log1p)
#   
#   real_size <- PCA %>% 
#     nest(data = -cell_type) %>% 
#     mutate(real_size = map_int(data, ~ n_distinct(.x$symbol)))
#   
#   area <- PCA %>%   
#     # remove non-numerical data to form a numerical data frame
#     select(cell_type, PC1, PC2) %>%
#     
#     # normalize principle component values
#     mutate(across(c("PC1", "PC2"), ~ .x %>% scale())) %>% 
#     
#     # nest by cell_type so as to calculate ellipse area for each cell type
#     nest(PC = -cell_type) %>% 
#     
#     # obtain covariance matrix for each cell type
#     mutate(cov = map(PC, ~ cov(.x))) %>% 
#     
#     # calculate the eigenvalues for the covariance matrix of each cell type
#     mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% 
#     
#     # transformation
#     mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>% # unnest(area)
#     
#     # below is the actual area for each ellipse
#     mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% 
#     
#     # collect size of each cluster as factors for weights
#     mutate(cluster_size = map_int(PC, ~ nrow(.x))) %>%
#     
#     # weight each area by the inverse of its cluster size
#     mutate(weighted_area = map2_dbl(area, cluster_size, ~ .x / .y))
#   
#   left_join(area, real_size)
#   
# }
# 
# 
# 
# 
# 
# 
# # NEW Functions for new data format ======================================================
# 
# # Functions for Hierarchical Analysis
# 
# ## 1 preprocess data
# 
# ### 1.1 string manipulation that converts level of interest (e.g "level_5") to its ancestor level (e.g "level_4")
# pre <- function(.level) {
#   .level %>% 
#     str_split("_") %>% 
#     {as.numeric(.[[1]][2])-1} %>% 
#     paste("level", ., sep = "_")
# }
# 
# ### 1.2 preprocess
# 
# preprocess <- function(.data, .level) {
#   
#   # load data
#   .data %>%
#     
#     tidybulk(sample, symbol, count) %>%
#     
#     # filter for the cell types of interest for gene marker selection
#     filter(is.na(!!as.symbol(.level))==F) %>%
#     
#     # Imputation of missing data within each level_5
#     # impute_missing_abundance(~ !!as.symbol(.level)) %>%
#     
#     # Group by ancestor
#     nest(data = - !!as.symbol(pre(.level))) %>%
#     
#     # Eliminate genes that are present in some but all cell types
#     # (can be still present in a subset of samples from each cell-type)
#     mutate(data = map(
#       data,
#       ~ .x %>%
#         nest(data = -c(symbol, !!as.symbol(.level))) %>%
#         add_count(symbol) %>%
#         filter(n == max(n)) %>%
#         unnest(data)
#     )) %>%
#     
#     # Imputation of missing data within each level_5
#     mutate(data = map(data, ~ .x %>% impute_missing_abundance(~ !!as.symbol(.level)))) %>%
#     
#     # scale count for further analysis
#     mutate(data=map(data, ~ .x %>%
#                       identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
#                       scale_abundance()
#     ))
# }
# 
# ## 2 contrast functions
# 
# ### 2.1 pairwise comparisons
# pairwise_contrast = function(.data, .level){
#   
#   .data %>% 
#     distinct(!!as.symbol(.level)) %>% 
#     mutate(!!as.symbol(.level) := paste0(.level, !!as.symbol(.level))) %>% 
#     
#     # Permute
#     mutate(cell_type2 := !!as.symbol(.level)) %>% 
#     expand(!!as.symbol(.level), cell_type2) %>% 
#     filter(!!as.symbol(.level) != cell_type2) %>% 
#     
#     # Create contrasts
#     mutate(contrast = sprintf("%s - %s", !!as.symbol(.level), cell_type2)) %>%
#     pull(contrast)
#   
# }
# 
# ### 2.2 create a contrast vector for limma::makeContrasts() or tidybulk::test_differential abundance()
# 
# mean_contrast <- function(.data, .level){
#   
#   # find all cell types
#   cell_types <- .data %>% 
#     distinct(!!as.symbol(.level)) %>% 
#     pull() %>% 
#     as.vector()
#   
#   # format cell_types with prefix
#   cell_types <- paste0(.level, cell_types)
#   
#   # initialise a vector called contrasts
#   contrasts <- 1: length(cell_types)
#   
#   # create all contrasts and store them in contrasts
#   for(i in 1:length(cell_types) ){
#     background = paste(cell_types[-i], collapse = "+")
#     divisor = length(cell_types[-i])
#     contrasts[i] <- sprintf("%s - (%s)/%s", cell_types[i], background, divisor)
#   }
#   
#   return(contrasts)
# }
# 
# ### 2.3 contrast with no average background (OBSOLETE)
# 
# mean_contrast2 <- function(.data, .level){
#   
#   # find all cell types
#   cell_types <- .data %>% 
#     distinct(!!as.symbol(.level)) %>% 
#     pull() %>% 
#     as.vector()
#   
#   # format cell_types with prefix
#   cell_types <- paste0(.level, cell_types)
#   
#   # initialise a vector called contrasts
#   contrasts <- 1: length(cell_types)
#   
#   # create all contrasts and store them in contrasts
#   for(i in 1: length(cell_types) ){
#     background = paste(cell_types[-i], collapse = "+")
#     contrasts[i] <- sprintf("%s-(%s)", cell_types[i], background)
#   }
#   
#   return(contrasts)
# }
# 
# ## 3 marker ranking & selection
# 
# select_markers_for_each_contrast = function(.markers, .sig_size){
#   .markers %>%
#     
#     # Group by contrast. Comparisons both ways.
#     pivot_longer(
#       cols = contains("___"),
#       names_to = c("stats", "contrast"), 
#       values_to = ".value", 
#       names_sep="___"
#     ) %>% 
#     
#     # Markers selection within each pair of contrast
#     nest(stat_df = -contrast) %>%
#     
#     # Reshape inside each contrast
#     mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>%
#     
#     # Rank
#     mutate(stat_df = map(stat_df, ~.x %>%
#                            filter(FDR < 0.05 & logFC > 2) %>%
#                            filter(logCPM > mean(logCPM)) %>%
#                            arrange(logFC %>% desc()) %>%
#                            slice(1: .sig_size)
#                          
#     )) %>%
#     
#     unnest(stat_df)
# }
# 
# ## 4 marker collection for each contrast
# 
# ### pairwise contrast method
# contrast_PW <- function(.tt, .level){
#   .tt %>%
#     
#     # Differential transcription
#     mutate(markers = map(
#       data,
#       ~ .x %>% 
#         test_differential_abundance(
#           ~ 0 + !!as.symbol(.level), 
#           .contrasts = pairwise_contrast(.x, .level),
#           action="only") 
#     ))
# }
# 
# ### mean contrast method
# contrast_MC <- function(.tt, .level){
#   .tt %>%
#     
#     # Differential transcription
#     mutate(markers = map(
#       data,
#       ~ test_differential_abundance(.x,
#                                     ~ 0 + !!as.symbol(.level), 
#                                     .contrasts = mean_contrast(.x, .level),
#                                     action="only") 
#     ))
# }
# 
# ### select signature genes and processing
# sig_select <- function(.contrast, .level, .sig_size) {
#   .contrast %>% 
#     
#     # Select markers from each contrast by rank of stats
#     mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, .sig_size))) %>%
#     
#     # Add original data info to the markers selected
#     mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
#     select(!!as.symbol(pre(.level)), markers) %>%
#     unnest(markers) %>%
#     
#     # make contrasts pretty
#     mutate(contrast_pretty = str_replace(contrast, .level, "") %>% str_replace(.level, ""))
# }
# 
# ## 5 calculate the area of confidence ellipses and the sum of their areas
# 
# ## 5.1 calculate the area of confidence ellipses and the sum of their areas
# ellipse <- function(.rdim, .level, .method) {
#   .rdim %>% 
#     
#     # remove non-numerical data to form a numerical data frame
#     select(!!as.symbol(.level), 
#            contains(str_sub(.method, end=-2L))) %>% 
#     
#     # normalize principle component values
#     mutate(across(contains(str_sub(.method, end=-2L)), scale)) %>% 
#     
#     # nest by cell_type so as to calculate ellipse area for each cell type
#     nest(dims = - !!as.symbol(.level)) %>% 
#     
#     # obtain covariance matrix for each cell type
#     mutate(cov = map(dims, ~ cov(.x))) %>% 
#     
#     # calculate the eigenvalues for the covariance matrix of each cell type
#     mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% 
#     
#     # transformation
#     mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>%
#     
#     # below is the actual area for each ellipse
#     mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% 
#     
#     # collect size of each cluster as factors for weights
#     mutate(cluster_size = map_int(dims, ~ nrow(.x))) %>%
#     
#     # weight each area by the inverse of its cluster size
#     mutate(weighted_area = map2_dbl(area, cluster_size, ~ .x / .y))
# }
# 
# ## 5.2 Ellipse area calculation 
# ellip_func <- function(.markers, .level, .method){
#   .markers %>% 
#     
#     # nest by ancestor cell types
#     nest(rdim = - !!as.symbol(pre(.level))) %>%
#     
#     # reduce dimension
#     mutate(rdim = map(rdim, ~ .x %>%
#                         distinct(sample, symbol, count_scaled, !!as.symbol(.level)))) %>%
#     mutate(rdim = map(rdim, ~ .x %>%
#                         reduce_dimensions(sample, symbol, count_scaled,
#                                           method = .method,
#                                           action = "add",
#                                           transform = log1p,
#                                           # check_duplicates is for Rtsne method
#                                           check_duplicates = FALSE) %>% 
#                         
#                         # save symbols for calculating real_size while reducing replicated rows resulted from symbol
#                         nest(data_symbol = c(symbol, count_scaled))
#     )) %>% 
#     
#     mutate(real_size = map_int(rdim, ~ .x$data_symbol %>% 
#                                  map_int(~ n_distinct(.x$symbol)) %>% 
#                                  unlist() %>% 
#                                  unique() )) %>% 
#     
#     mutate(area_df = map(rdim, ~ ellipse(.x, .level, .method) )
#     )
#   
# }
# 
# ## 5.3 Scale serialised ellip_func() output (a tibble called ellip_tb) for plotting
# ellip_scale <- function(.ellip_tb, .level) {
#   .ellip_tb %>% 
#     unnest(ellip) %>%
#     unnest(area_df) %>%
#     
#     # nest by ancestor cell type to rescale area for all sig_sizes
#     nest(cell_data = - !!as.symbol(pre(LEVEL))) %>%
#     mutate(cell_data = map(cell_data, ~ .x %>% 
#                              mutate(rescaled_area = area %>% 
#                                       scale(center = F))
#     )) %>%
#     
#     # nest by ancestor cell type to summarise areas for each real_size/sig_size
#     mutate(plot_data = map(cell_data, ~ .x %>%
#                              # sum all areas for each real_size for an ancestor node
#                              group_by(real_size) %>%
#                              summarise(sig_size,
#                                        stded_sum=sum(area, na.rm = T),
#                                        wted_sum = sum(weighted_area, na.rm = T),
#                                        rescaled_sum= sum(rescaled_area, na.rm = T)) %>%
#                              # remove duplicate rows
#                              distinct(real_size, sig_size, stded_sum, wted_sum, rescaled_sum) %>% 
#                              pivot_longer(ends_with("sum"), names_to='area_type', values_to="area_value")
#     ))
# }
# 
# ## 6 Silhoette function
# sil_func <- function(.markers, .level, .method){
#   .markers %>%
#     nest(reduced_dimensions = - !!as.symbol(pre(.level))) %>%
#     mutate(reduced_dimensions = map(reduced_dimensions, ~ .x %>%
#                                       distinct(sample, symbol, count_scaled, !!as.symbol(.level)))) %>%
#     
#     mutate(signature = map(reduced_dimensions, ~ .x %>% 
#                              pull(symbol) %>% 
#                              unique() )) %>% 
#     
#     mutate(real_size = map_int(signature, ~ length(.x))) %>% 
#     
#     mutate(reduced_dimensions = map(reduced_dimensions, ~ .x %>%
#                                       reduce_dimensions(sample, symbol, count_scaled,
#                                                         method = .method,
#                                                         transform = log1p,
#                                                         # check_duplicates is for Rtsne method
#                                                         check_duplicates = FALSE) %>% 
#                                       
#                                       # remove duplicates caused by symbol & counts to calculate the distance matrix
#                                       distinct(sample, !!as.symbol(.level), PC1, PC2)
#     )) %>%
#     
#     # calculate the dissimilarity matrix with PC values
#     mutate(distance = map(reduced_dimensions, ~ .x %>%
#                             select(contains(str_sub(.method, end = -2L))) %>%
#                             factoextra::get_dist(method = "euclidean")
#     )) %>%
#     
#     # calculate silhouette score
#     mutate(silhouette = map2(reduced_dimensions, distance,
#                              ~ silhouette(as.numeric(as.factor(`$`(.x, !!as.symbol(.level)))), .y)
#     )) %>%
#     mutate(silhouette = map(silhouette, ~ .x %>% summary())) %>%
#     mutate(silhouette = map(silhouette, ~ .x %>% 
#                               `$`(avg.width) ))%>% 
#     mutate(silhouette = unlist(silhouette))
#   
# }
# 
# ### for summary plot
# sil_tb <- function(.contrast, .level, .sig_size, .method) {
#   tibble(sig_size = 1: .sig_size) %>% 
#     
#     # select signature genes for each sig_size at each level
#     mutate(sil_df = map(sig_size, ~ sig_select(.contrast, .level, .x))) %>% 
#     
#     # calculate silhouette score for each ancestor cell type under sil_df
#     mutate(sil_df = map(sil_df, ~ sil_func(.x, .level, .method)))
# }
# 
# # signature size optimation functions:
# penalised_silhouette <- function(.plot_data, .penality_rate=0.2) {
#   .plot_data %>% 
#     
#     # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
#     # filter(real_size > 10) %>%
#     
#     # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
#     mutate(size_rescaled = rescale(real_size)) %>% 
#     
#     mutate(penalised_silhouette = silhouette - .penality_rate * size_rescaled) %>% 
#     
#     filter(penalised_silhouette==max(penalised_silhouette)) %>% 
#     
#     select(real_size, silhouette)
# }
# 
# ratio <- function(.plot_data) {
#   .plot_data %>% 
#     
#     # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
#     # filter(real_size > 10) %>%
#     
#     # calculate the difference between rescaled size of the max sil point with that of all other points
#     # the bigger different the better (even for negative numbers)
#     mutate(size_gap = real_size[which.max(silhouette)] - real_size) %>% 
#     
#     # calculate the difference between silhouette score of the max sil point and that of all other points
#     # the smaller the better
#     mutate(silhouette_gap = max(silhouette) - silhouette) %>% 
#     
#     # use min_max scaler to rescale size_difference to a scale (between 0 and 1)
#     mutate(size_gap_rescaled = rescale(size_gap)) %>%
#     
#     # use min_max scaler to rescale silhouette score difference to the same scale (between 0 and 1)
#     mutate(silhouette_gap_rescaled = rescale(silhouette_gap, c(1, 10))) %>% 
#     
#     # calculate the ratio between size_gap_rescaled/silhouette_gap, the bigger the better
#     mutate(ratio = map2_dbl(
#       silhouette_gap_rescaled, size_gap_rescaled, 
#       ~ if(.x == 0){0}else{.y / .x}
#     )) %>% 
#     
#     filter(ratio == max(ratio)) %>% 
#     
#     select(real_size, silhouette)
# }
# 
# # Functions for Non-hierarchical Analysis
# ## 2.1
# 
# mean_contrast0 <- function(.data){
#   
#   prefix <- "cell_type"
#   
#   # find all cell types
#   cell_types <- .data %>% 
#     distinct(cell_type) %>% 
#     pull() %>% 
#     as.vector()
#   
#   # format cell_types with prefix
#   cell_types <- paste0(prefix, cell_types)
#   
#   # initialise a vector called contrasts
#   contrasts <- 1: length(cell_types)
#   
#   # create all contrasts and store them in contrasts
#   for(i in 1:length(cell_types) ){
#     background = paste(cell_types[-i], collapse = "+")
#     divisor = length(cell_types[-i])
#     contrasts[i] <- sprintf("%s - (%s)/%s", cell_types[i], background, divisor)
#   }
#   
#   return(contrasts)
# }
# 
# get_contrasts_from_df0 = function(.data){
#   
#   .data %>% 
#     distinct(cell_type) %>% 
#     mutate(cell_type = paste0("cell_type",cell_type)) %>% 
#     
#     # Permute
#     mutate(cell_type2 = cell_type) %>% 
#     expand(cell_type, cell_type2) %>% 
#     filter(cell_type != cell_type2) %>% 
#     
#     # Create contrasts
#     mutate(contrast = sprintf("%s - %s", cell_type, cell_type2)) %>%
#     pull(contrast)
#   
# }
# 
# ## 2.3 marker collection for each contrast
# 
# # mean contrast method
# contrast_MC0 <- function(.tt){
#   .tt %>%
#     
#     # Differential transcription
#     mutate(markers = map(
#       data,
#       ~ test_differential_abundance(.x,
#                                     ~ 0 + cell_type, 
#                                     .contrasts = mean_contrast0(.x),
#                                     action="only")  ))
# }
# 
# # pairwise contrast method
# contrast_PW0 <- function(.tt){
#   .tt %>%
#     
#     # Differential transcription
#     mutate(markers = map(
#       data,
#       ~ test_differential_abundance(.x,
#                                     ~ 0 + cell_type, 
#                                     .contrasts = get_contrasts_from_df0(.x),
#                                     action="only")  ))
# }
# 
# # marker selection & processing
# sig_select0 <- function(.contrast, .sig_size) {
#   .contrast %>% 
#     
#     # Select markers from each contrast by rank of stats
#     mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, .sig_size))) %>%
#     
#     # Add original data info to the markers selected
#     mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
#     select(level_0, markers) %>%
#     unnest(markers) %>%
#     
#     # make contrasts pretty
#     mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
# }
# 
# ## 2.4 Ellipse area
# 
# ## 2.4.1 calculate the area of confidence ellipses and the sum of their areas
# ellipse0 <- function(.rdim, .method) {
#   .rdim %>% 
#     
#     # remove non-numerical data to form a numerical data frame
#     select(cell_type, 
#            contains(str_sub(.method, end=-2L))) %>% 
#     
#     # normalize principle component values
#     mutate(across(contains(str_sub(.method, end=-2L)), scale)) %>% 
#     
#     # nest by cell_type so as to calculate ellipse area for each cell type
#     nest(dims = - cell_type) %>% 
#     
#     # obtain covariance matrix for each cell type
#     mutate(cov = map(dims, ~ cov(.x))) %>% 
#     
#     # calculate the eigenvalues for the covariance matrix of each cell type
#     mutate(eigval = map(cov, ~ eigen(.x)$values)) %>% 
#     
#     # transformation
#     mutate(area = map(eigval, ~ sqrt(.x * qchisq(0.95, 2)))) %>%
#     
#     # below is the actual area for each ellipse
#     mutate(area = map_dbl(area, ~ prod(.x)*pi)) %>% 
#     
#     # collect size of each cluster as factors for weights
#     mutate(cluster_size = map_int(dims, ~ nrow(.x))) %>%
#     
#     # weight each area by the inverse of its cluster size
#     mutate(weighted_area = map2_dbl(area, cluster_size, ~ .x / .y))
# }
# 
# ## 2.4.2 Ellipse area calculation 
# ellip_func0 <- function(.markers, .method){
#   .markers %>% 
#     
#     # nest by ancestor cell types
#     nest(rdim = - level_0) %>%
#     
#     # reduce dimension
#     mutate(rdim = map(rdim, ~ .x %>%
#                         distinct(sample, symbol, count_scaled, cell_type))) %>%
#     mutate(rdim = map(rdim, ~ .x %>%
#                         reduce_dimensions(sample, symbol, count_scaled,
#                                           method = .method,
#                                           action = "add",
#                                           transform = log1p,
#                                           # check_duplicates is for Rtsne method
#                                           check_duplicates = FALSE) %>% 
#                         
#                         # save symbols for calculating real_size while reducing replicated rows resulted from symbol
#                         nest(data_symbol = c(symbol, count_scaled))
#     )) %>% 
#     
#     mutate(real_size = map_int(rdim, ~ .x$data_symbol %>% 
#                                  map_int(~ n_distinct(.x$symbol)) %>% 
#                                  unlist() %>% 
#                                  unique() )) %>% 
#     
#     mutate(area_df = map(rdim, ~ ellipse0(.x, .method) ))
#   
# }
# 
# ## 2.4.3 Scale serialised ellip_func() output (a tibble called ellip_tb) for plotting
# ellip_scale0 <- function(.ellip_tb) {
#   .ellip_tb %>% 
#     unnest(ellip) %>%
#     unnest(area_df) %>%
#     
#     # nest by ancestor cell type to rescale area for all sig_sizes
#     nest(cell_data = - level_0) %>%
#     mutate(cell_data = map(cell_data, ~ .x %>% 
#                              mutate(rescaled_area = area %>% 
#                                       scale(center = F))
#     )) %>%
#     
#     # nest by ancestor cell type to summarise areas for each real_size/sig_size
#     mutate(plot_data = map(cell_data, ~ .x %>%
#                              # sum all areas for each real_size for an ancestor node
#                              group_by(real_size) %>%
#                              summarise(sig_size,
#                                        stded_sum=sum(area, na.rm = T),
#                                        wted_sum = sum(weighted_area, na.rm = T),
#                                        rescaled_sum= sum(rescaled_area, na.rm = T)) %>%
#                              # remove duplicate rows
#                              distinct(real_size, sig_size, stded_sum, wted_sum, rescaled_sum) %>% 
#                              pivot_longer(ends_with("sum"), names_to='area_type', values_to="area_value")
#     ))
# }
# 
# ## 2.5 Silhouette function
# sil_func0 <- function(.markers, .method){
#   .markers %>% 
#     nest(rdim = - level_0) %>% 
#     mutate(rdim = map(rdim, ~ .x %>% 
#                         group_by(sample, symbol) %>% 
#                         summarise(count_scaled = mean(count_scaled), 
#                                   cell_type=unique(cell_type)) %>% 
#                         ungroup()
#     )) %>% 
#     mutate(signature = map(rdim, ~ .x %>% 
#                              pull(symbol) %>% 
#                              unique()
#     )) %>% 
#     mutate(real_size = map_int(signature, ~ length(.x))) %>% 
#     mutate(rdim = map(rdim, ~ .x %>%
#                         reduce_dimensions(sample, symbol, count_scaled,
#                                           method = .method,
#                                           transform = log1p,
#                                           # check_duplicates is for Rtsne method
#                                           check_duplicates = FALSE) %>% 
#                         
#                         # remove duplicates caused by symbol & counts to calculate the distance matrix
#                         distinct(sample, cell_type, PC1, PC2)
#     )) %>%
#     
#     # calculate the dissimilarity matrix with PC values
#     mutate(distance = map(rdim, ~ .x %>%
#                             select(contains(str_sub(.method, end = -2L))) %>%
#                             factoextra::get_dist(method = "euclidean")
#     )) %>%
#     
#     # calculate silhouette score
#     mutate(sil = map2(rdim, distance, 
#                       ~ silhouette(as.numeric(as.factor(`$`(.x, cell_type))), .y)
#     )) %>% 
#     mutate(sil = map(sil, ~ .x %>% summary())) %>%
#     mutate(sil = map(sil, ~ .x %>% `$`(avg.width) ))%>% 
#     mutate(sil = unlist(sil))
#   
# }
# 
# ### for summary plot
# sil_tb0 <- function(.contrast, .sig_size, .method) {
#   tibble(sig_size = 1: .sig_size) %>% 
#     
#     # select signature genes for each sig_size at each level
#     mutate(sil_df = map(sig_size, ~ sig_select0(.contrast, .x))) %>% 
#     
#     # calculate silhouette score for each ancestor cell type under sil_df
#     mutate(sil_df = map(sil_df, ~ sil_func0(.x, .method)))
# }
# 
# # Functions for Shiny App
# format_name <- function(.method) {
#   paste("markers", .method, sep = "_")
# }
# 
# cell_sig_select <- function(.markers) {
#   .markers %>% 
#     # obtain cell types in a node and nest by it to extract signatures for all cell types
#     mutate(cell_type = str_extract(contrast_pretty, "([a-z]|\\_)+(?=\\s)")) %>% 
#     nest(signature = - cell_type) %>% 
#     mutate(signature = map(signature, ~ .x %>% 
#                              pull(symbol) %>% 
#                              unique()
#     ))
# }
# 
# # NEW METHOD =======================================================
# # Hierarchical
# 
# # pairwise selection 1 marker at a time
# single_marker_pw_select <- function(.contrast, .level, .discard_num, .method) {
#   
#   # initialize variables
#   
#   contrast_copy <- .contrast %>% 
#     mutate(ancestor = !!as.symbol(pre(.level)))
#   
#   # initialise a signature tibble to store signature markers for each cell type in each iteration
#   signature <- tibble(
#     ancestor = contrast_copy %>% 
#       pull(ancestor)
#   ) %>% 
#     mutate(cumulative_signature = map(ancestor, ~ vector())) %>% 
#     mutate(last_silhouette = 0)
#   
#   # initialise an output tibble containing all results of interest
#   summary_tb <- tibble(
#     ancestor = character(),
#     new_challengers = list(),
#     winner = list(),
#     cumulative_signature = list(),
#     # reduced_dimensions = list(),
#     winning_silhouette = double()
#   )
#   
#   # set the base markers
#   contrast_pair_tb0 <-
#     
#     # contrast_copy contains all the statistics of all cell_type contrasts for each gene
#     contrast_copy %>%
#     
#     # select top 1 markers from each contrast, output is an unnested tibble
#     sig_select(.level, 1) %>%
#     
#     mutate(ancestor = !!as.symbol(pre(.level))) %>%
#     
#     nest(data = - ancestor) %>%
#     
#     mutate(data = map(data, ~ .x %>% 
#                         nest(new_challenger = - contrast_pretty) %>% 
#                         
#                         mutate(new_challenger = map(
#                           new_challenger, ~ .x %>% pull(symbol) %>% unique()))
#     )) %>% 
#     
#     mutate(new_challengers = map(data, ~.x %>% pull(new_challenger) %>% unlist())) %>% 
#     
#     mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
#     
#     mutate(cumulative_signature = winner) %>% 
#     
#     select(-data) %>% 
#     
#     mutate(silhouette = map2(
#       new_challengers, ancestor,
#       ~ silhouette_for_markers(.x, .y, .contrast, .level, .method) %>% 
#         select(
#           # reduced_dimensions, 
#           winning_silhouette = silhouette)
#     )) %>% 
#     
#     unnest(silhouette)
#   
#   
#   signature <- signature %>%
#     mutate(cumulative_signature = map2(
#       cumulative_signature, ancestor,
#       ~ .x %>% 
#         append(with(contrast_pair_tb0, cumulative_signature[ancestor==.y][[1]]))
#     )) %>% 
#     mutate(last_silhouette = map_dbl(
#       ancestor,
#       ~ with(contrast_pair_tb0, winning_silhouette[ancestor==.x])
#     ))
#   
#   summary_tb <- summary_tb %>% 
#     bind_rows(contrast_pair_tb0)
#   
#   # remove base markers from contrast_copy input before further selection
#   contrast_copy <- contrast_copy %>%
#     mutate(markers = map2(
#       markers, ancestor, 
#       ~ .x %>%
#         filter(!symbol %in% with(signature, cumulative_signature[ancestor==.y][[1]]))
#     ))
#   
#   # counter for number of genes discarded
#   j <- map_int(signature$cumulative_signature, ~ length(.x))
#   
#   # count the number of iterations
#   i <- 0
#   while (any(j < .discard_num) &
#          # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
#          all(map_int(contrast_copy$markers, 
#                      # hence the boundary should be the number of satisfactory genes selected
#                      ~ .x %>% select_markers_for_each_contrast(1) %>% nrow()) > 0)) {
#     
#     contrast_pair_tb <- 
#       
#       # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
#       contrast_copy %>% 
#       
#       # select top 1 markers from each contrast, output is an unnested tibble
#       sig_select(.level, 1) %>% 
#       
#       mutate(ancestor = !!as.symbol(pre(.level))) %>% 
#       
#       nest(data = - ancestor) %>% 
#       
#       mutate(data = map(data, ~ .x %>% 
#                           nest(new_challenger = - contrast_pretty) %>% 
#                           mutate(new_challenger = map_chr(
#                             new_challenger, 
#                             ~.x %>% pull(symbol) %>% unique()
#                           ))
#       )) %>% 
#       
#       unnest(data) %>% 
#       
#       mutate(challengers_for_silhouette = map2(
#         new_challenger, ancestor, 
#         ~ with(signature, cumulative_signature[ancestor==.y][[1]]) %>% 
#           append(.x)
#       )) %>% 
#       
#       mutate(silhouette = map2(
#         challengers_for_silhouette, ancestor, 
#         ~ silhouette_for_markers(.x, .y, .contrast, .level, .method) %>% 
#           select(-ancestor)
#       )) %>% 
#       
#       unnest(silhouette) %>% 
#       
#       group_by(ancestor) %>% 
#       
#       arrange(desc(silhouette), .by_group = TRUE) %>% 
#       
#       ungroup() %>% 
#       
#       mutate(is_greater = map2_lgl(
#         silhouette, ancestor, 
#         ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
#       )) %>% 
#       
#       nest(data = - ancestor) %>% 
#       
#       mutate(new_challengers = map(
#         data,
#         ~ .x %>% 
#           pull(new_challenger)
#       )) %>% 
#       
#       mutate(winner = map(
#         data,
#         ~ if(.x[1, ]$is_greater){
#           .x[1, ]$new_challenger
#         } else {NA}
#       )) %>% 
#       
#       mutate(cumulative_signature = pmap(
#         list(data, winner, ancestor),
#         ~ if(!is.na(..2)) {
#           with(..1[1, ], challengers_for_silhouette[[1]])
#         } else {
#           with(signature, cumulative_signature[ancestor==..3][[1]])
#         }
#       ))
#     
#     
#     # append the base + 1 markers that result in highest silhouette score
#     signature <- signature %>% 
#       
#       mutate(cumulative_signature = map(
#         ancestor,
#         ~ with(contrast_pair_tb, cumulative_signature[ancestor==.x][[1]])
#       )) %>% 
#       
#       mutate(last_silhouette = map2_dbl(
#         ancestor, last_silhouette,
#         ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x][[1]]))) {
#           with(contrast_pair_tb, data[ancestor==.x][[1]][[1, "silhouette"]])
#         } else {.y}
#       ))
#     
#     summary_tb <- summary_tb %>% 
#       bind_rows(
#         contrast_pair_tb %>% 
#           filter(!is.na(winner)) %>% 
#           # mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>% 
#           mutate(winning_silhouette = map_dbl(data, ~ .x$silhouette[1])) %>% 
#           select(-data)
#       )
#     
#     contrast_copy <- contrast_copy %>% 
#       mutate(markers = map2(
#         markers, ancestor, 
#         ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y][[1]]))){
#           .x %>% 
#             filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))
#         } else {
#           .x %>% 
#             filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]]))
#         }
#       ))
#     
#     # number of genes discarded for each node
#     j <- j + 
#       
#       # unsuccessful candidates
#       map_int(contrast_pair_tb$new_challengers, ~length(.x)) *
#       is.na(contrast_pair_tb$winner) +
#       
#       # winning candidates
#       map_int(contrast_pair_tb$winner, ~length(.x)) *
#       !is.na(contrast_pair_tb$winner)
#     
#     cat("genes discarded for each node: ", j, "\n")
#     cat("genes selected for each node: ", map_int(signature$cumulative_signature, ~ length(.x)),  "\n")
#     
#     i <- i + 1
#     cat("iteration: ", i, "\n")
#     
#   }
#   
#   output <- summary_tb %>% 
#     nest(signature_data = - ancestor) %>% 
#     mutate(level = .level, .before = ancestor)
#   
#   return(output)
# }
# 
# # calculates silhouette score for each set of signature (cumulative markers) at a signature size
# silhouette_for_markers <-function(.signature, .ancestor, .contrast, .level, .method) {
#   
#   .contrast %>%
#     
#     mutate(ancestor = !!as.symbol(pre(.level))) %>% 
#     
#     filter(ancestor == .ancestor) %>% 
#     
#     # filter markers that are in the signature
#     mutate(markers = map2(markers, ancestor, ~.x %>% 
#                             filter(symbol %in% .signature))) %>% 
#     
#     # format statistics from pairwise contrast
#     mutate(markers  = map(markers, 
#                           ~ .x %>% 
#                             # Group by contrast. Comparisons both ways.
#                             pivot_longer(
#                               cols = contains("___"),
#                               names_to = c("stats", "contrast"), 
#                               values_to = ".value", 
#                               names_sep="___"
#                             ) %>% 
#                             
#                             # Markers selection within each pair of contrast
#                             nest(stat_df = -contrast) %>%
#                             
#                             # Reshape inside each contrast
#                             mutate(stat_df = map(stat_df, ~.x %>% 
#                                                    pivot_wider(names_from = stats, 
#                                                                values_from = .value))) %>% 
#                             
#                             unnest(stat_df) )) %>% 
#     
#     # Add original data data to the markers selected
#     mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>% 
#     
#     # select only columns needed
#     select(-data) %>% 
#     
#     nest(silhouette_data = c(!!as.symbol(pre(.level)), markers)) %>% 
#     
#     mutate(silhouette_data = map(silhouette_data, ~ .x %>% unnest(markers))) %>% 
#     
#     # Calculate silhouette score for PCA plot resulted from the markers selected
#     mutate(silhouette_data = map(silhouette_data, ~ sil_func(.x, .level, .method))) %>% 
#     
#     mutate(reduced_dimensions = map(silhouette_data, ~ .x$rdim[[1]])) %>% 
#     mutate(silhouette = map_dbl(silhouette_data, ~ .x$sil)) %>% 
#     select(-silhouette_data)
#   
# }
# 
# # Non-hierarchical
# 
# single_marker_pw_select0 <- function(.contrast, .discard_num, .method) {
#   
#   # initialize variables
#   
#   contrast_copy <- .contrast %>% 
#     mutate(ancestor = level_0)
#   
#   # initialise a signature tibble to store signature markers for each cell type in each iteration
#   signature <- tibble(
#     ancestor = contrast_copy %>% 
#       pull(ancestor)
#   ) %>% 
#     mutate(cumulative_signature = map(ancestor, ~ vector())) %>% 
#     mutate(last_silhouette = 0)
#   
#   # initialise an output tibble containing all results of interest
#   summary_tb <- tibble(
#     ancestor = character(),
#     new_challengers = list(),
#     winner = list(),
#     cumulative_signature = list(),
#     # reduced_dimensions = list(),
#     winning_silhouette = double()
#   )
#   
#   # set the base markers
#   contrast_pair_tb0 <-
#     
#     # contrast_copy contains all the statistics of all cell_type contrasts for each gene
#     contrast_copy %>%
#     
#     # select top 1 markers from each contrast, output is an unnested tibble
#     sig_select0(1) %>%
#     
#     mutate(ancestor = level_0) %>%
#     
#     nest(data = - ancestor) %>%
#     
#     mutate(data = map(data, ~ .x %>% 
#                         nest(new_challenger = - contrast_pretty) %>% 
#                         
#                         mutate(new_challenger = map(
#                           new_challenger, ~ .x %>% pull(symbol) %>% unique()))
#     )) %>% 
#     
#     mutate(new_challengers = map(data, ~.x %>% pull(new_challenger) %>% unlist())) %>% 
#     
#     mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
#     
#     mutate(cumulative_signature = winner) %>% 
#     
#     select(-data) %>% 
#     
#     mutate(silhouette = map2(
#       new_challengers, ancestor,
#       ~ silhouette_for_markers0(.x, .y, .contrast, .method) %>% 
#         select(
#           # reduced_dimensions, 
#           winning_silhouette = silhouette)
#     )) %>% 
#     
#     unnest(silhouette)
#   
#   
#   signature <- signature %>%
#     mutate(cumulative_signature = map2(
#       cumulative_signature, ancestor,
#       ~ .x %>% 
#         append(with(contrast_pair_tb0, cumulative_signature[ancestor==.y][[1]]))
#     )) %>% 
#     mutate(last_silhouette = map_dbl(
#       ancestor,
#       ~ with(contrast_pair_tb0, winning_silhouette[ancestor==.x])
#     ))
#   
#   summary_tb <- summary_tb %>% 
#     bind_rows(contrast_pair_tb0)
#   
#   # remove base markers from contrast_copy input before further selection
#   contrast_copy <- contrast_copy %>%
#     mutate(markers = map2(
#       markers, ancestor, 
#       ~ .x %>%
#         filter(!symbol %in% with(signature, cumulative_signature[ancestor==.y][[1]]))
#     ))
#   
#   # counter for number of genes discarded
#   j <- map_int(signature$cumulative_signature, ~ length(.x))
#   
#   # count the number of iterations
#   i <- 0
#   while (any(j < .discard_num) &
#          # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
#          all(map_int(contrast_copy$markers, 
#                      # hence the boundary should be the number of satisfactory genes selected
#                      ~ .x %>% select_markers_for_each_contrast(1) %>% nrow()) > 0)) {
#     
#     contrast_pair_tb <- 
#       
#       # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
#       contrast_copy %>% 
#       
#       # select top 1 markers from each contrast, output is an unnested tibble
#       sig_select0(1) %>% 
#       
#       mutate(ancestor = level_0) %>% 
#       
#       nest(data = - ancestor) %>% 
#       
#       mutate(data = map(data, ~ .x %>% 
#                           nest(new_challenger = - contrast_pretty) %>% 
#                           mutate(new_challenger = map_chr(
#                             new_challenger, 
#                             ~.x %>% pull(symbol) %>% unique()
#                           ))
#       )) %>% 
#       
#       unnest(data) %>% 
#       
#       mutate(challengers_for_silhouette = map2(
#         new_challenger, ancestor, 
#         ~ with(signature, cumulative_signature[ancestor==.y][[1]]) %>% 
#           append(.x)
#       )) %>% 
#       
#       mutate(silhouette = map2(
#         challengers_for_silhouette, ancestor, 
#         ~ silhouette_for_markers0(.x, .y, .contrast, .method) %>% 
#           select(-ancestor)
#       )) %>% 
#       
#       unnest(silhouette) %>% 
#       
#       group_by(ancestor) %>% 
#       
#       arrange(desc(silhouette), .by_group = TRUE) %>% 
#       
#       ungroup() %>% 
#       
#       mutate(is_greater = map2_lgl(
#         silhouette, ancestor, 
#         ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
#       )) %>% 
#       
#       nest(data = - ancestor) %>% 
#       
#       mutate(new_challengers = map(
#         data,
#         ~ .x %>% 
#           pull(new_challenger)
#       )) %>% 
#       
#       mutate(winner = map(
#         data,
#         ~ if(.x[1, ]$is_greater){
#           .x[1, ]$new_challenger
#         } else {NA}
#       )) %>% 
#       
#       mutate(cumulative_signature = pmap(
#         list(data, winner, ancestor),
#         ~ if(!is.na(..2)) {
#           with(..1[1, ], challengers_for_silhouette[[1]])
#         } else {
#           with(signature, cumulative_signature[ancestor==..3][[1]])
#         }
#       ))
#     
#     
#     # append the base + 1 markers that result in highest silhouette score
#     signature <- signature %>% 
#       
#       mutate(cumulative_signature = map(
#         ancestor,
#         ~ with(contrast_pair_tb, cumulative_signature[ancestor==.x][[1]])
#       )) %>% 
#       
#       mutate(last_silhouette = map2_dbl(
#         ancestor, last_silhouette,
#         ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x][[1]]))) {
#           with(contrast_pair_tb, data[ancestor==.x][[1]][[1, "silhouette"]])
#         } else {.y}
#       ))
#     
#     summary_tb <- summary_tb %>% 
#       bind_rows(
#         contrast_pair_tb %>% 
#           filter(!is.na(winner)) %>% 
#           # mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>% 
#           mutate(winning_silhouette = map_dbl(data, ~ .x$silhouette[1])) %>% 
#           select(-data)
#       )
#     
#     contrast_copy <- contrast_copy %>% 
#       mutate(markers = map2(
#         markers, ancestor, 
#         ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y][[1]]))){
#           .x %>% 
#             filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))
#         } else {
#           .x %>% 
#             filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]]))
#         }
#       ))
#     
#     # number of genes discarded for each node
#     j <- j + 
#       
#       # unsuccessful candidates
#       map_int(contrast_pair_tb$new_challengers, ~length(.x)) *
#       is.na(contrast_pair_tb$winner) +
#       
#       # winning candidates
#       map_int(contrast_pair_tb$winner, ~length(.x)) *
#       !is.na(contrast_pair_tb$winner)
#     
#     cat("genes discarded for each node: ", j, "\n")
#     cat("genes selected for each node: ", map_int(signature$cumulative_signature, ~ length(.x)),  "\n")
#     
#     i <- i + 1
#     cat("iteration: ", i, "\n")
#     
#   }
#   
#   output <- summary_tb %>% 
#     nest(signature_data = - ancestor) %>% 
#     mutate(level = "non_hierarchy", .before = ancestor)
#   
#   return(output)
# }
# 
# # calculates silhouette score for each set of signature (cumulative markers) at a signature size
# silhouette_for_markers0 <-function(.signature, .contrast, .method) {
#   
#   .contrast %>%
#     
#     # filter markers that are in the signature
#     mutate(markers = map(markers, ~.x %>% 
#                            filter(symbol %in% .signature))) %>% 
#     
#     # format statistics from pairwise contrast
#     mutate(markers  = map(markers, 
#                           ~ .x %>% 
#                             # Group by contrast. Comparisons both ways.
#                             pivot_longer(
#                               cols = contains("___"),
#                               names_to = c("stats", "contrast"), 
#                               values_to = ".value", 
#                               names_sep="___"
#                             ) %>% 
#                             
#                             # Markers selection within each pair of contrast
#                             nest(stat_df = -contrast) %>%
#                             
#                             # Reshape inside each contrast
#                             mutate(stat_df = map(stat_df, ~.x %>% 
#                                                    pivot_wider(names_from = stats, 
#                                                                values_from = .value))) %>% 
#                             
#                             unnest(stat_df) )) %>% 
#     
#     # Add original data data to the markers selected
#     mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>% 
#     
#     # select only columns needed
#     select(-data) %>% 
#     
#     # create an input node (containing level_0 & markers) for the sil_func0 
#     nest(silhouette_data = c(level_0, markers)) %>% 
#     
#     mutate(silhouette_data = map(silhouette_data, ~ .x %>% unnest(markers))) %>% 
#     
#     # Calculate silhouette score for PCA plot resulted from the markers selected
#     mutate(silhouette_data = map(silhouette_data, ~ sil_func0(.x, .method))) %>% 
#     
#     mutate(reduced_dimensions = map(silhouette_data, ~ .x$rdim[[1]])) %>% 
#     mutate(silhouette = map_dbl(silhouette_data, ~ .x$sil)) %>% 
#     select(-silhouette_data)
#   
# }
# 
# 
# 
# # Modularised functions
# 

# Modularised functions ====================================================

# create input files for cibersortx cell signature selection using bulk RNA-seq data
produce_cibersortx_bulk_rnaseq_input <- function(.expression_df, .transcript, .sample, .cell_type, .count, 
                                                  .dir, .tree=NULL, .suffix=NULL){
  
  # Args: 
  # .expression_df is a tibble with shape: transcript | sample | cell_type | count
  # .transcript: tell the function which column in the data represents gene symbol
  # .sample: tell the function which column in the data represents sample
  # .cell_type: tell the function which column in the data represents cell_type
  # .count: tell the function which column in the data represents count
  # .tree: optional argument. if .tree = NULL, cell types will be sourced from those present in the .expression_df column
  #       if provided, the cell types used for cibersortx will be the leaves of the tree
  # .dir: the path to the directory where the two output files should be saved
  # .suffix: optional argument which adds suffix to the output filename: reference.txt and phenoclass.txt
  
  .transcript = enquo(.transcript)
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)
  
  .dir <- ifelse(grepl("\\/$", .dir), .dir, paste0(.dir, "/"))
  
  ## ref_names is produced to create header for reference file
  ref_names <- .expression_df %>% 
    {
     if (is.null(.tree)){
       (.)
     } else {
       (.) %>% 
         filter(!!.cell_type %in% as.phylo(.tree)$tip.label)
     }
       
    } %>% 
    mutate(!!.sample := str_replace_all(!!.sample, "\\.", "_")) %>%
    unite(cell_sample, c(!!.cell_type, !!.sample), sep = ".") %>% 
    pull(cell_sample) %>% 
    unique()
  
  
  # create reference file
  .expression_df %>% 
    
    # filter for leaf cell types in the tree
    {
      if (is.null(.tree)){
        (.)
      } else {
        (.) %>% 
          filter(!!.cell_type %in% as.phylo(.tree)$tip.label)
      }
      
    } %>% 
    select(!!.sample, !!.cell_type, !!.count, !!.transcript) %>% 
    mutate(!!.sample := str_replace_all(!!.sample, "\\.", "_")) %>%
    pivot_wider(names_from = c(!!.cell_type, !!.sample), values_from = !!.count, names_sep=".") %>% 
    `names<-`(ref_names %>% 
                str_extract(".*(?=\\.)") %>% 
                prepend(quo_name(.transcript))
    ) %>% 
    # save files as txt files (tab separated files)
    write_tsv(glue("{.dir}reference{.suffix}.txt"))
  
  # create phenoclass file
  tibble(cell_type = 
           {
             if (is.null(.tree)){
               .expression_df %>% 
                 distinct(!!.cell_type) %>% 
                 pull %>% 
                 as.character
             } else {
               as.phylo(.tree)$tip.label
             }
             
           }
         
         ) %>% 
    bind_cols(
      tibble(ref_names, value=1L) %>% 
        pivot_wider(names_from = ref_names, values_from = value)
    ) %>% 
    pivot_longer(-cell_type, names_to="ref_names", values_to="membership") %>% 
    
    # assign membership according to cibersortx tutorial 6: 
    # "1" indicates membership of the reference sample to the class as defined in that row,
    # "2" indicates the class that the sample will be compared against
    mutate(membership = if_else(str_detect(ref_names, cell_type), 1L, 2L)) %>% 
    
    pivot_wider(names_from = ref_names, values_from = membership) %>% 
    
    # set col_names to FALSE following cibersortx format for phenoclass
    write_tsv(glue("{.dir}phenoclass{.suffix}.txt"), col_names = FALSE)
  
}


# when stefano fixes missing count in as_summarisedExp(), remove the default NULL value for .count
main <- function(.input, .sample, .symbol, .count=NULL, .cell_type,
                 .is_hierarchy=TRUE, .level=NULL, 
                 .tree, .node=NULL,
                 .contrast_method, .ranking_method, .rank_stat=NULL, .bayes=NULL, 
                 .selection_method, .kmax=60, .discard_number=2000, .reduction_method = "PCA", .dims=2,
                 .optimisation_method, .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
                 .is_complete = TRUE) {
  
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  .cell_type = enquo(.cell_type)
  
  # subtree = tree_subset(.tree=.tree, .node=.node)
  
  if ( (!.is_hierarchy)|(.selection_method == "naive")) {
    
    .input %>%
      
      # adapt_tree(.tree = .tree, .node = .node) %>%
      # 
      # tree_and_signatures_to_database(tree=subtree, ., .sample=!!.sample, .cell_type=!!.cell_type,
      #                                .symbol=!!.symbol, .count=!!.count) %>%
      
      # # Remove redundant samples
      # remove_redundancy(.element=!!.sample, .feature=!!.symbol, .abundance=!!.count, correlation_threshold = 0.999, top = 500, method = "correlation") %>%
      # droplevels() %>%
      # 
      # # Eliminate suspicious samples
      # filter(!grepl("GSM3722278|GSM3722276|GSM3722277", !!.sample)) %>%
      
      # do_scaling(.sample = !!.sample, .symbol= !!.symbol , .count= !!.count, .cell_type=!!.cell_type) %>% 
      #   
      # do_imputation(.sample = !!.sample, .symbol=feature, .count= !!.count, .cell_type=!!.cell_type) %>% 
      

      # do_hierarchy(.sample=!!.sample,
      #              .symbol=!!.symbol,
      #              .cell_type = !!.cell_type,
      #              .tree = .tree,
      #              .is_hierarchy=.is_hierarchy,
      #              .level=.level) %>%
      
      # Input: data.frame columns_1 <int> | ...
      do_ranking(.sample=!!.sample, 
                 .symbol=!!.symbol,
                 .cell_type = !!.cell_type,
                 .ranking_method=.ranking_method, 
                 .contrast_method=.contrast_method, 
                 .rank_stat=.rank_stat, 
                 .bayes=.bayes,
                 .tree = .tree) %>%  
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      do_selection(.sample=!!.sample, 
                   .symbol=!!.symbol,
                   .selection_method=.selection_method, 
                   .reduction_method=.reduction_method, 
                   .discard_number=.discard_number, 
                   .kmax=.kmax,
                   .dims=.dims) %>% 
      
      do_optimisation(.optimisation_method=.optimisation_method, .symbol=!!.symbol) %>% 
      
      format_output(.is_complete = .is_complete)
    
  } else {
    
    .input %>%
      
      # adapt_tree(.tree = .tree, .node = .node) %>% 
      # 
      # tree_and_signatures_to_database(tree=subtree, ., .sample=!!.sample, .cell_type=!!.cell_type, 
      #                                .symbol=!!.symbol, .count=!!.count) %>% 
      
      # # Remove redundant samples
      # remove_redundancy(.element=!!.sample, .feature=!!.symbol, .abundance=!!.count, correlation_threshold = 0.999, top = 500, method = "correlation") %>%
      # droplevels() %>% 
      # 
      # # Eliminate suspicious samples
      # filter(!grepl("GSM3722278|GSM3722276|GSM3722277", !!.sample)) %>% 
      
      # do_scaling(.sample = !!.sample, .symbol= !!.symbol , .count= !!.count, .cell_type=!!.cell_type) %>% 
      #   
      # do_imputation(.sample = !!.sample, .symbol=feature, .count= !!.count, .cell_type=!!.cell_type) %>% 
      
 
      # do_hierarchy(.sample=!!.sample,
      #              .symbol=!!.symbol,
      #              .cell_type = !!.cell_type,
      #              .tree = .tree,
      #              .is_hierarchy=.is_hierarchy,
      #              .level=.level) %>%
      
      do_ranking(.sample=!!.sample, 
                 .symbol=!!.symbol,
                 .cell_type = !!.cell_type,
                 .ranking_method=.ranking_method, 
                 .contrast_method=.contrast_method, 
                 .rank_stat=.rank_stat, 
                 .bayes=.bayes,
                 .tree = .tree) %>% 
      
      mutate(level.copy = level) %>% 
      nest(data = -level.copy) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          
          do_selection(.sample=!!.sample, 
                       .symbol=!!.symbol,
                       .selection_method=.selection_method, 
                       .reduction_method=.reduction_method, 
                       .discard_number=.discard_number, 
                       .kmax=.kmax,
                       .dims=.dims) %>%
          
          do_optimisation(.optimisation_method = .optimisation_method, .symbol=!!.symbol) %>%
          
          format_output(.is_complete = .is_complete)
        
      )) %>% 
      
      unnest(data) %>% 
      select(-level.copy)
    
  }
}

# Add tree structure to raw data frame
tree_and_signatures_to_database = function(tree, signatures, .sample, .cell_type, .symbol, .count){
  
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .symbol = enquo(.symbol)
  .count = enquo(.count)

  signatures %>%

    # Add tree info
    left_join(
      tree %>%
        data.tree::Clone() %>%
        ToDataFrameTypeColFull(fill=NA) ,
      by = quo_name(.cell_type)
    ) %>%
    filter(level_1 %>% is.na %>% `!`) %>%

    # Reduce size
    mutate_if(is.character, as.factor) %>%
    droplevels %>%
    mutate(!!.count := !!.count %>% as.integer) %>%

    # Filter only symbol existing
    filter(!!.symbol %>% is.na %>% `!`) %>%

    # Aggregate
    aggregate_duplicates(!!.sample, !!.symbol, !!.count) %>%
    select(-one_of("merged_transcripts"))
}

# scale abundance
do_scaling <- function(.data_tree, .sample, .symbol, .count, .cell_type) {
  
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  .cell_type = enquo(.cell_type)
  
  .data_tree %>% 
    
    # ensure all genes are present to all level_1 cell types(missing genes for other level cell types can be imputed)
    nest(data = -c(level_1, !!.symbol)) %>%
    add_count(!!.symbol) %>%
    filter(n==n_distinct(.data_tree$level_1)) %>%
    select(-n) %>%
    unnest(data) %>%
    
    # Convert to SE
    as_SummarizedExperiment(!!.sample, !!.symbol, !!.count)  %>%
    
    # Scale with first degree imputation. 
    # This because there are no common genes to all samples. Note that quasiquotation doesn't work in a formula
    impute_missing_abundance(.formula = as.formula(sprintf("~ %s", quo_name(.cell_type)))) %>% 
    identify_abundant() %>%
    scale_abundance() %>%
    filter(!.imputed) %>% 
    
    select(-.imputed)  %>%
    
    # Just needed for the old version
    # select(-one_of("exposure_rate")) %>%
    
    # Calculate exposure for Bayes model
    mutate(exposure_rate = -log(multiplier))
}

# imputation
do_imputation <- function(.scaled_counts, .sample, .symbol, .count, .cell_type){
  
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  .cell_type = enquo(.cell_type)
  
  .scaled_counts %>% 
    # Convert to SE
    # as_SummarizedExperiment(.sample, .feature, count) %>%
    as_SummarizedExperiment(!!.sample, !!.symbol, .abundance = c(!!.count, count_scaled)) %>%
    
    # Hierarchical imputation. Suffix = "" equated to overwrite counts
    impute_missing_abundance(.formula = as.formula(sprintf("~ %s", quo_name(.cell_type))),
                             .abundance = c(!!.count, count_scaled)) %>%
    
    # AUTOMATE THIS!
    impute_missing_abundance(~ level_5, .abundance = c(!!.count, count_scaled)) %>%
    impute_missing_abundance(~ level_4, .abundance = c(!!.count, count_scaled)) %>%
    impute_missing_abundance(~ level_3, .abundance = c(!!.count, count_scaled)) %>%
    impute_missing_abundance(~ level_2, .abundance = c(!!.count, count_scaled)) %>%
    impute_missing_abundance(~ level_1, .abundance = c(!!.count, count_scaled)) %>% 
    
    # {
    #   for(level in (.) %>% colnames %>% str_extract("level\\_\\d") %>% .[!is.na(.)]) {
    #     (.) <- (.) %>% 
    #       impute_missing_abundance(~ !!as.symbol(level), .abundance = c(!!.count, count_scaled))
    #   }
    #   
    # } %>% 
    
    # Convert back to tibble
    as_tibble() %>%
    
    mutate(.imputed = if_any(contains("imputed"), ~ .x != 0)) %>% 
    
    select(-matches("imputed\\.\\d"))
}

# Input scale abundance

do_hierarchy <- function(.imputed_counts, .sample, .symbol, .cell_type, .tree, .is_hierarchy = TRUE, .level = NULL){
  
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)

  if (.is_hierarchy) { # scaling under each ancestor node at each level

    if (is.null(.level)) { # scale counts for all levels present

      tibble(level = names(.imputed_counts) %>% str_subset("level")) %>%

        mutate(tt = map(level, ~ .imputed_counts %>%

                          mutate(level_0 = "root") %>%

                          create_hierarchy_and_calculate_imputation_ratio(.level=.x, .sample=!!.sample, .symbol=!!.symbol))) %>%

        mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename(ancestor = pre(.y))))

    } else { # scale counts for the level specified by .level

      tibble(level = .level) %>%

        mutate(tt = map(level, ~ .imputed_counts %>%

                          mutate(level_0 = "root") %>%

                          create_hierarchy_and_calculate_imputation_ratio(.level=.x, .sample=!!.sample, .symbol=!!.symbol))) %>%

        mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename(ancestor = pre(.y))))
    }

  } else { # non-hierarchical: scaling all cell types under the root node

    .level <- "root"

    tibble(level = .level) %>%

      mutate(tt = map(level, ~ .imputed_counts %>%
                        
                        # non-hierarchical methods should only compare leaf cell types
                        filter(!!.cell_type %in% as.phylo(.tree)$tip.label) %>% 

                        # create a root column for pre(.level)
                        mutate(!!as.symbol(.level) := !!.cell_type) %>%
                        mutate(!!as.symbol(pre(.level)) := .level) %>%

                        create_hierarchy_and_calculate_imputation_ratio(.level=.x, .sample=!!.sample, .symbol=!!.symbol))) %>%

      mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename(ancestor = pre(.y))))

  }
}

create_hierarchy_and_calculate_imputation_ratio <- function(.imputed_counts, .level, .sample, .symbol) {
  # this preproces function ranged data in hierarchy(or non_hierarchy) and
  # calculates the imputation ratio for genes in each hierarchy
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)


  # load data
  .imputed_counts %>%

    # tidybulk(sample, symbol, count_scaled) %>% for imputed counts data
    # tidybulk(.sample = !!.sample, .transcript = !!.symbol, .abundance = !!.count) %>%

    # filter for cells at the level of interest. .level == level_1
    filter(!is.na(!!as.symbol(.level))) %>%

    # calculate the ratio of imputation for genes in a cell type
    nest(data = -c(!!.symbol, !!as.symbol(.level))) %>%

    # for a cell type some samples may miss genes in other samples: so for the same cell type genes may have different number of samples
    mutate(n_samples_per_gene = map_int(
      data,
      ~ .x %>%
        distinct(!!.sample) %>%
        nrow)) %>%

    unnest(data) %>%
    nest(data = -!!as.symbol(.level)) %>%

    mutate(n_samples= map_int(
      data,
      ~ .x %>%
        distinct(!!.sample) %>%
        nrow
    )) %>%

    unnest(data) %>%

    mutate(ratio_imputed_samples = 1 - n_samples_per_gene / n_samples) %>%

    # nest by ancestor
    nest(data = - !!as.symbol(pre(.level)))

}

pre <- function(.level) {

  if (.level == "root" ) {

    return("level_0")

  } else {

    ancestor <-

      .level %>%

      # str_split returns a list
      str_split("_") %>%

      {as.numeric(.[[1]][2])-1} %>%

      paste("level", ., sep = "_")

    return(ancestor)
  }
}

# Generate contrast for ranking

pairwise_contrast = function(.data, .level){

  .data %>%
    distinct(!!as.symbol(.level)) %>%
    mutate(!!as.symbol(.level) := paste0(.level, !!as.symbol(.level))) %>%

    # Permute
    mutate(cell_type2 = !!as.symbol(.level)) %>%
    tidyr::expand(!!as.symbol(.level), cell_type2) %>%
    filter(!!as.symbol(.level) != cell_type2) %>%

    # Create contrasts
    mutate(contrast = sprintf("%s - %s", !!as.symbol(.level), cell_type2)) %>%
    pull(contrast)

}

mean_contrast <- function(.data, .level){

  cell_types <- .data %>%
    distinct(!!as.symbol(.level)) %>%
    pull() %>%
    paste0(.level, .)

  .data %>%
    distinct(!!as.symbol(.level)) %>%
    mutate(!!as.symbol(.level) := paste0(.level, !!as.symbol(.level))) %>%
    mutate(background = map(
      !!as.symbol(.level), ~ cell_types[cell_types != .x])) %>%
    mutate(contrast = map2_chr(
      !!as.symbol(.level), background,
      ~ sprintf("%s - (%s)/%s", .x, paste(.y, collapse="+"), length(.y))
    )) %>%
    pull(contrast)
}

# Rank

do_ranking <- function(.hierarchical_counts, .sample, .symbol, .cell_type, 
                       .ranking_method, .contrast_method, 
                       .rank_stat=NULL, .bayes=NULL, .tree=NULL){

  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)

  # rank_stat takes either "Pvalue" or "logFC"
  .hierarchical_counts %>%

    # .bayes = .cellsig_theoretical_transcript_abundace_distribution
    .ranking_method(.sample=!!.sample,
                    .symbol=!!.symbol,
                    .cell_type = !!.cell_type,
                    .contrast_method=.contrast_method,
                    .rank_stat=.rank_stat,
                    .bayes=.bayes,
                    .tree = .tree) %>%

    # unnest(markers) %>%

    # filter out potential nodes in which no genes are considered significant by rank_by_stat
    filter(map_int(stat_df, ~nrow(.x)) != 0) %>%

    # assign ranks to genes
    mutate(stat_df = map(stat_df, ~ mutate(.x, rank = 1:nrow(.x)))) %>%

    nest(markers = -c(level, ancestor, data))

}

rank_edgR_quasi_likelihood <- function(.hierarchical_counts, .sample, .symbol, .cell_type, 
                                       .contrast_method, .rank_stat, 
                                       .bayes=NULL, .tree=NULL){

  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)

  .hierarchical_counts %>%
    unnest(tt) %>%

    # Differential transcription: generate contrast
    mutate(markers = map2(
      data, level,
      ~ .x %>%
        test_differential_abundance(
          .formula = as.formula(sprintf("~ 0 + %s", .y)),
          .sample = !!.sample,
          .transcript = !!.symbol,
          .abundance = count_scaled,
          .contrasts = .contrast_method(.x, .y),
          method = "edgeR_quasi_likelihood",
          action="only")
    )) %>%

    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ rank_by_stat(.x, .rank_stat) )) %>%

    unnest(markers) %>%

    # remove prefixes from contrast expressions
    mutate(contrast = map2_chr(contrast, level, ~ str_remove_all(.x, .y) ))

    # # only used for the user pipeline, not for benchmark
    # # filter for genes that have imputation rate less than 0.2
    # mutate(stat_df = pmap(
    #   list(data, contrast, stat_df),
    #
    #   ~..1 %>%
    #
    #     # filter for target cell type in data
    #     filter(!!.cell_type == str_extract(..2, ".*(?=\\s\\-)")) %>%
    #
    #     # select the symbols and ratio of imputed samples for that gene in the target cell type
    #     distinct(!!.symbol, ratio_imputed_samples) %>%
    #
    #     right_join(..3, by = quo_name(.symbol)) %>%
    #
    #     filter(ratio_imputed_samples < 0.2)
    # ))

}

rank_edgR_robust_likelihood_ratio <- function(.hierarchical_counts, .sample, .symbol, .cell_type, 
                                              .contrast_method, .rank_stat="PValue",
                                              .bayes=NULL, .tree=NULL){
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)

  .hierarchical_counts %>%
    unnest(tt) %>%

    # Differential transcription: generate contrast
    mutate(markers = map2(
      data, level,
      ~ .x %>%
        test_differential_abundance(
          .formula = as.formula(sprintf("~ 0 + %s", .y)),
          .sample = !!.sample,
          .transcript = !!.symbol,
          .abundance = count_scaled,
          .contrasts = .contrast_method(.x, .y),
          method = "edger_robust_likelihood_ratio",
          test_above_log2_fold_change = 1,
          action="only")
    )) %>%

    # Select markers from each contrast by rank of Pvalue
    mutate(markers = map(
      markers,
      ~ rank_by_stat(.x, "PValue")
    )) %>%

    unnest(markers) %>%

    # remove prefixes from contrast expressions
    mutate(contrast = map2_chr(contrast, level, ~ str_remove_all(.x, .y) ))

    # # only used for the user pipeline, not for benchmark
    # # filter for genes that have imputation rate less than 0.2
    # mutate(stat_df = pmap(
    #   list(data, contrast, stat_df),
    #
    #   ~..1 %>%
    #
    #     # filter for target cell type in data
    #     filter(!!.cell_type == str_extract(..2, ".*(?=\\s\\-)")) %>%
    #
    #     # select the symbols and ratio of imputed samples for that gene in the target cell type
    #     distinct(!!.symbol, ratio_imputed_samples) %>%
    #
    #     right_join(..3, by = quo_name(.symbol)) %>%
    #
    #     filter(ratio_imputed_samples < 0.2)
    #   ))
}

rank_by_stat <-  function(.markers, .rank_stat){

  .markers %>%

    # Group by contrast. Comparisons both ways.
    pivot_longer(
      cols = contains("___"),
      names_to = c("stats", "contrast"),
      values_to = ".value",
      names_sep="___"
    ) %>%

    # Markers selection within each pair of contrast
    nest(stat_df = -contrast) %>%

    # Reshape inside each contrast
    mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>%

    # Keep significant genes and rank the significant ones
    # THIS WILL HAVE TO CHANGE
    mutate(stat_df = map(
      stat_df,
      ~ .x %>%
        filter(FDR < 0.05 & logFC > 2) %>%
        filter(logCPM > mean(logCPM)) )) %>%

    mutate(stat_df = map(
      stat_df,
      ~ if(.rank_stat == "logFC"){
        .x %>% dplyr::arrange(desc(logFC))
      }else{
        .x %>% dplyr::arrange(PValue)
      }
    ))
}

rank_bayes <- function(.hierarchical_counts, .sample, .symbol, .cell_type, 
                       .contrast_method, .bayes, .tree, 
                       .rank_stat=NULL){

  # Args:
  # .bayes is the bayes data that have been imputed and obtained gene imputation ratio by the code:
  # counts_bayes_imputed %>%
  #   select(-level) %>%
  #   scale_input_counts(.is_hierarchy = TRUE)
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)

  .bayes %>%
    
    # force the column names of bayes data to be consistent with input expression data
    # dplyr::rename(!!.symbol := feature, !!.sample := sample, !!.cell_type := cell_type) %>%
    # do_hierarchy(.is_hierarchy = all(.hierarchical_counts$level != "root"), .tree = .tree,
    #              .sample=!!.sample, .symbol=!!.symbol, .cell_type= !!.cell_type) %>%

    unnest(tt) %>%

    mutate(contrast = map2(data, level, ~ .contrast_method(.x, .y))) %>%

    unnest(contrast) %>%
    mutate(contrast = str_remove_all(contrast, glue("{level}"))) %>%

    mutate(lower_quantile = map2(
      data, contrast,
      ~ .x %>%
        # filter for target cell type in the contrast
        filter(!!.cell_type == str_extract(.y, ".*(?=\\s\\-)")) %>%
        # filter out genes with imputation ratio greater than 0.2 (only used for user pipeline not benchmark)
        # filter(ratio_imputed_samples < 0.2) %>%
        select(!!.symbol, lower_quantile='25%')
    )) %>%

    mutate(mean_upper_quantile = map2(
      data, contrast,
      ~ {
        # obtain background cell type(s)
        background <- (.y) %>%
          str_extract("(?<=\\-\\s).*") %>%
          str_split("\\+") %>%
          unlist() %>%
          str_remove_all("(?<=\\/).*|\\W")

        (.x) %>%
          # calculate the mean 75% quantile of each gene over all background cell types
          filter(!!.cell_type %in% background) %>%
          group_by(!!.symbol) %>%
          summarise(!!.symbol, mean_upper_quantile = mean(`75%`)) %>%
          distinct() %>%
          ungroup()
      }
    )) %>%

    mutate(stat_df = map2(
      lower_quantile, mean_upper_quantile,
      ~ left_join(.x, .y, by= quo_name(.symbol))
    )) %>%
    select(-c(lower_quantile, mean_upper_quantile)) %>%

    mutate(stat_df = map(
      stat_df,
      ~ .x %>%
        mutate(difference = lower_quantile - mean_upper_quantile) %>%
        arrange(desc(difference))
    )) %>%

    select(-data) %>%

    nest(markers = -c(level, ancestor)) %>%

    right_join(.hierarchical_counts %>% unnest(tt), by = c("level", "ancestor")) %>%

    unnest(markers)

}



# Selection

do_selection <-
  function(.ranked, .sample, .symbol, .selection_method, .reduction_method, .kmax=NULL, .discard_number=NULL, .dims=2) {

    .sample = enquo(.sample)
    .symbol = enquo(.symbol)

    # .k is the number of genes selected from each cell_type contrast

    if (.selection_method == "naive") {

      .ranked %>%

        do_naive_selection(.sample=!!.sample,
                           .symbol=!!.symbol,
                           .kmax=.kmax,
                           .reduction_method=.reduction_method,
                           .dims=.dims)

    } else {

      .ranked %>%

        do_silhouette_selection(.sample=!!.sample,
                                .symbol=!!.symbol,
                                .discard_number=.discard_number,
                                .reduction_method=.reduction_method,
                                .dims=.dims)

    }

  }

## Naive selection

### calculate silhouette score for a series of sig_sizes

min_markers_per_contrast <- function(.markers, .dims, .symbol){
  
  .symbol = enquo(.symbol)
  
  k = 1L
  n_unique_markers = .markers %>% 
    mutate(top_k = map(stat_df, ~ .x %>% slice(1:k) %>% pull(!!.symbol))) %>% 
    pull(top_k) %>% 
    unlist %>% 
    n_distinct
  
  while (n_unique_markers < .dims) {
    k = k + 1L
    n_unique_markers = .markers %>% 
      mutate(top_k = map(stat_df, ~ .x %>% slice(1:k) %>% pull(!!.symbol))) %>% 
      pull(top_k) %>% 
      unlist %>% 
      n_distinct
  }
  
  return(k)
}

do_naive_selection <- function(.ranked, .sample, .symbol, .kmax, .reduction_method, .dims=2) {

  # Args:
  # .ranked: output from do_ranking
  # .kmax: maximum number of markers selected from each cell type contrast
  # .reduction_method: method used to reduce dimensions such as "PCA", "tSNE", "MSA", "UMAP"
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  
  .ranked %>% 
    
    # calculate the minimum number of genes need to be selected for feasible dimension reduction
    mutate(k0 = map_int(markers, ~ min_markers_per_contrast(.x, .dims=.dims, .symbol = !!.symbol))) %>% 
    
    # expand the number of markers selected to .kmax
    mutate(n_markers_from_each_contrast = map(k0, ~ .x: .kmax)) %>% 
    
    # clean unnecessary data
    select(-k0) %>%
    
    # expand the nuber of markers need to be selected
    unnest(n_markers_from_each_contrast) %>% 
    nest(data = -n_markers_from_each_contrast) %>% 
    
    # do naive selection according to the number of markers from each contrast and calculate silhouette score for the signature
    mutate(data = map2(
      data, n_markers_from_each_contrast,
      
      ~ naive_selection(.ranked=.x, .k=.y, .symbol=!!.symbol) %>%
        
        silhouette_function(.sample=!!.sample,
                            .symbol=!!.symbol,
                            .reduction_method=.reduction_method,
                            .dims=.dims) %>%
        
        # optional: remove reduced_dimensions matrix to save memory, if kept it can be used to plot PCA
        select(-reduced_dimensions)
    )) %>%
    
    # nest by ancestor nodes/cell types
    unnest(data) %>%
    nest(data = - c(level, ancestor))

}

naive_selection <- function(.ranked, .k, .symbol) {

  # Args:
  # .ranked: output from do_ranking()
  # .k: the number of genes selected from each cell_type contrast
  .symbol = enquo(.symbol)

  .ranked %>%

    # selection markers from each contrast
    mutate(markers = map(
      markers,
      ~ .x %>%
        mutate(stat_df = map(stat_df, ~ .x %>% dplyr::slice(1: .k))) %>%
        unnest(stat_df)
    )) %>%

    # collect from which contrasts signature genes are extracted
    mutate(children = map(markers, ~ .x %>%
                            select(contrast, !!.symbol, rank) %>%
                            nest(enriched = - contrast)
    )) %>%

    # Add original expression data info to the markers selected for dimension reduction,
    # use inner_join to ensure symbols are present in both markers and data
    mutate(markers = map2(markers, data, ~ inner_join(.x, .y, by=quo_name(.symbol)))) %>%

    # remove unnecessary column
    select(-data) %>%

    # collect signature genes selected
    mutate(signature = map(markers, ~ .x %>% pull(!!.symbol) %>% unlist %>% unique)) %>%

    # number of the signature genes
    mutate(real_size = map_int(signature, ~ length(.x)))

}

### calculate silhouette score

silhouette_function <- function(.selected, .sample, .symbol, .reduction_method, .dims=2){

  .sample = enquo(.sample)
  .symbol = enquo(.symbol)

  .selected %>%

    # reduce dimensions
    mutate(reduced_dimensions = map2(
      markers, level,
      ~ dimension_reduction(.markers=.x, .level=.y,
                            .sample = !!.sample,
                            .symbol = !!.symbol,
                            .reduction_method=.reduction_method,
                            .dims=.dims)
    )) %>%

    # calculate distance matrix using PC1 & PC2
    mutate(distance = map(
      reduced_dimensions,
      ~ distance_matrix(.x, .reduction_method=.reduction_method)
    )) %>%

    # calculate silhouette score
    mutate(silhouette = pmap_dbl(
      list(reduced_dimensions, distance, level),
      ~ silhouette_score(..1, ..2, ..3)
    )) %>%

    # remove unnecessary columns
    select(-c(markers, distance))

}

dimension_reduction <- function(.markers, .level, .sample, .symbol, .reduction_method, .dims=2) {

  .sample = enquo(.sample)
  .symbol = enquo(.symbol)

  .markers %>%

    distinct(!!.sample, !!.symbol, count_scaled, !!as.symbol(.level)) %>%

    reduce_dimensions(.element = !!.sample,
                      .feature = !!.symbol,
                      .abundance = count_scaled,
                      action = "get",
                      method = .reduction_method,
                      .dims = .dims,
                      log_transform = TRUE,
                      top = Inf,
                      scale = FALSE,
                      check_duplicates = FALSE)
}

distance_matrix <- function(.reduced_dimensions, .reduction_method){

  .reduced_dimensions %>%

    select(contains(str_sub(.reduction_method, end = -2L))) %>%

    factoextra::get_dist(method = "euclidean")
}

silhouette_score <- function(.reduced_dimensions, .distance, .level){

  .reduced_dimensions %>%

    pull(!!as.symbol(.level)) %>%

    as.factor() %>%

    as.numeric() %>%

    silhouette(.distance) %>%

    summary() %>%

    .$avg.width

}

## Silhouette selection

do_silhouette_selection <-
  function(.ranked, .sample, .symbol, .discard_number, .reduction_method, .dims=2) {

    .sample = enquo(.sample)
    .symbol = enquo(.symbol)

    # initialize variables

    # ranked_copy is created as a pool of markers for selection,
    # which continuously decrease with each iterative selection,
    # input .rank is used for calculating silhouette score for the selected markers
    ranked_copy <- .ranked

    # initialise a signature tibble to store signature markers for each cell type in each iteration
    signature <- .ranked %>%
      select(level, ancestor) %>%
      mutate(signature = map(ancestor, ~ vector())) %>%
      mutate(last_silhouette = 0)

    # initialise an output tibble containing all results of interest
    summary_tb <- tibble(
      level = character(),
      ancestor = character(),
      new_challengers = list(),
      winner = list(),
      children = list(),
      signature = list(),
      # reduced_dimensions = list(),
      silhouette = double()
    )

    # set the base markers
    contrast_pair_tb0 <-

      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      .ranked %>%
      
      # calculate the minimum number of genes need to be selected for feasible dimension reduction
      mutate(k0 = map_int(markers, ~ min_markers_per_contrast(.x, .dims=.dims, .symbol = !!.symbol))) %>% 
      
      # nest by k0 because naive selection takes the whole ranked data frame
      nest(data = -k0) %>% 
      
      # select the minimum number of base markers for each ancestor node
      mutate(data = map2(data, k0, ~ naive_selection(.ranked=.x, .k=.y, .symbol=!!.symbol))) %>% 
      
      unnest(data) %>% 
      
      # clean up data
      select(-k0) %>% 

      dplyr::rename(new_challengers = signature) %>%

      mutate(winner = map(new_challengers, ~ unique(.x))) %>%

      mutate(signature = winner) %>%

      silhouette_function(.sample=!!.sample, .symbol=!!.symbol, .reduction_method=.reduction_method, .dims=.dims) %>%

      select(-c(reduced_dimensions, real_size))


    signature <- signature %>%

      # append cumulative markers
      mutate(signature = map2(
        signature, ancestor,
        ~ .x %>%
          append(with(contrast_pair_tb0, signature[ancestor==.y][[1]]))
      )) %>%

      # append silhouette scores for these markers
      mutate(last_silhouette = map_dbl(
        ancestor,
        ~ with(contrast_pair_tb0, silhouette[ancestor==.x])
      ))

    summary_tb <- summary_tb %>%
      bind_rows(contrast_pair_tb0)

    # remove base markers from contrast_copy input before further selection
    ranked_copy <- ranked_copy %>%
      mutate(markers = map2(
        markers, ancestor,
        ~ .x %>%
          unnest(stat_df) %>%
          filter(!(!!.symbol %in% with(signature, signature[ancestor==.y][[1]]))) %>%
          nest(stat_df = - contrast)
      ))

    # counter for number of genes discarded
    j <- map_int(signature$signature, ~ length(.x))

    # count the number of iterations
    i <- 0L
    while (any(j < .discard_number) &
           # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
           all(map_int(ranked_copy$markers,
                       # hence the boundary should be the number of satisfactory genes selected
                       ~ .x %>% unnest(stat_df) %>% nrow()) > 0)) {

      contrast_pair_tb <-

        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        ranked_copy %>%

        # select top 1 markers from each contrast, ignore the signature output
        naive_selection(.k=1, .symbol=!!.symbol) %>%

        # pick the one new challenger from each contrast
        select(-c(markers, signature, real_size)) %>%
        unnest(children) %>%
        unnest(enriched) %>%
        dplyr::rename(new_challenger = !!.symbol) %>%

        # append the new challenger from each contrast to the base markers for that ancestor node
        mutate(challengers_for_silhouette = map2(
          new_challenger, ancestor,
          ~ with(signature, signature[ancestor==.y][[1]]) %>%
            append(.x)
        )) %>%

        # calculate silhouette score for the challengers from each contrast
        mutate(silhouette = map2_dbl(
          challengers_for_silhouette, ancestor,
          ~ silhouette_for_markers(.ranked=.ranked, .sample=!!.sample, .symbol=!!.symbol,
                                   .signature=.x, .ancestor=.y, 
                                   .reduction_method=.reduction_method, .dims=.dims) %>%
            pull(silhouette)
        )) %>%

        # arrange silhouette score in a descending manner within each ancestor node
        group_by(ancestor) %>%
        arrange(desc(silhouette), .by_group = TRUE) %>%
        ungroup() %>%

        # check if the silhouette score for the challengers is greater than previous silhouette score
        mutate(is_greater = map2_lgl(
          silhouette, ancestor,
          ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
        )) %>%

        # nest under ancestor node to select markers that is TRUE for is_greater
        nest(data = - c(level, ancestor)) %>%

        # record new_challengers
        mutate(new_challengers = map(data, ~ .x %>% pull(new_challenger))) %>%

        # check if the biggest silhouette score is greater than previous score, if true we have a winner, else no winner
        mutate(winner = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ]$new_challenger
        } else {NA}
        )) %>%

        # record which contrast the winner comes from
        mutate(children = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ] %>%
            select(contrast, !!.symbol := new_challenger, rank) %>%
            nest(enriched = -contrast)
        } else {NA}
        )) %>%


        # cummulative signature: winner + previously selected
        mutate(signature = pmap(
          list(data, winner, ancestor),
          ~ if(!is.na(..2)) {
            with(..1[1, ], challengers_for_silhouette[[1]])
          } else {
            with(signature, signature[ancestor==..3][[1]])
          }
        )) %>%

        # silhouette score
        mutate(silhouette = map_dbl(data, ~ .x[[1, "silhouette"]]))


      # append the base + 1 markers that result in highest silhouette score
      signature <- signature %>%

        mutate(signature = map(
          ancestor,
          ~ with(contrast_pair_tb, signature[ancestor==.x][[1]])
        )) %>%

        mutate(last_silhouette = map2_dbl(
          ancestor, last_silhouette,
          ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x]))) {
            with(contrast_pair_tb, silhouette[ancestor==.x])
          } else {.y}
        ))

      # append the winning signatures into the output summary table
      summary_tb <- summary_tb %>%
        bind_rows(
          contrast_pair_tb %>%
            filter(!is.na(winner)) %>%
            # mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>%
            select(-data)
        )

      # remove the signatures and unsuccessful genes from the selection list(ranked_copy)
      ranked_copy <- ranked_copy %>%
        mutate(markers = map2(
          markers, ancestor,
          ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y]))){
            .x %>%
              unnest(stat_df) %>%
              filter(!(!!.symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))) %>%
              nest(stat_df = -contrast)
          } else {
            .x %>%
              unnest(stat_df) %>%
              filter(!!.symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]])) %>%
              nest(stat_df = -contrast)
          }
        ))

      # number of genes discarded for each node
      j <- j +

        # unsuccessful candidates
        map_int(contrast_pair_tb$new_challengers, ~length(.x)) *
        is.na(contrast_pair_tb$winner) +

        # winning candidates
        map_int(contrast_pair_tb$winner, ~length(.x)) *
        !is.na(contrast_pair_tb$winner)

      cat("genes discarded for each node: ", j, "\n")
      cat("genes selected for each node: ", map_int(signature$signature, ~ length(.x)),  "\n")

      i <- i + 1L
      cat("iteration: ", i, "\n")

    }

    # format output for optimisation
    output <- summary_tb %>%
      mutate(real_size = map_int(signature, ~ length(.x))) %>%
      nest(data = - c(level, ancestor))

    return(output)
  }

silhouette_for_markers <-function(.ranked, .sample, .symbol, .signature, .ancestor, .reduction_method, .dims=2) {

  .sample=enquo(.sample)
  .symbol=enquo(.symbol)

  .ranked %>%

    filter(ancestor == .ancestor) %>%

    select(-markers) %>%

    # filter markers that are in the signature
    mutate(data = map(data, ~.x %>%
                        filter(!!.symbol %in% .signature))) %>%

    # format input for silhouette_function
    dplyr::rename(markers = data) %>%

    silhouette_function(.sample=!!.sample,
                        .symbol=!!.symbol,
                        .reduction_method = .reduction_method,
                        .dims=.dims)

}

# Optimisation

do_optimisation <- function(.selected,
                            .symbol,
                            .optimisation_method,
                            .penalty_rate=0.2,
                            .kernel = "normal",
                            .bandwidth = 0.05,
                            .gridsize = 100){
  .symbol = enquo(.symbol)

  {
    if (.optimisation_method == "penalty"){

      .selected %>%
        
        mutate(optimal_size = map_int(
          data, 
          ~ if(nrow(.x)==1L){
            .x %>% pull(real_size)
          }else{
            .x %>% penalised_silhouette(.penalty_rate=.penalty_rate)
          }
          ))

    } else if(.optimisation_method == "curvature") {

      .selected %>%
        
        mutate(optimal_size = map_int(
          data,
          ~ if(nrow(.x) == 1L){
            .x %>% pull(real_size)
          }else{
            .x %>% 
              curvature_of_kernel_smoothed_trend(.kernel=.kernel, .bandwidth=.bandwidth, .gridsize=.gridsize)
            }
        ))
    }

  } %>%

    { # this if statement below is for non_hierarchical method when the optimal size selected is less than
      # the number of cell types and deconvolution wouldn't be possible

      n_cell_type <- map_int((.)$data,
                              ~ .x$children %>%
                                .[[1]] %>%
                                .$contrast %>%
                                str_extract(".*(?=\\s\\-)") %>%
                                n_distinct)

      if ((.)$level == "root" &&
          (.)$optimal_size < n_cell_type){

        (.) %>%
          mutate(optimal_size = map_int(
            data,
            ~ with(.x, real_size[which(real_size >= n_cell_type) %>% min])
            ))

      } else {(.)}

    } %>%

    unnest(data) %>%

    filter(real_size <= optimal_size) %>%

    nest(children = -c(level, ancestor)) %>%

    mutate(signature = map(children, ~ .x %>% tail(1) %>% pull(signature) %>% unlist() )) %>%

    mutate(silhouette = map_dbl(children, ~ .x %>% tail(1) %>% pull(silhouette) %>% unlist() )) %>%

    mutate(children = map(
      children,
      ~ .x %>%
        unnest(children) %>%
        unnest(enriched) %>%
        distinct(contrast, !!.symbol, rank) %>%
        nest(enriched = -contrast)
    ))

}


curvature <- function(.drv1, .drv2){
  abs(.drv2) / (1 + .drv1^2)^(3/2)
}

curvature_of_kernel_smoothed_trend <- function(.plot_data,
                                               .kernel = "normal",
                                               .bandwidth = 0.05,
                                               .gridsize = 100){
  .plot_data %>%
    
    mutate(head = "head") %>% 
    nest(data = -head) %>% 

    mutate(data = map(
      data,
      ~ .x %>%
        mutate(size.rescaled = rescale(real_size))
    )) %>%

    mutate(smoothed.estimate = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette,
                drv = 0L, degree=2, kernel = .kernel,
                bandwidth = .bandwidth, gridsize = .gridsize) %>%
        as_tibble() %>%
        `colnames<-`(c("grid", "estimate"))
    )) %>%

    mutate(first.derivative = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette,
                drv = 1L, degree=2, kernel = .kernel,
                bandwidth = .bandwidth, gridsize = .gridsize) %>%
        as_tibble() %>%
        `colnames<-`(c("grid", "deriv1"))
    )) %>%

    mutate(second.derivative = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette,
                drv = 2L, degree=2, kernel = .kernel,
                bandwidth = .bandwidth, gridsize = .gridsize) %>%
        as_tibble() %>%
        `colnames<-`(c("grid", "deriv2"))
    )) %>%

    mutate(smoothed = pmap(
      list(smoothed.estimate, first.derivative, second.derivative),
      ~..1 %>%
        left_join(..2, by = "grid") %>%
        left_join(..3, by = "grid")
    )) %>%

    mutate(smoothed = map(
      smoothed,
      ~ .x %>%
        mutate(curvature = map2_dbl(
          deriv1, deriv2,
          ~ curvature(.drv1=.x, .drv2=.y)
        ))
    )) %>%

    select(-c(smoothed.estimate, first.derivative, second.derivative)) %>%

    # optimal_size
    mutate(optimal_size = map_dbl(
      smoothed,
      ~ with(.x, grid[which(peaks(curvature))[which.max(curvature[peaks(curvature)])]])
    )) %>%

    mutate(optimal_size = map2_dbl(
      data, optimal_size,
      ~ with(.x, size.rescaled[which.min(abs(size.rescaled - .y))])
    )) %>%

    mutate(optimal_size = map2_int(
      data, optimal_size,
      ~ .x %>%
        with(real_size[size.rescaled == .y])
    )) %>%
    
    # set the minimum signature size recommended to be 10.
    # note that some ancestor nodes might have fewer than 10 marker select and
    # the signature size will be the maximum number of markers selected for that node
    mutate(optimal_size = ifelse(optimal_size<10L, 10L, optimal_size)) %>% 
    
    pull(optimal_size)

}

penalised_silhouette <- function(.plot_data, .penalty_rate=0.2) {

  .plot_data %>%

    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    # filter(real_size > 10) %>%

    # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
    mutate(size_rescaled = rescale(real_size)) %>%

    mutate(penalised_silhouette = silhouette - .penalty_rate * size_rescaled) %>%

    filter(penalised_silhouette == max(penalised_silhouette)) %>%

    pull(real_size) %>% 
    
    # set the minimum signature size recommended to be 10.
    # note that some ancestor nodes might have fewer than 10 marker select and
    # the signature size will be the maximum number of markers selected for that node
    { if( (.) < 10L ){10L} else{(.)}  }
}

# Format output

format_output <- function(.optimised, .is_complete=TRUE){

  if (!.is_complete) {

    .optimised %>%

      select(node = ancestor, signature)

  } else {

    .optimised

  }


}

# Benchmark evaluation
evaluation <- function(.signature_df, .stream, .markers,
                       .imputed_counts, .sample, .symbol, .cell_type,
                       .mixture, 
                       .reduction_method, .dims=2, .tree){
  
  .stream = enquo(.stream)
  .markers = enquo(.markers)
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)

  
  .signature_df %>% 
    
    # silhouette evaluation
    mutate(silhouette = map(
      !!.markers, 
      ~ silhouette_evaluation(
        .markers = .x,
        .reduction_method = .reduction_method,
        .tree = .tree,
        .imputed_counts = .imputed_counts,
        .sample = !!.sample,
        .symbol = !!.symbol,
        .cell_type = !!.cell_type,
        .dims = .dims)
    )) %>% 
    
    mutate(avg_silhouette = map_dbl(silhouette, ~ mean(.x$sil_width))) %>%
    
    mutate(silhouette = map(
      silhouette, 
      ~ .x %>% 
        group_by(!!.cell_type) %>% 
        summarise(cluster_silhouette = mean(sil_width), cluster_size = n()) %>% 
        distinct() %>% 
        ungroup()
    )) %>% 
    
    # deconvolution evaluation
    # for each mixture, combine with the signatures from all methods
    expand_grid(.mixture, .) %>% 
    
    mutate(deconvolution = map2(
      !!.markers, mix, 
      ~ deconvolution_evaluation(
        .markers = .x, 
        .mixture = .y, 
        .tree = .tree,
        .imputed_counts = .imputed_counts,
        .sample = !!.sample,
        .symbol = !!.symbol,
        .cell_type = !!.cell_type)
    )) %>% 
    
    # mse by method
    mutate(MSE = map_dbl(
      deconvolution,
      ~ mean((.x$estimated_proportion - .x$proportion)^2)
    )) %>% 
    nest(data=-!!.stream) %>% 
    # mutate(median_MSE_over_mixes = map_dbl(data, ~ median(.x$MSE))) %>%
    mutate(mean_MSE_over_mixes = map_dbl(data, ~ mean(.x$MSE))) %>%
    unnest(data) %>% 
    
    # mse by cell type
    unnest(deconvolution) %>% 
    mutate(squared_error = (estimated_proportion - proportion)^2) %>% 
    nest(data = -c(!!.stream, !!.cell_type)) %>% 
    mutate(mean_MSE_for_cell_type = map_dbl(data, ~ mean(.x$squared_error))) %>% 
    unnest(data) %>% 
    
    select(-c(!!.markers, mixture_ID, mix, replicate, estimated_proportion, proportion, squared_error))
  
}

silhouette_evaluation <- function(.markers, .reduction_method, .dims=2, .tree, 
                                  .imputed_counts, .sample, .symbol, .cell_type){

  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)

  # modified silhouette_function and silhouette_score for evaluation

  silhouette_info <- function(.reduced_dimensions, .distance, .level, .cell_type){
    
    .cell_type = enquo(.cell_type)

    silhouette <- .reduced_dimensions %>%

      pull(!!as.symbol(.level)) %>%

      as.factor() %>%

      as.numeric() %>%

      silhouette(.distance)
    

    tibble(cluster = silhouette[, "cluster"],
           neighbor = silhouette[, "neighbor"],
           sil_width= silhouette[, "sil_width"]
           ) %>%

      mutate(!!.cell_type := .reduced_dimensions %>%
               pull(!!as.symbol(.level)) %>%
               as.factor())

  }

  silhouette_function <- function(.selected, .sample, .symbol, .reduction_method, .dims, .cell_type){

    .sample = enquo(.sample)
    .symbol = enquo(.symbol)
    .cell_type = enquo(.cell_type)

    .selected %>%

      # reduce dimensions
      mutate(reduced_dimensions = map2(
        markers, level,
        ~ dimension_reduction(.markers=.x, .level=.y,
                              .sample=!!.sample,
                              .symbol=!!.symbol,
                              .reduction_method=.reduction_method,
                              .dims=.dims)
      )) %>%

      # calculate distance matrix using PC1 & PC2
      mutate(distance = map(
        reduced_dimensions,
        ~ distance_matrix(.x, .reduction_method=.reduction_method)
      )) %>%

      # calculate silhouette score
      mutate(silhouette = pmap(
        list(reduced_dimensions, distance, level),
        ~ silhouette_info(.reduced_dimensions=..1, .distance=..2,  .level=..3, .cell_type=!!.cell_type)
      )) %>%

      # remove unnecessary columns
      select(-c(markers, distance))

  }

  .imputed_counts %>%
    mutate(level = "root") %>% 
    mutate(root = !!.cell_type) %>% 
    filter(!!.cell_type %in% as.phylo(.tree)$tip.label) %>% 
    filter(!!.symbol %in% .markers) %>%
    nest(markers = - level) %>%
    # calculate silhouette score for all signatures combined in each method
    silhouette_function(.sample=!!.sample, .symbol=!!.symbol, .cell_type=!!.cell_type,
                        .reduction_method = .reduction_method, .dims=.dims) %>%
    select(silhouette) %>%
    unnest(silhouette)
}

deconvolution_evaluation <- function(.markers, .mixture, .tree, 
                                     .imputed_counts, .sample, .symbol, .cell_type){

  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)

  # filter out data for signature genes
  reference <- .imputed_counts %>%
    
    filter(!!.cell_type %in% as.phylo(.tree)$tip.label) %>% 
    filter(!!.symbol %in% .markers) %>%

    # reshape the input matrix for deconvolve_cellularity():
    select(!!.symbol, !!.cell_type, !!.sample, count_scaled) %>%
    group_by(!!.symbol, !!.cell_type) %>%
    summarise(count_scaled_median = median(count_scaled)) %>%
    ungroup() %>%
    pivot_wider(id_cols = !!.symbol, names_from = !!.cell_type, values_from = count_scaled_median) %>%
    # must be a matrix
    tidybulk::as_matrix(rownames = !!.symbol)

  tidybulk::deconvolve_cellularity(
    .data = .mixture,
    .sample = replicate,
    .transcript = symbol, # the column in the mixture is called symbol not provided by the pipeline
    .abundance = count_mix,
    reference = reference,
    method = "llsr",
    prefix = "llsr_",
    action = "get",
    intercept = FALSE) %>%

    pivot_longer(cols=starts_with("llsr_"),
                 names_prefix ="llsr_",
                 names_to=quo_name(.cell_type),
                 values_to="estimated_proportion") %>%
    
    # join by the true proportion
    left_join(.mixture %>%
                unnest(data_samples) %>%
                distinct(replicate, !!.cell_type, proportion),
              by = c(quo_name(.cell_type), "replicate")
              )

}

# subset tree

# adapt_tree modifies the cell_type in the data according to the tree structure specified
adapt_tree <- function(.data, .tree, .node=NULL){
  
  # if no .node of interest is given, the original dataframe will be returned
  if (is.null(.node)){.data}

  # if subtree of interest is given, cell types in data will be modified accordingly
  else{

    subtree <- tree_subset(.tree = .tree, .node = .node)

    .tree %>%
      ToDataFrameTypeColFull(fill=NA) %>%

      # find ancestors of each cell type in the full tree
      unite("ancestors", contains("level"), sep=".", na.rm = TRUE) %>%
      mutate(ancestors = map2(
        ancestors, cell_type,
        ~.x %>%
          paste("Tissue", ., sep=".") %>%
          str_split("\\.") %>%
          unlist() %>%
          .[. != .y] %>%
          # arrange the ancestors from most recent to the furthest
          rev()
      )) %>%

      # get all the cell types in the subtree of interest
      mutate(cell_type_subtree = subtree$Get("level") %>% names() %>% list()) %>%

      # if the cell type in the full tree is present in the subtree, leave it as it is;
      # otherwise find its most recent ancestor present in the subtree
      mutate(cell_type_to_be = pmap_chr(
        list(cell_type, ancestors, cell_type_subtree),
        ~ if (! ..1 %in% ..3){
          ..2[which(..2 %in% ..3)[1]]
        } else {..1}
      )) %>%
      select(cell_type, cell_type_to_be) %>%

      left_join(.data, ., by = "cell_type") %>%

      # remove the old cell type
      select(-cell_type) %>%

      # use cell type to be as the new cell type
      dplyr::rename(cell_type = cell_type_to_be)

  }

}

tree_subset <- function(.tree, .node=NULL) {

  if(is.null(.node)){.tree}

  else{

    tree_subset_tibble(.tree = .tree, .node = .node) %>%

      # reverse the row order of tibble to make the root node on top
      # .[rev(1: nrow(.)), ] %>%

      update_tree_edge_matrix() %>%

      as.phylo() %>%

      as.Node(replaceUnderscores=FALSE)

  }

}

tree_subset_tibble <- function(.tree, .node){

  # Args:
  # .tree: a Node R6 object from data.tree
  # .node: a character specifying the cell type of interest

  tree_tbl <- .tree %>% as.phylo %>% as_tibble

  if (.node == "Tissue") { return( tree_tbl %>% filter(label == "Tissue") )

  } else {

    return(sibling(tree_tbl, .node) %>%
             bind_rows(
               tree_subset_tibble(.tree, parent(tree_tbl, .node) %>% .$label),
               .)
    )
  }

}

update_tree_edge_matrix <- function(.tree_tbl) {

  # Args:
  # .tree_tbl: a tree_tibble object from tidytree

  tree_phylo <- .tree_tbl %>% as.phylo

  .tree_tbl %>%

    mutate(node = map_int(label,
                          ~ if(.x %in% tree_phylo$tip.label){
                            which(.x == tree_phylo$tip.label)
                          } else {
                            which(.x == tree_phylo$node.label) + Ntip(tree_phylo)
                          }
    )) %>%

    mutate(parent = rep(parent %>% unique() %>% rank(),
                        times = table(parent)) + Ntip(tree_phylo)
    )

}

# Shiny App functions
# get leaf nodes at a specific level
get_leaf_nodes_at_a_level <- function(.tree, .level){
  
  if(.level == "root"){
    .tree %>% as.phylo %>% .$tip.label
  }else{
    L  = .level %>% str_split("_") %>% {as.numeric(.[[1]][2])+1}
    
    Clone(.tree, pruneFun = function(node) node$level <= L) %>% 
      as.phylo %>% 
      .$tip.label
  }
}

get_target_cell_markers <- function(.data, .k=6){
  
  .data %>%
    unnest(data) %>% 
    unnest(children) %>% 
    mutate(contrast = str_extract(contrast, ".*(?=\\s\\-)")) %>% 
    rename(target = contrast) %>% 
    select(target, enriched) %>% 
    unnest(enriched) %>% 
    distinct(target, symbol, rank) %>% 
    nest(enriched = -target) %>% 
    mutate(enriched = map_chr(enriched, ~ .x %>% 
                                arrange(rank) %>% 
                                pull(symbol) %>% 
                                vector_to_formatted_text(.k=.k)
    ))
}

vector_to_formatted_text <- function(.vector, .k=6){
  
  L <- length(.vector)
  
  if(L <= .k){text = paste(.vector, collapse = " ")}
  else{
    
    text=c()
    
    for (i in seq(1, L-.k, by = .k)) {
      text = text %>% append(paste(.vector[i:(i+.k-1)], collapse = " "))
    }
    
    text = text %>% 
      append(paste(.vector[(L%/%.k*.k+1):L], collapse = " ")) %>% 
      paste(collapse = "\n")
  }
  
  return(text)
}


# old preprocessing where data rectangularisation(same set of genes for all cell types), 
# imputation, scale_abundance was done
# preprocess <- function(.transcript, .level) {
#   
#   # load data
#   .transcript %>%
#     
#     tidybulk(sample, symbol, count) %>%
#     
#     # aggregate duplicate sample/gene pairs in the data
#     # aggregate_duplicates(sample, symbol, count) %>%
# 
#     # rectangularise data
#     nest(data = -c(symbol, cell_type)) %>%
#     add_count(symbol) %>%
#     filter(n == max(n)) %>%
#     unnest(data) %>%
# 
#     # Imputation of missing data
#     impute_missing_abundance(~ cell_type) %>%
# 
#     # scale counts
#     identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
#     scale_abundance() %>%
#     
#     # filter for cells at the level of interest. .level == level_1
#     filter(is.na(!!as.symbol(.level))==FALSE) %>%
#     
#     # nest by ancestor
#     nest(data = - !!as.symbol(pre(.level)))
#   
# }

# old rank_bayes function that uses a bayes dataframe with only symbols and quantiles without imputation
# rank_bayes <- function(.preprocessed, .contrast_method, .bayes, .rank_stat=NULL){
#  
#   .preprocessed %>% 
#     
#     unnest(tt) %>% 
#     
#     mutate(contrast = map2(data, level, ~ .contrast_method(.x, .y))) %>% 
#     
#     unnest(contrast) %>% 
#     mutate(contrast = str_remove_all(contrast, glue("{level}"))) %>% 
#     
#     mutate(lower_quantile = map(
#       contrast,
#       ~ .bayes %>% 
#         filter(cell_type == str_extract(.x, ".*(?=\\s\\-)")) %>% 
#         select(symbol, lower_quantile='25%')
#       # arrange(symbol)
#     )) %>% 
#     
#     mutate(mean_upper_quantile = map(
#       contrast,
#       ~ {
#         background <- (.x) %>% 
#           str_extract("(?<=\\-\\s).*") %>% 
#           str_split("\\+") %>% 
#           unlist() %>% 
#           str_remove_all("(?<=\\/).*|\\W")
#         
#         .bayes %>% 
#           # calculate the mean 75% quantile of each gene over all background cell types
#           filter(cell_type %in% background) %>% 
#           group_by(symbol) %>% 
#           summarise(symbol, mean_upper_quantile = mean(`75%`)) %>% 
#           distinct() %>% 
#           ungroup()
#       }
#     )) %>% 
#     
#     mutate(stat_df = map2(
#       lower_quantile, mean_upper_quantile,
#       ~ inner_join(.x, .y, by= "symbol")
#     )) %>% 
#     select(-c(lower_quantile, mean_upper_quantile)) %>% 
#     
#     mutate(stat_df = map(
#       stat_df,
#       ~ .x %>% 
#         mutate(difference = lower_quantile - mean_upper_quantile) %>% 
#         arrange(desc(difference))
#     ))
#     
#     # nest(markers = -c(level, ancestor, data))
#   
# }