# devtools::install_github("stemangiola/nanny@convert-to-S3", force = TRUE)
# devtools::install_github("stemangiola/tidybulk@dev", force = TRUE)



# library(plotly)
# library(future)
# library(furrr)
# library(proxy)
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

main <- function(.transcript, .is_hierarchy=TRUE, .level=NULL, 
                 .contrast_method, .ranking_method, .rank_stat=NULL, .bayes=NULL, 
                 .selection_method, .kmax=60, .discard_number=2000, .reduction_method = "PCA",
                 .optimisation_method, .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
                 .is_complete = FALSE) {
  
  
  if ( (!.is_hierarchy)|(.selection_method == "naive")) {
    
    .transcript %>%
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      # scale_input_counts(.is_hierarchy=.is_hierarchy, .level=.level) %>% 
      
      # Input: data.frame columns_1 <int> | ...
      do_ranking(.ranking_method=.ranking_method, .contrast_method=.contrast_method, .rank_stat=.rank_stat, .bayes=.bayes) %>%  
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      do_selection(.selection_method=.selection_method, .reduction_method=.reduction_method, .discard_number=.discard_number, .kmax=.kmax) %>% 
      
      do_optimisation(.optimisation_method=.optimisation_method) %>% 
      
      format_output(.is_complete = .is_complete)
    
  } else {
    
    .transcript %>%
      
      # Input: data.frame columns_1 <int> | ...
      # Output: 
      # scale_input_counts(.is_hierarchy=.is_hierarchy, .level=.level) %>%
      
      do_ranking(.ranking_method=.ranking_method, .contrast_method=.contrast_method, .rank_stat=.rank_stat, .bayes=.bayes) %>% 
      
      mutate(level.copy = level) %>% 
      nest(data = -level.copy) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          
          do_selection(.selection_method=.selection_method, .reduction_method=.reduction_method, .discard_number=.discard_number, .kmax=.kmax) %>%
          
          do_optimisation(.optimisation_method = .optimisation_method) %>%
          
          format_output(.is_complete = .is_complete)
        
      )) %>% 
      
      unnest(data) %>% 
      select(-level.copy)
    
  }
}

# Input scale abundance

scale_input_counts <- function(.transcript, .is_hierarchy = TRUE, .level = NULL){
  
  if (.is_hierarchy) { # scaling under each ancestor node at each level
    
    if (is.null(.level) == TRUE) { # scale counts for all levels present
      
      tt_hierarchy <- 
        
        tibble(level = names(.transcript) %>% str_subset("level")) %>% 
        
        mutate(tt = map(level, ~ .transcript %>% 
                          
                          mutate(level_0 = "root") %>% 
                          
                          preprocess(.x))) %>% 
        
        mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename("ancestor" = pre(.y))))
      
    } else { # scale counts for the level specified by .level
      
      tt_hierarchy <- 
        
        tibble(level = .level) %>% 
        
        mutate(tt = map(level, ~ .transcript %>% 
                          
                          mutate(level_0 = "root") %>% 
                          
                          preprocess(.x))) %>% 
        
        mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename("ancestor" = pre(.y))))
    }
    
    return(tt_hierarchy)
    
  } else { # non-hierarchical: scaling all cell types under the root node
    
    .level <- "root"
    
    tt_non_hierarchy <- 
      
      tibble(level = .level) %>% 
      
      mutate(tt = map(level, ~ .transcript %>%
                        
                        # create a root column for pre(.level)
                        mutate(!!as.symbol(.level) := cell_type) %>% 
                        mutate(!!as.symbol(pre(.level)) := .level) %>% 
                        
                        preprocess(.x))) %>% 
      
      mutate(tt = map2(tt, level, ~ .x %>% dplyr::rename("ancestor" = pre(.y))))
    
    return(tt_non_hierarchy)
  }
}

preprocess <- function(.transcript, .level) {
  
  # load data
  .transcript %>%
    
    tidybulk(sample, symbol, count) %>%
    
    # aggregate duplicate sample/gene pairs in the data
    # aggregate_duplicates(sample, symbol, count) %>% 
    
    # rectangularise data
    nest(data = -c(symbol, cell_type)) %>%
    add_count(symbol) %>%
    filter(n == max(n)) %>%
    unnest(data) %>% 
    
    # Imputation of missing data
    impute_missing_abundance(~ cell_type) %>%
    
    # scale counts
    identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
    scale_abundance() %>% 
    
    # filter for cells at the level of interest
    filter(is.na(!!as.symbol(.level))==FALSE) %>%
    
    # nest by ancestor
    nest(data = - !!as.symbol(pre(.level)))
  
  # # scale count for further analysis
  # mutate(data=map(data, ~ .x %>%
  #                   identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
  #                   scale_abundance()
  # ))
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

do_ranking <- function(.preprocessed, .ranking_method, .contrast_method, .rank_stat, 
                       .bayes){
  
  # rank_stat takes either "Pvalue" or "logFC" 
  .preprocessed %>%
    
    # .bayes = .cellsig_theoretical_transcript_abundace_distribution
    .ranking_method(.contrast_method=.contrast_method, .rank_stat=.rank_stat, .bayes=.bayes) %>% 
    
    # filter out potential nodes in which no genes are considered significant by rank_by_logFC
    filter(map_int(markers, ~ .x %>% unnest(stat_df) %>% nrow()) != 0)
  
}

rank_edgR_quasi_likelihood <- function(.preprocessed, .contrast_method, .rank_stat, .bayes=NULL){
  
  .preprocessed %>%
    unnest(tt) %>% 
    
    # Differential transcription: generate contrast
    mutate(markers = map2(
      data, level,
      ~ .x %>% 
        test_differential_abundance(
          as.formula(sprintf("~ 0 + %s", .y)),
          .contrasts = .contrast_method(.x, .y),
          method = "edgeR_quasi_likelihood",
          action="only") 
    )) %>% 
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ rank_by_stat(.x, .rank_stat) )) %>% 
    
    # remove prefixes from contrast expressions
    mutate(markers = map2(
      markers, level,
      ~ .x %>% 
        mutate(contrast = map2_chr(contrast, .y, ~ str_replace_all(.x, .y, "")))
    ))
}

rank_edgR_robust_likelihood_ratio <- function(.preprocessed, .contrast_method, .rank_stat="PValue",
                                              .bayes=NULL){
  
  .preprocessed %>%
    unnest(tt) %>% 
    
    # Differential transcription: generate contrast
    mutate(markers = map2(
      data, level,
      ~ .x %>% 
        test_differential_abundance(
          as.formula(sprintf("~ 0 + %s", .y)),
          .contrasts = .contrast_method(.x, .y),
          method = "edger_robust_likelihood_ratio",
          test_above_log2_fold_change = 1,
          action="only") 
    )) %>% 
    
    # Select markers from each contrast by rank of Pvalue
    mutate(markers = map(markers, ~ rank_by_stat(.x, "PValue") )) %>% 
    
    # remove prefixes from contrast expressions
    mutate(markers = map2(
      markers, level,
      ~ .x %>% 
        mutate(contrast = map2_chr(contrast, .y, ~ str_replace_all(.x, .y, "")))
    ))
}

rank_by_stat <-  function(.markers, .rank_stat){
  
  # .rank_stat = enquo(.rank_stat)
  
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
    
    # Filter out insignificant genes and rank the significant ones
    
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


rank_bayes <- function(.preprocessed, .contrast_method, .bayes, .rank_stat=NULL){
  
  .preprocessed %>% 
    
    unnest(tt) %>% 
    
    mutate(descendants = map2(data, level, ~ .contrast_method(.x, .y))) %>% 
    
    unnest(descendants) %>% 
    mutate(descendants = str_remove_all(descendants, "level_\\d")) %>% 
    separate(descendants, into = c("target", "background"), sep = " - ") %>% 
    
    mutate(background = str_split(background, "\\+")) %>% 
    mutate(background = map(background, ~ .x %>% str_remove_all("\\W|\\d$"))) %>% 
    
    mutate(lower_quantile = map(
      target,
      ~ .bayes %>% 
        filter(cell_type == .x) %>% 
        select(symbol, lower_quantile='25%') %>% 
        arrange(symbol)
    )) %>% 
    
    mutate(mean_upper_quantile = map(
      background,
      ~ .bayes %>% 
        # calculate the mean 75% quantile of each gene over all background cell types
        filter(cell_type %in% .x) %>% 
        group_by(symbol) %>% 
        summarise(symbol, mean_upper_quantile = mean(`75%`)) %>% 
        distinct() %>% 
        ungroup()
    )) %>% 
    
    mutate(stat_df = map2(
      lower_quantile, mean_upper_quantile,
      ~ inner_join(.x, .y, by= "symbol")
    )) %>% 
    select(-c(lower_quantile, mean_upper_quantile)) %>% 
    
    mutate(stat_df = map(
      stat_df,
      ~ .x %>% 
        mutate(difference = lower_quantile - mean_upper_quantile) %>% 
        arrange(desc(difference))
    )) %>% 
    
    nest(markers = -c(level, ancestor, data))
  
}

# Selection

do_selection <- 
  function(.ranked, .selection_method, .reduction_method, .kmax, .discard_number) {
    
    # .k is the number of genes selected from each cell_type contrast
    
    if (.selection_method == "naive") {
      
      .ranked %>% 
        
        do_naive_selection(.kmax=.kmax, .reduction_method=.reduction_method)
      
    } else {
      
      .ranked %>% 
        
        single_marker_pw_selection_using_silhouette(.discard_number=.discard_number, .reduction_method=.reduction_method)
      
    }
    
  }

## Naive selection

### calculate silhouette score for a series of sig_sizes

do_naive_selection <- function(.ranked, .kmax, .reduction_method) {
  
  # Args:
  # .ranked: output from do_ranking
  # .kmax: maximum number of markers selected from each cell type contrast
  # .reduction_method: method used to reduce dimensions such as "PCA", "tSNE", "MSA"
  
  tibble(number_of_markers_from_each_contrast = 1: .kmax) %>% 
    
    # select signature and calculate silhouette score 
    mutate(data = map(
      number_of_markers_from_each_contrast,
      ~ naive_selection(.ranked=.ranked, .x) %>% 
        silhouette_function(.reduction_method=.reduction_method)
    )) %>% 
    
    # nest by ancestor nodes/cell types
    unnest(data) %>%
    nest(data = - c(level, ancestor))
  
}

naive_selection <- function(.ranked, .k) {
  
  # Args:
  # .ranked: output from do_ranking()
  # .k: the number of genes selected from each cell_type contrast
  
  .ranked %>% 
    
    # selection markers from each contrast
    mutate(markers = map(
      markers,
      ~ .x %>% 
        mutate(stat_df = map(stat_df, ~ .x %>% dplyr::slice(1: .k))) %>% 
        unnest(stat_df)
    )) %>% 
    
    # Add original data info to the markers selected, 
    # use inner_join to ensure symbols are present in both markers and data
    mutate(markers = map2(markers, data, ~ inner_join(.x, .y, by="symbol"))) %>%
    
    # remove unnecessary column
    select(-data) %>% 
    
    # collect from which contrasts signature genes are extracted
    # mutate(contrast = map(markers, ~ .x %>% distinct(contrast, symbol))) %>% 
    
    # collect signature genes selected
    mutate(signature = map(markers, ~ .x$symbol %>% unique())) %>% 
    
    # number of the signature genes
    mutate(real_size = map_int(signature, ~ length(.x)))
  
}

### calculate silhouette score

silhouette_function <- function(.selected, .reduction_method){
  
  .selected %>% 
    
    # reduce dimensions
    mutate(reduced_dimensions = map2(
      markers, level, 
      ~ dimension_reduction(.x, .y, .reduction_method=.reduction_method)
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

dimension_reduction <- function(.markers, .level, .reduction_method) {
  
  .markers %>% 
    
    distinct(sample, symbol, count_scaled, !!as.symbol(.level)) %>% 
    
    reduce_dimensions(sample, symbol, count_scaled, 
                      action = "get",
                      method = .reduction_method,
                      # .dims = 2,
                      transform = log1p,
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

single_marker_pw_selection_using_silhouette <- 
  function(.ranked, .discard_number, .reduction_method) {
    
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
      winning_contrast = list(),
      signature = list(),
      # reduced_dimensions = list(),
      silhouette = double()
    )
    
    # set the base markers
    contrast_pair_tb0 <- 
      
      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      .ranked %>%
      
      # select top 1 markers from each contrast
      naive_selection(1) %>%
      
      dplyr::rename("new_challengers" = "signature") %>% 
      
      mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
      
      mutate(winning_contrast = map(
        markers,
        ~ .x %>% pull(contrast) %>% unique()
        # distinct(contrast, symbol) %>% 
        # mutate(contrast = contrast %>% str_extract(".*(?=\\s\\-)")) %>% 
        # mutate(contrast_symbol = map2_chr(contrast, symbol, ~ paste(.x, .y, sep = "."))) %>% 
        # pull(contrast_symbol)
      )) %>% 
      
      mutate(signature = winner) %>% 
      
      silhouette_function(.reduction_method=.reduction_method) %>% 
      
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
          filter(!symbol %in% with(signature, signature[ancestor==.y][[1]])) %>% 
          nest(stat_df = - contrast)
      ))
    
    # counter for number of genes discarded
    j <- map_int(signature$signature, ~ length(.x))
    
    # count the number of iterations
    i <- 0
    while (any(j < .discard_number) &
           # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
           all(map_int(ranked_copy$markers, 
                       # hence the boundary should be the number of satisfactory genes selected
                       ~ .x %>% unnest(stat_df) %>% nrow()) > 0)) {
      
      contrast_pair_tb <- 
        
        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        ranked_copy %>% 
        
        # select top 1 markers from each contrast, ignore the signature output
        naive_selection(1) %>% 
        
        # pick the one new challenger from each contrast
        mutate(markers = map(
          markers,
          ~ .x %>% 
            nest(new_challenger = - contrast) %>% 
            mutate(new_challenger = map_chr(new_challenger, ~.x %>% distinct(symbol) %>% pull()))
        )) %>% 
        unnest(markers) %>% 
        select(-c(signature, real_size)) %>% 
        
        # append the new challenger from each contrast to the base markers for that ancestor node
        mutate(challengers_for_silhouette = map2(
          new_challenger, ancestor, 
          ~ with(signature, signature[ancestor==.y][[1]]) %>% 
            append(.x)
        )) %>% 
        
        # calculate silhouette score for the challengers from each contrast
        mutate(silhouette = map2_dbl(
          challengers_for_silhouette, ancestor, 
          ~ silhouette_for_markers(.ranked=.ranked, .signature=.x, .ancestor=.y, .reduction_method=.reduction_method) %>% 
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
        mutate(winning_contrast = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ]$contrast # %>% str_extract(".*(?=\\s\\-)")
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
              filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]])) %>% 
              nest(stat_df = -contrast)
          } else {
            .x %>% 
              unnest(stat_df) %>% 
              filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]])) %>% 
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
      
      i <- i + 1
      cat("iteration: ", i, "\n")
      
    }
    
    # format output for optimisation
    output <- summary_tb %>% 
      mutate(real_size = map_int(signature, ~ length(.x))) %>% 
      nest(data = - c(level, ancestor))
    
    return(output)
  }

silhouette_for_markers <-function(.ranked, .signature, .ancestor, .reduction_method) {
  
  .ranked %>%
    
    filter(ancestor == .ancestor) %>% 
    
    select(-markers) %>%
    
    # filter markers that are in the signature
    mutate(data = map(data, ~.x %>% 
                        filter(symbol %in% .signature))) %>% 
    
    # format input
    dplyr::rename("markers" = "data") %>% 
    
    silhouette_function(.reduction_method = .reduction_method)
  
}

# Optimisation

do_optimisation <- function(.selected, 
                            .optimisation_method, 
                            .penalty_rate=0.2,
                            .kernel = "normal", 
                            .bandwidth = 0.05, 
                            .gridsize = 100){
  
  if(.optimisation_method == "penalty") {
    
    .selected %>% 
      
      mutate(optimal_size = map_int(data, ~ penalised_silhouette(.x, .penalty_rate=.penalty_rate))) %>% 
      
      unnest(data) %>% 
      
      filter(real_size == optimal_size)
    
  } else if (.optimisation_method == "curvature") {
    
    .selected %>% 
      
      curvature_of_kernel_smoothed_trend(.kernel=.kernel, .bandwidth=.bandwidth, .gridsize=.gridsize) %>% 
      
      unnest(data) %>% 
      
      filter(real_size <= optimal_size) %>% 
      
      select(-c(size.rescaled, smoothed)) %>% 
      
      nest(data = -c(level, ancestor, signature)) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          unnest(winner, winning_contrast) %>% 
          nest(enriched = -winning_contrast) %>% 
          mutate(enriched = map(enriched, ~ .x %>% pull(winner) %>% unique()))
      ))
  }
  
}

curvature <- function(.drv1, .drv2){
  abs(.drv2) / (1 + .drv1^2)^(3/2)
}

curvature_of_kernel_smoothed_trend <- function(.plot_data, 
                                               .kernel = "normal", 
                                               .bandwidth = 0.05, 
                                               .gridsize = 100){
  .plot_data %>% 
    
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
    
    mutate(optimal_size = ifelse(optimal_size<10, 10, optimal_size))
  
}

penalised_silhouette <- function(.plot_data, .penalty_rate=0.2) {
  
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    # filter(real_size > 10) %>%
    
    # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
    mutate(size_rescaled = rescale(real_size)) %>% 
    
    mutate(penalised_silhouette = silhouette - .penalty_rate * size_rescaled) %>% 
    
    filter(penalised_silhouette == max(penalised_silhouette)) %>% 
    
    pull(real_size)
}

# Format output

format_output <- function(.optimised, .is_complete=FALSE){
  
  if (!.is_complete) {
    
    .optimised %>% 
      
      select(node = ancestor, signature)
    
  } else {
    
    .optimised
    
  }
  
  
}

