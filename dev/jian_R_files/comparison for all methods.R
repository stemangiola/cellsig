
library(tidyverse)
library(ggplot2)
library(tidybulk)
library(cluster)
library(proxy)
library(factoextra)
library(stringr)
library(scales)
library(cellsig)
library(tidyHeatmap)

# Load data
load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")

# All functions

## Functions for hierarchy analysis

## 1 preprocess data

### 1.1 string manipulation that converts level of interest (e.g "level_5") to its ancestor level (e.g "level_4")
pre <- function(.level) {
  .level %>% 
    str_split("_") %>% 
    {as.numeric(.[[1]][2])-1} %>% 
    paste("level", ., sep = "_")
}

### 1.2 preprocess

preprocess <- function(.data, .level) {
  
  # load data
  .data %>%
    
    tidybulk(sample, symbol, count) %>%
    
    # filter for the cell types of interest for gene marker selection
    filter(is.na(!!as.symbol(.level))==F) %>%
    
    # Imputation of missing data within each level_5
    # impute_missing_abundance(~ !!as.symbol(.level)) %>%
    
    # Group by ancestor
    nest(data = - !!as.symbol(pre(.level))) %>%
    
    # Eliminate genes that are present in some but all cell types
    # (can be still present in a subset of samples from each cell-type)
    mutate(data = map(
      data,
      ~ .x %>%
        nest(data = -c(symbol, !!as.symbol(.level))) %>%
        add_count(symbol) %>%
        filter(n == max(n)) %>%
        unnest(data)
    )) %>%
    
    # Imputation of missing data within each level_5
    mutate(data = map(data, ~ .x %>% impute_missing_abundance(~ !!as.symbol(.level)))) %>%
    
    # scale count for further analysis
    mutate(data=map(data, ~ .x %>%
                      identify_abundant(factor_of_interest = !!as.symbol(.level)) %>%
                      scale_abundance()
    ))
}

## 2 contrast functions

### 2.1 pairwise comparisons
get_contrasts_from_df = function(.data, .level){
  
  .data %>% 
    distinct(!!as.symbol(.level)) %>% 
    mutate(!!as.symbol(.level) := paste0(.level, !!as.symbol(.level))) %>% 
    
    # Permute
    mutate(cell_type2 := !!as.symbol(.level)) %>% 
    expand(!!as.symbol(.level), cell_type2) %>% 
    filter(!!as.symbol(.level) != cell_type2) %>% 
    
    # Create contrasts
    mutate(contrast = sprintf("%s - %s", !!as.symbol(.level), cell_type2)) %>%
    pull(contrast)
  
}

### 2.2 create a contrast vector for limma::makeContrasts() or tidybulk::test_differential abundance()

mean_contrast <- function(.data, .level){
  
  # find all cell types
  cell_types <- .data %>% 
    distinct(!!as.symbol(.level)) %>% 
    pull() %>% 
    as.vector()
  
  # format cell_types with prefix
  cell_types <- paste0(.level, cell_types)
  
  # initialise a vector called contrasts
  contrasts <- 1: length(cell_types)
  
  # create all contrasts and store them in contrasts
  for(i in 1:length(cell_types) ){
    background = paste(cell_types[-i], collapse = "+")
    divisor = length(cell_types[-i])
    contrasts[i] <- sprintf("%s - (%s)/%s", cell_types[i], background, divisor)
  }
  
  return(contrasts)
}

## 3 marker ranking & selection

select_markers_for_each_contrast = function(.markers, .sig_size){
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
    
    # Rank
    mutate(stat_df = map(stat_df, ~.x %>%
                           filter(FDR < 0.05 & logFC > 2) %>%
                           filter(logCPM > mean(logCPM)) %>%
                           arrange(logFC %>% desc()) %>%
                           slice(1: .sig_size)
                         
    )) %>%
    
    unnest(stat_df)
}

## 4 marker collection for each contrast

contrast_PW <- function(.tt, .level){
  .tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + !!as.symbol(.level), 
                                    .contrasts = get_contrasts_from_df(.x, .level),
                                    action="only") 
    ))
}

contrast_MC <- function(.tt, .level){
  .tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + !!as.symbol(.level), 
                                    .contrasts = mean_contrast(.x, .level),
                                    action="only") 
    ))
}

sig_select <- function(.contrast, .level, .sig_size) {
  .contrast %>% 
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, .sig_size))) %>%
    
    # Add original data info to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(!!as.symbol(pre(.level)), markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, .level, "") %>% str_replace(.level, ""))
}


## 6 Silhouette score calculation for nodes in a level using .method of dimension reduction
sil_func <- function(.markers, .level, .method){
  .markers %>%
    nest(rdim = - !!as.symbol(pre(.level))) %>%
    mutate(rdim = map(rdim, ~ .x %>%
                        distinct(sample, symbol, count_scaled, !!as.symbol(.level)))) %>%
    mutate(rdim = map(rdim, ~ .x %>%
                        reduce_dimensions(sample, symbol, count_scaled,
                                          method = .method,
                                          action = "add",
                                          transform = log1p,
                                          # check_duplicates is for Rtsne method
                                          check_duplicates = FALSE) %>% 
                        
                        # save symbols for calculating real_size while reducing replicated rows resulted from symbol
                        nest(data_symbol = c(symbol, count_scaled))
    )) %>%
    
    # calculate the dissimilarity matrix with PC values
    mutate(distance = map(rdim, ~ .x %>%
                            select(contains(str_sub(.method, end = -2L))) %>%
                            factoextra::get_dist(method = "euclidean")
    )) %>%
    
    # calculate silhouette score
    mutate(sil = map2(rdim, distance,
                      ~ silhouette(as.numeric(as.factor(`$`(.x, !!as.symbol(.level)))), .y)
    )) %>%
    mutate(sil = map(sil, ~ .x %>% summary())) %>%
    mutate(sil = map(sil, ~ .x %>% 
                       `$`(avg.width) ))%>% 
    mutate(sil = unlist(sil)) %>% 
    
    # obtain the actual number of signature genes
    mutate(real_size = map_int(rdim, ~ .x$data_symbol %>% 
                                 map_int(~ n_distinct(.x$symbol)) %>% 
                                 unlist() %>% 
                                 unique() ))
  
}

sil_tb <- function(.contrast, .level, .sig_size, .method) {
  tibble(sig_size = 1: .sig_size) %>% 
    
    # select signature genes for each sig_size at each level
    mutate(sil_df = map(sig_size, ~ sig_select(.contrast, .level, .x))) %>% 
    
    # calculate silhouette score for each ancestor cell type under sil_df
    mutate(sil_df = map(sil_df, ~ sil_func(.x, .level, .method)))
}


# Functions for no hierarchy analysis

## 2.1

mean_contrast0 <- function(.data){
  
  prefix <- "cell_type"
  
  # find all cell types
  cell_types <- .data %>% 
    distinct(cell_type) %>% 
    pull() %>% 
    as.vector()
  
  # format cell_types with prefix
  cell_types <- paste0(prefix, cell_types)
  
  # initialise a vector called contrasts
  contrasts <- 1: length(cell_types)
  
  # create all contrasts and store them in contrasts
  for(i in 1:length(cell_types) ){
    background = paste(cell_types[-i], collapse = "+")
    divisor = length(cell_types[-i])
    contrasts[i] <- sprintf("%s - (%s)/%s", cell_types[i], background, divisor)
  }
  
  return(contrasts)
}

get_contrasts_from_df0 = function(.data){
  
  .data %>% 
    distinct(cell_type) %>% 
    mutate(cell_type = paste0("cell_type",cell_type)) %>% 
    
    # Permute
    mutate(cell_type2 = cell_type) %>% 
    expand(cell_type, cell_type2) %>% 
    filter(cell_type != cell_type2) %>% 
    
    # Create contrasts
    mutate(contrast = sprintf("%s - %s", cell_type, cell_type2)) %>%
    pull(contrast)
  
}

## 2.3 marker collection for each contrast

# mean contrast method
contrast_MC0 <- function(.tt){
  .tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = mean_contrast0(.x),
                                    action="only")  ))
}

# pairwise contrast method
contrast_PW0 <- function(.tt){
  .tt %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = get_contrasts_from_df0(.x),
                                    action="only")  ))
}

# marker selection & processing
sig_select0 <- function(.contrast, .sig_size) {
  .contrast %>% 
    
    # Select markers from each contrast by rank of stats
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x, .sig_size))) %>%
    
    # Add original data info to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(level_0, markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
}

## 2.5 Silhouette function
sil_func0 <- function(.markers, .method){
  .markers %>% 
    nest(rdim = - level_0) %>% 
    mutate(rdim = map(rdim, ~ .x %>% 
                        distinct(sample, symbol, count_scaled, cell_type))) %>% 
    mutate(rdim = map(rdim, ~ .x %>%
                        reduce_dimensions(sample, symbol, count_scaled,
                                          method = .method,
                                          action = "add",
                                          transform = log1p,
                                          # check_duplicates is for Rtsne method
                                          check_duplicates = FALSE) %>% 
                        
                        # save symbols for calculating real_size while reducing replicated rows resulted from symbol
                        nest(data_symbol = c(symbol, count_scaled))
    )) %>%
    
    # calculate the dissimilarity matrix with PC values
    mutate(distance = map(rdim, ~ .x %>%
                            select(contains(str_sub(.method, end = -2L))) %>%
                            factoextra::get_dist(method = "euclidean")
    )) %>%
    
    # calculate silhouette score
    mutate(sil = map2(rdim, distance, 
                      ~ silhouette(as.numeric(as.factor(`$`(.x, cell_type))), .y)
    )) %>% 
    mutate(sil = map(sil, ~ .x %>% summary())) %>%
    mutate(sil = map(sil, ~ .x %>% `$`(avg.width) ))%>% 
    mutate(sil = unlist(sil)) %>% 
    
    # obtain the actual number of signature genes
    mutate(real_size = map_int(rdim, ~ .x$data_symbol %>% 
                                 map_int(~ n_distinct(.x$symbol)) %>% 
                                 unlist() %>% 
                                 unique() ))
  
}

sil_tb0 <- function(.contrast, .sig_size, .method) {
  tibble(sig_size = 1: .sig_size) %>% 
    
    # select signature genes for each sig_size at each level
    mutate(sil_df = map(sig_size, ~ sig_select0(.contrast, .x))) %>% 
    
    # calculate silhouette score for each ancestor cell type under sil_df
    mutate(sil_df = map(sil_df, ~ sil_func0(.x, .method)))
}

# Hierarchical analysis

## Preprocess data at ALL LEVELS EXCEPT FOR no hierarchy========================

contrast_all <- 
  # tibble(level = 1:5) %>% 
  # mutate(level = paste("level", level, sep = "_")) %>% 
  # 
  # # preprocess data
  # mutate(tt = map(level, ~ counts %>% 
  #                   mutate(level_0 = "cell") %>% 
  #                   preprocess(.x))) %>% 
  tt_all %>% 
  
  # generate contrast by pairwise comparison
  mutate(contrast_PW = map2(tt, level, ~ contrast_PW(.x, .y) )) %>% 
  
  # generate contrast by mean contrast method
  mutate(contrast_MC = map2(tt, level, ~ contrast_MC(.x, .y) ))

saveRDS(contrast_all, "contrast_all.rds")

contrast_all <- readRDS("dev/intermediate_data/contrast_all.rds")

# Hierarchy + Pairwise ============

sig_size <- 20

H_PW <- contrast_all %>%
  # slice(3) %>%
  
  # calculate silhouette score for all signature sizes at each level of interest
  mutate(sil = map2(contrast_PW, level, ~ sil_tb(.x, .y, sig_size, "tSNE"))) %>%
  
  # select what is needed to plot sil score for each ancestor node(type) changing over sig_size
  mutate(sil = map(sil, ~ .x %>% 
                     unnest(sil_df) %>% 
                     mutate(ancestor = select(., contains("level_")) %>% 
                              as_vector() ) %>% 
                     select(sig_size, real_size, sil, ancestor)
  )) %>% 
  select(sil) %>% 
  unnest(sil) %>% 
  mutate(analysis = "PW+H")


# Hierarchy + Mean contrast ============

sig_size <- 20

H_MC <- contrast_all %>%
  # slice(3) %>% 
  mutate(sil = map2(contrast_MC, level, ~ sil_tb(.x, .y, sig_size, "tSNE"))) %>% 
  mutate(sil = map(sil, ~ .x %>% 
                     unnest(sil_df) %>% 
                     mutate(ancestor = select(., contains("level_")) %>% 
                              as_vector() ) %>% 
                     select(sig_size, real_size, sil, ancestor)
  )) %>% 
  select(sil) %>% 
  unnest(sil) %>% 
  mutate(analysis = "MC+H")


# No hierarchy Analysis

## Preprocess for non-hierarchical analysis ====================================
tt_simple <- counts %>%
  
  mutate(level_0 = "cell") %>% 
  
  tidybulk(sample, symbol, count) %>%
  
  # filter for the cell types of interest for gene marker selection
  filter(is.na(cell_type)==F) %>%
  
  # Group by ancestor
  nest(data = - level_0) %>%
  
  # Eliminate genes that are present in some but all cell types
  # (can be still present in a subset of samples from each cell-type)
  mutate(data = map(
    data,
    ~ .x %>%
      nest(data = -c(symbol, cell_type)) %>%
      add_count(symbol) %>%
      filter(n == max(n)) %>%
      unnest(data)
  )) %>%
  
  # Imputation of missing data within each level_5
  mutate(data = map(data, ~ .x %>% impute_missing_abundance(~ cell_type))) %>%
  
  # scale count for further analysis
  mutate(data=map(data, ~ .x %>%
                    identify_abundant(factor_of_interest = cell_type) %>%
                    scale_abundance()
  ))

saveRDS(tt_simple, "tt_simple.rds")

tt_simple <- readRDS("dev/intermediate_data/tt_simple.rds")

contrast_NH_MC <- tt_simple %>% 
  contrast_MC0()

saveRDS(contrast_NH_MC, "contrast_NH_MC.rds")

contrast_NH_MC <- readRDS("dev/intermediate_data/contrast_NH.rds")

contrast_NH_PW <- tt_simple %>% 
  contrast_PW0()

saveRDS(contrast_NH_PW, "contrast_NH_PW.rds")

contrast_NH_PW <- readRDS("dev/intermediate_data/contrast_NH_pairwise.rds")

## No hierarchy + Mean contrast ===============
sig_size <- 20

NH_MC <- 
  sil_tb0(contrast_NH_MC, sig_size, "tSNE") %>% 
  unnest(sil_df) %>% 
  mutate(ancestor_type = level_0) %>% 
  select(sig_size, real_size, sil, ancestor_type) %>% 
  mutate(analysis = "MC+NH")

## No hierarchy + Pairwise comparison ===============
sig_size <- 20

NH_PW <- 
  sil_tb0(contrast_NH_PW, sig_size, "tSNE") %>% 
  unnest(sil_df) %>% 
  mutate(ancestor_type = level_0) %>% 
  select(sig_size, real_size, sil, ancestor_type) %>% 
  mutate(analysis = "PW+NH")


# Combine all sil_score from all three types of analysis

final_PCA <- bind_rows(H_PW, H_MC, NH_MC, NH_PW)

summary_plot_PCA <- final_PCA %>% 
  ggplot(aes(sig_size, sil, 
             group=interaction(ancestor_type, analysis), 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5))

summary_plot_PCA

ggsave("summary_plot_PCA.png", summary_plot_PCA)

final_tSNE <- bind_rows(H_PW, H_MC, NH_MC, NH_PW)

summary_plot_tSNE <- final_tSNE %>% 
  ggplot(aes(sig_size, sil, 
             group=interaction(ancestor_type, analysis), 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5))

summary_plot_tSNE

ggsave("summary_plot_tSNE.png", summary_plot_tSNE)

# compare 4 marker selection strategie by sil score in different cell types (faceted by nodal cell type)
final_tSNE %>% 
  ggplot(aes(sig_size, sil, 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5)) +
  facet_wrap(~ ancestor_type) +
  ggtitle("Compare four marker selection strategies by silhouette score over 20 signature sizes") +
  theme(plot.title = element_text(hjust = 0.5))

# Below is a summary bar plot comparing silhouette score for four selection methods + CIBERSORT=========
# It requires code from these files: cibersortx.R, hierarchy.R, no_hierarchy.R

sil_all_methods <- bind_rows(ciber_sil, hierarchy_PW_sil, hierarchy_MC_sil, NH_PW_sil, NH_MC_sil) %>% 
  mutate(method = c("CIBERSORTx", "H + PW", "H + MC", "NH + PW", "NH + MC")) %>% 
  select(method, sil, real_size) %>% 
  arrange(desc(sil))

sil_all_methods_plot <- sil_all_methods %>% 
  ggplot(aes(method, sil, fill = method)) +
  geom_col() +
  geom_text(aes(label = signif(sil, 3)), vjust = 1.5)

ggsave("sil_all_methods_plot.png", sil_all_methods_plot)

# Summary TREND PLOT of silhouette score vs real size for all nodes and all methods
# hierarchy data only ======================
naive_methods_H <- list.files("naive methods silhouette score data/", pattern = "[^(NH)]\\.rds")
new_methods_H <- list.files("intermediate_data/", pattern = "[pm].*L[0-9]\\.rds")

new_H <- map_dfr(new_methods_H, ~ readRDS(paste("intermediate_data", .x, sep = "/")))
naive_H <- map_dfr(naive_methods_H, ~ readRDS(paste("naive methods silhouette score data", .x, sep = "/")))

naive_H <- naive_H %>% 
  mutate(method = rep(
    c("mean_contrast.naive.hierarchy", 
      "pairwise.naive.hierarchy"), each=14)
  ) %>% 
  mutate(plot_data = map(plot_data, ~.x %>% select(real_size, silhouette=sil))) %>% 
  unnest(plot_data)

new_H <- new_H %>% 
  mutate(method = rep(
    c("mean_contrast.silhouette.hierarchy", 
      "pairwise.silhouette.hierarchy"), each=14)
  ) %>% 
  unnest(signature_data) %>% 
  mutate(real_size = map_int(cumulative_signature, ~length(.x))) %>% 
  rename(silhouette = winning_silhouette) %>% 
  select(ancestor, real_size, silhouette, method)

full_data_hierarchy <- rbind(new_H, naive_H)

full_data_hierarchy %>% 
  ggplot(aes(real_size, silhouette, color = method))+
  geom_point(size=0.1)+
  geom_line()+
  xlim(0, 100)+
  facet_wrap(~ ancestor) +
  guides(color = guide_legend(
    title.position = "left",
    nrow = 2,
    byrow = TRUE))+
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.title.align = 0.5
  )

ggsave("x.png", x)

# Non-hieararchical data only ==========================

naive_methods_NH <- list.files("naive methods silhouette score data/", pattern = "NH")
new_methods_NH <- list.files("intermediate_data/", pattern = "[pm].*NH")

new_NH <- map_dfr(new_methods_NH, ~ readRDS(paste("intermediate_data", .x, sep = "/")))
naive_NH <- map_dfr(naive_methods_NH, ~ readRDS(paste("naive methods silhouette score data", .x, sep = "/")))

naive_NH <- naive_NH %>% 
  mutate(method = c("mean_contrast.naive.hierarchy", "pairwise.naive.hierarchy")) %>% 
  mutate(plot_data = map(plot_data, ~.x %>% select(real_size, silhouette=sil))) %>% 
  unnest(plot_data)

new_NH <- new_NH %>% 
  mutate(method = "mean_contrast.silhouette.hierarchy") %>% 
  unnest(signature_data) %>% 
  mutate(real_size = map_int(cumulative_signature, ~length(.x))) %>% 
  rename(silhouette = winning_silhouette) %>% 
  select(ancestor, real_size, silhouette, method)

full_data_NH <- rbind(new_NH, naive_NH)

full_data_NH %>% 
  ggplot(aes(real_size, silhouette, color = method))+
  geom_point(size=0.1)+
  geom_line()+
  xlim(0, 1000)+
  theme(legend.position = "bottom",
        legend.title = element_text(size=10)
  )

# Hierarchical and non-hierarchical in the same plot =================
naive_methods <- list.files("dev/topInf_scaleFALSE/unoptimised/", pattern = ".*naive\\..*")
new_methods <- list.files("dev/topInf_scaleFALSE/unoptimised/", pattern = ".*silhouette\\..*")

new <- map_dfr(new_methods, ~ readRDS(paste0("dev/topInf_scaleFALSE/unoptimised/", .x)))
naive <- map_dfr(naive_methods, ~ readRDS(paste0("dev/topInf_scaleFALSE/unoptimised/", .x)))

naive <- naive %>% 
  mutate(method = rep(
    c("mean_contrast.naive.hierarchy",
      "mean_contrast.naive.non_hierarchy",
      "pairwise.naive.hierarchy",
      "pairwise.naive.non_hierarchy"
    ), times=c(14, 1, 14, 1))
  ) %>% 
  mutate(plot_data = map(plot_data, ~.x %>% select(real_size, silhouette))) %>% 
  unnest(plot_data)

new <- new %>% 
  mutate(method = rep(
    c("mean_contrast.silhouette.hierarchy",
      "mean_contrast.silhouette.non_hierarchy",
      "pairwise.silhouette.hierarchy"), times=c(14, 1, 14))
  ) %>% 
  unnest(signature_data) %>% 
  mutate(real_size = map_int(cumulative_signature, ~length(.x))) %>% 
  rename(silhouette = winning_silhouette) %>% 
  select(ancestor, real_size, silhouette, method)

full_df <- rbind(new, naive)

full_df %>% 
  filter(str_detect(method, "silhouette")) %>% # change for naive or silhouette methods
  unnest(data) %>% 
  ggplot(aes(real_size, silhouette, color = method))+
  geom_point(alpha= 0.5)+
  geom_line()+
  xlim(0, 100)+
  facet_wrap(~ ancestor) +
  guides(color = guide_legend(
    title.position = "left",
    nrow = 2,
    byrow = TRUE))+
  labs(title = "silhouette selection", y = "mean silhouette score", x="signature size") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size=12),
    legend.title.align = 0.5
  )

# 20/05/2021
# PCA plots for CD4, CD8, Macrophage using each of the 8 methods ===================================

all_methods_filtered_signature <- 
  all_methods_comparison %>% 
  select(-reduced_dimensions) %>% 
  mutate(tt_filtered = map(
    signature,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x)
  ))

cell_types <- tt_non_hierarchy %>% 
  unnest(tt) %>% 
  unnest(data) %>% 
  pull(cell_type) %>% 
  unique()

# Macrophage ============================================================

macrophage <- all_methods_filtered_signature %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~.x %>% 
      filter(cell_type %in% c("macrophage_M1", "macrophage_M2")) %>% 
      
      # filter genes that have zero variance
      group_by(symbol) %>% 
      mutate(standard_deviation = sd(count_scaled)) %>% 
      filter(standard_deviation != 0) %>% 
      ungroup() %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        top = Inf,
                        scale = FALSE,
                        action = "get",
                        .dims = 10)
  ))

# 1
## PCA
macrophage %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.non_hierarchy")

## tSNE
macrophage %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(tSNE1, tSNE2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.non_hierarchy")

## Heatmap
macrophage %>% 
  pluck("tt_filtered", 1) %>% 
  select(sample, cell_type, contains("PC")) %>% 
  pivot_longer(cols = contains("PC"), names_to = "PC", values_to = "value") %>% 
  heatmap(.row = PC, .column = sample, .value = value) %>% 
  add_tile(cell_type)


# 8
## PCA
macrophage %>% 
  pluck("tt_filtered", 8) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("cibersortx")

## Heatmap
macrophage %>% 
  pluck("tt_filtered", 8) %>% 
  select(sample, cell_type, contains("PC")) %>% 
  pivot_longer(cols = contains("PC"), names_to = "PC", values_to = "value") %>% 
  heatmap(.row = sample, .column = PC, .value = value) %>% 
  add_tile(cell_type)


# CD4 ==========================================
CD4 <- all_methods_filtered_signature %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~.x %>% 
      filter(cell_type %in% c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg",  
                              "t_helper_h1", "t_helper_h2", "t_helper_h17")) %>% 
      
      # filter genes that have zero variance
      group_by(symbol) %>% 
      mutate(standard_deviation = sd(count_scaled)) %>% 
      filter(standard_deviation != 0) %>% 
      ungroup() %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 10)
  ))


x <- all_methods_filtered_signature %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~.x %>% 
      filter(cell_type %in% c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg",  
                              "t_helper_h1", "t_helper_h2", "t_helper_h17")) %>% 
      
      # filter genes that have zero variance
      group_by(symbol) %>% 
      mutate(standard_deviation = sd(count_scaled)) %>% 
      filter(standard_deviation != 0) %>% 
      ungroup()
  ))

y <- x %>% pluck("tt_filtered", 1)

pca <- y %>% 
  select(symbol, sample, count_scaled) %>% 
  pivot_wider(names_from = symbol, values_from = count_scaled, id_cols = sample) %>% 
  select(-sample) %>% 
  data.frame() %>% 
  prcomp(center=TRUE, scale=TRUE)
# number_of_pcs_needed(0.95)

summary(pca)

pca$pc_num
plot(pca$scree)

number_of_pcs_needed <- function(.df, .percent_var){
  sd_pca <- .df %>% 
    data.frame() %>% 
    prcomp(center=TRUE, scale=TRUE) %>% 
    .$sdev
  
  num <- sd_pca %>% 
    (function(x){which(cumsum(x^2/sum(x^2)) >= .percent_var) %>% min()})
  
  scree_plot <- list(x = 1:length(sd_pca),
                     y = sd_pca) 
  
  return(list(pc_num = num, scree = scree_plot))
}


# 1
## PCA
CD4 %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.non_hierarchy; scale=FALSE")

## tSNE
CD4 %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(tSNE1, tSNE2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.non_hierarchy")

## Heatmap
CD4 %>% 
  pluck("tt_filtered", 1) %>% 
  select(sample, cell_type, contains("PC")) %>% 
  pivot_longer(cols = contains("PC"), names_to = "PC", values_to = "value") %>% 
  heatmap(.row = PC, .column = sample, .value = value) %>% 
  add_tile(cell_type)

# 8
## PCA
CD4 %>% 
  pluck("tt_filtered", 8) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("cibersortx")

## Heatmap
CD4 %>% 
  pluck("tt_filtered", 8) %>% 
  select(sample, cell_type, contains("PC")) %>% 
  pivot_longer(cols = contains("PC"), names_to = "PC", values_to = "value") %>% 
  heatmap(.row = sample, .column = PC, .value = value) %>% 
  add_tile(cell_type)

# CD8 ==================================================================
CD8 <- all_methods_filtered_signature %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~.x %>% 
      filter(cell_type %in% c("t_CD8_memory_central", "t_CD8_memory_effector",  "t_CD8_naive"))%>% 
      
      # filter genes that have zero variance
      group_by(symbol) %>% 
      mutate(standard_deviation = sd(count_scaled)) %>% 
      filter(standard_deviation != 0) %>% 
      ungroup() %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        top = Inf,
                        scale = FALSE,
                        action = "get",
                        .dims = 10)
  ))

# 1
#PCA
CD8 %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.non_hierarchy")

## tSNE
CD8 %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(tSNE1, tSNE2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.non_hierarchy")

## Heatmap
CD8 %>% 
  pluck("tt_filtered", 1) %>% 
  select(sample, cell_type, contains("PC")) %>% 
  pivot_longer(cols = contains("PC"), names_to = "PC", values_to = "value") %>% 
  heatmap(.row = PC, .column = sample, .value = value) %>% 
  add_tile(cell_type)

# 8
## PCA
CD8 %>% 
  pluck("tt_filtered", 8) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("cibersortx")

## Heatmap
CD8 %>% 
  pluck("tt_filtered", 8) %>% 
  select(sample, cell_type, contains("PC")) %>% 
  pivot_longer(cols = contains("PC"), names_to = "PC", values_to = "value") %>% 
  heatmap(.row = sample, .column = PC, .value = value) %>% 
  add_tile(cell_type)


# node-signature pairwise =================
pairwise.silhouette.hierarchy <- readRDS("intermediate_data/pairwise.silhouette.hierarchy.rds")

macrophage_signature <- pairwise.silhouette.hierarchy %>% 
  slice(9) %>% 
  pull(cumulative_signature) %>% 
  unlist() %>% 
  unique()

CD4_signature <- pairwise.silhouette.hierarchy %>% 
  slice(c(8, 12, 14)) %>% 
  pull(cumulative_signature) %>% 
  unlist() %>% 
  unique()

CD8_signature <- pairwise.silhouette.hierarchy %>% 
  slice(c(10, 13)) %>% 
  pull(cumulative_signature) %>% 
  unlist() %>% 
  unique()

nodal_filtered <- 
  tibble(cell = c("macrophage", "CD4", "CD8"),
         signature = list(macrophage_signature, CD4_signature, CD8_signature)
  ) %>% 
  mutate(tt_filtered = map(
    signature,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x)
  ))

# macrophage ==========================================

macrophage_node <- nodal_filtered %>% 
  slice(1) %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~ .x %>% 
      filter(cell_type %in% c("macrophage_M1", "macrophage_M2")) %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 6 # macrophage has only 6 signature genes
      )
  ))


## PCA
macrophage_node %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("pairwise.silhouette.hierarchy; scale=FALSE")

# CD4 ==========================================

CD4_node <- nodal_filtered %>% 
  slice(2) %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~ .x %>% 
      filter(cell_type %in% c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg",  
                              "t_helper_h1", "t_helper_h2", "t_helper_h17")) %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 10)
  ))


## PCA
CD4_node %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("pairwise.silhouette.hierarchy; scale=FALSE")

# CD8 ==========================================

CD8_node <- nodal_filtered %>% 
  slice(3) %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~ .x %>% 
      filter(cell_type %in% c("t_CD8_memory_central", "t_CD8_memory_effector",  "t_CD8_naive")) %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 10)
  ))


## PCA
CD8_node %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("pairwise.silhouette.hierarchy; scale=FALSE")





# node-signature mean_contrast =================
mean_contrast.silhouette.hierarchy <- readRDS("intermediate_data/mean_contrast.silhouette.hierarchy.rds")

macrophage_signature_mc <- mean_contrast.silhouette.hierarchy %>% 
  slice(9) %>% 
  pull(cumulative_signature) %>% 
  unlist() %>% 
  unique()

CD4_signature_mc <- mean_contrast.silhouette.hierarchy %>% 
  slice(c(8, 12, 14)) %>% 
  pull(cumulative_signature) %>% 
  unlist() %>% 
  unique()

CD8_signature_mc <- mean_contrast.silhouette.hierarchy %>% 
  slice(c(10, 13)) %>% 
  pull(cumulative_signature) %>% 
  unlist() %>% 
  unique()

nodal_filtered_mc <- 
  tibble(cell = c("macrophage", "CD4", "CD8"),
         signature = list(macrophage_signature_mc, CD4_signature_mc, CD8_signature_mc)
  ) %>% 
  mutate(tt_filtered = map(
    signature,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x)
  ))

# macrophage ==========================================

macrophage_node_mc <- nodal_filtered_mc %>% 
  slice(1) %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~ .x %>% 
      filter(cell_type %in% c("macrophage_M1", "macrophage_M2")) %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 6 # macrophage has only 6 signature genes
      )
  ))


## PCA
macrophage_node_mc %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy; scale=FALSE")

# CD4 ==========================================

CD4_node_mc <- nodal_filtered_mc %>% 
  slice(2) %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~ .x %>% 
      filter(cell_type %in% c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg",  
                              "t_helper_h1", "t_helper_h2", "t_helper_h17")) %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  ))


## PCA
CD4_node_mc %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy; pure 14 CD4 signature")

# CD8 ==========================================

CD8_node_mc <- nodal_filtered_mc %>% 
  slice(3) %>% 
  mutate(tt_filtered = map(
    tt_filtered,
    ~ .x %>% 
      filter(cell_type %in% c("t_CD8_memory_central", "t_CD8_memory_effector",  "t_CD8_naive")) %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 10)
  ))


## PCA
CD8_node_mc %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy; scale=FALSE")





markers_level <- 
  mean_contrast.silhouette.hierarchy %>% 
  nest(number_of_markers = - level) %>% 
  mutate(number_of_markers = map(
    number_of_markers,
    ~ .x %>% 
      pull(cumulative_signature) %>% 
      unlist() %>% 
      unique()
  ))

CD4_and_level_noise <- markers_level %>% 
  mutate(markers_not_in_CD4 = map(number_of_markers,
                                  ~ .x[!.x %in% CD4_signature_mc])) %>% 
  mutate(level_markers_and_CD4 = map(markers_not_in_CD4,
                                     ~.x %>% append(CD4_signature_mc))) %>% 
  mutate(tt_filtered = map(
    level_markers_and_CD4,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      filter(cell_type %in% c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg",  
                              "t_helper_h1", "t_helper_h2", "t_helper_h17")) %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  ))

# Level 1 noise
CD4_and_level_noise %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + Level 1 noise")

# Level 2 noise
CD4_and_level_noise %>% 
  pluck("tt_filtered", 2) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + Level 2 noise")

# Level 3 noise
CD4_and_level_noise %>% 
  pluck("tt_filtered", 3) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + Level 3 noise")

# Level 4 noise
CD4_and_level_noise %>% 
  pluck("tt_filtered", 4) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + Level 4 noise")

# Level 5 noise
CD4_and_level_noise %>% 
  pluck("tt_filtered", 5) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + Level 5 noise")


# Same number of noise + CD4 markers
CD4_and_same_number_noise <- markers_level %>% 
  
  mutate(markers_not_in_CD4 = map(
    number_of_markers,
    ~ .x[!.x %in% CD4_signature_mc])) %>% 
  
  slice(1:4) %>% 
  
  mutate(sampled_noise = map(
    markers_not_in_CD4,
    ~.x %>% 
      sample(21)
  )) %>% 
  
  mutate(sampled_noise_and_CD4 = map(
    sampled_noise,
    ~.x %>% 
      append(CD4_signature_mc)
  )) %>% 
  
  mutate(tt_filtered = map(
    sampled_noise_and_CD4,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      filter(cell_type %in% c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg",  
                              "t_helper_h1", "t_helper_h2", "t_helper_h17")) %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  ))

# 21 Level 1 noise
CD4_and_same_number_noise %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + 21 Level 1 noise")

# 21 Level 2 noise
CD4_and_same_number_noise %>% 
  pluck("tt_filtered", 2) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + 21 Level 2 noise")

# 21 Level 3 noise
CD4_and_same_number_noise %>% 
  pluck("tt_filtered", 3) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + 21 Level 3 noise")

# 21 Level 4 noise
CD4_and_same_number_noise %>% 
  pluck("tt_filtered", 4) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4 + 21 Level 4 noise")


xx <- all_methods_comparison %>% 
  slice(2) %>% 
  select(-reduced_dimensions) %>% 
  mutate(signature = map(signature, ~ .x[!.x %in% markers_L1_L2])) %>% 
  mutate(tt_filtered = map(
    signature,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      filter(cell_type %in% c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg",  
                              "t_helper_h1", "t_helper_h2", "t_helper_h17")) %>% 
      # filter genes that have zero variance
      group_by(symbol) %>% 
      mutate(standard_deviation = sd(count_scaled)) %>% 
      filter(standard_deviation != 0) %>% 
      ungroup() %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  ))

xx %>% 
  pluck("tt_filtered", 1) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.hierarchy;CD4;removing signature from L1 & L2")

non_CD4_markers <- all_methods_comparison %>% 
  pluck("signature", 2) %>% 
  .[! . %in% CD4_signature_mc]

ten <- sample(non_CD4_markers, 24)
twenty <- sample(non_CD4_markers, 48)
thirty <- sample(non_CD4_markers, 72)

ten <- sample(markers_L1_L2, 14)
twenty <- sample(markers_L1_L2, 29)
thirty <- sample(markers_L1_L2, 43)

xx <- tibble(
  mix = c("10% non_CD4 + CD4", "20% non CD4 + CD4", "30% non CD4 + CD4"),
  signature = list(append(CD4_signature_mc, ten), 
                   append(CD4_signature_mc, twenty),
                   append(CD4_signature_mc, thirty))
) %>% 
  mutate(tt_filtered = map(
    signature,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      filter(cell_type %in% c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg",  
                              "t_helper_h1", "t_helper_h2", "t_helper_h17")) %>% 
      # filter genes that have zero variance
      group_by(symbol) %>% 
      mutate(standard_deviation = sd(count_scaled)) %>% 
      filter(standard_deviation != 0) %>% 
      ungroup() %>% 
      
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  ))

# 10% non_CD4 + CD4
xx %>% 
  pluck("tt_filtered", 3) %>% 
  ggplot(aes(PC1, PC2, color = cell_type), label=sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("43 (30%) L1&L2 markers + 44 CD4 markers")







# comparisons for all methods using curvature optimisation

# import unoptimised data from all methods
naive <- list.files("topInf_scaleFALSE/unoptimised/", pattern = ".*naive\\..*\\..*")
silhouette <- list.files("topInf_scaleFALSE/unoptimised/", pattern = ".*silhouette\\..*\\..*")

naive_df <- map_dfr(naive, ~ readRDS(paste0("topInf_scaleFALSE/unoptimised/", .x))) %>%
  mutate(method = rep(str_replace_all(naive, '\\.rds', ''), c(14, 1, 14, 1)))

o <- rep(str_replace_all(silhouette, '\\.rds', ''), c(14, 1, 14))
silhouette_df <- map_dfr(silhouette, ~ readRDS(paste0("topInf_scaleFALSE/unoptimised/", .x))) %>% 
  mutate(method = o) %>% 
  mutate(data = map(data, ~ .x %>% dplyr::rename(signature = cumulative_signature)))
rm(o)

full_df <- bind_rows(naive_df, silhouette_df) %>% 
  nest(signature = -method)


full_df <- full_df %>% 
  mutate(signature = map(
    signature,
    ~.x %>% 
      do_optimisation(.optimisation_method = "curvature")
  ))

silhouette_score <- function(.reduced_dimensions, .distance, .level){
  
  .reduced_dimensions %>% 
    
    pull(!!as.symbol(.level)) %>% 
    
    as.factor() %>% 
    
    as.numeric() %>% 
    
    silhouette(.distance) %>% 
    
    summary()
  
}

silhouette_function <- function(.selected, .reduction_method){
  
  .selected %>% 
    
    # reduce dimensions
    mutate(reduced_dimensions = map2(
      markers, level, 
      ~ dimension_reduction(.x, .y, .reduction_method)
    )) %>% 
    
    # calculate distance matrix using PC1 & PC2
    mutate(distance = map(
      reduced_dimensions,
      ~ distance_matrix(.x, .reduction_method)
    )) %>% 
    
    # calculate silhouette score
    mutate(silhouette = pmap(
      list(reduced_dimensions, distance, level),
      ~ silhouette_score(..1, ..2, ..3)
    )) %>% 
    
    # remove unnecessary columns
    select(-c(markers, distance))
  
}

all_methods_silhouette <- full_df %>% 
  mutate(signature = map(signature, ~.x %>% pull(signature) %>% unlist() %>% unique())) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score for all signatures combined in each method
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

cibersortx <- readRDS("topInf_scaleFALSE/cibersortx.new.rds")
cibersort_signature <- cibersortx$signature[[1]]

cibersortx <- tibble(method = "cibersortx") %>% 
  mutate(signature = list(cibersort_signature)) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

# silhouette
all_methods_comparison <- all_methods_silhouette %>% 
  bind_rows(cibersortx) %>% 
  mutate(cluster.silhouette = map(silhouette, ~ .x$clus.avg.widths)) %>% 
  mutate(avg.silhouette = map_dbl(silhouette, ~ .x$avg.width)) %>% 
  select(-c(reduced_dimensions, silhouette)) %>% 
  arrange(desc(avg.silhouette)) %>% 
  mutate(method = str_remove_all(method, "\\.new") %>% 
           str_remove_all("\\.unOP"))

all_methods_comparison %>%
  unnest(cluster.silhouette) %>% 
  ggplot(aes(x=reorder(method, avg.silhouette), y=cluster.silhouette, colour=method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  theme(axis.text.x = element_blank())


mix100 <- readRDS("intermediate_data/mix100.rds")

deconvolution_all_methods <- mix100 %>% 
  tidyr::expand(
    nesting(method=all_methods_comparison$method, 
            signature=all_methods_comparison$signature),
    nesting(mixture_ID, mix)
  ) %>% 
  
  # reference
  mutate(reference = map(
    signature,
    
    # filter out data for signature genes
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      
      # reshape the input matrix for deconvolve_cellularity():
      select(symbol, cell_type, sample, count_scaled) %>% 
      group_by(symbol, cell_type) %>% 
      summarise(count_scaled_median = median(count_scaled)) %>% 
      ungroup() %>% 
      pivot_wider(id_cols = symbol, names_from = cell_type, values_from = count_scaled_median) %>% 
      tidybulk::as_matrix(rownames = symbol) # must be a matrix
  )) %>% 
  
  # deconvolution
  mutate(deconvolution = map2(
    mix, reference,
    ~ tidybulk::deconvolve_cellularity(
      .x,
      replicate, symbol, count_mix,
      reference = .y,
      method = "llsr", 
      prefix = "llsr_", 
      action = "get") %>% 
      
      pivot_longer(cols=starts_with("llsr_"), 
                   names_prefix ="llsr_", 
                   names_to="cell_type", 
                   values_to="estimated_proportion") %>%
      
      left_join(.x %>% 
                  unnest(data_samples) %>% 
                  distinct(replicate, cell_type, proportion))
  ))


deconvolution_all_methods %>% 
  select(method, mixture_ID, deconvolution) %>% 
  mutate(method = str_remove_all(method, "\\.new") %>% str_remove_all("\\.unOP") ) %>% 
  
  # macro.mse for ordering
  mutate(MSE = map_dbl(
    deconvolution,
    ~ mean((.x$estimated_proportion - .x$proportion)^2)
  )) %>% 
  # nest(data=-method) %>%
  # mutate(macro.mse = map_dbl(data, ~ mean(.x$MSE))) %>%
  # unnest(data) %>%
  
  # mse by cell type
  unnest(deconvolution) %>% 
  mutate(squared.error = (estimated_proportion - proportion)^2) %>% 
  nest(data = -c(method, cell_type)) %>% 
  mutate(mse.cell = map_dbl(data, ~ mean(.x$squared.error))) %>%
  mutate(macro.mse = map_dbl(data, ~ mean(.x$MSE))) %>%
  group_by(method) %>% 
  summarise(method, cell_type, data, mse.cell, macro.mse) %>%
  # summarise(method, cell_type, data, mse.cell, median = median(mse.cell)) %>%
  ungroup() %>% 
  
  ggplot(aes(x=reorder(method, -macro.mse), y=log(mse.cell))) +
  # ggplot(aes(x=reorder(method, -median), y=log(mse.cell))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = cell_type), 
              position=position_jitter(0.2)) +
  # theme(axis.text.x = element_blank())
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))

# boxplot of macro.mse by method
deconvolution_all_methods %>% 
  select(method, mixture_ID, deconvolution) %>% 
  mutate(method = str_remove_all(method, "\\.new") %>% str_remove_all("\\.unOP") ) %>% 
  mutate(MSE = map_dbl(
    deconvolution,
    ~ mean((.x$estimated_proportion - .x$proportion)^2)
  )) %>% 
  group_by(method) %>% 
  summarise(MSE, macro.mse=mean(MSE)) %>% 
  ungroup() %>% 
  ggplot(aes(x= reorder(method, -macro.mse), y=log(MSE), colour = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  # theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1)) +
  theme(axis.text.x = element_blank())

