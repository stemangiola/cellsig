library(tidyverse)
library(ggplot2)
library(cluster)
library(proxy)
library(factoextra)
library(stringr)
library(scales)
# library(tidySummarizedExperiment)
library(tidybulk)


# load signatures from all 8 methods
all_methods_comparison <- readRDS("topInf_scaleFALSE/all_methods_comparison.new.rds")
all_methods_comparison

counts_non_hierarchy <- readRDS("./dev/intermediate_data/counts_non_hierarchy.rds")


mix_base = readRDS("./dev/counts_imputed.rds") %>% 
  
  # rename columns and calculate count_scaled
  dplyr::rename(symbol = .feature, sample = .sample) %>% 
  # mutate(count_scaled = count / exp(exposure_rate)) %>% 
  
  # # keep genes present in all
  # filter(count_scaled %>% is.na %>% `!`) %>%
  # nest(data = -c(cell_type, symbol)) %>%
  # add_count(symbol) %>%
  # filter(n == max(n)) %>%
  # unnest(data) %>% 
  # select(-n)

saveRDS(mix_base, "./dev/intermediate_data/mix_base.rds", compress = "xz")


# have a 100 mixtures future map
# EMPY YOUR SESSION FROM BIG OBJECTS FIRST! 
# LOOK AT NEXT TO "Import Dataset". DO YOU SEE > 10Gb. 
# IF YES, DELETE OBJECTS.
library(future)
library(furrr)
plan(multisession, workers = 15)
options(future.globals.maxSize= +Inf)

mix100 <- 
  
  tibble(mixture_ID = 1:100) %>%
  
  # mix
  mutate(mix = map(mixture_ID, ~ {
    proportions = 
      gtools::rdirichlet(1, rep(1, length(unique(mix_base$cell_type)))) %>%
      as.data.frame() %>%
      setNames(unique(mix_base$cell_type))
    
    cellsig::generate_mixture_from_proportion_matrix(mix_base, proportions)
  }))
    
saveRDS(mix100, "./dev/intermediate_data/mix100.rds", compress = "xz")

deconvolution_all_methods <-
  
  mix100 %>% 
  
  # for each mixture, combine with the signatures from 8 methods
  tidyr::expand(
    nesting(method=all_methods_comparison$method, 
            signature=all_methods_comparison$signature),
    nesting(mixture_ID, mix)
  )

deconvolution_all_methods <- deconvolution_all_methods %>% 
  
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

# mse <- deconvolution_all_methods %>% 
#   # calculate deconvolution MSE for each method
#   mutate(MSE = map_dbl(
#     deconvolution,
#     ~ mean((.x$estimated_proportion - .x$proportion)^2)
#   )) %>% 
#   nest(MSE.macro = -c(method, signature)) %>% 
#   mutate(MSE.macro = map_dbl(MSE.macro, ~ mean(.x$MSE))) %>% 
#   mutate(method = map_chr(method, ~ .x %>% str_replace_all("\\.new", "") )) %>% 
#   arrange(MSE.macro)
# 
# 
# deconvolution_all_methods %>% 
#   select(method, mixture_ID, deconvolution) %>% 
#   mutate(method = map_chr(method, ~ .x %>% str_replace_all("\\.new", "") )) %>% 
#   mutate(MSE = map_dbl(
#     deconvolution,
#     ~ mean((.x$estimated_proportion - .x$proportion)^2)
#   )) %>% 
#   nest(data=-method) %>% 
#   mutate(macro.mse = map_dbl(data, ~ mean(.x$MSE))) %>% 
#   unnest(data) %>% 
#   
#   unnest(deconvolution) %>% 
#   mutate(squared.error = (estimated_proportion - proportion)^2) %>% 
#   nest(data = -c(method, cell_type)) %>% 
#   mutate(mse.cell = map_dbl(data, ~ mean(.x$squared.error))) %>%
#   mutate(macro.mse = map_dbl(data, ~ mean(.x$MSE))) %>% 
#   
#   ggplot(aes(x=reorder(method, -macro.mse), y=sqrt(mse.cell))) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(color = cell_type), 
#               position=position_jitter(0.2)) +
#   # theme(axis.text.x = element_blank())
#   theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
# 
# deconvolution_all_methods %>% 
#   select(method, mixture_ID, deconvolution) %>% 
#   mutate(method = map_chr(method, ~ .x %>% str_replace_all("\\.new", "") )) %>% 
#   mutate(MSE = map_dbl(
#     deconvolution,
#     ~ mean((.x$estimated_proportion - .x$proportion)^2)
#   )) %>% 
#   group_by(method) %>% 
#   summarise(MSE, macro.mse=mean(MSE)) %>% 
#   ungroup() %>% 
#   ggplot(aes(x= reorder(method, -macro.mse), y=sqrt(MSE), colour = method)) +
#   geom_boxplot() +
#   geom_jitter(position=position_jitter(0.2), alpha=0.5) +
#   theme(axis.text.x = element_blank())
# 
# # =======================================================
# deconvolution_all_methods <- all_methods_comparison %>% 
#   select(-reduced_dimensions) %>% 
#   mutate(deconvolution = map(
#     signature,
#     
#     # filter out data for signature genes
#     ~ tt_non_hierarchy %>% 
#       unnest(tt) %>% 
#       unnest(data) %>% 
#       filter(symbol %in% .x) %>% 
#       
#       # reshape the input matrix for deconvolve_cellularity():
#       select(symbol, cell_type, sample, count_scaled) %>% 
#       group_by(symbol, cell_type) %>% 
#       summarise(count_scaled_median = median(count_scaled)) %>% 
#       ungroup() %>% 
#       pivot_wider(id_cols = symbol, names_from = cell_type, values_from = count_scaled_median) %>% 
#       tidybulk::as_matrix(rownames = symbol) # must be a matrix
#     
#   )) %>% 
#   
#   # deconvolution
#   mutate(deconvolution = map(
#     deconvolution,
#     ~ tidybulk::deconvolve_cellularity(
#       mix,
#       replicate, symbol, count_mix,
#       reference = .x %>% tidybulk::as_matrix(rownames = symbol), # must be a matrix
#       method = "llsr", 
#       prefix = "llsr_", 
#       action = "get"
#     ) %>% 
#       
#       pivot_longer(cols=starts_with("llsr_"), 
#                    names_prefix ="llsr_", 
#                    names_to="cell_type", 
#                    values_to="estimated_proportion") %>%
#       left_join(mix %>% 
#                   unnest(data_samples) %>% 
#                   distinct(replicate, cell_type, proportion)
#       )
#     
#   )) %>% 
#   
#   # calculate mean squared error between estimated proportion and true proportion(from mix)
#   mutate(MSE = map_dbl(
#     deconvolution,
#     ~ mean((.x$estimated_proportion - .x$proportion)^2)
#   ))
# 
# 
# 
# # MSE plot ================================================
# mse %>% 
#   ggplot(aes(x = reorder(method, -MSE.macro), y = MSE.macro, fill = method)) +
#   geom_bar(stat = "identity") +
#   theme(axis.text.x = element_blank()) + 
#   geom_text(aes(x = reorder(method, -MSE.macro), 
#                 y = MSE.macro, 
#                 label = signif(MSE.macro, 2)), 
#             vjust = 0) +
#   ggtitle("All methods comparison using deconvolution mean squared error")
# 
# # silhouette
# deconvolution_all_methods %>% 
#   ggplot(aes(x = reorder(method, silhouette), y = silhouette, fill = method)) +
#   geom_bar(stat = "identity") +
#   theme(axis.text.x = element_blank()) +
#   geom_text(aes(label = signif(silhouette, 2)), vjust = 1)
# 
# #====================================
# 
# # stefano's original code
# tidybulk::deconvolve_cellularity(
#   mix,
#   replicate, symbol, count_mix,
#   #reference = ...,
#   method = "llsr", 
#   prefix = "llrs_", 
#   action = "get"
# ) %>%
#   pivot_longer(cols=starts_with("llrs_"), names_prefix ="llrs_", names_to="cell_type", values_to="estimated_proportion") %>%
#   left_join(
#     mix %>% 
#       unnest(data_samples) %>% 
#       distinct(replicate, cell_type, proportion)
#   )
# 
# # debug technique
# deconvolution_all_methods <- deconvolution_all_methods %>% 
#   
#   mutate(deconvolution = map(
#     deconvolution,
#     # debugging technique
#     ~ {browser(); (.)} %>%  tidybulk::deconvolve_cellularity(
#       mix,
#       replicate, symbol, count_mix,
#       reference = .x,
#       method = "llsr", 
#       prefix = "llsr_", 
#       action = "get"
#     ) 
#   ))



 