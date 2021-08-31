source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")
library(patchwork)
library(tidyverse)
library(tidybulk)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
input_directory = args[1]
output_directory = args[2]

# import data
mix100 <- readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/mix100.rds")
counts_non_hierarchy <- 
  readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_non_hierarchy.rds")

# TO BE DELETED!
indices <- sample(1:length(dir(input_directory)), 7)
input_directory = "dev/benchmark_results/"
output_directory = "dev/benchmark_results/benchmark_plot.pdf"


plot_data <- dir(input_directory) %>% 
  `names<-`(dir(input_directory)) %>% 
  
  # TO BE DELETED!
  .[indices] %>% 
  
  map_dfr(~ readRDS(glue("{input_directory}{.x}")), .id = "stream") %>% 
  mutate(stream = str_remove(stream, "\\.rds")) %>% 
  nest(signature = - stream) %>% 
  mutate(signature = map(signature, ~ .x$signature %>% unlist() %>% unique())) %>% 
  
  # bind cibersortx
  bind_rows(
    read_delim("dev/jian_R_files/cibersortx/CIBERSORTx_Job21_phenoclass_1.CIBERSORTx_Job21_reference_1.bm.K999.txt", 
               "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
      pull(NAME) %>% 
      list() %>% 
      tibble(stream = "cibersortx", signature = .)
  ) %>% 
  
  # silhouette evaluation
  mutate(silhouette = map(
    signature, 
    ~ silhouette_evaluation(
      .signature = .x,
      .reduction_method = "PCA",
      .preprocessed_non_hierarchy = counts_non_hierarchy)
    )) %>% 
  
  mutate(avg_silhouette = map_dbl(silhouette, ~ mean(.x$sil_width))) %>% 
  
  mutate(silhouette = map(
    silhouette, 
    ~ .x %>% 
      group_by(cell_type) %>% 
      summarise(cluster_silhouette = mean(sil_width)) %>% 
      distinct() %>% 
      ungroup()
    )) %>% 
  
  # deconvolution evaluation
  # for each mixture, combine with the signatures from all methods
  expand_grid(mix100 %>% dplyr::slice(1:10), .) %>% 
  
  mutate(deconvolution = map2(
    signature, mix, 
    ~ deconvolution_evaluation(
      .signature = .x, 
      .mix=.y, 
      .preprocessed_non_hierarchy = counts_non_hierarchy)
    )) %>% 
  
  # mse by method
  mutate(MSE = map_dbl(
    deconvolution,
    ~ mean((.x$estimated_proportion - .x$proportion)^2)
  )) %>% 
  nest(data=-stream) %>% 
  # mutate(median_MSE_over_mixes = map_dbl(data, ~ median(.x$MSE))) %>% 
  mutate(mean_MSE_over_mixes = map_dbl(data, ~ mean(.x$MSE))) %>% 
  unnest(data) %>% 
  
  # mse by cell type
  unnest(deconvolution) %>% 
  mutate(squared_error = (estimated_proportion - proportion)^2) %>% 
  nest(data = -c(stream, cell_type)) %>% 
  # mutate(median_MSE_for_cell_type = map_dbl(data, ~ median(.x$squared_error))) %>% 
  mutate(mean_MSE_for_cell_type = map_dbl(data, ~ mean(.x$squared_error))) %>% 
  unnest(data) %>% 
  
  select(-c(signature, mixture_ID, mix, replicate, estimated_proportion, proportion, squared_error))


# silhouette plot

boxplot_silhouette <- plot_data %>%

  nest(data = -c(stream, silhouette, avg_silhouette)) %>%
  select(-data) %>%
  unnest(silhouette) %>%
  
  ggplot(aes(x=reorder(stream, avg_silhouette), y=cluster_silhouette)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = cell_type), 
              position=position_jitter(0.2)) +
  labs(y = "silhouette score of cell type clusters",
       title = "benchmark by silhouette score",
       tag = "A",
       caption = "(streams are arranged by mean silhoutte score ascendingly.)"
       ) +
  
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1, face = "bold"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
        )

  
# deconvolution plot
boxplot_deconvolution_by_cell_type <- plot_data %>% 
  nest(data = -c(stream, cell_type, mean_MSE_for_cell_type, mean_MSE_over_mixes)) %>% 
  select(-data) %>% 
  
  ggplot(aes(x=reorder(stream, -mean_MSE_over_mixes), y=log(mean_MSE_for_cell_type))) +
  # ggplot(aes(x=reorder(method, -median), y=log(mse.cell))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = cell_type), 
              position=position_jitter(0.2)) +
  
  labs(title = "benchmark by mean deconvolution MSE over 100 mixes for cell types",
       tag = "B",
       caption = "(streams are arranged by mean deconvolution MSE over 100 mixes descendingly.)"
  ) +
  
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1, face = "bold"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
  )


boxplot_deconvolution_by_method <- plot_data %>% 
  nest(data = -c(stream, MSE, mean_MSE_over_mixes)) %>% 
  select(-data) %>% 
  
  ggplot(aes(x=reorder(stream, -mean_MSE_over_mixes), y=log(MSE), colour = stream)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  labs(title = "benchmark by deconvolution MSE over 100 mixes",
       tag = "C",
       caption = "(streams are arranged by mean deconvolution MSE over 100 mixes descendingly.)"
  ) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )


patch1 <- boxplot_silhouette + boxplot_deconvolution_by_cell_type + theme(legend.position = "right")


pdf(file = output_directory,   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

# Step 2: Create the plot with R code
patch1 / boxplot_deconvolution_by_method

# Step 3: Run dev.off() to create the file!
dev.off()



