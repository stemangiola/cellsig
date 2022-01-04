# devtools::install_github("stemangiola/tidybulk@sparse-counts", force = TRUE)
source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
input_directory = args[1]
output_directory = args[2]
dim = args[3] %>% as.integer()

# import data
mix100 <- readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/mix100.rds")

counts_imputed <- 
  readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_imputed.rds") %>% 
  rename(symbol = feature)

new_tree <- 
  read_yaml("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/new_tree.yaml") %>% as.Node

# # TO BE DELETED!
# indices <- sample(1:length(dir(input_directory)), 7)
# input_directory = "dev/benchmark_results_multiPC_NH/"
# output_directory = "dev/benchmark_results_multiPC_NH/"
# dim = 4L


plot_data <- dir(input_directory, pattern = glue(".*{dim}\\.rds")) %>%
  `names<-`(dir(input_directory, pattern = glue(".*{dim}\\.rds"))) %>% 
  
  # # TO BE DELETED!
  # .[1] %>%
  
  map_dfr(~ readRDS(glue("{input_directory}{.x}")), .id = "stream") %>% 
  mutate(stream = str_remove(stream, "\\.rds")) %>% 
  nest(signature = - stream) %>% 
  mutate(signature = map(signature, ~ .x$signature %>% unlist() %>% unique())) %>% 
  
  # bind cibersortx
  bind_rows(
    read_delim("dev/jian_R_files/cibersortx/CIBERSORTx_Job22_phenoclass_new_tree.CIBERSORTx_Job22_reference_new_tree.bm.K999.txt", 
               "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
      pull(NAME) %>% 
      list() %>% 
      tibble(stream = "cibersortx", signature = .)
  ) %>% 
  
  # for differentiating hierarchical and non_hierarchical methods
  mutate(is_hierarchy = case_when(
    str_detect(stream, "^hierarchical") ~ "hierarchical",
    str_detect(stream, "non") ~ "non_hierarchical",
    TRUE ~ "cibersortx"
  ), .before = signature) %>% 
  
  evaluation(.stream = stream, .markers = signature, 
             .mixture = mix100, 
             .imputed_counts = counts_imputed, .sample = sample, .symbol = symbol, .cell_type = cell_type,
             .tree = new_tree, .reduction_method = "PCA", .dims = 2)



# silhouette plot

boxplot_silhouette <- plot_data %>%

  nest(data = -c(stream, is_hierarchy, silhouette, avg_silhouette)) %>%
  select(-data) %>%
  unnest(silhouette) %>%
  
  # ggplot(aes(x=reorder(stream, avg_silhouette), y=cluster_silhouette)) +
  ggplot(aes(x=reorder(stream, cluster_silhouette, median), y=cluster_silhouette)) +
  geom_boxplot(aes(fill = is_hierarchy),
               alpha = 0.4, 
               outlier.shape = NA) +
  geom_jitter(aes(colour = cell_type,
                  size = cluster_size
                  ), 
              alpha = 0.8,
              position=position_jitter(0.2)) +
  labs(y = "mean silhouette score for clusters",
       title = "benchmark by mean silhouette score for cell type clusters",
       # tag = "A",
       caption = "(streams are arranged by medians of mean silhouette score ascendingly.)"
       ) +
  
  guides(
    color = guide_legend(title.position = "left", ncol = 7, byrow = FALSE, order=3),
    
    fill = guide_legend(title = NULL, ncol = 1, byrow = FALSE, order =1),
    
    size = guide_legend(title.position = "left", ncol = 1, byrow = FALSE, order=2)
  ) +
  
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust = 1, face = "bold", size = 7),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(angle = 90),
        legend.title.align = 0.5,
        legend.position = "bottom",
        legend.spacing = unit(0, "cm"),
        legend.box.spacing = unit(0, "cm"),
        plot.margin=unit(c(0, 2.5, 0, 2.5), "cm")
        )

  
# # deconvolution plot
# boxplot_deconvolution_by_cell_type <- plot_data %>% 
#   nest(data = -c(stream, is_hierarchy, cell_type, MSE_over_cell_type, mean_MSE_over_mixes)) %>% 
#   select(-data) %>% 
#   
#   # ggplot(aes(x=reorder(stream, -mean_MSE_over_mixes), y=log10(MSE_over_cell_type))) +
#   ggplot(aes(x=reorder(stream, -MSE_over_cell_type, median), y=log10(MSE_over_cell_type))) +
#   geom_boxplot(aes(fill = is_hierarchy), 
#                alpha = 0.4, 
#                outlier.shape = NA) +
#   geom_jitter(aes(color = cell_type), 
#               alpha = 0.8,
#               position=position_jitter(0.2)) +
#   
#   labs(title = "benchmark by cell type-specific deconvolution MSE over 100 mixes",
#        # tag = "B",
#        caption = "(streams are arranged by median of cell-specific deconvolution MSE descendingly.)"
#   ) +
#   
#   guides(
#     color = guide_legend(title.position = "left", ncol = 7, byrow = FALSE, order=2),
#     
#     fill = guide_legend(title = NULL, ncol = 1, byrow = FALSE, order=1)
#     ) +
#   
#   theme(axis.text.x = element_text(angle=60, vjust=1, hjust = 1, face = "bold", size = 7),
#         axis.title.x = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_text(angle = 90),
#         legend.title.align = 0.5,
#         legend.position = "bottom",
#         legend.spacing = unit(0, "cm"),
#         legend.box.spacing = unit(0, "lines"),
#         plot.margin=unit(c(0, 2.5, 0, 2.5), "cm")
#   )
# 
# 
# boxplot_deconvolution_by_method <- plot_data %>% 
#   nest(data = -c(stream, is_hierarchy, MSE, mean_MSE_over_mixes)) %>% 
#   select(-data) %>% 
#   
#   # ggplot(aes(x=reorder(stream, -mean_MSE_over_mixes), y=log10(MSE), colour = stream)) +
#   ggplot(aes(x=reorder(stream, -MSE, median), y=log10(MSE))) +
#   geom_boxplot(aes(fill = is_hierarchy), 
#                alpha = 0.4, 
#                outlier.shape = NA) +
#   geom_jitter(position=position_jitter(0.2), alpha=0.5) +
#   labs(title = "benchmark by deconvolution MSE over 100 mixes",
#        # tag = "C",
#        caption = "(streams are arranged by median deconvolution MSE descendingly.)"
#   ) +
#   theme(axis.text.x = element_text(angle=60, vjust=1, hjust = 1, face = "bold", size = 7),
#         axis.title.x = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         plot.margin=unit(c(0, 2.5, 0, 2.5), "cm")
#   )


# patch1 <- boxplot_silhouette + 
#   boxplot_deconvolution_by_cell_type + 
#   theme(legend.position = "bottom",
#         legend.spacing.x = unit(0, "lines")
#           )

saveRDS(plot_data, 
        file = glue("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/plot_data_dim{dim}.rds"), 
        compress = "xz" )

ggsave(glue("{output_directory}boxplot_silhouette_dim{dim}.png"),
       boxplot_silhouette,
       width=35, height=20, unit="cm")

# ggsave(glue("{output_directory}boxplot_deconvolution_by_cell_type_dim{dim}.png"),
#        boxplot_deconvolution_by_cell_type,
#        width=35, height=20, unit="cm")
# 
# ggsave(glue("{output_directory}boxplot_deconvolution_by_method_dim{dim}.png"),
#        boxplot_deconvolution_by_method,
#        width=35, height=20, unit="cm")
# 
# pdf(file = glue("{output_directory}benchmark_plot_dim{dim}.pdf"), 
#     paper = "a4r", height = 8.3, width = 11.7) # The height of the plot in inches
# 
# # Step 2: Create the plot with R code
# boxplot_silhouette
# 
# boxplot_deconvolution_by_cell_type
# 
# boxplot_deconvolution_by_method
# 
# # Step 3: Run dev.off() to create the file!
# dev.off()



