# How to call this file
# From bash
# Rscript database_plus_tree_as_input_benchmark_plot_as_output.R
#
# Output -> a PDF file
# Input nothing
#
# What does it do
# Produce a pdf running the whole pipeline, including all methods, and testing with sihuette and deconvolution

source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")
# library(furrr)
# plan(multisession, workers = 10)

cellsig_theoretical_transcript_abundance_distribution <- 
  readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/cellsig_theoretical_transcript_abundance_distribution.rds") %>% 
  # take the first row of the duplicated data so that each gene in a cell type has only one set of quantiles
  nest(data = - c(symbol, cell_type)) %>% 
  mutate(data = map(data, ~.x[1, ])) %>% 
  # ensure all cell types have the same set of genes
  add_count(symbol) %>%
  filter(n == max(n)) %>%
  unnest(data)

counts_hierarchy <- readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/counts.rds") %>% 
  scale_input_counts(.is_hierarchy = TRUE)

saveRDS(counts_hierarchy, file = "/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_hierarchy.rds", compress = "xz")

# HIERARCHICAL RESULTS
# Create stream through argument combination
tibble(
  contrast = c(mean_contrast, pairwise_contrast, NA),
  contrast_name = c("mean_contrast", "pairwise_contrast", NA),
  rank = c(rank_edgR_quasi_likelihood, rank_edgR_robust_likelihood_ratio, rank_bayes),
  rank_name = c("edgR", "edgR_robust", "bayes"),
  rank_stat = c("logFC", "PValue", NA),
  selection = c("silhouette", "naive", NA),
  optimisation = c("penalty", "curvature", NA)
  ) %>% 
  tidyr::expand(nesting(contrast, contrast_name), nesting(rank, rank_name), rank_stat, selection, optimisation) %>%
  
  # Drop arguments for some methods
  filter(!( rank_name == "edgR_robust" & rank_stat == "logFC")) %>%
  filter(!(selection == "naive" & optimisation == "curvature")) %>% 
  mutate(rank_stat = map2(
    rank_stat, rank_name,
    ~ if (.y == "bayes"){.x = NULL} else {.x}
    )) %>%
  
  mutate(bayes = map(rank_name, 
                     ~ if(.x == "bayes"){cellsig_theoretical_transcript_abundance_distribution
                       }else(NULL))) %>% 
  
  filter(!(is.na(contrast) | is.na(rank_stat) | is.na(selection) | is.na(optimisation)) ) %>% 
  distinct() %>% 
  
  mutate(benchmark = pmap(
    list(contrast, rank, rank_stat, bayes, selection, optimisation),
    ~ main(counts_hierarchy, .is_hierarchy=TRUE,
            .contrast_method = ..1, .ranking_method = ..2, .rank_stat = ..3, .bayes = ..4,
            .selection_method = ..5, .kmax = 60, .discard_number = 2000, .reduction_method = "PCA",
           .optimisation_method= ..6,
           .is_complete = FALSE)
  )) %>% 
  
  select(-bayes) %>% 
  saveRDS(file = "./dev/signature_hierarchial_methods.rds", compress = "xz")

  # unite("stream", c(contrast, rank_name, rank_stat, selection, optimisation),  sep="_", remove = FALSE)


counts_non_hierarchy <- readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/counts.rds") %>% 
  scale_input_counts(.is_hierarchy = FALSE)

saveRDS(counts_non_hierarchy, 
        file = "/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_non_hierarchy.rds", compress = "xz")

# NON-HIERARCHICAL RESULTS
# Create stream through argument combination
tibble(
  contrast = c(mean_contrast, pairwise_contrast, NA),
  contrast_name = c("mean_contrast", "pairwise_contrast", NA),
  rank = c(rank_edgR_quasi_likelihood, rank_edgR_robust_likelihood_ratio, rank_bayes),
  rank_name = c("edgR", "edgR_robust", "bayes"),
  rank_stat = c("logFC", "PValue", NA),
  selection = c("silhouette", "naive", NA),
  optimisation = c("penalty", "curvature", NA)
) %>% 
  tidyr::expand(nesting(contrast, contrast_name), nesting(rank, rank_name), rank_stat, selection, optimisation) %>%
  
  # Drop arguments for some methods
  filter(!( rank_name == "edgR_robust" & rank_stat == "logFC")) %>%
  filter(!(selection == "naive" & optimisation == "curvature")) %>% 
  mutate(rank_stat = map2(
    rank_stat, rank_name,
    ~ if (.y == "bayes"){.x = NULL} else {.x}
  )) %>%
  
  mutate(bayes = map(rank_name, 
                     ~ if(.x == "bayes"){cellsig_theoretical_transcript_abundance_distribution
                     }else(NULL))) %>% 
  
  filter(!(is.na(contrast) | is.na(rank_stat) | is.na(selection) | is.na(optimisation)) ) %>% 
  distinct() %>% 
  
  mutate(benchmark = pmap(
    list(contrast, rank, rank_stat, bayes, selection, optimisation),
    ~ main(counts_non_hierarchy, .is_hierarchy=FALSE,
           .contrast_method = ..1, .ranking_method = ..2, .rank_stat = ..3, .bayes = ..4,
           .selection_method = ..5, .kmax = 60, .discard_number = 10000, .reduction_method = "PCA",
           .optimisation_method= ..6,
           .is_complete = FALSE)
  )) %>% 
  
  select(-bayes) %>% 
  saveRDS(file = "./dev/signature_non_hierarchial_methods.rds", compress = "xz")

  
 
