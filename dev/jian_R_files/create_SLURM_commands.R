# Save preliminary files


# Create SLURM commands

# Create commands combinations
# Create stream through argument combination
tibble(
  is_hierarchy = c("hierarchical", "non_hierarchical", NA),
  #contrast = c(mean_contrast, pairwise_contrast, NA),
  contrast_name = c("mean_contrast", "pairwise_contrast", NA),
  #rank = c(rank_edgR_quasi_likelihood, rank_edgR_robust_likelihood_ratio, rank_bayes),
  rank_name = c("edgR", "edgR_robust", "bayes"),
  rank_stat = c("logFC", "PValue", NA),
  selection = c("silhouette", "naive", NA),
  optimisation = c("penalty", "curvature", NA)
) %>% 
  tidyr::expand(is_hierarchy, contrast_name, rank_name, rank_stat, selection, optimisation) %>%
  
  # Drop arguments for some methods
  filter(!( rank_name == "edgR_robust" & rank_stat == "logFC")) %>%
  filter(!(selection == "naive" & optimisation == "curvature")) %>% 
  mutate(rank_stat = map2_chr(
    rank_stat, rank_name,
    ~ if (.y == "bayes"){.x = "_"} else {.x}
  )) %>%
  
  mutate(bayes = map_int(rank_name,  ~ if(.x == "bayes"){1L}else(0L))) %>% 
  
  filter(!(is.na(is_hierarchy) | is.na(contrast_name) | is.na(rank_stat) | is.na(selection) | is.na(optimisation)) ) %>% 
  distinct() %>%
  
  # Add hierarchical argument
  mutate(R_script = if_else(is_hierarchy=="hierarchical", 
                            "produce_benchmark_plot_from_database_plus_tree_hierarchy.R", 
                            "produce_benchmark_plot_from_database_plus_tree_non_hierarchy.R")) %>% 
  unite("R_command", 
        c(R_script, is_hierarchy, contrast_name, rank_name, rank_stat, bayes, selection, optimisation), 
        remove=FALSE, sep=" ") %>% 

  mutate(SLURM_command = glue::glue("sbatch Rscript {R_command}")) %>% 
  pull(SLURM_command) %>%
  write_lines("./dev/SLURM_benchmark_pipeline.sh")
  
  # unite("R_command", 
  #       c("RScript", R_script, is_hierarchy, contrast_name, rank_name, rank_stat, selection, optimisation), 
  #       remove=FALSE, sep=" ")

# Then, you will call
# >sh SLURM_benchmark_pipeline.sh