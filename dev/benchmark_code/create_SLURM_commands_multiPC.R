# Save preliminary files
library(glue)
library(tidyverse)

input_directory = "dev/benchmark_code/"
output_directory = "dev/benchmark_results_multiPC/"
tab = "\t"
# Create SLURM commands

# Create commands combinations
# Create stream through argument combination
table_of_commands = 
  tibble(
    is_hierarchy = c("hierarchical", "non_hierarchical", NA),
    contrast = c("mean_contrast", "pairwise_contrast", NA),
    rank = c("edgR", "edgR_robust", "bayes"),
    rank_stat = c("logFC", "PValue", NA),
    selection = c("silhouette", "naive", NA),
    dims = c(2, 4, 10),
    optimisation = c("penalty", "curvature", NA)
  ) %>% 
  tidyr::expand(is_hierarchy, contrast, rank, rank_stat, selection, optimisation, dims) %>%
  
  # Drop arguments for some methods
  # 1. edgR_robust should only use PValue as rank_stat
  filter(!( rank == "edgR_robust" & rank_stat == "logFC")) %>%
  
  # 2. curvature optimisation is not meaningful for naive selection methods
  filter(!(selection == "naive" & optimisation == "curvature")) %>% 
  
  # 3. the computational burden is too high for non_hierarchical using pairwise_contrast and silhouette selection  
  filter(!(is_hierarchy == "non_hierarchical" & contrast == "pairwise_contrast" & selection == "silhouette")) %>% 
  
  # if rank method is bayes then rank_stat should be left as an absent argument but not filered by !is.na()
  mutate(rank_stat = map2_chr(
    rank_stat, rank,
    ~ if (.y == "bayes"){.x = "_"} else {.x}
  )) %>%
  
  # create an argument for whether the bayes dataset is provided
  mutate(bayes = map_int(rank,  ~ if(.x == "bayes"){1L}else(0L))) %>% 
  
  # remove rows containing na and duplicate streams
  drop_na() %>% 
  distinct() %>% 

  # Define output file
  unite("options", everything(), sep = " ", remove = FALSE) %>% 
  unite("output_file", -c(bayes, options), sep = "_", remove = FALSE) %>% 
  mutate(output_file = glue("{output_directory}{output_file}.rds")) %>% 
  
  # Choose the R_script to run: hierarchical or non-hierarchical
  mutate(R_script = ifelse(is_hierarchy == "hierarchical", 
                           glue("Rscript {input_directory}produce_benchmark_multiPC_hierarchy.R"),
                           glue("Rscript {input_directory}produce_benchmark_multiPC_non_hierarchy.R")
                           )) %>% 
  
  # Add the command line arguments and output file to the R_script
  mutate(command = sprintf("%s:\n\t%s %s %s", output_file, R_script, options, output_file))
  
  # This below produces an error file.
  # mutate(command = sprintf("%s:\n\t%s %s %s > dev/AAA_err.stderr  2>&1", output_file, R_script, options, output_file))
  # alternatively glue works the same way but \n\t cannot be evaluated consecutively
  # mutate(command = glue("{output_file}:\n{tab}{R_script} {options} {output_file}"))

output = sprintf("%s%s", output_directory, dir("dev/benchmark_results_multiPC/"))
  
# pull command and write it to a makeflow file with necessary heading format
table_of_commands %>%
  
  filter(! output_file %in% output) %>% 
  # pull(output_file)
  
  pull(command) %>% 
  
  # Add SLURM requirements
  purrr::prepend("CATEGORY=yes_no_hierarchy\nMEMORY=80000\nCORES=2\nWALL_TIME=172800") %>%
  # purrr::prepend("CATEGORY=yes_no_hierarchy\nMEMORY=30000\nCORES=2\nWALL_TIME=86400") %>% 

  write_lines("dev/benchmark_code/benchmark_multiPC.makeflow")

