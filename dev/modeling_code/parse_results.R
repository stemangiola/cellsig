library(tidyverse)
library(cellsig)
library(tidybulk)
library(tidySummarizedExperiment)
library(magrittr)
library(furrr)
library(data.tree)
library(yaml)
library(glue)
plan(multisession, workers = 10)

# Read arguments
args = commandArgs(trailingOnly=TRUE)
directory_in = args[1]
tree_in_path = args[2]

counts_bayes = 
  dir(directory_in, pattern = "result.rds", full.names = T) %>%
  future_map_dfr(~ {
    x = readRDS(.x) %>% mutate(file = .x)
    attr(x, "fit") = NULL
    attr(x, "rng") = NULL
    x
  }) 


job::job({ counts_bayes %>% saveRDS(glue("{directory}/counts_bayes_{prefix}.rds"), compress = "xz") })


counts_bayes_parsed = 
  counts_bayes %>%
  
  #unite( "sample", c(cell_type, level), remove = FALSE) %>%
  tree_and_signatures_to_database(

    read_yaml("dev/tree_kamran.yaml") %>% as.Node, 

    ., 
    cell_type, 
    cell_type, 
    .feature,  
    `50%`
  )  


job::job({ counts_bayes_parsed %>% saveRDS(glue("{directory}/counts_bayes_{prefix}_parsed.rds"), compress = "xz") })


counts_bayes_imputed = 
  counts_bayes_parsed %>%

  # Convert to SE
  as_SummarizedExperiment(cell_type, .feature, c(`10%`, `50%`,  `90%`, log_mean, log_sd) ) %>%
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts

  when("level_4" %in% colnames(.) ~ impute_missing_abundance(.,~ level_4, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd), ~ (.))) %>%
  when("level_3" %in% colnames(.) ~ impute_missing_abundance(.,~ level_3, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd), ~ (.))) %>%
  when("level_2" %in% colnames(.) ~ impute_missing_abundance(.,~ level_2, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd), ~ (.))) %>%
  when("level_1" %in% colnames(.) ~ impute_missing_abundance(.,~ level_1, .abundance = c(`10%`, `50%`,  `90%`, log_mean, log_sd), ~ (.))) %>% 

  
  # Convert back to tibble
  as_tibble() %>%
  
  # Merge the imputed column

  mutate(.imputed = .imputed | .imputed.1 | .imputed.2 | .imputed.3  ) %>%
  dplyr::select(-c(  .imputed.1 , .imputed.2 ,.imputed.3 )) 

job::job({ counts_bayes_imputed %>% saveRDS(glue("{directory}/counts_bayes_{prefix}_parsed_imputed.rds"), compress = "xz") })

