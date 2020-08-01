# devtools::github_install("stemangiola/nanny@convert-to-S3")
# devtools::github_install("stemangiola/tidybulk@dev")

library(tidyverse)
library(plotly)
library(nanny)
library(ggrepel)
library(GGally)
library(tidyHeatmap)
library(furrr)
plan(multiprocess, workers=5)

# To be loaded after all libraries
library(tidybulk)

load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")

tt =
  counts %>%
  tidybulk(sample, symbol, count) %>%

  # Group by level because otherwise samples are duplicated
  nest(data = -level) %>%

  # Redefine factors inside each level
  mutate(data = future_map(data, ~ droplevels(.x))) %>%

  # Fill missing data. There are many genes that
  # are not shared by the majority of samples
  mutate(data = future_map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%

  # Scale and reduce dimensions
  mutate(data = future_map(
    data,
    ~ .x %>%
      scale_abundance() %>%
      reduce_dimensions(method="MDS") %>%
      reduce_dimensions(method="PCA") %>%
      reduce_dimensions(method="tSNE")
  )) %>%

  # Cluster
  mutate(data = future_map(data,~ cluster_elements(.x, method="SNN")))

counts_imm_epi_de <-
  tt %>%

  # Investigate one level
  filter(level==1) %>%
  unnest(data) %>%

  # investigate one cell type pair
  filter(cell_type %in% c("immune_cell", "epithelial")) %>%
  mutate(cell_type = as.character(cell_type) ) %>%

  #test
  test_differential_abundance( ~ cell_type, sample, symbol, count_scaled) %>%

  # Nest
  nest_subset(data = -symbol) %>%

  # Select markers
  filter(FDR < 0.05 & abs(logFC) > 2) %>%
  filter(logCPM > mean(logCPM)) %>%
  arrange(logFC %>% desc()) %>%
  slice(1:10) %>%
  unnest()

