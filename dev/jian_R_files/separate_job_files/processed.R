library(tidyverse)
library(tidybulk)
library(ggplot2)

load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")

processed <- counts %>% 
  pivot_longer(contains("level_"), names_to = "level", values_to="cell_type2") %>% 
  mutate(cell_type2 = ifelse(is.na(cell_type2), cell_type, cell_type2)) %>% 
  filter(is.na(cell_type2)==F & is.na(symbol)==F) %>%
  tidybulk(sample, symbol, count) %>%
  nest(data = -level) %>%
  mutate(data = map(data, ~ droplevels(.x))) %>%
  # mutate(data = map(data, ~aggregate_duplicates(.x))) %>%
  # mutate(data = map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
  mutate(data = map(data, ~ .x %>%
                      identify_abundant(factor_of_interest = cell_type2) %>%
                      scale_abundance(sample, symbol, count)
  ))


geneNames <- counts %>%
  select(symbol) %>%
  filter(is.na(symbol)==F) %>%
  distinct() %>%
  pull()