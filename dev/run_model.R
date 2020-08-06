library(tidyverse)
library(magrittr)
library(cellsig)


# counts =
#   readRDS(file="~/PhD/deconvolution/ARMET/dev/counts_infer_NB.rds") %>%
#   #left_join( (.) %>% distinct(`symbol original`) %>% mutate(run = sample(1:5, n(), replace = T))) %>%
#   #filter(`Cell type category` == "house_keeping" | run == my_run) %>%
#   mutate(symbol = `symbol original`) %>%
#
#   # Replace the category house_keeping
#   mutate(`house keeping` = `Cell type category` == "house_keeping") %>%
#   rename(temp = `Cell type category`) %>%
#   left_join(
#     (.) %>%
#       filter(temp != "house_keeping") %>%
#       distinct(sample, temp, level) %>%
#       rename(`Cell type category` = temp)
#   ) %>%
#   select(-temp) %>%
#   filter(`Cell type category` %>% is.null %>% `!`) %>%
#
#   # Adapt it as input
#   select(sample, symbol, count, `Cell type category`, level, `count scaled`, `house keeping`)

my_level = 2



fit =
  cellsig::counts %>%

  # # subset
  # nest(data = -c(symbol, house_keeping)) %>%
  # nest(hk = -house_keeping) %>%
  # arrange(house_keeping ) %>%
  # mutate(how_many = c(5000, 500)) %>%
  # mutate(hk = map2(hk, how_many, ~ .x %>% sample_n(.y))) %>%
  # unnest(hk) %>%
  # unnest(data) %>%

  # Infer
  ref_intercept_only(my_level, cores = 10, approximate_posterior = T)


