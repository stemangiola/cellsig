# The genes selected for each of the 20 signature sizes are stored in all_contrasts_L2

all_contrasts_L2 <- readRDS("dev/all_contrasts_L2.rds")

# The error occurs from sig_size = 3:
# "x Problem with `mutate()` input `distance`. x negative length vectors are not allowed."
# because the distance matrix starts to exceeds the memory for the maximum length of a vector in R.

all_contrasts_L2 %>% 
  slice(1:3) %>% 
  mutate(sil_df = map(sil_df, ~ .x %>% sil_func(LEVEL)))

# When sig_size = 2 (NO ERROR example):
s2 <- all_contrasts_L2 %>% 
  
  # slice the contrast() output at sig_size 2
  slice(2) %>% 
  unnest(sil_df) %>% 
  
  # reduce dimension to get PC values
  nest(pca = - !!as.symbol(pre(LEVEL))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     reduce_dimensions(sample, symbol, count_scaled,
                                       method = "PCA",
                                       action = "add",
                                       transform = log1p) %>%
                     nest(data_symbol = c(symbol, count_scaled))))

y <- s2 %>% 
  mutate(real_size = map_int(pca, ~ .x$data_symbol %>% 
                           map_int(~ n_distinct(.x$symbol)) %>% 
                           unlist() %>% 
                           unique()
                         ))


s2 %>% unnest(pca) %>% pluck("data_symbol", 1)
  select(contains("PC")) %>% 
  str()

a <- s2 %>% 
  mutate(distance = map(pca, ~ .x %>%
                          select(contains("PC")) %>%
                          dist() ))

# when sig_size = 3 (ERROR Example):

s3 <- all_contrasts_L2 %>% 
  
  # slice the contrast() output at sig_size 2
  slice(3) %>% 
  unnest(sil_df) %>% 
  
  # reduce dimension to get PC values
  nest(pca = - !!as.symbol(pre(LEVEL))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     reduce_dimensions(sample, symbol, count_scaled,
                                       method = "PCA",
                                       action = "add",
                                       transform = log1p) ))

s3 %>% unnest(pca) %>% 
  select(contains("PC")) %>% 
  str()


b <- s3 %>% 
  mutate(distance = map(pca, ~ .x %>%
                          select(contains("PC")) %>%
                          dist()))

# when sig_size = 4 (ERROR Example):

s4 <- all_contrasts_L2 %>% 
  
  # slice the contrast() output at sig_size 2
  slice(4) %>% 
  unnest(sil_df) %>% 
  
  # reduce dimension to get PC values
  nest(pca = - !!as.symbol(pre(LEVEL))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     reduce_dimensions(sample, symbol, count_scaled,
                                       method = "PCA",
                                       action = "add",
                                       transform = log1p) ))

s4 %>% unnest(pca) %>% 
  select(contains("PC")) %>% 
  str()

c <- s4 %>% 
  mutate(distance = map(pca, ~ .x %>%
                          select(contains("PC")) %>%
                          factoextra::get_dist(method = "euclidean") ))

# when sig_size = 6 (ERROR Example):

s6 <- all_contrasts_L2 %>% 
  
  # slice the contrast() output at sig_size 2
  slice(6) %>% 
  unnest(sil_df) %>% 
  
  # reduce dimension to get PC values
  nest(pca = - !!as.symbol(pre(LEVEL))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     reduce_dimensions(sample, symbol, count_scaled,
                                       method = "PCA",
                                       action = "add",
                                       transform = log1p) ))

s6 %>% unnest(pca) %>% 
  select(contains("PC")) %>% 
  str()

d <- s6 %>% 
  mutate(distance = map(pca, ~ .x %>%
                          select(contains("PC")) %>%
                          factoextra::get_dist(method = "euclidean") ))

# calculate silhouette score
d %>%  
  mutate(sil = map2(pca, distance,
                    ~ silhouette(as.numeric(as.factor(`$`(.x, !!as.symbol(LEVEL)))), .y) ))
  