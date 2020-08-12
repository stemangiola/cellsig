raw_df = readRDS("~/PhD/deconvolution/ARMET/dev/counts_infer_NB.rds")

sample_ct = 
  raw_df %>%
  distinct(sample, `Cell type category`, level) %>%
  filter(`Cell type category` != "house_keeping")

raw_df_ct = 
  raw_df %>%
  mutate(house_keeping = `Cell type category` == "house_keeping" ) %>%
  select(-`Cell type category`) %>%
  left_join(sample_ct) %>%
  mutate(symbol = gsub("house_keeping_", "", symbol))



counts = 
  raw_df_ct %>%
  select(
    sample, 
    level, 
    cell_type = `Cell type category`, 
    cell_type_formatted =  `Cell type formatted`,
    cell_tyoe_original = `Cell type`,
    symbol,
    count,
    house_keeping
  ) %>%
  mutate_if(is.character, as.factor)


save(counts, file="data/counts.rda", compress = "xz")



