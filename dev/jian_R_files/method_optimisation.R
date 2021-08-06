# Optimisation

# Finalised Functions =========================

penalised_silhouette <- function(.plot_data, .penality_rate=0.2) {
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    # filter(real_size > 10) %>%
    
    # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
    mutate(size_rescaled = rescale(real_size)) %>% 
    
    mutate(penalised_silhouette = silhouette - .penality_rate * size_rescaled) %>% 
    
    filter(penalised_silhouette==max(penalised_silhouette)) %>% 
    
    select(real_size, silhouette)
}

ratio <- function(.plot_data) {
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    # filter(real_size > 10) %>%
    
    # calculate the difference between rescaled size of the max sil point with that of all other points
    # the bigger different the better (even for negative numbers)
    mutate(size_gap = real_size[which.max(silhouette)] - real_size) %>% 
    
    # calculate the difference between silhouette score of the max sil point and that of all other points
    # the smaller the better
    mutate(silhouette_gap = max(silhouette) - silhouette) %>% 
    
    # use min_max scaler to rescale size_difference to a scale (between 0 and 1)
    mutate(size_gap_rescaled = rescale(size_gap)) %>%
    
    # use min_max scaler to rescale silhouette score difference to the same scale (between 0 and 1)
    mutate(silhouette_gap_rescaled = rescale(silhouette_gap, c(1, 10))) %>% 
    
    # calculate the ratio between size_gap_rescaled/silhouette_gap, the bigger the better
    mutate(ratio = map2_dbl(
      silhouette_gap_rescaled, size_gap_rescaled, 
      ~ if(.x == 0){0}else{.y / .x}
    )) %>% 
    
    filter(ratio == max(ratio)) %>% 
    
    select(real_size, silhouette)
}


# optimal_size is obsolete!
optimal_size <- function(.plot_data) {
  
  # if the highest silhouette score is the first point then select this point (no left points)
  if (which.max(.plot_data$sil) == 1) {
    op_size <- c(sig_size = .plot_data$sig_size[1], 
                 real_size = .plot_data$real_size[1])
    
    # if the highest silhouette score is the last point then select the best left point (no right points)
  } else if (which.max(.plot_data$sil) == nrow(.plot_data)) {
    
    # find the index of the point that gives the best ratio to the left of max silhouette score point
    lop_index <- which.max(.plot_data$ratio[.plot_data$ratio > 0])
    
    op_size <- c(sig_size = .plot_data$sig_size[lop_index],
                 real_size = .plot_data$real_size[lop_index])
    
  } else {
    
    # find the index that gives the optimal point to the left of max silhouette score point
    lop_index <- which.max(.plot_data$ratio[.plot_data$ratio > 0])
    
    # find the index that gives the optimal point to the right of max silhouette score point
    rop_index <- which.max(.plot_data$sil) +
      which.max(.plot_data$ratio[.plot_data$ratio < 0])
    
    # choose the sizes of the optimal left point if its silhouette score is bigger than the optimal right point
    if (.plot_data$sil[lop_index] >= .plot_data$sil[rop_index]) {
      op_size <- c(sig_size = .plot_data$sig_size[lop_index],
                   real_size = .plot_data$real_size[lop_index])
    } else {
      # choose the max silhouette point if silhouette score of the optim right point is bigger than that of the left
      op_size <- c(sig_size = .plot_data$sig_size[which.max(.plot_data$sil)],
                   real_size = .plot_data$real_size[which.max(.plot_data$sil)])
    }
  }
  return(op_size)
}


# Hierarchical and non-hierarchical in the same plot
# Read in and process data ======================================================

naive_methods <- list.files("naive methods silhouette score data/")
new_methods <- list.files("intermediate_data/", pattern = "^[pm].*(NH|L\\d)\\.rds$")

new <- map_dfr(new_methods, ~ readRDS(paste("intermediate_data", .x, sep = "/")))
naive <- map_dfr(naive_methods, ~ readRDS(paste("naive methods silhouette score data", .x, sep = "/")))

naive <- naive %>% 
  # print(n=30) 
  mutate(method = rep(
    c("mean_contrast.naive.hierarchy",
      "mean_contrast.naive.non_hierarchy",
      "pairwise.naive.hierarchy",
      "pairwise.naive.non_hierarchy"
    ), times=c(14, 1, 14, 1)),
    .before = plot_data
  ) %>% 
  mutate(plot_data = map(
    plot_data, 
    ~.x %>% select(real_size, silhouette)))

new <- new %>% 
  # print(n=30)
  mutate(method = rep(
    c("mean_contrast.silhouette.hierarchy",
      "mean_contrast.silhouette.non_hierarchy",
      "pairwise.silhouette.hierarchy"), times=c(14, 1, 14))
  ) %>% 
  unnest(signature_data) %>% 
  mutate(real_size = map_int(cumulative_signature, ~length(.x))) %>% 
  rename(silhouette = winning_silhouette) %>% 
  select(level, ancestor, real_size, silhouette, method) %>% 
  nest(plot_data = -c(level, ancestor, method))

full_data <- bind_rows(new, naive)

# 
naive <- list.files("topInf_scaleFALSE/unoptimised/", pattern = ".*naive\\..*\\..*")
silhouette <- list.files("topInf_scaleFALSE/unoptimised/", pattern = ".*silhouette\\..*\\..*")

naive_df <- map_dfr(naive, ~ readRDS(paste0("topInf_scaleFALSE/unoptimised/", .x))) %>% 
  mutate(method = rep(str_replace_all(naive, '\\.unOP\\.new\\.rds', ''), c(14, 1, 14, 1)))

o <- rep(str_replace_all(silhouette, '\\.unOP\\.new\\.rds', ''), c(14, 1, 14))
silhouette_df <- map_dfr(silhouette, ~ readRDS(paste0("topInf_scaleFALSE/unoptimised/", .x))) %>% 
  mutate(method = o)
rm(o)

full_df <- silhouette_df %>% 
  bind_rows(naive_df)

# 1 Optimisation by Penalised silhouette ==========================
## 1.1 choose penalty rate===============

penalty_inspect <- function(.plot_data, .penality_rate) {
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    # filter(real_size > 10) %>%
    
    # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
    mutate(size_rescaled = rescale(real_size)) %>% 
    
    mutate(penalised_silhouette = silhouette - .penality_rate * size_rescaled)
}

plot_data <- pairwise.silhouette.hierarchy.unOP %>% 
  pluck("data", 1) # change the number to see the effect of penalty rate for other nodes

penalty <- tibble(penalty_rate = seq(0, 1, 0.1)) %>% 
  mutate(data = map(penalty_rate, ~ penalty_inspect(plot_data, .x))) %>% 
  mutate(optimal_size = map_int(data, ~ with(.x, real_size[which.max(penalised_silhouette)]) )) %>% 
  mutate(optimal_silhouette = map_dbl(data, ~with(.x, max(penalised_silhouette))))

penalty %>% 
  unnest(data) %>% 
  ggplot(aes(real_size, penalised_silhouette)) +
  geom_point(size=0.2) +
  geom_point(
    aes(x = optimal_size, y = optimal_silhouette),
    colour = "red",
    size = 1.5
  ) +
  geom_text(
    aes(x = optimal_size, 
        y = optimal_silhouette,
        label = sprintf("%i , %.3f", optimal_size, optimal_silhouette)
    ),
    position = position_nudge(y = -0.1)
  ) +
  geom_line() +
  facet_wrap(~ penalty_rate)

## 1.2 Output optimal sizes and their silhouette ===============
op <- full_df %>% 
  mutate(optimal_size = map_int(data, ~ penalised_silhouette(.x, 0.8))) %>% 
  mutate(optimal_silhouette = map2_dbl(data, optimal_size, ~ with(.x, silhouette[real_size == .y])))

saveRDS(op, "optimal_size_penalty.rds")

## 1.3 Global Evaluation ======================

full_df %>% 
  unnest(data) %>% 
  ggplot(aes(real_size, silhouette, color = method))+
  geom_point(size=0.1)+
  geom_point(
    data = op,
    aes(x = optimal_size, y = optimal_silhouette),
    colour = "black",
    size=1
  ) +
  geom_line()+
  xlim(0, 100)+
  facet_wrap(~ ancestor) +
  guides(color = guide_legend(
    title.position = "left",
    nrow = 2,
    byrow = TRUE))+
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.title.align = 0.5,
        plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("penalty optimisation, penalty rate = 0.2")

# 2 Optimisation by Ratio ==========================

## 2.1 choose scaling interval===============

ratio_inspect <- function(.plot_data) {
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    # filter(real_size > 10) %>%
    
    # calculate the difference between rescaled size of the max sil point with that of all other points
    # the bigger different the better (even for negative numbers)
    mutate(size_gap = real_size[which.max(silhouette)] - real_size) %>% 
    
    # calculate the difference between silhouette score of the max sil point and that of all other points
    # the smaller the better
    mutate(silhouette_gap = max(silhouette) - silhouette) %>% 
    
    # use min_max scaler to rescale size_difference to a scale (between 0 and 1)
    mutate(size_gap_rescaled = rescale(size_gap)) %>%
    
    # use min_max scaler to rescale silhouette score difference to the same scale (between 0 and 1)
    mutate(silhouette_gap_rescaled = rescale(silhouette_gap, c(1, 10))) %>% 
    
    # calculate the ratio between size_gap_rescaled/silhouette_gap, the bigger the better
    mutate(ratio = map2_dbl(
      silhouette_gap_rescaled, size_gap_rescaled, 
      ~ if(.x == 0){0}else{.y / .x}
    ))
}


ratio <- function(.plot_data) {
  .plot_data %>% 
    
    # too few markers won't be able to resolve cell types in a large mixed cohort hence remove them
    # filter(real_size > 10) %>%
    
    # use min_max scaler to rescale real_size to the same scale as silhouette score (between 0 and 1)
    mutate(size_rescaled = rescale(real_size)) %>% 
    
    mutate(ratio = silhouette / size_rescaled) %>% 
    
    filter(ratio != Inf) %>% 
    
    filter(ratio == max(ratio)) %>% 
    
    select(real_size, silhouette)
}

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 1) %>% 
  ratio()

## 2.2 Output optimal sizes and their silhouette ===============
op2 <- full_data %>% 
  mutate(optimisation = map(plot_data, ~ ratio(.x))) %>% 
  mutate(optimal_size = map_int(optimisation, ~ .x$real_size)) %>% 
  mutate(optimal_silhouette = map_dbl(optimisation, ~ .x$silhouette)) %>% 
  select(-optimisation, -plot_data)

## 2.3 Global Evaluation ======================

full_data %>% 
  unnest(plot_data) %>% 
  ggplot(aes(real_size, silhouette, color = method))+
  geom_point(size=0.1)+
  geom_point(
    data = op2,
    aes(x = optimal_size, y = optimal_silhouette),
    colour = "black",
    size=1
  ) +
  geom_line()+
  xlim(0, 100)+
  facet_wrap(~ ancestor) +
  guides(color = guide_legend(
    title.position = "left",
    nrow = 2,
    byrow = TRUE))+
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.title.align = 0.5,
        plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("ratio optimisation, scale = c(1, 10)")



saveRDS(op2, "optimal_size_ratio.rds")



# Find the gene signature for optimal sizes ========================

# Signatures from optimisation by penalty

tt_non_hierarchy <- scale_input_counts(counts, .is_hierarchy = FALSE)

## new method data

new <- map_dfr(
  new_methods, ~ readRDS(paste("intermediate_data", .x, sep = "/"))
) %>% 
  mutate(method = rep(
    c("mean_contrast.silhouette.hierarchy",
      "mean_contrast.silhouette.non_hierarchy",
      "pairwise.silhouette.hierarchy"), times=c(14, 1, 14)),
    .before = signature_data
  )

## naive method data

naive <- map_dfr(
  naive_methods, 
  ~ readRDS(paste("naive methods silhouette score data", .x, sep = "/"))
) %>% 
  mutate(method = rep(
    c("mean_contrast.naive.hierarchy",
      "mean_contrast.naive.non_hierarchy",
      "pairwise.naive.hierarchy",
      "pairwise.naive.non_hierarchy"
    ), times=c(14, 1, 14, 1)),
    .before = plot_data
  ) %>% 
  rename(signature_data = plot_data)


full_data <- bind_rows(new, naive)


signature_all_methods <-
  
  full_data %>% 
  
  # add the optimal signature size to the original data frame with signatures
  left_join(op, by = c("level", "ancestor", "method")) %>% 
  
  # create a column of signature size to match with the optimal signature size
  mutate(signature_data = map(
    signature_data,
    ~ .x %>% 
      mutate(signature_size = map_int(cumulative_signature, ~length(.x)))
  )) %>% 
  
  # select the signature set having the optimal size
  mutate(signature_data = map2(
    signature_data, optimal_size,
    ~ .x %>% 
      filter(signature_size == .y) %>% 
      pull(cumulative_signature) %>% 
      unlist()
  )) %>% 
  
  # nest by method because we want to combine all the signatures collected from each hierarchy to compare with that from the Nonh-hierarchical approach
  nest(data = -method) %>% 
  
  # combine all the signatures for each method
  mutate(signature = map(
    data,
    ~ .x %>%
      pull(signature_data) %>%
      unlist() %>%
      unique()
  )) %>% 
  
  # calculate silhouette score for signatures collected
  mutate(silhouette = map(
    signature,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>%
      sil_func("root", METHOD) %>%
      select(reduced_dimensions, silhouette, real_size)
  )) %>% 
  
  unnest(silhouette) %>% 
  
  # remove unnecessary column
  select(-data)


# add signatures selected from cibersort
signature_all_methods <- bind_rows(signature_all_methods, cibersort_silhouette) %>% 
  arrange(desc(silhouette))

# summary bar plot comparing all methods using silhouette score

signature_all_methods %>% 
  ggplot(aes(method, silhouette, fill = method)) +
  geom_col() +
  geom_text(aes(label = round(silhouette, 3)), vjust = 1.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("All methods comparison using silhouette score")

# PCA plot 
signature_all_methods %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(PC1, PC2, color = root), label=sample) +
  geom_point() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mean_contrast.silhouette.non_hierarchy")

saveRDS(tt_non_hierarchy, "tt_non_hierarchy.rds", compress = "xz")



# x <- tt_hierarchy %>%
#   filter(symbol %in% sig1) %>%
#   count(sample, symbol)
#   sil_func0(METHOD)
# sig1 <- signature_all_methods %>% pluck("signature", 1)
# sig2 <- signature_all_methods %>% pluck("signature", 2)

# # Combine preprocessed data(tt_L1-5) from all hierarchies for filtering signatures
# 
# tt_L_filenames <- list.files("intermediate_data/", pattern = "^t.*L")
# 
# tt_hierarchy <-
#   
#   # load all tt_L1-5 into a data frame (takes long time)
#   map_dfr(
#     tt_L_filenames,
#     ~ readRDS(paste0("intermediate_data/", .x))) %>% 
#   
#   # unify columns
#   mutate(level_0 = "cell") %>% 
#   select(level_0, data) %>% 
#   
#   # select only distinct symbol, sample, count_scaled needed for reduce_dimensions() and silhouette score
#   unnest(data) %>% 
#   distinct(symbol, sample, count_scaled, .keep_all = TRUE) %>% 
#   
#   # select only columns needed
#   select(level_0, symbol, sample, cell_type, count_scaled)
# 
# saveRDS(tt_hierarchy, "tt_hierarchy.rds")


# Manual optimisation ==================
mean_contrast.silhouette.hierarchy.unOP <- 
  readRDS("topInf_scaleFALSE/unoptimised/mean_contrast.silhouette.hierarchy.unOP.new.rds")

# level 1; ancestor: root
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 1) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 1) %>% 
               filter(real_size == 30),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: root, size  = 30")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 1) %>% 
  filter(real_size == 30) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 1) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_1), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("root, size = 30")


# level_2; ancestor: immune cell
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 2) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 2) %>% 
               filter(real_size == 54),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: immune cell, size  = 54")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 2) %>% 
  filter(real_size == 54) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 2) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_2), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("immune, size = 54")


# level_3; ancestor: mono_derived
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 3) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 3) %>% 
               filter(real_size == 12),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: mono_derived, size  = 12")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 3) %>% 
  filter(real_size == 12) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 3) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_3), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("mono_derived, size = 12")

# level_3; ancestor: b_cell
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 4) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 4) %>% 
               filter(real_size == 13),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: t cell, size  = 13")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 4) %>% 
  filter(real_size == 13) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 4) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_3), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("b_cell, size = 13")


# level_3; ancestor: b_cells
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 5) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 5) %>% 
               filter(real_size == 26),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: b_cell, size  = 26")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 5) %>% 
  filter(real_size == 26) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 5) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_3), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("b_cell, size = 26")


# level_3; ancestor: b_cell
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 6) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 6) %>% 
               filter(real_size == 5),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: b_cell, size  = 5")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 6) %>% 
  filter(real_size == 5) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 6) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_3), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("b_cell, size = 5")


# level_3; ancestor: NK
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 7) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 7) %>% 
               filter(real_size == 10),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: NK, size  = 10")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 7) %>% 
  filter(real_size == 10) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 7) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_3), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("NK, size = 10")


# level_4; ancestor: t_CD4
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 8) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 8) %>% 
               filter(real_size == 30),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: t_CD4, size  = 30")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 8) %>% 
  filter(real_size == 30) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 8) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_4), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("t_CD4, size = 30")


# level_4; ancestor: macrophage
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 9) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 9) %>% 
               filter(real_size == 9),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: macrophage, size  = 9")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 9) %>% 
  filter(real_size == 9) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 9) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_4), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("macrophage, size = 9")


# level_4; ancestor: t_CD8
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 10) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 10) %>% 
               filter(real_size == 13),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: t_CD8, size  = 13")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 10) %>% 
  filter(real_size == 13) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 10) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_4), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("t_CD8, size = 13")


# level_4; ancestor: nk_primed
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 11) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 11) %>% 
               filter(real_size == 12),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: nk_primed, size  = 12")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 11) %>% 
  filter(real_size == 12) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 12) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_4), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("nk_primed, size = 12")


# level_5; ancestor: t_CD4_memory
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 12) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 12) %>% 
               filter(real_size == 3),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: t_CD4_memory, size  = 3")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 12) %>% 
  filter(real_size == 3) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 13) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_5), label = sample) +
  geom_point() +
  # stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("t_CD4_memory, size = 3")


# level_5; ancestor: t_helper
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 13) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 13) %>% 
               filter(real_size == 10),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: t_helper, size  = 10")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 13) %>% 
  filter(real_size == 10) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 14) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_5), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("t_helper, size = 10")

# level_5; ancestor: t_helper
mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 14) %>% 
  ggplot(aes(real_size, silhouette)) +
  geom_point() +
  geom_line() +
  geom_point(data = mean_contrast.silhouette.hierarchy.unOP %>% 
               pluck("data", 14) %>% 
               filter(real_size == 7),
             aes(x = real_size, y = silhouette),
             color = "red") +
  ggtitle("mean contrast silhouette hierarchy; ancestor: t_helper, size  = 7")

mean_contrast.silhouette.hierarchy.unOP %>% 
  pluck("data", 14) %>% 
  filter(real_size == 7) %>% 
  mutate(pca = map(
    cumulative_signature,
    ~ contrast_MC_H %>% 
      pluck("data", 15) %>% 
      filter(symbol %in% .x) %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  )) %>% 
  pluck("pca", 1) %>% 
  ggplot(aes(PC1, PC2, colour = level_5), label = sample) +
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("t_helper, size = 7")

optimal_size.man <- c(30, 54, 12, 13, 26, 5, 10, 30, 9, 13, 12, 3, 10, 7)

op.man <- mean_contrast.silhouette.hierarchy.unOP %>% 
  mutate(optimal_size = optimal_size.man) %>% 
  unnest(data) %>% 
  filter(real_size == optimal_size)

nodal_signature <- list()

for (i in 1:14){
  nodal_signature[[i]] <- op.man$cumulative_signature[[i]]
}

names(nodal_signature) <- op.man$ancestor

nodal_signature$root <- nodal_signature %>% flatten_chr() %>% unique()
nodal_signature$immune_cell <- nodal_signature[-1] %>% flatten_chr() %>% unique()
nodal_signature$mono_derived <- c(nodal_signature$mono_derived, nodal_signature$macrophage) %>% unique()
nodal_signature$t_cell <- c(nodal_signature$t_cell, nodal_signature$t_CD4, 
                            nodal_signature$t_CD8, nodal_signature$t_CD4_memory, 
                            nodal_signature$t_CD8_memory, nodal_signature$t_helper) %>% unique()
nodal_signature$natural_killer <- c(nodal_signature$natural_killer, nodal_signature$nk_primed) %>% unique()
nodal_signature$t_CD4 <- c(nodal_signature$t_CD4, 
                           nodal_signature$t_CD4_memory, 
                           nodal_signature$t_helper) %>% unique()
nodal_signature$t_CD8 <- c(nodal_signature$t_CD8, nodal_signature$t_CD8_memory) %>% unique()

cell_types <- tt_non_hierarchy %>% 
  unnest(tt) %>% 
  unnest(data) %>% 
  pull(cell_type) %>% 
  unique()

nodal_cell_types <- list(
  root = c("epithelial", "endothelial", "fibroblast", "monocyte", "macrophage_M2", 
           "neutrophil", "b_naive", "t_CD4_memory_central", "t_CD8_memory_central", 
           "t_CD4_memory_effector", "t_CD8_memory_effector", "macrophage_M1", 
           "eosinophil", "b_memory", "t_reg", "t_CD8_naive", "t_gamma_delta", 
           "t_helper_h1", "t_helper_h17", "t_helper_h2", "nk_resting", "nk_primed_IL2_PDGFD",
           "nk_primed_IL2", "mast_cell"),
  immune_cell = c("monocyte", "macrophage_M2", "neutrophil", "b_naive", 
                  "t_CD4_memory_central", "t_CD8_memory_central", 
                  "t_CD4_memory_effector", "t_CD8_memory_effector", "macrophage_M1", 
                  "eosinophil", "b_memory", "t_reg", "t_CD8_naive", "t_gamma_delta", 
                  "t_helper_h1", "t_helper_h17", "t_helper_h2", "nk_resting", "nk_primed_IL2_PDGFD",
                  "nk_primed_IL2", "mast_cell"),
  mono_derived = c("macrophage_M1", "macrophage_M2", "monocyte"),
  t_cell = c("t_CD4_memory_central", "t_CD4_memory_effector", 
             "t_CD8_memory_central", "t_CD8_memory_effector", "t_CD8_naive",   
             "t_gamma_delta",  "t_reg", "t_helper_h1", "t_helper_h17", "t_helper_h2"),
  granulocyte = c("neutrophil", "eosinophil",  "mast_cell"),
  b_cell = c("b_naive", "b_memory"),
  natural_killer = c("nk_resting", "nk_primed_IL2_PDGFD", "nk_primed_IL2"),
  t_CD4 = c("t_CD4_memory_central", "t_CD4_memory_effector", "t_reg", "t_helper_h1",
            "t_helper_h17", "t_helper_h2"),
  macrophage = c("macrophage_M1", "macrophage_M2"),
  t_CD8 = c("t_CD8_memory_central", "t_CD8_memory_effector", "t_CD8_naive"),
  nk_primed = c("nk_primed_IL2_PDGFD", "nk_primed_IL2"),
  t_CD4_memory = c("t_CD4_memory_central", "t_CD4_memory_effector"),
  t_CD8_memory = c("t_CD8_memory_central", "t_CD8_memory_effector"),
  t_helper = c("t_helper_h1", "t_helper_h17", "t_helper_h2", "t_reg")
)

nodal_signature_all <- mean_contrast.silhouette.hierarchy.unOP %>% 
  mutate(optim_size = optimal_size.man) %>% 
  unnest(data) %>% 
  filter(real_size == optim_size) %>% 
  select(-c(new_challengers, winner, optim_size)) %>% 
  mutate(nodal_signature = nodal_signature) %>% 
  mutate(nodal_cell_type = nodal_cell_types) %>% 
  mutate(reduced_dimensions = map2(
    nodal_signature, nodal_cell_type,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      filter(cell_type %in% .y)
  )) %>% 
  mutate(reduced_dimensions = map(
    reduced_dimensions,
    ~ .x %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  ))


nodal_signature_all %>% 
  pluck("reduced_dimensions", 1) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("root, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 2) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("immune_cell, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 3) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("mono_derived, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 4) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("t_cell, 0.2")


nodal_signature_all %>% 
  pluck("reduced_dimensions", 5) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("granulocyte, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 6) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("b_cell, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 7) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("natural_killer, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 8) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("t_CD4, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 9) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("macrophage, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 10) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("t_CD8, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 11) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("nk_primed, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 12) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("t_CD4_memory, 0.2")

nodal_signature_all %>% 
  pluck("reduced_dimensions", 13) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("t_CD8_memory, 0.2")


nodal_signature_all %>% 
  pluck("reduced_dimensions", 14) %>% 
  ggplot(aes(PC1, PC2, colour = cell_type), label = sample) +
  geom_point() +
  stat_ellipse(type = 't') +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("t_helper, 0.2")

# contrast optim_size selected by setting penalty rate to 0.2
mean_contrast.silhouette.hierarchy <- mean_contrast.silhouette.hierarchy.unOP %>% 
  do_optimisation("penalised")

nodal_signature <- list()

for (i in 1:14){
  nodal_signature[[i]] <- mean_contrast.silhouette.hierarchy$cumulative_signature[[i]]
}

names(nodal_signature) <- mean_contrast.silhouette.hierarchy$ancestor

nodal_signature$root <- nodal_signature %>% flatten_chr() %>% unique()
nodal_signature$immune_cell <- nodal_signature[-1] %>% flatten_chr() %>% unique()
nodal_signature$mono_derived <- c(nodal_signature$mono_derived, nodal_signature$macrophage) %>% unique()
nodal_signature$t_cell <- c(nodal_signature$t_cell, nodal_signature$t_CD4, 
                            nodal_signature$t_CD8, nodal_signature$t_CD4_memory, 
                            nodal_signature$t_CD8_memory, nodal_signature$t_helper) %>% unique()
nodal_signature$natural_killer <- c(nodal_signature$natural_killer, nodal_signature$nk_primed) %>% unique()
nodal_signature$t_CD4 <- c(nodal_signature$t_CD4, 
                           nodal_signature$t_CD4_memory, 
                           nodal_signature$t_helper) %>% unique()
nodal_signature$t_CD8 <- c(nodal_signature$t_CD8, nodal_signature$t_CD8_memory) %>% unique()


nodal_signature_all <- mean_contrast.silhouette.hierarchy %>% 
  select(-c(new_challengers, winner, optimal_size)) %>% 
  mutate(nodal_signature = nodal_signature) %>% 
  mutate(nodal_cell_type = nodal_cell_types) %>% 
  mutate(reduced_dimensions = map2(
    nodal_signature, nodal_cell_type,
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      filter(cell_type %in% .y)
  )) %>% 
  mutate(reduced_dimensions = map(
    reduced_dimensions,
    ~ .x %>% 
      reduce_dimensions(sample, symbol, count_scaled,
                        method = METHOD,
                        action = "get",
                        top = Inf,
                        scale = FALSE,
                        .dims = 2)
  ))

# boxplot comparing manual signature for mean_contrast.silhouette.hierarchy vs 0.2 signature for all other methods

naive <- list.files("topInf_scaleFALSE/", pattern = ".*naive\\..*\\..*")
silhouette <- list.files("topInf_scaleFALSE/", pattern = ".*silhouette\\..*\\..*")

silhouette[1] <- "op.man.rds"
silhouette

silhouette_score <- function(.reduced_dimensions, .distance, .level){
  
  .reduced_dimensions %>% 
    
    pull(!!as.symbol(.level)) %>% 
    
    as.factor() %>% 
    
    as.numeric() %>% 
    
    silhouette(.distance) %>% 
    
    summary()
  
}

silhouette_function <- function(.selected, .reduction_method){
  
  .selected %>% 
    
    # reduce dimensions
    mutate(reduced_dimensions = map2(
      markers, level, 
      ~ dimension_reduction(.x, .y, .reduction_method)
    )) %>% 
    
    # calculate distance matrix using PC1 & PC2
    mutate(distance = map(
      reduced_dimensions,
      ~ distance_matrix(.x, .reduction_method)
    )) %>% 
    
    # calculate silhouette score
    mutate(silhouette = pmap(
      list(reduced_dimensions, distance, level),
      ~ silhouette_score(..1, ..2, ..3)
    )) %>% 
    
    # remove unnecessary columns
    select(-c(markers, distance))
  
}

yy <- full_df %>% 
  nest(signature = -method) %>% 
  mutate(signature = map(signature, ~.x %>% pull(signature) %>% unlist() %>% unique())) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score for all signatures combined in each method
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

cibersortx <- readRDS("topInf_scaleFALSE/cibersortx.new.rds")
cibersort_signature <- cibersortx$signature[[1]]

cibersortx <- tibble(method = "cibersortx") %>% 
  mutate(signature = list(cibersort_signature)) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

yy <- yy %>% bind_rows(cibersortx)

yy %>% 
  mutate(cluster.silhouette = map(silhouette, ~ .x$clus.avg.widths)) %>% 
  mutate(avg.silhouette = map_dbl(silhouette, ~ .x$avg.width)) %>% 
  select(-c(reduced_dimensions, silhouette)) %>% 
  unnest(cluster.silhouette) %>% 
  ggplot(aes(x=reorder(method, avg.silhouette), y=cluster.silhouette, colour=method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  theme(axis.text.x = element_blank())

# Kernel smooth optimisation =============
library(KernSmooth)
library(splus2R)

mean_contrast.silhouette.hierarchy.unOP <- 
  readRDS("topInf_scaleFALSE/unoptimised/mean_contrast.silhouette.hierarchy.unOP.new.rds")

pairwise.silhouette.hierarchy.unOP <- 
  readRDS("topInf_scaleFALSE/unoptimised/pairwise.silhouette.hierarchy.unOP.new.rds")

mean_contrast.naive.hierarchy.unOP <- 
  readRDS("topInf_scaleFALSE/unoptimised/mean_contrast.naive.hierarchy.unOP.new.rds")

pairwise.naive.hierarchy.unOP <- 
  readRDS("topInf_scaleFALSE/unoptimised/pairwise.naive.hierarchy.unOP.new.rds")

curvature <- function(.drv1, .drv2){
  abs(.drv2) / (1 + .drv1^2)^(3/2)
}

curvature_of_kernel_smoothed_trend <- function(.plot_data, 
                                               .kernel = "normal", 
                                               .bandwidth = 0.05, 
                                               .gridsize = 100){
  .plot_data %>% 
    
    mutate(data = map(
      data,
      ~ .x %>% 
        mutate(size.rescaled = rescale(real_size))
    )) %>% 
    
    mutate(smoothed.estimate = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette, 
                drv = 0L, degree=2, kernel = .kernel, 
                bandwidth = .bandwidth, gridsize = .gridsize) %>% 
        as_tibble() %>% 
        `colnames<-`(c("grid", "estimate"))
    )) %>% 
    
    mutate(first.derivative = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette, 
                drv = 1L, degree=2, kernel = .kernel, 
                bandwidth = .bandwidth, gridsize = .gridsize) %>% 
        as_tibble() %>% 
        `colnames<-`(c("grid", "deriv1"))
    )) %>% 
    
    mutate(second.derivative = map(
      data,
      ~ locpoly(.x$size.rescaled, .x$silhouette, 
                drv = 2L, degree=2, kernel = .kernel, 
                bandwidth = .bandwidth, gridsize = .gridsize) %>% 
        as_tibble() %>% 
        `colnames<-`(c("grid", "deriv2"))
    )) %>% 
    
    mutate(smoothed = pmap(
      list(smoothed.estimate, first.derivative, second.derivative),
      ~..1 %>% 
        left_join(..2, by = "grid") %>% 
        left_join(..3, by = "grid")
    )) %>% 
    
    mutate(smoothed = map(
      smoothed,
      ~ .x %>% 
        mutate(curvature = map2_dbl(
          deriv1, deriv2,
          ~ curvature(.x, .y)
        ))
    )) %>% 
    
    select(-c(smoothed.estimate, first.derivative, second.derivative)) %>% 
    
    # optimal_size
    mutate(optimal_size = map_dbl(
      smoothed,
      ~ with(.x, grid[which(peaks(curvature))[which.max(curvature[peaks(curvature)])]]) 
    )) %>% 
    
    mutate(optimal_size = map2_dbl(
      data, optimal_size,
      ~ with(.x, size.rescaled[which.min(abs(size.rescaled - .y))])
    )) %>% 
    
    mutate(optimal_size = map2_int(
      data, optimal_size,
      ~ .x %>% 
        with(real_size[size.rescaled == .y])
    )) %>% 
    
    mutate(optimal_size = ifelse(optimal_size<10, 10, optimal_size))
  
}

x <- mean_contrast.silhouette.hierarchy.unOP %>% 
  
  mutate(data = map(
    data,
    ~ .x %>% 
      mutate(size.rescaled = rescale(real_size))
  )) %>% 
  
  mutate(smoothed.estimate = map(
    data,
    ~ locpoly(.x$size.rescaled, .x$silhouette, 
              drv = 0L, degree=2, kernel = "normal", 
              bandwidth = 0.05, gridsize = 100) %>% 
      as_tibble() %>% 
      `colnames<-`(c("grid", "estimate"))
  )) %>% 
  
  mutate(first.derivative = map(
    data,
    ~ locpoly(.x$size.rescaled, .x$silhouette, 
              drv = 1L, degree=2, kernel = "normal", 
              bandwidth = 0.05, gridsize = 100) %>% 
      as_tibble() %>% 
      `colnames<-`(c("grid", "deriv1"))
  )) %>% 
  
  mutate(second.derivative = map(
    data,
    ~ locpoly(.x$size.rescaled, .x$silhouette, 
              drv = 2L, degree=2, kernel = "normal", 
              bandwidth = 0.05, gridsize = 100) %>% 
      as_tibble() %>% 
      `colnames<-`(c("grid", "deriv2"))
  )) %>% 
  
  mutate(smoothed = pmap(
    list(smoothed.estimate, first.derivative, second.derivative),
    ~..1 %>% 
      left_join(..2, by = "grid") %>% 
      left_join(..3, by = "grid")
  )) %>% 
  
  mutate(smoothed = map(
    smoothed,
    ~ .x %>% 
      mutate(curvature = map2_dbl(
        deriv1, deriv2,
        ~ curvature(.x, .y)
      ))
  )) %>% 
  
  select(-c(smoothed.estimate, first.derivative, second.derivative))

p1 <- x %>% 
  pluck("data", 1) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 1), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 1),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 1),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 1) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 1) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.14, y=0.7, label="40", color="tomato") +
  ggtitle("root node, manual.size=30")


p2 <- x %>% 
  pluck("data", 2) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 2), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 2),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 2),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 2) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 2) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.18, y=0.65, label="72", color="tomato") +
  ggtitle("immune node, manual.size=54")

p3 <- x %>% 
  pluck("data", 3) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 3), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 3),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 3),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 3) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 3) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.13, y=0.62, label="13", color="tomato") +
  ggtitle("mono_derived, manual.size=12")

p4 <- x %>% 
  pluck("data", 4) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 4), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 4),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 4),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 4) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 4) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.2, y=0.55, label="14", color="tomato") +
  ggtitle("t_cell, manual.size=13")


p5 <- x %>% 
  pluck("data", 5) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 5), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 5),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 5),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 5) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 5) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.05, y=0.33, label="5", color="tomato") +
  ggtitle("granulocyte, manual.size=26")

p6 <- x %>% 
  pluck("data", 6) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 6), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 6),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 6),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 6) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 6) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.08, y=0.82, label="6", color="tomato") +
  ggtitle("b_cell, manual.size=5")


p7 <- x %>% 
  pluck("data", 7) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 7), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 7),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 7),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 7) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 7) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.03, y=0.85, label="3", color="tomato") +
  ggtitle("natural_killer, manual.size=10")


p8 <- x %>% 
  pluck("data", 8) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 8), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 8),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 8),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 8) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 8) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.04, y=0.46, label="5", color="tomato") +
  ggtitle("t_CD4, manual.size=30")


p9 <- x %>% 
  pluck("data", 9) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 9), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 9),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 9),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 9) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 9) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.13, y=0.82, label="9", color="tomato") +
  ggtitle("macrophage, manual.size=9")

p10 <- x %>% 
  pluck("data", 10) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 10), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 10),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 10),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 10) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 10) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.28, y=0.83, label="14", color="tomato") +
  ggtitle("t_CD8, manual.size=13")


p11 <- x %>% 
  pluck("data", 11) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 11), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 11),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 11),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 11) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 11) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.07, y=0.65, label="10", color="tomato") +
  ggtitle("nk_primed, manual.size=12")


p12 <- x %>% 
  pluck("data", 12) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 12), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 12),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 12),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 12) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 12) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.1, y=4, label="3", color="tomato") +
  ggtitle("t_CD4_memory, manual.size=3")


p13 <- x %>% 
  pluck("data", 13) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 13), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 13),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 13),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 13) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 13) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.08, y=0.42, label="3", color="tomato") +
  ggtitle("t_CD8_memory, manual.size=10")


p14 <- x %>% 
  pluck("data", 14) %>% 
  ggplot(aes(size.rescaled, silhouette)) +
  geom_line() +
  # geom_point(colour="blue", alpha=0.5) + 
  geom_point(aes(grid, estimate), data=pluck(x, "smoothed", 14), colour="black") +
  geom_point(aes(grid, deriv1/10), data=pluck(x, "smoothed", 14),  colour="dodgerblue") +
  geom_point(aes(grid, rescale(curvature)), data=pluck(x, "smoothed", 14),  colour="tomato")+
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 14) %>% filter(curvature==max(curvature[peaks(curvature)])),
             colour = "tomato"
  ) +
  geom_point(aes(grid, estimate), 
             data=pluck(x, "smoothed", 14) %>% 
               with(.[which.min(abs(deriv1-1)), ]),
             colour = "dodgerblue"
  ) +
  # annotate(geom="text", x=0.1, y=0.6, label="3", color="tomato") +
  ggtitle("t_helper, manual.size=5")


p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+
  plot_layout(ncol=4, nrow=4)+
  plot_annotation(
    title = "mean_contrast.silhouette.hierarchy, bandwidth=0.05"
  )

x %>% 
  # manual.op.size
  mutate(manual.op.size = optimal_size.man) %>%
  
  # deriv1.op.size
  mutate(deriv1.op.size = map_dbl(
    smoothed,
    ~ with(.x, grid[which.min(abs(deriv1-1))])
  )) %>%
  
  mutate(
    deriv1.op.size = map2_dbl(
      data, deriv1.op.size,
      ~ with(.x, size.rescaled[which.min(abs(size.rescaled - .y))])
    )
  ) %>%
  mutate(deriv1.op.size = map2_int(
    data, deriv1.op.size,
    ~ with(.x, real_size[size.rescaled==.y])
  )) %>%
  
  # curvature.op.size
  mutate(curvature.op.size = map_dbl(
    smoothed,
    ~ with(.x, grid[which(peaks(curvature))[which.max(curvature[peaks(curvature)])]]) 
  )) %>% 
  
  mutate(curvature.op.size = map2_dbl(
    data, curvature.op.size,
    ~ with(.x, size.rescaled[which.min(abs(size.rescaled - .y))])
  )) %>% 
  
  mutate(curvature.op.size = map2_int(
    data, curvature.op.size,
    ~ .x %>% 
      with(real_size[size.rescaled == .y])
  )) %>% 
  
  # op.size when penalty rate = 0.8
  mutate(penalty0.8.op.size = map_int(
    data,
    ~.x %>% 
      penalised_silhouette(.penalty_rate = 0.8)
  ))

# comparison ==========================

# import summary data from all methods
naive <- list.files("topInf_scaleFALSE/", pattern = ".*naive\\..*\\..*")
silhouette <- list.files("topInf_scaleFALSE/", pattern = ".*silhouette\\..*\\..*")

naive_df <- map_dfr(naive, ~ readRDS(paste0("topInf_scaleFALSE/", .x))) %>% 
  select(level, ancestor, real_size, signature, silhouette) %>% 
  mutate(method = rep(str_replace_all(naive, '\\.rds', ''), c(14, 1, 14, 1)))

o <- rep(str_replace_all(silhouette, '\\.rds', ''), c(14, 1, 14))
silhouette_df <- map_dfr(silhouette, ~ readRDS(paste0("topInf_scaleFALSE/", .x))) %>% 
  select(level, ancestor, real_size, signature=cumulative_signature, silhouette) %>% 
  mutate(method = o)
rm(o)

mean_contrast.silhouette.hierarchy.curvature <- mean_contrast.silhouette.hierarchy.unOP %>% 
  curvature_of_kernel_smoothed_trend() %>% 
  unnest(data) %>% 
  filter(real_size == curvature.op.size) %>% 
  rename(optimal_size = curvature.op.size) %>% 
  select(level, ancestor, real_size, signature=cumulative_signature, silhouette) %>% 
  mutate(method = "mean_contrast.silhouette.hierarchy.curvature")

op.man <- mean_contrast.silhouette.hierarchy.unOP %>% 
  mutate(optimal_size = optimal_size.man) %>% 
  unnest(data) %>% 
  filter(real_size == optimal_size) %>% 
  select(level, ancestor, real_size, signature=cumulative_signature, silhouette) %>% 
  mutate(method = "mean_contrast.silhouette.hierarchy.manual")

full_df <- silhouette_df %>% 
  bind_rows(mean_contrast.silhouette.hierarchy.curvature, op.man, naive_df)


all_methods_silhouette <- full_df %>% 
  nest(signature = -method) %>% 
  mutate(signature = map(signature, ~.x %>% pull(signature) %>% unlist() %>% unique())) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score for all signatures combined in each method
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

cibersortx <- readRDS("topInf_scaleFALSE/cibersortx.new.rds")
cibersort_signature <- cibersortx$signature[[1]]

cibersortx <- tibble(method = "cibersortx") %>% 
  mutate(signature = list(cibersort_signature)) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

all_methods_comparison <- all_methods_silhouette %>% 
  bind_rows(cibersortx) %>% 
  mutate(cluster.silhouette = map(silhouette, ~ .x$clus.avg.widths)) %>% 
  mutate(avg.silhouette = map_dbl(silhouette, ~ .x$avg.width)) %>% 
  select(-c(reduced_dimensions, silhouette)) %>% 
  arrange(desc(avg.silhouette)) %>% 
  mutate(method = str_remove_all(method, "\\.new"))


all_methods_comparison %>%
  unnest(cluster.silhouette) %>% 
  ggplot(aes(x=reorder(method, avg.silhouette), y=cluster.silhouette, colour=method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  theme(axis.text.x = element_blank())

# silhouette

mix100 <- readRDS("intermediate_data/mix100.rds")

# mix100 <- tibble(mixture_ID = 1:100) %>%
#   
#   # mix
#   mutate(mix = map(mixture_ID, ~ {
#     proportions = 
#       gtools::rdirichlet(1, rep(1, length(unique(mix_base$cell_type)))) %>%
#       as.data.frame() %>%
#       setNames(unique(mix_base$cell_type))
#     
#     cellsig::generate_mixture_from_proportion_matrix(mix_base, proportions)
#   }))

deconvolution_all_methods <- mix100 %>% 
  tidyr::expand(
    nesting(method=all_methods_comparison$method, 
            signature=all_methods_comparison$signature),
    nesting(mixture_ID, mix)
  ) %>% 
  
  # reference
  mutate(reference = map(
    signature,
    
    # filter out data for signature genes
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      
      # reshape the input matrix for deconvolve_cellularity():
      select(symbol, cell_type, sample, count_scaled) %>% 
      group_by(symbol, cell_type) %>% 
      summarise(count_scaled_median = median(count_scaled)) %>% 
      ungroup() %>% 
      pivot_wider(id_cols = symbol, names_from = cell_type, values_from = count_scaled_median) %>% 
      tidybulk::as_matrix(rownames = symbol) # must be a matrix
  )) %>% 
  
  # deconvolution
  mutate(deconvolution = map2(
    mix, reference,
    ~ tidybulk::deconvolve_cellularity(
      .x,
      replicate, symbol, count_mix,
      reference = .y,
      method = "llsr", 
      prefix = "llsr_", 
      action = "get") %>% 
      
      pivot_longer(cols=starts_with("llsr_"), 
                   names_prefix ="llsr_", 
                   names_to="cell_type", 
                   values_to="estimated_proportion") %>%
      
      left_join(.x %>% 
                  unnest(data_samples) %>% 
                  distinct(replicate, cell_type, proportion))
  ))

# mse <- deconvolution_all_methods %>% 
#   # calculate deconvolution MSE for each method
#   mutate(MSE = map_dbl(
#     deconvolution,
#     ~ mean((.x$estimated_proportion - .x$proportion)^2)
#   )) %>% 
#   nest(MSE.macro = -c(method, signature)) %>% 
#   mutate(MSE.macro = map_dbl(MSE.macro, ~ mean(.x$MSE))) %>% 
#   mutate(method = map_chr(method, ~ .x %>% str_replace_all("\\.new", "") )) %>% 
#   arrange(MSE.macro)


deconvolution_all_methods %>% 
  select(method, mixture_ID, deconvolution) %>% 
  mutate(method = str_remove_all(method, "\\.new") ) %>% 
  
  # macro.mse for ordering
  mutate(MSE = map_dbl(
    deconvolution,
    ~ mean((.x$estimated_proportion - .x$proportion)^2)
  )) %>% 
  # nest(data=-method) %>% 
  # mutate(macro.mse = map_dbl(data, ~ mean(.x$MSE))) %>%
  # unnest(data) %>% 
  
  # mse by cell type
  unnest(deconvolution) %>% 
  mutate(squared.error = (estimated_proportion - proportion)^2) %>% 
  nest(data = -c(method, cell_type)) %>% 
  mutate(mse.cell = map_dbl(data, ~ mean(.x$squared.error))) %>%
  # mutate(macro.mse = map_dbl(data, ~ mean(.x$MSE))) %>% 
  group_by(method) %>% 
  summarise(method, cell_type, data, mse.cell, median = median(mse.cell)) %>% 
  ungroup() %>% 
  
  # ggplot(aes(x=reorder(method, -macro.mse), y=log(mse.cell))) +
  ggplot(aes(x=reorder(method, -median), y=log(mse.cell))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = cell_type), 
              position=position_jitter(0.2)) +
  # theme(axis.text.x = element_blank())
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))

# boxplot of macro.mse by method
deconvolution_all_methods %>% 
  select(method, mixture_ID, deconvolution) %>% 
  mutate(method = map_chr(method, ~ .x %>% str_replace_all("\\.new", "") )) %>% 
  mutate(MSE = map_dbl(
    deconvolution,
    ~ mean((.x$estimated_proportion - .x$proportion)^2)
  )) %>% 
  group_by(method) %>% 
  summarise(MSE, macro.mse=mean(MSE)) %>% 
  ungroup() %>% 
  ggplot(aes(x= reorder(method, -macro.mse), y=log(MSE), colour = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
theme(axis.text.x = element_blank())

