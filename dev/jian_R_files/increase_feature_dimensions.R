library(tidyverse)
library(tidybulk)
library(umap)

tt_non_hierarchy <- 
  readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/tt_non_hierarchy.rds")

pca <- tt_non_hierarchy %>% 
  unnest(tt) %>% 
  unnest(data) %>% 
  select(sample, symbol, count_scaled) %>% 
  pivot_wider(names_from = symbol, values_from = count_scaled) %>% 
  as_matrix(rownames = sample) %>% dim()
  prcomp()

#scree plot
plot(x = 1:30, y = pca$sdev[1:30], xlab = "number of PC", ylab = "standard devidation explained",
     main = "Scree plot")

# proportion of variance explained by individual PCs
barplot(100*(pca$sdev^2 / sum(pca$sdev^2))[1:30],
        ylab = "proportion of variance explained (%)",
        xlab = "number of PC",
        names.arg = 1:30,
        axisnames = TRUE)

# cumulative proportion of variance explained by PCs
cumsum(pca$sdev^2 / sum(pca$sdev^2))[1:30]

barplot(100*cumsum(pca$sdev^2 / sum(pca$sdev^2))[1:30], 
        ylab = "cumulative porportion of variance (%)",
        xlab = "number of PC",
        ylim = c(0, 100),
        names.arg = 1:30,
        axisnames = TRUE)
abline(h=80, col = "red", lty=2)


counts_imputed_non_hierarchy <- readRDS("dev/intermediate_data/counts_imputed_non_hierarchy.rds")

pca <- counts_imputed_non_hierarchy %>% 
  unnest(tt) %>% 
  unnest(data) %>% 
  select(sample, symbol, count_scaled) %>% 
  pivot_wider(names_from = symbol, values_from = count_scaled) %>% 
  as_matrix(rownames = sample) %>%
  prcomp()

#scree plot
png("dev/jian_R_files/pca_var_indiv.png", height = 800, width = 800)
# par(mfrow=c(1, 2))
# plot(x = 1:30, y = pca$sdev[1:30], xlab = "number of PC", ylab = "standard devidation explained",
#      main = "Scree plot")

# proportion of variance explained by individual PCs
barplot(100*(pca$sdev^2 / sum(pca$sdev^2))[1:30],
        ylab = "proportion of variance explained (%)",
        xlab = "number of PC",
        names.arg = 1:30,
        axisnames = TRUE,
        font.lab = 2,
        cex.lab = 2,
        cex.axis = 1.5)

dev.off()


# cumulative proportion of variance explained by PCs
# cumsum(pca$sdev^2 / sum(pca$sdev^2))[1:30]

png("dev/jian_R_files/pca_var_cumu.png", height = 800, width = 800)

barplot(100*cumsum(pca$sdev^2 / sum(pca$sdev^2))[1:30], 
        ylab = "cumulative porportion of variance (%)",
        xlab = "number of PC",
        ylim = c(0, 100),
        names.arg = 1:30,
        axisnames = TRUE,
        font.lab = 2,
        cex.lab = 2,
        cex.axis = 1.5)
abline(h=90, col = "red", lty=2, lwd=2)

# par(mfrow=c(1, 1))
dev.off()

# produce umap plot using 2 components(default) from signatures of best stream by deconvolution:
# hierarchical_mean_contrast_bayes___silhouette_curvature

# cibersortx signature
cibersortx <- tibble(stream = "cibersortx") %>% 
  mutate(signature = 
           read_delim(
             "dev/jian_R_files/cibersortx/CIBERSORTx_Job21_phenoclass_1.CIBERSORTx_Job21_reference_1.bm.K999.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
           pull(NAME) %>% list
  )

x <- tibble(stream = "hierarchical_mean_contrast_bayes___silhouette_curvature") %>% 
  mutate(signature = hierarchical_mean_contrast_bayes___silhouette_curvature %>% 
           pull(signature) %>% unlist %>% unique %>% list) %>% 
  bind_rows(cibersortx) %>% 
  mutate(data = counts_imputed_non_hierarchy %>% unnest(tt) %>% unnest(data) %>% list) %>% 
  
  mutate(data = map2(data, signature, ~.x %>% filter(symbol %in% .y))) %>% 
  mutate(data = map(data, ~ .x %>%
                      pivot_longer(contains("level"), names_to = "level", values_to = "cell_at_level") %>% 
                      mutate(across(contains("cell"), ~ as.character(.x))) %>% 
                      drop_na() %>% 
                      filter(cell_type == cell_at_level) %>% 
                      select(sample, symbol, count_scaled, cell_type, level) %>% 
                      nest(cell_level = -c(sample, symbol, count_scaled))
                    )) %>% 
  mutate(reduced_dimension = map(data, ~ .x %>% 
                                   distinct(sample, symbol, count_scaled) %>% 
                                   reduce_dimensions(sample, symbol, count_scaled, 
                                                     method = "UMAP", 
                                                     .dims=2, 
                                                     top=Inf, scale=FALSE, 
                                                     action = "get")
                                   )) %>% 

  mutate(reduced_dimension = map2(
    reduced_dimension, data,
    ~.x %>% 
      left_join(.y %>% distinct(sample, cell_level), by="sample")
  ))
  

library(patchwork)

plot1 <- x %>% 
  pluck("reduced_dimension", 1) %>%
  unnest(cell_level) %>% 
  ggplot(aes(UMAP1, UMAP2, color=cell_type)) +
  geom_point()+
  ggtitle("hierarchical.meanContrast.bayes.silhouette.curvature, by cell_type")

plot2 <- x %>% 
  pluck("reduced_dimension", 1) %>%
  unnest(cell_level) %>% 
  ggplot(aes(UMAP1, UMAP2, color=level)) +
  geom_point() +
  ggtitle("hierarchical.meanContrast.bayes.silhouette.curvature, by level")

plot3 <- x %>% 
  pluck("reduced_dimension", 2) %>%
  unnest(cell_level) %>% 
  ggplot(aes(UMAP1, UMAP2, color=cell_type)) +
  geom_point()+
  ggtitle("cibersortx, by cell_type")

plot4 <- x %>% 
  pluck("reduced_dimension", 2) %>%
  unnest(cell_level) %>% 
  ggplot(aes(UMAP1, UMAP2, color=level)) +
  geom_point() +
  ggtitle("cibersortx, by level")


(plot1 | plot2) / (plot3 | plot4) +
  plot_layout(guides = "collect") +
  plot_annotation(title="UMAP dims=2")

# produce UMAP plot with .dim=10

y <- tibble(stream = "hierarchical_mean_contrast_bayes___silhouette_curvature") %>% 
  mutate(signature = hierarchical_mean_contrast_bayes___silhouette_curvature %>% 
           pull(signature) %>% unlist %>% unique %>% list) %>% 
  bind_rows(cibersortx) %>% 
  mutate(data = counts_imputed_non_hierarchy %>% unnest(tt) %>% unnest(data) %>% list) %>% 
  
  mutate(data = map2(data, signature, ~.x %>% filter(symbol %in% .y))) %>% 
  mutate(data = map(data, ~ .x %>%
                      pivot_longer(contains("level"), names_to = "level", values_to = "cell_at_level") %>% 
                      mutate(across(contains("cell"), ~ as.character(.x))) %>% 
                      drop_na() %>% 
                      filter(cell_type == cell_at_level) %>% 
                      select(sample, symbol, count_scaled, cell_type, level) %>% 
                      nest(cell_level = -c(sample, symbol, count_scaled))
  )) %>% 
  
  mutate(reduced_dimension = map(data, ~ .x %>% 
                                   distinct(sample, symbol, count_scaled) %>% 
                                   reduce_dimensions(sample, symbol, count_scaled, 
                                                     method = "PCA", 
                                                     .dims=10, 
                                                     top=Inf, scale=FALSE, 
                                                     action = "get")
                                 )) %>% 
  
  mutate(reduced_dimension = map2(
    reduced_dimension, data,
    ~.x %>% 
      left_join(.y %>% distinct(sample, cell_level), by="sample")
  ))


library(patchwork)

plot1 <- y %>% 
  pluck("reduced_dimension", 1) %>%
  unnest(cell_level) %>% 
  ggplot(aes(UMAP1, UMAP2, color=cell_type)) +
  geom_point()+
  ggtitle("hierarchical.meanContrast.bayes.silhouette.curvature, by cell_type")

plot2 <- y %>% 
  pluck("reduced_dimension", 1) %>%
  unnest(cell_level) %>% 
  ggplot(aes(UMAP1, UMAP2, color=level)) +
  geom_point() +
  ggtitle("hierarchical.meanContrast.bayes.silhouette.curvature, by level")

plot3 <- y %>% 
  pluck("reduced_dimension", 2) %>%
  unnest(cell_level) %>% 
  ggplot(aes(UMAP1, UMAP2, color=cell_type)) +
  geom_point()+
  ggtitle("cibersortx, by cell_type")

plot4 <- y %>% 
  pluck("reduced_dimension", 2) %>%
  unnest(cell_level) %>% 
  ggplot(aes(UMAP1, UMAP2, color=level)) +
  geom_point() +
  ggtitle("cibersortx, by level")


(plot1 | plot2) / (plot3 | plot4) +
  plot_layout(guides = "collect") +
  plot_annotation(title="UMAP dims=10")
