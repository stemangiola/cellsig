library(tidyverse)
library(tidybulk)

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

  
