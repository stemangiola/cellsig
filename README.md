cellsig: a Bayesian sparse multilevel modelling approach
================

### Explore database through the web user interface

[Shiny app](https://shiny.wehi.edu.au/khan.k/cellsigdb)

### Installation

**Github**

``` {r}
install.packages(devtools)
devtools::install_github("stemangiola/cellsig")

```

### Example

#### Here, we'll demonstrate an example of how to utilize the Bayesian multilevel noise-modelling of a transcriptome dataset.


```{r}
library(tidyverse)
library(magrittr)
library(tidybulk)
library(tidySummarizedExperiment)
library(glue)
library(cellsig)
library(yaml)
library(data.tree)
library(rstan)


## Load the example dataset and the tree file with cell type hierarchy

dataset <- readRDS("dev/test_data/count_dataset.rds")

dataset

##  A tibble: 1,934,920 × 5
##    cell_type  symbol database    sample      count
##    <chr>      <fct>  <fct>       <fct>       <int>
##  1 epithelial TSPAN6 ENCSR822SUG ENCFF556AGI  1330
##  2 epithelial TSPAN6 ENCSR822SUG ENCFF106YFW  1378
##  3 epithelial TSPAN6 ENCSR000AAD ENCFF377UBC  1388
##  4 epithelial TSPAN6 ENCSR000AAD ENCFF193YHT  1216
##  5 epithelial TSPAN6 ENCSR118TVR ENCFF380AIN  2750
##  6 epithelial TSPAN6 ENCSR118TVR ENCFF129BNX  4444
##  7 epithelial TSPAN6 ENCSR373BDG ENCFF108PJF  2598
##  8 epithelial TSPAN6 ENCSR000COX ENCFF300LEY  2287
##  9 epithelial TSPAN6 ENCSR000CUN ENCFF374KZN  1571
## 10 epithelial TSPAN6 ENCSR000AAL ENCFF826JBV  2883
##  ℹ 1,934,910 more rows
# ℹ Use `print(n = ...)` to see more rows

tree <- read_yaml("dev/test_data/tree.yaml") %>% as.Node()

tree 

##                   levelName
## 1 Tissue                   
## 2  ¦--epithelial           
## 3  ¦--fibroblast           
## 4  ¦--endothelial          
## 5  °--immune_cell          
## 6      °--natural_killer   
## 7          ¦--nk_cd56bright
## 8          °--nk_cd56dim   

```


#### Create a hierarchical signature data frame from a tree and and signature database

```{r}
# Firstly, we have prepare an input dataset for the modelling

dataset_input <- dataset %>%
  
  # Imputation
  as_SummarizedExperiment(sample, symbol, count) %>% 
  impute_missing_abundance(~ cell_type, force_scaling = TRUE) %>% # imputing the missing cell type-transcript abundances
  as_tibble() %>% 
  filter(count %>% is.na %>% `!`) %>% 
  
  # Parsing
  tree_and_signatures_to_database(  ## Assigning the cellular hierarchy levels in the dataset 
    tree, ## and removing any unwanted cell type (which is not present in the tree file) from the dataset
    .,
    .sample,
    cell_type,
    .feature,
    count
  ) %>%
  identify_abundant(.sample, .feature, count) %>%
  scale_abundance(.sample, .feature, count) %>%
  dplyr::select(-count_scaled) %>%
  filter(!.imputed) %>% 
  select(-.imputed) %>% 

  # CREATE INPUTS
  select(-cell_type) %>%
  pivot_longer(
    contains("level_"), names_prefix="level_", 
    names_to = "level", values_to="cell_type",
    names_transform=list(level=as.integer)
  ) %>%
  filter(cell_type %>% is.na %>% `!`) %>%
  mutate(count = as.integer(count))

## No group or design set. Assuming all samples belong to one group.
## tidybulk says: the sample with largest library size ENCFF467LOW was chosen as reference for scaling
## Warning: tidybulk says: There are < 100 features/genes that are present in all you samples. Because edgeR::calcNormFactors does not allow NAs, the scaling is performed on that limited set of features.genes. The scaling could not be accurate, it is adivasble to perform impute_missing_abundance() before scaling. It is possible to filter the imputed counts after scaling.

<<<<<<< HEAD
```


#### Now, we'll perform the modelling on the prepared input dataset

``` {r}
=======


>>>>>>> 19cc4f1861def579fb7781b60e5c193c8b25fcd2
modelled_dataset <- dataset_input |>
    
    cellsig_multilevel_varing_intercept(
      .sample, 
      .feature,
      count, 
      cell_type,
      multiplier, 
      database
    )

## Warning: There were 6 divergent transitions after warmup. See
## https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
## to find out why this is a problem and how to eliminate them.Warning: Examine the pairs() plot to diagnose sampling problems

## Warning: The largest R-hat is 1.17, indicating chains have not mixed.
## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#r-hatWarning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.

## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#bulk-essWarning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles #may be unreliable.

## Running the chains for more iterations may help. See
## https://mc-stan.org/misc/warnings.html#tail-ess

###############################################################################################################
## You can also follow-up the progress from the viewer window
###############################################################################################################

## Click the Refresh button to see progress of the chains
## starting worker pid=8968 on localhost:11465 at 13:13:21.944
## starting worker pid=5644 on localhost:11465 at 13:13:22.110
## starting worker pid=16376 on localhost:11465 at 13:13:22.253

## SAMPLING FOR MODEL 'mixed_effect' NOW (CHAIN 1).
## Chain 1: 
## Chain 1: Gradient evaluation took 0.002212 seconds
## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 22.12 seconds.
## Chain 1: Adjust your expectations accordingly!
## Chain 1: 
## Chain 1: 
## Chain 1: Iteration:   1 / 350 [  0%]  (Warmup)

## SAMPLING FOR MODEL 'mixed_effect' NOW (CHAIN 2).
## Chain 2: 
## Chain 2: Gradient evaluation took 0.002853 seconds
## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 28.53 seconds.
## Chain 2: Adjust your expectations accordingly!
## Chain 2: 
## Chain 2: 
## Chain 2: Iteration:   1 / 350 [  0%]  (Warmup)

## SAMPLING FOR MODEL 'mixed_effect' NOW (CHAIN 3).
## Chain 3: 
## Chain 3: Gradient evaluation took 0.004355 seconds
## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 43.55 seconds.
## Chain 3: Adjust your expectations accordingly!
## Chain 3: 
## Chain 3: 
## Chain 3: Iteration:   1 / 350 [  0%]  (Warmup)
## Chain 3: Iteration:  35 / 350 [ 10%]  (Warmup)
## Chain 2: Iteration:  35 / 350 [ 10%]  (Warmup)
## Chain 1: Iteration:  35 / 350 [ 10%]  (Warmup)
## Chain 3: Iteration:  70 / 350 [ 20%]  (Warmup)
## Chain 2: Iteration:  70 / 350 [ 20%]  (Warmup)
## Chain 1: Iteration:  70 / 350 [ 20%]  (Warmup)
## Chain 3: Iteration: 105 / 350 [ 30%]  (Warmup)
## Chain 1: Iteration: 105 / 350 [ 30%]  (Warmup)
## Chain 2: Iteration: 105 / 350 [ 30%]  (Warmup)
## Chain 3: Iteration: 140 / 350 [ 40%]  (Warmup)
## Chain 1: Iteration: 140 / 350 [ 40%]  (Warmup)
## Chain 2: Iteration: 140 / 350 [ 40%]  (Warmup)
## Chain 3: Iteration: 175 / 350 [ 50%]  (Warmup)
## Chain 2: Iteration: 175 / 350 [ 50%]  (Warmup)
## Chain 1: Iteration: 175 / 350 [ 50%]  (Warmup)
## Chain 3: Iteration: 210 / 350 [ 60%]  (Warmup)
## Chain 2: Iteration: 210 / 350 [ 60%]  (Warmup)
## Chain 1: Iteration: 210 / 350 [ 60%]  (Warmup)
## Chain 3: Iteration: 245 / 350 [ 70%]  (Warmup)
## Chain 3: Iteration: 251 / 350 [ 71%]  (Sampling)
## Chain 2: Iteration: 245 / 350 [ 70%]  (Warmup)
## Chain 2: Iteration: 251 / 350 [ 71%]  (Sampling)
## Chain 1: Iteration: 245 / 350 [ 70%]  (Warmup)
## Chain 1: Iteration: 251 / 350 [ 71%]  (Sampling)
## Chain 3: Iteration: 285 / 350 [ 81%]  (Sampling)
## Chain 2: Iteration: 285 / 350 [ 81%]  (Sampling)
## Chain 1: Iteration: 285 / 350 [ 81%]  (Sampling)
## Chain 3: Iteration: 320 / 350 [ 91%]  (Sampling)
## Chain 2: Iteration: 320 / 350 [ 91%]  (Sampling)

modelled_dataset

## # A tibble: 700 × 13
##    .feature_idx .feature cell_type          mean   se_mean        sd `10%`  `50%`   `90%` n_eff  Rhat log_mean log_sd
##           <int> <chr>    <chr>             <dbl>     <dbl>     <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl>    <dbl>  <dbl>
##  1            1 AADAT    endothelial      434.     196.      3281.    53.7  174    506.    279. 1.00     5.11   1.01 
##  2            2 AADAT    epithelial       836.     117.      2079.    62.9  289   1650.    317. 0.997    5.68   1.43 
##  3            3 AADAT    fibroblast     24696.   21908.    381473.   298.  1219   5160.    303. 1.00     7.14   1.40 
##  4            4 AADAT    immune_cell       17.0      2.83      51.3    0      3.5   37     329. 0.998    1.64   1.43 
##  5            5 AADAT    natural_killer    13.4      2.62      43.9    0      4     25.1   280. 0.997    1.59   1.31 
##  6            6 AADAT    nk_cd56bright      1.46     0.225      4.08   0      0      4.10  328. 0.997    0.387  0.794
##  7            7 AADAT    nk_cd56dim        88.9     33.0      585.     1     12     90.3   314. 1.00     2.60   1.60 
##  8            8 ADAD2    endothelial       10.5      2.57      45.6    0      2     18.1   315. 1.00     1.21   1.25 
##  9            9 ADAD2    epithelial        13.0      8.86     154.     0      0      3     303. 1.00     0.408  1.05 
## 10           10 ADAD2    fibroblast         2.19     0.732     12.9    0      0      2     311. 0.998    0.312  0.796
# ℹ 690 more rows

```

