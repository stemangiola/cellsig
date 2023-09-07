cellsig: a Bayesian sparse multilevel modelling approach for modelling
variability of celltype specific transcriptional profiles across studies
================

### Explore database through the web user interface

[Shiny app](https://shiny.wehi.edu.au/khan.k/cellsigdb)

### Installation

**Github**

``` r
install.packages("devtools")
devtools::install_github("stemangiola/cellsig")
```

### Example

#### Here, we’ll demonstrate an example of how to utilize the Bayesian multilevel noise-modelling of a transcriptome dataset.

##### Load the required packages

``` r
library(tidyverse)
library(magrittr)
library(tidybulk)
library(tidySummarizedExperiment)
library(cellsig)
library(yaml)
library(data.tree)
library(rstan)
```

##### Load the example count dataset

``` r
dataset <- readRDS("dev/test_data/count_dataset.rds")

dataset
```

    ## # A tibble: 1,000 × 5
    ##    cell_type   sample      symbol  database    count
    ##    <chr>       <fct>       <fct>   <fct>       <int>
    ##  1 endothelial ENCFF117AWL CFLAR   ENCSR000AAB 36711
    ##  2 endothelial ENCFF117AWL MYLIP   ENCSR000AAB  2554
    ##  3 endothelial ENCFF117AWL BTK     ENCSR000AAB     2
    ##  4 endothelial ENCFF117AWL CD22    ENCSR000AAB    60
    ##  5 endothelial ENCFF117AWL SEMA3B  ENCSR000AAB    62
    ##  6 endothelial ENCFF117AWL GRAMD1B ENCSR000AAB   168
    ##  7 endothelial ENCFF117AWL LCP2    ENCSR000AAB    50
    ##  8 endothelial ENCFF117AWL PRSS8   ENCSR000AAB     0
    ##  9 endothelial ENCFF117AWL CS      ENCSR000AAB 10642
    ## 10 endothelial ENCFF117AWL CTSA    ENCSR000AAB  4014
    ## # ℹ 990 more rows

##### Load the exmaple tree file with cell type hierarchy

``` r
tree <- read_yaml("dev/test_data/tree.yaml") %>% as.Node() # tree file was generated using data.tree package

tree 
```

    ##                   levelName
    ## 1 Tissue                   
    ## 2  ¦--epithelial           
    ## 3  ¦--fibroblast           
    ## 4  ¦--endothelial          
    ## 5  °--immune_cell          
    ## 6      °--natural_killer   
    ## 7          ¦--nk_cd56bright
    ## 8          °--nk_cd56dim

##### Create a hierarchical signature data frame from a tree and and signature database

``` r
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
```

    ## No group or design set. Assuming all samples belong to one group.

    ## tidybulk says: the sample with largest library size ENCFF467LOW was chosen as reference for scaling

    ## Warning in add_scaled_counts_bulk.calcNormFactor(df, reference, .sample =
    ## !!.sample, : tidybulk says: There are < 100 features/genes that are present in
    ## all you samples. Because edgeR::calcNormFactors does not allow NAs, the scaling
    ## is performed on that limited set of features.genes. The scaling could not be
    ## accurate, it is adivasble to perform impute_missing_abundance() before scaling.
    ## It is possible to filter the imputed counts after scaling.

``` r
dataset_input
```

    ## # A tibble: 1,800 × 11
    ##    .feature .sample count sample database  symbol  .abundant   TMM multiplier
    ##    <chr>    <chr>   <int> <chr>  <chr>     <chr>   <lgl>     <dbl>      <dbl>
    ##  1 AADAT    1405_0      0 1405_0 GSE152571 AADAT   FALSE      1.14       4.22
    ##  2 AADAT    1405_0      0 1405_0 GSE152571 AADAT   FALSE      1.14       4.22
    ##  3 AADAT    1405_0      0 1405_0 GSE152571 AADAT   FALSE      1.14       4.22
    ##  4 ADAD2    1405_0      1 1405_0 GSE152571 ADAD2   FALSE      1.14       4.22
    ##  5 ADAD2    1405_0      1 1405_0 GSE152571 ADAD2   FALSE      1.14       4.22
    ##  6 ADAD2    1405_0      1 1405_0 GSE152571 ADAD2   FALSE      1.14       4.22
    ##  7 AGAP3    1405_0    951 1405_0 GSE152571 AGAP3   TRUE       1.14       4.22
    ##  8 AGAP3    1405_0    951 1405_0 GSE152571 AGAP3   TRUE       1.14       4.22
    ##  9 AGAP3    1405_0    951 1405_0 GSE152571 AGAP3   TRUE       1.14       4.22
    ## 10 ANAPC15  1405_0    436 1405_0 GSE152571 ANAPC15 TRUE       1.14       4.22
    ## # ℹ 1,790 more rows
    ## # ℹ 2 more variables: level <int>, cell_type <chr>

##### Now, we’ll perform the modelling on the prepared input dataset

``` r
modelled_dataset <- dataset_input |>
    
    cellsig_multilevel_varing_intercept(
      .sample, 
      .feature,
      count, 
      cell_type,
      multiplier, 
      database
    )
```

``` r
##################################################################################################
## You can also follow-up the progress from the viewer window
##################################################################################################

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
```

``` r
modelled_dataset
```

    ## # A tibble: 700 × 13
    ##    .feature_idx .feature cell_type        mean se_mean     sd `10%` `50%`  `90%`
    ##           <int> <chr>    <chr>           <dbl>   <dbl>  <dbl> <dbl> <dbl>  <dbl>
    ##  1            1 AADAT    endothelial    2.78e2 2.45e+1 4.42e2  37    156. 5.41e2
    ##  2            2 AADAT    epithelial     1.27e3 2.97e+2 5.10e3  65.8  305  1.45e3
    ##  3            3 AADAT    fibroblast     4.97e3 2.10e+3 3.65e4 300.   953  4.62e3
    ##  4            4 AADAT    immune_cell    1.84e1 4.64e+0 8.34e1   0      3  3.41e1
    ##  5            5 AADAT    natural_killer 2.86e1 1.02e+1 1.80e2   0      4  3.5 e1
    ##  6            6 AADAT    nk_cd56bright  4.68e0 1.82e+0 3.23e1   0      0  5.10e0
    ##  7            7 AADAT    nk_cd56dim     7.46e1 2.40e+1 3.29e2   1     14  9.92e1
    ##  8            8 ADAD2    endothelial    6.14e0 7.34e-1 1.06e1   0      2  1.7 e1
    ##  9            9 ADAD2    epithelial     5.27e0 2.61e+0 4.51e1   0      0  2.10e0
    ## 10           10 ADAD2    fibroblast     1.61e0 3.51e-1 6.23e0   0      0  3.10e0
    ## # ℹ 690 more rows
    ## # ℹ 4 more variables: n_eff <dbl>, Rhat <dbl>, log_mean <dbl>, log_sd <dbl>
