cellsig: a Bayesian sparse multilevel modelling approach for modelling
variability of celltype specific transcriptional profiles across studies
================

## Explore database through the web user interface

[Shiny app](https://shiny.wehi.edu.au/khan.k/cellsigdb)

## Installation

**Github**

``` r
install.packages("devtools")
devtools::install_github("stemangiola/cellsig")
```

## Example

### Here, we’ll demonstrate an example of how to utilize the Bayesian multilevel noise-modelling of a transcriptome dataset.

### Load the required packages

``` r
library(tidyverse)
library(tidybulk)
library(tidySummarizedExperiment)
library(cellsig)
library(yaml)
library(data.tree)
library(rstan)
```

### Load the example count dataset

``` r
dataset <- readRDS("dev/test_data/count_dataset.rds")
```

The count dataset is in tibble format (as the data is not rectangular),
where the `cell_type` column refers to the cell types included in the
dataset, `sample` and `database` columns list the RNA-seq studies and
samples within a study for a particular cell type. `symbol` column lists
the gene transcripts and their abundances are indexed in the `count`
column.

``` r
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

### Load the exmaple tree file with cell type hierarchy

To utilise the cellular hierarchy information in the modelling, we
generated a `tree` file using `data.tree` package.

``` r
tree <- read_yaml("dev/test_data/tree.yaml") %>% as.Node()
```

In our example `tree` file, we have classified our cell type of
interests into 4 broad first level categories, namely- `epithelial`,
`endothelial`, `fibroblast` and `immune cell`. Branching from the
`immune cell` ancestor node, we put `natural killer` cells on the second
level of this hierarchy, which further branched into the two subsets of
NK cells `nk_cd56bright` and `nk_cd56dim` on the third and final level.

``` r
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

### Create a hierarchical signature data frame from a tree and and signature database

Now, we have to prepare the input dataset with the cellular hierarchy
embedded within it, so that the modelling step can utilise the cellular
hierarchy information. The dataset has to be rectangular (equal amount
of cell type-transcript pairs for every included cell types) for the
integration of cellular hierarchy levels, therefore, we impute the
missing transcript abundances for each cell type with
`impute_missing_abundance` function of `tidybulk` package. Next, we used
`tree_and_signatures_to_database` function of `cellsig` to assign the
respective cellular hierarchy levels for each cell type. Once we have
merged the cellular hierarchy, we dropped the imputed transcripts from
the dataset.

``` r
dataset_input <- dataset %>%
  
  # Imputation
  as_SummarizedExperiment(sample, symbol, count) %>% 
  impute_missing_abundance(~ cell_type, force_scaling = TRUE) %>% 
  as_tibble() %>% 
  filter(!is.na(count)) %>% 
  
  # Parsing
  tree_and_signatures_to_database(
    tree,
    .,
    .sample,
    cell_type,
    .feature,
    count
  ) %>%
  identify_abundant(.sample, .feature, count) %>% ## we're identifying the most abundant transcripts for each sample
  scale_abundance(.sample, .feature, count) %>% ## scaling and normalisation by "trimmed mean of M values (TMM)" approach
  dplyr::select(-count_scaled) %>% ## We 
  filter(!.imputed) %>%  ## we're dropping the transcripts which were imputed 
  select(-c(.imputed, .sample, .feature)) %>%  ## removing redundant columns

  # Adjusting the columns to prepare the input dataset for modelling 
  select(-cell_type) %>%
  pivot_longer(
    contains("level_"), names_prefix="level_", 
    names_to = "level", values_to="cell_type",
    names_transform=list(level=as.integer)
  ) %>%
  filter(!is.na(cell_type)) %>%
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

    ## # A tibble: 1,800 × 9
    ##    count sample database  symbol  .abundant   TMM multiplier level cell_type    
    ##    <int> <chr>  <chr>     <chr>   <lgl>     <dbl>      <dbl> <int> <chr>        
    ##  1     0 1405_0 GSE152571 AADAT   FALSE      1.14       4.22     1 immune_cell  
    ##  2     0 1405_0 GSE152571 AADAT   FALSE      1.14       4.22     2 natural_kill…
    ##  3     0 1405_0 GSE152571 AADAT   FALSE      1.14       4.22     3 nk_cd56bright
    ##  4     1 1405_0 GSE152571 ADAD2   FALSE      1.14       4.22     1 immune_cell  
    ##  5     1 1405_0 GSE152571 ADAD2   FALSE      1.14       4.22     2 natural_kill…
    ##  6     1 1405_0 GSE152571 ADAD2   FALSE      1.14       4.22     3 nk_cd56bright
    ##  7   951 1405_0 GSE152571 AGAP3   TRUE       1.14       4.22     1 immune_cell  
    ##  8   951 1405_0 GSE152571 AGAP3   TRUE       1.14       4.22     2 natural_kill…
    ##  9   951 1405_0 GSE152571 AGAP3   TRUE       1.14       4.22     3 nk_cd56bright
    ## 10   436 1405_0 GSE152571 ANAPC15 TRUE       1.14       4.22     1 immune_cell  
    ## # ℹ 1,790 more rows

In our prepared input dataset, we have generated some extra information,
such as `TMM` and `multiplier` columns refer to the estimated library
size factors and trimmed mean of M values, respectively for each
sample-transcript pairs. `.abundant` column indicates whether one
transcript is one of the most abundant transcripts within a sample.
`level` column highlights the cellular hierarchy levels of the cell
types listed in `cell_type`.

### Now, we’ll perform the modelling on the prepared input dataset

``` r
modelled_dataset <- dataset_input %>%
    
    cellsig_multilevel_varing_intercept(
      sample, 
      symbol,
      count, 
      cell_type,
      multiplier, 
      database
    )
```

We can also follow-up the progress from the `Viewer` window on Rstudio.

``` r
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
    ##    .feature_idx symbol cell_type    mean se_mean     sd `10%` `50%`  `90%` n_eff
    ##           <int> <chr>  <chr>       <dbl>   <dbl>  <dbl> <dbl> <dbl>  <dbl> <dbl>
    ##  1            1 AADAT  endotheli… 2.74e2  50.3   6.32e2  49.8  148.  455.   158.
    ##  2            2 AADAT  epithelial 1.60e3 690.    1.21e4  40.9  268. 1709.   307.
    ##  3            3 AADAT  fibroblast 1.66e3 155.    2.49e3 258.   992. 2952    257.
    ##  4            4 AADAT  immune_ce… 2.10e1   5.16  8.85e1   0      3    38    294.
    ##  5            5 AADAT  natural_k… 4.78e1  31.9   5.54e2   0      3    28.3  302.
    ##  6            6 AADAT  nk_cd56br… 2.29e0   0.537 9.75e0   0      0     5    329.
    ##  7            7 AADAT  nk_cd56dim 5.48e1  13.8   1.73e2   1     13   117    158.
    ##  8            8 ADAD2  endotheli… 7.94e0   1.31  2.34e1   0      1    19.2  318.
    ##  9            9 ADAD2  epithelial 1.43e0   0.465 6.67e0   0      0     2    206.
    ## 10           10 ADAD2  fibroblast 2.85e0   0.945 1.63e1   0      0     2    296.
    ## # ℹ 690 more rows
    ## # ℹ 3 more variables: Rhat <dbl>, log_mean <dbl>, log_sd <dbl>

The resulting dataframe contains the modelled variability information
such as - `mean` indicating the mean abundance of a cell type-transcript
pair, `se_mean` and `sd` refer to the standard error and standard
deviation of the mean values, respectively. `10%`, `50%`, and `90%`
values highlight the quantile distributions of abundances for each cell
type-transcript pair. `n_eff` values points to the effective sample
size, and `Rhat` refers to the Rhat statistic which is a measure of
chain equilibrium (should be within 0.05 of 1.0) estimated by the `stan`
program. `log_mean` and `log_sd`values represents the log-transformed
`mean` and `sd` values, respectively.

#### Finally, the modelled abundance variability information can further be used in obtaining downstream transcriptomic applications such as cell type specific marker identification approaches.
