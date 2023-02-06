README
================

# Explore database through the web user interface

[Shiny app](https://shiny.wehi.edu.au/khan.k/cellsigdb)

# Example

Create a hierarchical signature data frame from a tree and and signature
database

``` r
library(here)
```

    ## here() starts at /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/cellsig

``` r
library(tibble)
library(readr)
library(cellsig)
```

    ## Warning: replacing previous import 'Rcpp::LdFlags' by 'RcppParallel::LdFlags'
    ## when loading 'cellsig'

    ## Warning: replacing previous import 'magrittr::set_names' by 'purrr::set_names'
    ## when loading 'cellsig'

    ## Warning: replacing previous import 'magrittr::extract' by 'tidyr::extract' when
    ## loading 'cellsig'

``` r
test_random_intercept |>
    
    #mutate(multiplier = exp(exposure_rate)) |> 
    cellsig_multilevel_varing_intercept(
      .sample, 
      .feature,
      count, 
      cell_type,
      multiplier, 
      database, 
      pass_fit = TRUE
    )
```

    ## Warning: There were 4 divergent transitions after warmup. See
    ## https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
    ## to find out why this is a problem and how to eliminate them.

    ## Warning: Examine the pairs() plot to diagnose sampling problems

    ## Warning: The largest R-hat is 1.14, indicating chains have not mixed.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#r-hat

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#bulk-ess

    ## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#tail-ess

    ## # A tibble: 50 × 13
    ##    .feature…¹ .feat…² cell_…³    mean se_mean      sd  `10%`  `50%`  `90%` n_eff
    ##         <int> <chr>   <fct>     <dbl>   <dbl>   <dbl>  <dbl>  <dbl>  <dbl> <dbl>
    ##  1          1 AGBL5-… endoth… 2.4 e-1 5.02e-2 8.78e-1 0         0   1   e0  306.
    ##  2          2 ASGR1   endoth… 9.65e+2 6.53e+1 1.07e+3 1.65e2  692.  2.02e3  271.
    ##  3          3 BCL2L2… endoth… 4.15e+2 4.92e+1 7.39e+2 3.38e1  220   8.69e2  226.
    ##  4          4 BTG2    endoth… 1.51e+4 2.83e+3 5.05e+4 1.79e3 7692.  2.55e4  318.
    ##  5          5 C2orf66 endoth… 5.82e+1 7.48e+0 1.24e+2 1.90e0   22.5 1.28e2  277.
    ##  6          6 CACNB4  endoth… 8.53e+1 1.19e+1 1.63e+2 4   e0   35   2.13e2  190.
    ##  7          7 CCL21   endoth… 6.18e+1 1.28e+1 2.31e+2 0         4   1.67e2  323.
    ##  8          8 CHURC1… endoth… 1.23e+2 1.09e+1 1.93e+2 5.9 e0   59   2.93e2  310.
    ##  9          9 CLCN3P1 endoth… 9.29e+0 1.21e+0 1.94e+1 0         1   2.71e1  258.
    ## 10         10 COL14A1 endoth… 5.36e+1 6.60e+0 9.53e+1 1   e0   19.5 1.27e2  208.
    ## # … with 40 more rows, 3 more variables: Rhat <dbl>, log_mean <dbl>,
    ## #   log_sd <dbl>, and abbreviated variable names ¹​.feature_idx, ²​.feature,
    ## #   ³​cell_type
