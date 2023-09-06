README
================

# Explore database through the web user interface

[Shiny app](https://shiny.wehi.edu.au/khan.k/cellsigdb)

# Example

Create a hierarchical signature data frame from a tree and and signature
database

``` r
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
    
    cellsig_multilevel_varing_intercept(
      .sample, 
      .feature,
      count, 
      cell_type,
      multiplier, 
      database
    )
```

    ## Warning: There were 3 divergent transitions after warmup. See
    ## https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
    ## to find out why this is a problem and how to eliminate them.

    ## Warning: Examine the pairs() plot to diagnose sampling problems

    ## Warning: The largest R-hat is 1.08, indicating chains have not mixed.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#r-hat

    ## Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#bulk-ess

    ## Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
    ## Running the chains for more iterations may help. See
    ## https://mc-stan.org/misc/warnings.html#tail-ess

    ## # A tibble: 50 × 13
    ##    .feature_…¹ .feat…² cell_…³    mean se_mean     sd  `10%`  `50%`  `90%` n_eff
    ##          <int> <chr>   <fct>     <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>
    ##  1           1 AGBL5-… endoth… 3.57e-1 1.35e-1 2.40e0    0      0   1   e0  316.
    ##  2           2 ASGR1   endoth… 8.55e+2 4.55e+1 7.88e2  152.   648   1.73e3  300.
    ##  3           3 BCL2L2… endoth… 4.21e+2 5.22e+1 6.24e2   23.9  212.  9.26e2  143.
    ##  4           4 BTG2    endoth… 1.64e+4 5.32e+3 9.36e4 1871.  7344.  2.21e4  309.
    ##  5           5 C2orf66 endoth… 7.28e+1 1.63e+1 2.78e2    2     25.5 1.29e2  290.
    ##  6           6 CACNB4  endoth… 1.36e+2 3.19e+1 5.33e2    5     47.5 2.14e2  279.
    ##  7           7 CCL21   endoth… 5.88e+1 1.91e+1 3.39e2    0      5   1.01e2  315.
    ##  8           8 CHURC1… endoth… 1.47e+2 3.38e+1 5.94e2    4     49   2.89e2  309.
    ##  9           9 CLCN3P1 endoth… 1.21e+1 3.05e+0 5.13e1    0      2   1.81e1  283.
    ## 10          10 COL14A1 endoth… 9.87e+1 3.15e+1 5.37e2    1     20.5 1.20e2  290.
    ## # … with 40 more rows, 3 more variables: Rhat <dbl>, log_mean <dbl>,
    ## #   log_sd <dbl>, and abbreviated variable names ¹​.feature_idx, ²​.feature,
    ## #   ³​cell_type
