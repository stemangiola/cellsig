

#' ref_intercept_only
#'
#' @description This function calls the stan model.
#'
#'
#' @importFrom tibble tibble
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr mutate_if
#' @import dplyr
#'
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom tidyr drop_na
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes median_qi
#' @import tidybayes
#'
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#'
#' @import stringr
#'
#' @import data.tree
#'
#' @param mix A matrix
#' @param my_design A matrix
#' @param cov_to_test A character string
#' @param fully_bayesian A boolean
#' @param is_mix_microarray A boolean
#' @param verbose A boolean
#' @param save_report A boolean
#' @param custom_ref A matrix
#' @param multithread A boolean
#' @param do_debug A boolean
#' @param cell_type_root A character string
#' @param choose_internal_ref A design matrix
#' @param omit_regression A boolean
#' @param save_fit A boolean
#' @param seed An integer
#'
#' @return An ARMET object
#'
#' @export
#'
ref_intercept_only = function(reference,
                              level,
                              cores = 8,
                              approximate_posterior = F
) {
  # Former parameters
  X = matrix(rep(1, 1))
  do_regression = F
  full_bayesian = T
  omit_regression =                   T
  save_fit =                          T
  seed =                              NULL
  iterations = 250
  sampling_iterations = 100
  levels = 1:4

  # Add fake lambda, sigma
  reference =
    reference %>%
    mutate(lambda = 1, sigma_raw = 1)

  shards = cores #* 2
  is_level_in = shards %>% `>` (0) %>% as.integer

  # Global properties - derived by previous analyses of the whole reference dataset
  sigma_intercept = 1.3420415
  sigma_slope = -0.3386389
  sigma_sigma = 1.1720851
  lambda_mu_mu = 5.612671
  lambda_sigma = 7.131593

  # Non centered
  lambda_mu_prior = c(6.2, 1)
  lambda_sigma_prior =  c(log(3.3) , 1)
  lambda_skew_prior =  c(-2.7, 1)
  sigma_intercept_prior = c(1.9 , 0.1)

  # Set up tree structure
  levels_in_the_tree = 1:4

  tree =
    tree %>%
    data.tree::Clone() %>%	{
      # Filter selected levels
      data.tree::Prune(., function(x)
        x$level <= max(levels_in_the_tree) + 1)
      .
    }

  ct_in_nodes =
    tree %>%
    data.tree::ToDataFrameTree("name", "level", "C", "count", "isLeaf") %>%
    as_tibble %>%
    arrange(level, C) %>%
    filter(!isLeaf) %>%
    pull(count)

  # Get the number of leafs for every level
  ct_in_levels = foreach(l = levels_in_the_tree + 1, .combine = c) %do% {
    data.tree::Clone(tree) %>%
      when((.) %>% data.tree::ToDataFrameTree("level") %>% pull(2) %>% max %>% `>` (l)  ~ {
                    (.)
                    data.tree::Prune(., function(x)
                      x$level <= l)
                    (.)
                  }, ~(.)
        )  %>%
      data.tree::Traverse(., filterFun = isLeaf) %>%
      length()
  }

  n_nodes = ct_in_nodes %>% length
  n_levels = ct_in_levels %>% length

  singles_lv2 = tree$Get("C1", filterFun = isLeaf) %>% na.omit %>% as.array
  SLV2 = length(singles_lv2)
  parents_lv2 = tree$Get("C1", filterFun = isNotLeaf) %>% na.omit %>% as.array
  PLV2 = length(parents_lv2)

  singles_lv3 = tree$Get("C2", filterFun = isLeaf) %>% na.omit %>% as.array
  SLV3 = length(singles_lv3)
  parents_lv3 = tree$Get("C2", filterFun = isNotLeaf) %>% na.omit %>% as.array
  PLV3 = length(parents_lv3)

  singles_lv4 = tree$Get("C3", filterFun = isLeaf) %>% na.omit %>% as.array
  SLV4 = length(singles_lv4)
  parents_lv4 = tree$Get("C3", filterFun = isNotLeaf) %>% na.omit %>% as.array
  PLV4 = length(parents_lv4)


  # Prepare data frames -
  # For G house keeing first

  reference_filtered =
    reference %>%
    select(
      level,
      sample,
      symbol,
      `cell_type`,
      `count`,
      lambda,
      sigma_raw,
      `house_keeping`
    ) %>%

    # Bug after I deleted FANTOM5 I have to rerun infer NB. Some genes are not in all cell types anynore
    # Other bug
    filter(`cell_type` %>% is.na %>% `!`) %>%

    # Check if this is still important
    group_by(level) %>%
    do((.) %>% inner_join(
      (.) %>%
        distinct(symbol, `cell_type`) %>%
        count(symbol) %>%
        filter(n == max(n))
    )) %>%
    ungroup() %>%


    # left_join(n_markers, by=c("ct1", "ct2")) %>%
    # filter_reference(mix) %>%
    # select(-ct1, -ct2, -rank, -`n markers`) %>%
    # distinct %>%

    # Select cell types in hierarchy
    inner_join(
      tree %>%
        data.tree::ToDataFrameTree("cell_type", "C", "C1", "C2", "C3", "C4") %>%
        as_tibble %>%
        select(-1)

    )

  res1 = run_model_ref(
    tree,
    reference_filtered,
    shards,
    level,
    T,
    approximate_posterior,
    iterations = iterations,
    sampling_iterations = sampling_iterations
  )

  res1[[1]] %>% filter(!query) %>% distinct(symbol, `cell_type`, G, S, sample) %>%

    # Attach lambda sigma
    left_join(
      res1[[2]] %>% rstan::summary(c("lambda_log", "sigma_inv_log")) %$% summary %>%
        as_tibble(rownames = "par") %>%
        separate(par, c("par", "G"), sep = "\\[|\\]", extra = "drop") %>%
        mutate(G = G %>% as.integer) %>%
        select(par, G, "50%") %>%
        spread(par, `50%`),
      by = c("G")
    ) %>%

    # Attach exposure
    left_join(
      res1[[2]] %>%	rstan::summary("exposure_rate") %$% summary %>%
        as_tibble(rownames="par") %>%
        separate(par, c("par", "S"), sep="\\[|\\]", extra = "drop") %>%
        mutate(S = S %>% as.integer) %>%
        select(par, S, "50%") %>%
        rename(exposure = `50%`) %>%
        select(-par)
    ) %>%

    # Replace the category house_keeping
    mutate(`house_keeping` = `cell_type` == "house_keeping") %>%
    rename(temp = `cell_type`) %>%
    left_join(
      (.) %>%
        filter(temp != "house_keeping") %>%
        distinct(sample, S, temp) %>%
        rename(`cell_type` = temp)
    ) %>%
    select(-temp) %>%

    # attach fit
    add_attr(res1[[2]], "fit")
}

#' Add attribute to abject
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}


#' @importFrom rstan vb
#' @importFrom rstan sampling
#'
#' @export
#'
run_model_ref = function(tree,
                         reference_filtered,
                         shards,
                         lv,
                         full_bayesian,
                         approximate_posterior,
                         exposure_posterior = tibble(.mean = 0, .sd = 0)[0,],
                         iterations = 250,
                         sampling_iterations = 100) {


  reference_filtered =
    reference_filtered %>%
    inner_join(
      # Filter on level considered
      tree %>%
        data.tree::Clone() %>%
        {
          # Filter selected levels
          data.tree::Prune(., function(x)
            x$level <= lv + 1)
          .
        } %>%
        {
          .$Get("name", filterFun  = isLeaf)
        } %>%
        as_tibble() %>%
        setNames("cell_type")
    )

  # Check if there are not house_keeping
  if (reference_filtered %>% filter(`house_keeping`) %>% nrow %>% equals(0))
    stop("No house_keeping genes in your reference data frame")

  df = ref_format(reference_filtered) %>% distinct()

  G = df %>% filter(!`query`) %>% distinct(G) %>% nrow()
  GM = df %>% filter(!`house_keeping`) %>% distinct(symbol) %>% nrow()

  # For  reference MPI inference
  counts_baseline =
    df %>%

    # Eliminate the query part, not the house_keeping of the query
    filter(!`query` | `house_keeping`)  %>%

    format_for_MPI(shards)

  S = counts_baseline %>% distinct(sample) %>% nrow()
  N = counts_baseline %>% distinct(idx_MPI, `count`, `count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
  M = counts_baseline %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max

  lambda_log = 	  counts_baseline %>% filter(!query) %>% distinct(G, lambda) %>% arrange(G) %>% pull(lambda)
  sigma_inv_log = counts_baseline %>% filter(!query) %>% distinct(G, sigma_raw) %>% arrange(G) %>% pull(sigma_raw)

  counts_baseline_to_linear =
    counts_baseline %>%
    filter_house_keeping_query_if_fixed(full_bayesian) %>%
    arrange(G, S) %>%
    mutate(counts_idx = 1:n()) %>%
    mutate(S = S %>% as.factor %>% as.integer)

  counts_linear = counts_baseline_to_linear %>%  pull(`count`)
  G_to_counts_linear = counts_baseline_to_linear %>% pull(G)
  G_linear = G_to_counts_linear
  S_linear = counts_baseline_to_linear %>% pull(S)

  CL = length(counts_linear)
  S = counts_baseline_to_linear %>% distinct(S) %>% nrow


  MPI_data = get_MPI_df_ref(counts_baseline_to_linear,
                            counts_baseline,
                            shards,
                            lv)

  # library(rstan)
  # fileConn<-file("~/.R/Makevars")
  # writeLines(c( "CXX14FLAGS += -O2","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
  # close(fileConn)
  # ARMET_tc_model = rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_ref.stan", auto_write = F)

  Sys.setenv("STAN_NUM_THREADS" = shards)

  list(df,
       switch(
         approximate_posterior %>% sum(1),

         # HMC
         sampling(
           stanmodels$ARMET_ref,
           #rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan"),
           chains = 3,
           cores = 3,
           iter = iterations,
           warmup = iterations - sampling_iterations,
           data = MPI_data,
           #pars=
           # c("prop_1", "prop_2", "prop_3", sprintf("prop_%s", letters[1:9])) %>%
           # c("alpha_1", sprintf("alpha_%s", letters[1:9])) %>%
           # c("exposure_rate") %>%
           # c("lambda_UFO") %>%
           # c("prop_UFO") %>%
           # c(additional_par_to_save),
           save_warmup = FALSE
         ) %>%
           {
             (.)  %>% rstan::summary() %$% summary %>% as_tibble(rownames = "par") %>% arrange(Rhat %>% desc) %>% print
             (.)
           },

         vb_iterative(
           stanmodels$ARMET_ref,
           #rstan::stan_model("~/PhD/deconvolution/ARMET/inst/stan/ARMET_tc_fix_hierarchical.stan"),
           output_samples = 500,
           iter = 50000,
           tol_rel_obj = 0.005,
           algorithm = "meanfield",
           data = MPI_data,
           # pars = c(
           # 	"prop_1",
           # 	"prop_2",
           # 	"prop_3",
           # 	"prop_4",
           # 	"exposure_rate",
           # 	"lambda_log",
           # 	"sigma_inv_log",
           # 	"sigma_intercept_dec"
           # ),
           # #,
           init = function ()
             list(lambda_log = lambda_log, sigma_inv_log = sigma_inv_log) # runif(G,  lambda_log - 1, lambda_log + 1)	)

         )
       ))

}

#' @importFrom nanny as_matrix
#'
#' @export
#'
get_MPI_df_ref = function(counts_baseline_to_linear,
                          counts_baseline,
                          shards_in_levels,
                          lv) {
  list(
    counts_idx_lv_MPI =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv)  %>%
      distinct(idx_MPI, counts_idx, `count MPI row`) %>%
      spread(idx_MPI,  counts_idx) %>%
      select(-`count MPI row`) %>%
      replace(is.na(.),-999 %>% as.integer) %>%
      as_matrix() %>% t %>% 		as.data.frame,

    size_counts_idx_lv_MPI =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv)   %>%
      distinct(idx_MPI, counts_idx, `count MPI row`) %>%
      count(idx_MPI) %>%
      pull(n) %>%
      when(length(.) == 0 ~ 0, ~ (.)) %>%
      as.array,

    # Count indexes
    counts_G_lv_MPI =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv)   %>%
      distinct(idx_MPI, G, `count MPI row`)  %>%
      spread(idx_MPI,  G) %>%
      select(-`count MPI row`) %>%
      replace(is.na(.),-999 %>% as.integer) %>%
      as_matrix() %>% t %>% 		as.data.frame,

    size_counts_G_lv_MPI =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv)   %>%
      distinct(idx_MPI, G, `count MPI row`)  %>%
      count(idx_MPI) %>%
      pull(n) %>%
      when(length(.) == 0 ~ 0, ~ (.)) %>%
      as.array,

    counts_G_lv_MPI_non_redundant =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv)   %>%
      distinct(idx_MPI, G)  %>%
      group_by(idx_MPI) %>% do((.) %>% rowid_to_column("count MPI row")) %>% ungroup() %>%
      spread(idx_MPI,  G) %>%
      select(-`count MPI row`) %>%
      replace(is.na(.),-999 %>% as.integer) %>%
      as_matrix() %>% t %>% 		as.data.frame,

    size_counts_G_lv_MPI_non_redundant =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv)   %>%
      distinct(idx_MPI, G)  %>%
      count(idx_MPI) %>%
      pull(n) %>%
      when(length(.) == 0 ~ 0, ~ (.)) %>%
      as.array,

    counts_G_lv_MPI_non_redundant_reps =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv)   %>%
      distinct(idx_MPI, G, `count MPI row`)  %>%
      left_join((.) %>% count(idx_MPI, G)) %>%
      distinct(idx_MPI, G, n) %>%
      group_by(idx_MPI) %>% do((.) %>% rowid_to_column("count MPI row")) %>% ungroup() %>%
      distinct(idx_MPI, n, `count MPI row`) %>%
      spread(idx_MPI,  n) %>%
      select(-`count MPI row`) %>%
      replace(is.na(.),-999 %>% as.integer) %>%
      as_matrix() %>% t %>% 		as.data.frame,

    # Count indexes
    counts_S_lv_MPI =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv) %>%
      distinct(idx_MPI, S, `count MPI row`)  %>%
      spread(idx_MPI,  S) %>%
      select(-`count MPI row`) %>%
      replace(is.na(.),-999 %>% as.integer) %>%
      as_matrix() %>% t %>% 		as.data.frame,

    size_counts_S_lv_MPI =
      counts_baseline_to_linear %>%
      parse_baseline(shards_in_levels, lv) %>%
      distinct(idx_MPI, S, `count MPI row`)   %>%
      count(idx_MPI) %>%
      pull(n) %>%
      when(length(.) == 0 ~ 0, ~ (.)) %>%
      as.array,


    G_linear_MPI =
      counts_baseline %>% filter(level == lv) %>%

      # I have fixed this for right order
      select(level, G, GM, sprintf("C%s", lv)) %>%
      distinct() %>%
      arrange(GM,!!as.symbol(sprintf("C%s", lv))) %>%

      #distinct(G, GM, C, level) %>%
      left_join(tibble(level = lv, shards = shards_in_levels)) %>%
      format_for_MPI_from_linear_dec(lv) %>%
      distinct(idx_MPI, G, `count MPI row`) %>%
      spread(idx_MPI,  G) %>%
      select(-`count MPI row`) %>%
      replace(is.na(.),-999 %>% as.integer) %>%
      as_matrix() %>% t %>% 		as.data.frame,

    size_G_linear_MPI =
      counts_baseline %>% filter(level == lv) %>%
      # I have fixed this for right order
      select(level, G, GM, sprintf("C%s", lv)) %>%
      distinct() %>%
      arrange(GM,!!as.symbol(sprintf("C%s", lv))) %>%

      #distinct(G, GM, C, level) %>%
      left_join(tibble(level = lv, shards = shards_in_levels)) %>%
      format_for_MPI_from_linear_dec(lv) %>%
      distinct(idx_MPI, G, `count MPI row`)   %>%
      count(idx_MPI) %>%
      pull(n) %>%
      when(length(.) == 0 ~ 0, ~ (.)) %>%
      as.array
  )
}

#' @export
#'
ref_format = function(ref) {
  # Get reference based on mix genes
  ref %>% mutate(`query` = FALSE)	%>%

    # Add marker symbol indeces
    left_join((.) %>%
                filter(!`house_keeping`) %>%
                distinct(`symbol`) %>%
                mutate(M = 1:n())) %>%

    # Add sample indeces
    arrange(!`query`) %>% # query first
    mutate(S = factor(sample, levels = .$sample %>% unique) %>% as.integer) %>%

    # Add house_keeping into Cell type label
    mutate(`cell_type` = ifelse(`house_keeping`, "house_keeping", `cell_type`)) %>%

    # # Still needed? NOT because I have sample, ct unique, no redundancy
    # anti_join(
    # 	(.) %>%
    # 		filter(`house_keeping` & !`query`) %>%
    # 		distinct(symbol, level) %>%
    # 		group_by(symbol) %>%
    # 		arrange(level) %>%
    # 		slice(2:max(n(), 2)) %>% # take away house_keeping from level 2 above
    # 		ungroup()
    # ) %>%

  # If house_keeping delete level infomation
  mutate(level = ifelse(`house_keeping`, NA, level)) %>%

    # Create unique symbol ID
    unite(ct_symbol, c("cell_type", "symbol"), remove = F) %>%

    # Add gene idx
    left_join(
      (.) %>%
        filter(!`query`) %>%
        distinct(`cell_type`, ct_symbol, `house_keeping`) %>%
        arrange(!`house_keeping`, ct_symbol) %>% # house_keeping first
        mutate(G = 1:n())
    ) %>%
    left_join(
      (.) %>%
        filter(!`house_keeping` & !`query`) %>%
        distinct(level, symbol) %>%
        arrange(level, symbol) %>%
        mutate(GM = 1:n()) %>%
        select(-level)
    )

}


