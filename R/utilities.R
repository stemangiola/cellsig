#
##' create_tree_object
##'
##' @description create tree object that is in data directory
##'
##' @export
#create_tree_object = function(my_ref = cellsig::counts) {
# 
# 
# 
#   #yaml:: yaml.load_file("~/PhD/deconvolution/ARMET/data/tree.yaml") %>%
#   tree =
#     yaml::yaml.load_file("dev/tree.yaml") %>%
#     data.tree::as.Node() %>%
# 
#     {
# 
#       my_ct = my_ref %>% distinct(`cell_type`) %>% pull(1) %>% as.character
# 
#       # Filter if not in reference
#       data.tree::Prune(., pruneFun = function(x) ( x$name %in% my_ct ))
# 
#       # Sort tree by name
#       data.tree::Sort(., "name")
# 
#       # Add C indexes
#       .$Set(C =
#               tibble(
#                 name = .$Get('name'),
#                 level = .$Get('level')
#               ) %>%
#               left_join((.) %>% arrange(level, name) %>%	 	mutate(C = 0:(n(
# 
#               ) - 1)))	%>%
#               pull(C))
#       .$Set(C1 = get_idx_level(., 1))
#       .$Set(C2 = get_idx_level(., 2))
#       .$Set(C3 = get_idx_level(., 3))
#       .$Set(C4 = get_idx_level(., 4))
#       #		if(max(levels)>1) for(l in 2:max(levels)) { my_c = sprintf("C%s", l); .$Set(	my_c = get_idx_level(.,2)	); . }
# 
#       # Set cell_type label
#       .$Set("cell_type" = .$Get("name"))
# 
#       .
# 
#     }
# 
#   save(tree, file="data/tree.rda", compress = "xz")
# 
#   # ancestor_child = tree %>% get_ancestor_child
#   #
#   # save(ancestor_child, file="data/ancestor_child.rda", compress = "gzip")
# }

#' get_idx_level
#'
get_idx_level = function(tree, my_level) {
  left_join(
    tree %>% data.tree::ToDataFrameTree("name") %>% as_tibble,
    data.tree::Clone(tree) %>%
      {
        data.tree::Prune(., function(x)
          x$level <= my_level + 1)
        .
      } %>%
      data.tree::ToDataFrameTree("level", "C", "isLeaf", "name") %>%
      as_tibble %>%
      filter(isLeaf) %>%
      left_join(
        tree %>%
          data.tree::ToDataFrameTree("name", "level", "isLeaf") %>%
          as_tibble %>%
          filter(level <= my_level & isLeaf) %>%
          mutate(isAncestorLeaf = !isLeaf) %>%
          select(name, isAncestorLeaf)

      ) %>%
      arrange(isAncestorLeaf) %>%
      mutate(my_C = 1:n()) %>%
      select(name, my_C)
  ) %>%
    pull(my_C)
}

#' format_for_MPI
#'
#' @importFrom tibble rowid_to_column
#'
#' @description Format reference data frame for MPI
format_for_MPI = function(df, shards) {
  df %>%

    left_join((.) %>%
                distinct(G) %>%
                arrange(G) %>%
                mutate(idx_MPI = head(
                  rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
                )), by = "G") %>%
    arrange(idx_MPI, G) %>%

    # Decide start - end location
    group_by(idx_MPI) %>%
    do((.) %>%
         left_join(
           (.) %>%
             distinct(sample, G) %>%
             arrange(G) %>%
             count(G) %>%
             mutate(end = cumsum(n)) %>%
             mutate(start = c(
               1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1)
             ))
           , by = "G")) %>%
    ungroup() %>%

    # Add ct_symbol MPI rows indexes - otherwise spread below gives error
    left_join(
      (.) %>%
        group_by(idx_MPI) %>%
        distinct(G) %>%
        arrange(G) %>%
        mutate(`symbol MPI row` = 1:n()) %>%
        ungroup,
      by = c("G", "idx_MPI")
    ) %>%

    # Add counts MPI rows indexes
    group_by(idx_MPI) %>%
    arrange(G) %>%
    mutate(`count MPI row` = 1:n()) %>%
    # do( (.) %>% arrange(G) %>% rowid_to_column("count MPI row") ) %>%
    ungroup

}

filter_house_keeping_query_if_fixed =  function(.data, full_bayesian) {
  .data %>%

    # If full Bayesian false just keep house_keeping
    when(!full_bayesian ~ (.) %>% filter(`house_keeping` & `query`), ~ (.))
}

parse_baseline = function(.data, shards_in_levels, lv) {
  .data %>%
    filter(level == lv) %>%
    distinct(
      sample,
      symbol,
      `cell_type`,
      level,
      count,
      counts_idx,
      G,
      GM,
      S,
      `house_keeping`
    ) %>%
    left_join(tibble(level = lv, shards = shards_in_levels), by = "level") %>%
    format_for_MPI_from_linear()
}

#' @import magrittr
#' @importFrom tibble rowid_to_column
format_for_MPI_from_linear = function(df) {
  shards = df %>% arrange(shards %>% desc) %>% slice(1) %>% pull(shards)
  if (shards %>% length %>% equals(0))
    shards = 1


  df %>%

    left_join((.) %>%
                distinct(GM) %>%
                arrange(GM) %>%
                mutate(idx_MPI = head(
                  rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
                )), by = "GM") %>%
    arrange(idx_MPI, GM, G) %>%

    # Add counts MPI rows indexes
    group_by(idx_MPI) %>%
    arrange(GM, G) %>%
    do((.) %>% rowid_to_column("count MPI row")) %>%
    ungroup

}

get_ancestor_child = function(tree){
  tree %>% ToDataFrameTypeColFull %>% distinct(level_1, level_2) %>% setNames(c("ancestor", "cell_type")) %>%
    bind_rows(
    tree %>% ToDataFrameTypeColFull %>% distinct(level_2, level_3) %>% setNames(c("ancestor", "cell_type"))
  ) %>%
    bind_rows(
      tree %>% ToDataFrameTypeColFull %>% distinct(level_3, level_4) %>% setNames(c("ancestor", "cell_type"))
    ) %>%
    bind_rows(
      tree %>% ToDataFrameTypeColFull %>% distinct(level_4, level_5) %>% setNames(c("ancestor", "cell_type"))
    ) %>%
    filter(ancestor != `cell_type`)
}

#' @importFrom tibble rowid_to_column
format_for_MPI_from_linear_dec = function(df, lv) {
  shards = df %>% arrange(shards %>% desc) %>% slice(1) %>% pull(shards)
  if (shards %>% length %>% equals(0))
    shards = 1

  df %>%

    left_join((.) %>%
                distinct(GM) %>%
                arrange(GM) %>%
                mutate(idx_MPI = head(
                  rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
                ))) %>%
    arrange(idx_MPI, GM, !!as.symbol(sprintf("C%s", lv))) %>%

    # Add counts MPI rows indexes
    group_by(idx_MPI) %>%
    arrange(GM, !!as.symbol(sprintf("C%s", lv))) %>%
    do((.) %>% rowid_to_column("count MPI row")) %>%
    ungroup

}

#' vb_iterative
#'
#' @description Runs iteratively variational bayes until it suceeds
#'
#' @importFrom rstan vb
#'
#' @param model A Stan model
#' @param output_samples An integer of how many samples from posteriors
#' @param iter An integer of how many max iterations
#' @param tol_rel_obj A real
#'
#' @return A Stan fit object
#'
vb_iterative = function(model,
                        output_samples = 1000,
                        iter,
                        tol_rel_obj,
                        algorithm = "fullrank",
                        ...) {
  res = NULL
  i = 0
  while (res %>% is.null | i > 5) {
    res = tryCatch({
      my_res = vb(
        model,
        output_samples = output_samples,
        iter = iter,
        tol_rel_obj = tol_rel_obj,
        algorithm = algorithm,
        ...
        #, pars=c("counts_rng", "exposure_rate", additional_parameters_to_save)
      )
      boolFalse <- T
      return(my_res)
    },
    error = function(e) {
      i = i + 1
      writeLines(sprintf("Further attempt with Variational Bayes: %s", e))
      return(NULL)
    },
    finally = {

    })
  }

  return(res)
}

draws_to_tibble = function(fit, par, x) {
  
  par_names = names(fit) %>% grep(sprintf("%s", par), ., value = T)
  
  fit %>%
    rstan::extract(par_names, permuted=F) %>% 
    as.data.frame %>% 
    as_tibble() %>%
    mutate(.iteration = 1:n()) %>% 
    pivot_longer(
      names_to = c("dummy", ".chain", ".variable", x),  
      cols = contains(par), 
      names_sep = "\\.|\\[|,|\\]|:",
      values_to = ".value"
    ) %>%
    mutate(.chain = as.integer(.chain), !!as.symbol(x) := as.integer(!!as.symbol(x))) %>%
    select(-dummy) %>%
    arrange(.variable, !!as.symbol(x),  .chain) %>%
    group_by(.variable, !!as.symbol(x)) %>%
    mutate(.draw = 1:n()) %>%
    ungroup() %>%
    select(!!as.symbol(x), .chain, .iteration, .draw ,.variable ,     .value)
  
}

summary_to_tibble = function(fit, par, x) {
  
  par_names = names(fit) %>% grep(sprintf("%s", par), ., value = T)
  
  fit %>%
    rstan::summary(par_names) %$%
    summary %>%
    as_tibble(rownames = ".variable") %>% tidyr::extract(col = .variable, into = c(".variable", x), "(.+)\\[(.+)\\]", convert = T) 
  
  
}

#' Get matrix from tibble
#'
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @examples
#'
#' library(dplyr)
#'
#' tidybulk::se_mini |> tidybulk() |> select(feature, count) |> head() |> as_matrix(rownames=feature)
#'
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  rownames = enquo(rownames)
  tbl %>%
    
    # Through warning if data frame is not numerical beside the rownames column (if present)
    when(
      do_check &&
        tbl %>%
        # If rownames defined eliminate it from the data frame
        when(!quo_is_null(rownames) ~ (.)[,-1], ~ (.)) %>%
        dplyr::summarise_all(class) %>%
        tidyr::gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% not() %>% any() ~ {
        warning("tidybulk says: there are NON-numerical columns, the matrix will NOT be numerical")
        (.)
      },
      ~ (.)
    ) %>%
    as.data.frame() %>%
    
    # Deal with rownames column if present
    when(
      !quo_is_null(rownames) ~ magrittr::set_rownames(., tbl %>% pull(!!rownames)) %>%
        select(-1),
      ~ (.)
    ) %>%
    
    # Convert to matrix
    as.matrix()
}
