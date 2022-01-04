
# load all necessary functions YOU DO NOT WANT TO CHANGE DEFAULT parameter values for your own run ================
source("dev/jian_R_files/function_jian.R")

produce_cibersortx_bulk_rnaseq_input <- function(.expression_df, .transcript, .sample, .cell_type, .count, 
                                                 .dir, .tree=NULL, .suffix=NULL){
  
  # Args: 
  # .expression_df is a tibble with shape: transcript | sample | cell_type | count
  # .transcript: tell the function which column in the data represents gene symbol
  # .sample: tell the function which column in the data represents sample
  # .cell_type: tell the function which column in the data represents cell_type
  # .count: tell the function which column in the data represents count
  # .tree: optional argument. if .tree = NULL, cell types will be sourced from those present in the .expression_df column
  #       if provided, the cell types used for cibersortx will be the leaves of the tree
  # .dir: the path to the directory where the two output files should be saved e.g. cellsig/dev/
  # .suffix: optional argument which adds suffix to the output filename: reference.txt and phenoclass.txt
  
  .transcript = enquo(.transcript)
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)
  
  .dir <- ifelse(grepl("\\/$", .dir), .dir, paste0(.dir, "/"))
  
  ## ref_names is produced to create header for reference file
  ref_names <- .expression_df %>% 
    {
      if (is.null(.tree)){
        (.)
      } else {
        (.) %>% 
          filter(!!.cell_type %in% as.phylo(.tree)$tip.label)
      }
      
    } %>% 
    mutate(!!.sample := str_replace_all(!!.sample, "\\.", "_")) %>%
    unite(cell_sample, c(!!.cell_type, !!.sample), sep = ".") %>% 
    pull(cell_sample) %>% 
    unique()
  
  
  # create reference file
  .expression_df %>% 
    
    # filter for leaf cell types in the tree
    {
      if (is.null(.tree)){
        (.)
      } else {
        (.) %>% 
          filter(!!.cell_type %in% as.phylo(.tree)$tip.label)
      }
      
    } %>% 
    select(!!.sample, !!.cell_type, !!.count, !!.transcript) %>% 
    mutate(!!.sample := str_replace_all(!!.sample, "\\.", "_")) %>%
    pivot_wider(names_from = c(!!.cell_type, !!.sample), values_from = !!.count, names_sep=".") %>% 
    `names<-`(ref_names %>% 
                str_extract(".*(?=\\.)") %>% 
                prepend(quo_name(.transcript))
    ) %>% 
    # save files as txt files (tab separated files)
    write_tsv(glue("{.dir}reference{.suffix}.txt"))
  
  # create phenoclass file
  tibble(cell_type = 
           {
             if (is.null(.tree)){
               .expression_df %>% 
                 distinct(!!.cell_type) %>% 
                 pull %>% 
                 as.character
             } else {
               as.phylo(.tree)$tip.label
             }
             
           }
         
  ) %>% 
    bind_cols(
      tibble(ref_names, value=1L) %>% 
        pivot_wider(names_from = ref_names, values_from = value)
    ) %>% 
    pivot_longer(-cell_type, names_to="ref_names", values_to="membership") %>% 
    
    # assign membership according to cibersortx tutorial 6: 
    # "1" indicates membership of the reference sample to the class as defined in that row,
    # "2" indicates the class that the sample will be compared against
    mutate(membership = if_else(str_detect(ref_names, cell_type), 1L, 2L)) %>% 
    
    pivot_wider(names_from = ref_names, values_from = membership) %>% 
    
    # set col_names to FALSE following cibersortx format for phenoclass
    write_tsv(glue("{.dir}phenoclass{.suffix}.txt"), col_names = FALSE)
  
}


do_imputation <- function(.scaled_counts, .sample, .symbol, .count, .cell_type){
  
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  .cell_type = enquo(.cell_type)
  
  .scaled_counts %>% 
    # Convert to SE
    # as_SummarizedExperiment(.sample, .feature, count) %>%
    as_SummarizedExperiment(!!.sample, !!.symbol, .abundance = c(!!.count, count_scaled)) %>%
    
    # Hierarchical imputation. Suffix = "" equated to overwrite counts
    impute_missing_abundance(.formula = as.formula(sprintf("~ %s", quo_name(.cell_type))),
                             .abundance = c(!!.count, count_scaled)) %>%
    
    # AUTOMATE THIS!
    # impute_missing_abundance(~ level_5, .abundance = c(!!.count, count_scaled)) %>%
    impute_missing_abundance(~ level_4, .abundance = c(!!.count, count_scaled)) %>%
    impute_missing_abundance(~ level_3, .abundance = c(!!.count, count_scaled)) %>%
    impute_missing_abundance(~ level_2, .abundance = c(!!.count, count_scaled)) %>%
    impute_missing_abundance(~ level_1, .abundance = c(!!.count, count_scaled)) %>% 
    
    # Convert back to tibble
    as_tibble() %>%
    
    mutate(.imputed = if_any(contains("imputed"), ~ .x != 0)) %>% 
    
    select(-matches("imputed\\.\\d"))
}

create_hierarchy_and_calculate_imputation_ratio <- function(.imputed_counts, .level, .sample, .symbol) {
  # this preproces function ranged data in hierarchy(or non_hierarchy) and
  # calculates the imputation ratio for genes in each hierarchy
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  
  
  # load data
  .imputed_counts %>%

    # filter for cells at the level of interest. .level == level_1
    filter(!is.na(!!as.symbol(.level))) %>%
    
    # nest by ancestor
    nest(data = - !!as.symbol(pre(.level)))
  
}


counts_tree_to_gene_markers = function(.input, .sample, .symbol, .count, .cell_type,
                                       .is_hierarchy=TRUE, .level=NULL, 
                                       .tree, .node=NULL,
                                       .contrast_method, .ranking_method, .rank_stat=NULL, .bayes=NULL, 
                                       .selection_method, .kmax=60, .discard_number=2000, .reduction_method = "PCA", .dims=2,
                                       .optimisation_method, .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
                                       .is_complete = TRUE) {
  
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  .cell_type = enquo(.cell_type)
  
  subtree = tree_subset(.tree=.tree, .node=.node)
  
  if ( (!.is_hierarchy)|(.selection_method == "naive")) {
    
    .input %>%
      
      # adapt_tree(.tree = .tree, .node = .node) %>%
      # 
      # tree_and_signatures_to_database(tree=subtree, ., .sample=!!.sample, .cell_type=!!.cell_type,
      #                                 .symbol=!!.symbol, .count=!!.count) %>%
      # 
      # 
      # 
      # # comment out these six lines when using your single cell data, keep when using bulk data
      # # Remove redundant samples
      # remove_redundancy(.element=!!.sample, .feature=!!.symbol, .abundance=!!.count, correlation_threshold = 0.999, top = 500, method = "correlation") %>%
      # droplevels() %>%
      # 
      # # Eliminate suspicious samples
      # filter(!grepl("GSM3722278|GSM3722276|GSM3722277", !!.sample)) %>%
      # 
      # do_scaling(.sample = !!.sample, .symbol= !!.symbol , .count= !!.count, .cell_type=!!.cell_type) %>%
      # 
      # do_imputation(.sample = !!.sample, .symbol=feature, .count = !!.count, .cell_type=!!.cell_type) %>%
      # 
      # rename(symbol = feature) %>% 
      # ##
      

      do_hierarchy(.sample=!!.sample,
                   .symbol=!!.symbol,
                   .cell_type = !!.cell_type,
                   .tree = .tree,
                   .is_hierarchy=.is_hierarchy,
                   .level=.level) %>%
      
      # Input: data.frame columns_1 <int> | ...
      do_ranking(.sample=!!.sample, 
                 .symbol=!!.symbol,
                 .cell_type = !!.cell_type,
                 .ranking_method=.ranking_method, 
                 .contrast_method=.contrast_method, 
                 .rank_stat=.rank_stat, 
                 .bayes=.bayes,
                 .tree = .tree) 
      
      # Input: data.frame columns_1 <int> | ...
      # Output:
      do_selection(.sample=!!.sample,
                   .symbol=!!.symbol,
                   .selection_method=.selection_method,
                   .reduction_method=.reduction_method,
                   .discard_number=.discard_number,
                   .kmax=.kmax,
                   .dims=.dims) %>%

      do_optimisation(.optimisation_method=.optimisation_method, .symbol=!!.symbol) %>%

      format_output(.is_complete = .is_complete)
    
  } else {
    
    .input %>%
      
      # adapt_tree(.tree = .tree, .node = .node) %>%
      # 
      # tree_and_signatures_to_database(tree=subtree, ., .sample=!!.sample, .cell_type=!!.cell_type,
      #                                 .symbol=!!.symbol, .count=!!.count) %>%
      # 
      # 
      # 
      # # comment out these six lines when using your single cell data, keep when using bulk data
      # # Remove redundant samples
      # remove_redundancy(.element=!!.sample, .feature=!!.symbol, .abundance=!!.count, correlation_threshold = 0.999, top = 500, method = "correlation") %>%
      # droplevels() %>%
      # 
      # # Eliminate suspicious samples
      # filter(!grepl("GSM3722278|GSM3722276|GSM3722277", !!.sample)) %>%
      # 
      # do_scaling(.sample = !!.sample, .symbol= !!.symbol , .count= !!.count, .cell_type=!!.cell_type) %>%
      # 
      # do_imputation(.sample = !!.sample, .symbol=feature, .count = !!.count, .cell_type=!!.cell_type) %>%
      # 
      # rename(symbol = feature) %>% 
      # ##
      
      
      # do_hierarchy(.sample=!!.sample,
      #              .symbol=!!.symbol,
      #              .cell_type = !!.cell_type,
      #              .tree = .tree,
      #              .is_hierarchy=.is_hierarchy,
      #              .level=.level) %>%
      # 
      # do_ranking(.sample=!!.sample,
      #            .symbol=!!.symbol,
      #            .cell_type = !!.cell_type,
      #            .ranking_method=.ranking_method,
      #            .contrast_method=.contrast_method,
      #            .rank_stat=.rank_stat,
      #            .bayes=.bayes,
      #            .tree = .tree) %>%
      # 
      # mutate(level.copy = level) %>%
      # nest(data = -level.copy)
      
      mutate(data = map(
        data,
        ~ .x %>%

          # do_selection(.sample=!!.sample,
          #              .symbol=!!.symbol,
          #              .selection_method=.selection_method,
          #              .reduction_method=.reduction_method,
          #              .discard_number=.discard_number,
          #              .kmax=.kmax,
          #              .dims=.dims) %>% 

          do_optimisation(.optimisation_method = .optimisation_method, .symbol=!!.symbol) %>%

          format_output(.is_complete = .is_complete)

      )) %>% 

      unnest(data) %>%
      select(-level.copy)
    
  }
}


# load input data =============================================================================

# MODIFY THIS: load your own tree
zijie_tree = read_yaml("dev/intermediate_data/tree_zijie.yaml") %>% as.Node
tree_zijie_level5 = read_yaml("dev/intermediate_data/tree_zijie_level5.yaml") %>% as.Node

# MODIFY THIS: load your own input data
zijie_input_expression_data

# FOR SINGLE CELL DATA:
# To select markers for your own data, run the code below:
zijie_input_expression_data %>% 
  
  # rename your abundance column to count_scaled 
  # because count_scaled should not be provided by the user but generated from the count column within the function after do_scaling()
  # however you do not need to scale count (i.e. abundance in your data) therefore I have disabled do_scaling and 
  # just rename abundance to count_scaled
  rename(count_scaled = abundance) %>% 
  
  counts_tree_to_gene_markers(
    
    .sample = sample, .symbol = transcript, .cell_type = cell_type_formatted, .count = count_scaled,
    
    .is_hierarchy=TRUE,

    .tree = tree, 

    .contrast_method = pairwise_contrast, 

    .ranking_method = rank_edgR_quasi_likelihood, .rank_stat = "PValue",
    
    .selection_method = "silhouette", .dims = 4,
    
    .optimisation_method = "penalty") %>% 
  
  # save the output signature of the function
  # PLEASE MODIFY THIS LINE OF CODE
  saveRDS("your_directory/name_of_the_output_file.rds", compress = "xz")

# FOR BULK-RNA-seq DATA:
# preprocess data
counts_imputed_zijie <- counts %>%
  
  adapt_tree(.tree = tree_zijie_level5) %>%
  
  tree_and_signatures_to_database(tree=tree_zijie_level5, ., .sample=sample, .cell_type=cell_type,
                                  .symbol=symbol, .count=count) %>% 
  
  # comment out these six lines when using your single cell data, keep when using bulk data
  # Remove redundant samples
  remove_redundancy(.element=sample, .feature=symbol, .abundance=count, 
                    correlation_threshold = 0.999, top = 500, method = "correlation") %>%
  droplevels() %>%
  
  # Eliminate suspicious samples
  filter(!grepl("GSM3722278|GSM3722276|GSM3722277", sample)) %>% 
  
  do_scaling(.sample = sample, .symbol= symbol , .count= count, .cell_type=cell_type)

counts_imputed_zijie <- counts_imputed_zijie %>% 
  
  do_imputation(.sample = sample, .symbol=feature, .count = count, .cell_type= cell_type) %>%
  
  rename(symbol = feature) %>% 
  
  saveRDS("dev/intermediate_data/counts_imputed_zijie.rds", compress = "xz")

# signature selection
# counts_imputed_zijie %>% 
counts_imputed_zijie %>%  
  
  counts_tree_to_gene_markers(
    
    .sample = sample, .symbol = symbol, .cell_type = cell_type, .count = NULL,
    
    .is_hierarchy=TRUE,
    
    .tree = tree_zijie_level5, # change to your tree!!!
    
    .contrast_method = pairwise_contrast, 
    
    .ranking_method = rank_edgR_quasi_likelihood, .rank_stat = "PValue",
    
    .selection_method = "silhouette", .dims=4,
    
    .optimisation_method = "penalty") %>% 
  
  # save the output signature of the function
  # PLEASE MODIFY THIS LINE OF CODE
  saveRDS("dev/intermediate_data/zijie_bulk_signature_before_selection.rds", compress = "xz")

zijie_bulk_signature_final <- zijie_bulk_signature %>% 
  
  counts_tree_to_gene_markers(
    
    .sample = sample, .symbol = symbol, .cell_type = cell_type, .count = NULL,
    
    .is_hierarchy=TRUE,
    
    .tree = tree_zijie_level5, # change to your tree!!!
    
    .contrast_method = pairwise_contrast, 
    
    .ranking_method = rank_edgR_quasi_likelihood, .rank_stat = "PValue",
    
    .selection_method = "silhouette", .dims=4,
    
    .optimisation_method = "penalty") %>% 
  
  # save the output signature of the function
  # PLEASE MODIFY THIS LINE OF CODE
  saveRDS(zijie_bulk_signature_final, "dev/intermediate_data/zijie_bulk_signature.rds", compress = "xz")

# x3 <- 
# zijie_bulk_signature_levels %>% 
#   pluck("data", 3) %>% 
x3 %>% 
  do_silhouette_selection(
    .sample = sample,
    .symbol = symbol,
    .reduction_method = "PCA",
    .discard_number = 10,
    .dims = 4
  )


# Jian's example ====================================================================

# this is my example data
toy_data <- readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/toy_data.rds")

# create a small tree that fits the cell types in the toy data
# cell types on the leaves of the tree must be present in your input data: such as epithelial t_cell, b_cell here
# internal cell types such are immune_cell and Tissue here do not need to be in your input data
jian_tree <- Node$new("Tissue")
jian_tree$AddChild("epithelial")
jian_tree$AddChild("immune_cell")
jian_tree$immune_cell$AddChild("t_cell")
jian_tree$immune_cell$AddChild("b_cell")


jian_example_output = 
  
  # input
  toy_data %>% 
  
  # feed input into the function to obtain cell signatures
  # I have set the parameters to the combination that gives good markers in terms of deconvolution score
  counts_tree_to_gene_markers(
    
    # first specify the names of the essential variables(columns) in the input data frame so that the function knows which is which
    .sample = sample, .symbol = feature, .cell_type = cell_type, .count = count_scaled,
    
    # choose hierarchical method: there are 2 options: hierarchical or non_hierarchical by setting .is_hierarchy to TRUE or FALSE
    # hierarchical method is better for your purpose
    .is_hierarchy=TRUE,
    
    # specify the tree structure
    .tree = jian_tree, 
    
    # specify the contrast method:
    # there are 2 options: mean_contrast and pairwise_contrast
    # choose pairwise_contrast for your purpose
    # essentially to select markers we need to compare the differential expression between genes in different cell types
    # say we have 2 cell types: t_CD4, t_CD8, b_cell
    # pairwise_contrast generates 3*2=6 comparisons in the form t_CD4 - t_CD8, t_CD4 - b_cell, t_CD8 - t_CD4, t_CD8 - b_cell, b_cell - t_CD4, b_cell - t_CD8
    # whereas mean_contrast generates 3 comparisons in the form of t_CD4 - (t_CD8 + b_cell)/2, t_CD8 - (t_CD4 + b_cell)/2, b_cell - (t_CD4 + t_CD8)/2
    # the minus sign - doesn't mean subtraction here but only indicates which compares with which
    # cell type before "-" is called the target cell whereas cell type(s) after the "-" is/are termed background cell type(s)
    # say we have gene Myb, cell type contrast "t_CD4 - t_CD8" means that we compare Myb expression in t_CD4 with that in t_CD8
    # whereas cell type contrast "t_CD4 - (t_CD8 + b_cell)/2" means we compare Myb expression in t_CD4 with the mean of Myb expression in t_CD8 and b_cell
    .contrast_method = pairwise_contrast, 
    
    # specify ranking method:
    # there are 3 options: rank_edgR_quasi_likelihood, rank_edgR_robust_likelihood_ratio, rank_bayes
    # choose rank_edgR_quasi_likelihood
    # for the ranking method rank_edgR_quasi_likelihood, there are 2 ranking statistic: "PValue" and "logFC"
    # choose "PValue"
    # we rank the genes for each cell type contrast by their level of differential expression
    # the ones with lowest PValue or the highest logFC are ranked top 1 (most differentially expressed or most enriched in the target cell type)
    .ranking_method = rank_edgR_quasi_likelihood, .rank_stat = "PValue",
    
    # specify marker selection method:
    # there are 2 methods: "naive" and "silhouette"
    # naive method selects the top k genes from the ranked list for each cell type contrast so that the markers are enriched 
    # in the target cell type of each contrast
    # I set the default value for k to be 60 which means (kmax set to 10 in this example so as not to waste your time)
    # the top k genes (k belongs to [1, 60]) will be iteratively selected from each cell type contrast to make up the signature set for the target cell type
    .selection_method = "naive", .kmax=10,
    
    # specify optimisation method:
    # there are 2 options: "curvature" and "penalty"
    # out of 60 sets of signatures for each target cell type in a contrast which set performs the best?
    # optimisation by "penalty" selects the best signature size based on a metric called penalised silhouette score
    .optimisation_method = "penalty")
