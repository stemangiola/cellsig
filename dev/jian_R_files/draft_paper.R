single_marker_pw_select <- 
  
  function(.contrast, .level, .target_size=NULL, .discard_num=NULL, .method="PCA") {
    
    # initialize variables
    
    contrast_copy <- .contrast %>% 
      mutate(ancestor = !!as.symbol(pre(.level)))
    
    # initialise a signature tibble to store signature markers for each cell type in each iteration
    signature <- tibble(
      ancestor = contrast_copy %>% 
        pull(ancestor)
    ) %>% 
      mutate(markers_cumu = map(ancestor, ~ vector()))
    
    # initialise an output tibble containing all results of interest
    contrast_summary_tb <- vector()
    
    # set the base markers
    contrast_pair_tb0 <- 
      
      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(.level, 1) %>% 
      
      mutate(ancestor = !!as.symbol(pre(.level))) %>% 
      
      nest(data = - ancestor) %>% 
      
      # enables markers from each contrast pair to be entitled to every other contrast so that 
      # all the base genes rather than one base gene per contrast can be incorporated in each iteration
      mutate(data = map(data, ~ .x %>% 
                          nest(markers = - contrast_pretty) %>% 
                          expand(contrast_pretty, markers) %>% 
                          unnest(markers)
      )) %>% 
      # base markers
      mutate(base_markers = map(data, ~ .x %>% pull(symbol) %>% unique() ))
    
    # remove base markers from contrast_copy input before further selection
    contrast_copy <- contrast_copy %>%
      left_join(contrast_pair_tb0, by= "ancestor", suffix=c("", ".y")) %>% 
      select(-ends_with(".y")) %>% 
      mutate(markers = map2(markers, base_markers, ~.x %>% 
                              filter(!symbol %in% .y) ))
    
    sil_score <- tibble(ancestor = contrast_copy %>% 
                          pull(ancestor)) %>% 
      mutate(sil_pre = map_dbl(ancestor, ~ 0))
    
    # counter for number of genes discarded
    j <- map_int(contrast_pair_tb0$base_markers, ~ length(.x))
    
    condition <- 
      purrr::when(
        !is.null(.target_size) ~ any(map_int(signature$markers_cumu, ~ length(.x)) < .target_size),
        !is.null(.discard_num) ~ any(j < .discard_num),
        ~ stop("please supply either target signature size or the numer of genes to be discarded")
      )
    
    while (condition & 
           (!identical(map_int(contrast_copy$markers, ~ nrow(.x)), integer(0))) ) {
      
      contrast_pair_tb <- 
        
        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        contrast_copy %>% 
        
        # select top 1 markers from each contrast, output is an unnested tibble
        sig_select(.level, 1) %>% 
        
        mutate(ancestor = !!as.symbol(pre(.level))) %>% 
        
        # combine markers at the base level before iteration
        bind_rows(contrast_pair_tb0 %>% 
                    select(- base_markers) %>% 
                    unnest(data) ) %>% 
        
        nest(data = - ancestor) %>% 
        
        mutate(data = map(data, ~ .x %>% sil_for_each_contrast(.level, .method) )) %>% 
        
        mutate(is_sil_greater = 
                 map2_lgl(data, ancestor, 
                          ~ if (.x$sil_pair[1]> with(sil_score, sil_pre[ancestor==.y])) 
                          {TRUE} else {FALSE} )) %>% 
        
        mutate(sil_current = map_dbl(data, ~ .x[[1, "sil_pair"]] )) %>% 
        
        mutate(markers_to_filter = map2(data, is_sil_greater, ~
                                          if(.y == TRUE) {
                                            .x %>% 
                                              # select the first marker (with highest silhouette score)
                                              slice(1) %>% 
                                              pull(markers_new) %>% 
                                              unlist()
                                          } else {
                                            .x %>% 
                                              # select all the unique markers from all contrasts for removal
                                              pull(markers_new) %>% 
                                              unlist() %>% 
                                              unique()
                                          } )) %>% 
        # if we want the marker, keep contrast information, else, we don't care if the contrast is correct
        mutate(contrast = map_chr(data, ~ .x$contrast_pretty[[1]])) %>% 
        
        mutate(markers_cumu = map(ancestor, ~ vector())) %>% 
        
        select(-data)
      
      sil_score <- sil_score %>% 
        mutate(sil_pre = map2_dbl(sil_pre, ancestor, ~ 
                                    if (with(contrast_pair_tb, is_sil_greater[ancestor == .y]) )
                                    {with(contrast_pair_tb, sil_current[ancestor == .y])} else{.x} ))
      
      # append the base + 1 markers that result in highest silhouette score
      signature <- signature %>% 
        mutate(markers_cumu = map2(markers_cumu, ancestor, ~ 
                                     if (with(contrast_pair_tb, is_sil_greater[ancestor == .y]) ) 
                                     {.x %>% 
                                         append(
                                           with(contrast_pair_tb, markers_to_filter[ancestor == .y][[1]])
                                         ) %>% 
                                         unique()} else {.x} ))
      
      contrast_summary_tb <- contrast_summary_tb %>% 
        append(
          contrast_pair_tb %>% 
            mutate(markers_cumu = map2(markers_cumu, ancestor,
                                       ~ .x %>% 
                                         append(
                                           with(signature, markers_cumu[ancestor==.y][[1]])
                                         ))) %>% 
            nest(data = - is_sil_greater) %>% 
            with(data[is_sil_greater==TRUE])
        )
      
      contrast_copy <- contrast_copy %>% 
        mutate(markers = map2(markers, ancestor,
                              ~ if(with(contrast_pair_tb, is_sil_greater[ancestor==.y])) {
                                .x %>% 
                                  filter(!symbol %in% with(signature, markers_cumu[ancestor==.y][[1]]))
                              } else {
                                .x %>% 
                                  filter(!symbol %in% with(contrast_pair_tb, markers_to_filter[ancestor==.y][[1]]))
                              } 
        ))
      
      # number of genes discarded for each node
      j <- j + 
        (map_int(contrast_pair_tb$markers_to_filter, ~ length(.x)) -
           map_int(contrast_pair_tb0$base_markers, ~ length(.x)) ) * 
        (!contrast_pair_tb$is_sil_greater)
      
      cat("genes discarded for each node: ", j, "\n")
      cat("genes selected for each node: ", map_int(signature$markers_cumu, ~ length(.x)),  "\n")
    }
    
    return(contrast_summary_tb %>% 
             bind_rows() %>% 
             nest(sig_df = - ancestor))
  }

dat <- structure(list(x = c(1.158362492, 1.1430148, 1.11058971, 1.120573931, 
                            1.149219113, 1.123851641, 1.096910013), 
                      y = c(1.322219295, 1.267171728, 1.252853031, 1.260071388, 
                            1.278753601, 1.276461804, 1.222716471
                      )), .Names = c("x", "y"), class = "data.frame", row.names = c(NA, -7L))
dat %>% 
  cov() %>% 
  eigen() %>% 
  mutate(vec = map_dbl(values, ~ sqrt(.x * qchisq(0.95, 2))))

vec <- sqrt(qchisq(0.95, 2) * eig$values)

pi * vec[1] * vec[2]


x <- list(
  list(-1, x = 1, y = c(2), z = "a"),
  list(-2, x = 4, y = c(5, 6), z = "b"),
  list(-3, x = 8, y = c(9, 10, 11))
  

)
x
product(c(1, 2))
t <- c(1, 2)
product(t)
prod(t)
prod(c(3, 2))

A <-  matrix(c(1:3), nrow = 3, ncol = 3)
A
mean(A)

x <- c(1, 2, 3)
y <- c(4, 5, 6)
df <- tibble(x = x, y = y)
dist(df, method = 'euclidean', upper = T)


tt %>% 
  filter(level == 1) %>% 
  select(data) %>% 
  
  distinct(cell_type) %>% 
  
  # Permute
  mutate(cell_type2 = cell_type) %>% 
  expand(cell_type, cell_type2) %>% 
  filter(cell_type != cell_type2) %>% 
  
  # Create contrasts
  mutate(contrast = sprintf("cell_type%s - cell_type%s", cell_type, cell_type2)) %>%
  pull(contrast)

strings <- rep(c("cell_typeendothelial - cell_typeepithelial", "cell_typeendothelial - cell_typefibroblast"), 5)
random <- rnorm(10)
index <- 1:10
df <- tibble(index, random, strings)
df

df %>% dplyr::distinct(strings, .keep_all = T)
d

df %>% 
  separate(col = strings , into = c('strings', 'strings2'), sep = "-", remove = T) %>% 
  mutate(across(c("strings", "strings2"), ~ trimws(.x)))



df <- tibble()

x <- list(a = list(foo = 1:2, bar = 3:4), b = list(baz = 5:6))
x
map_depth(x, 2, paste, collapse = "/")

cell_types <- tt %>% 
  filter(level == LEVEL) %>% 
  unnest(data) %>% 
  distinct(cell_type) %>% 
  pull() 

prefix = "cell_type"

cell_types <- paste(prefix, cell_types, sep="")

comparison <- 1: length(cell_types)

for(i in 1: length(cell_types) ){
  backgrounds = cell_types[-i]
  divisor = length(backgrounds)
  background = paste(backgrounds, collapse = "+")
  comparison[i] <- sprintf("%s-(%s)/%s", cell_types[i], background, divisor)
}

comparison




select_markers_for_each_contrast = function(.data){
  .data %>%
    
    # Group by contrast. Comparisons both ways.
    pivot_longer(
      cols = contains("___"),
      names_to = c("stats", "contrast"), 
      values_to = ".value", 
      names_sep="___"
    ) %>% 
    
    # Markers selection
    nest(stat_df = -contrast) %>%
    
    # Reshape inside each contrast
    mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value))) %>%
    
    # Rank
    mutate(stat_df = map(stat_df, ~.x %>%
                           filter(FDR < 0.05 & logFC > 2) %>%
                           filter(logCPM > mean(logCPM)) %>%
                           arrange(logFC %>% desc())
                         
    )) %>%
    unnest(stat_df)
  
}

contrast <- function(tt, LEVEL){
  tt %>%
    
    # Investigate one level
    filter(level == LEVEL) %>%
    
    # Differential transcription
    mutate(markers = map(
      data,
      ~ test_differential_abundance(.x,
                                    ~ 0 + cell_type, 
                                    .contrasts = make_contrasts(tt, LEVEL),
                                    action="only") 
    )) %>%
    
    # Select rank from each contrast
    mutate(markers = map(markers, ~ select_markers_for_each_contrast(.x))) %>%
    
    # Select markers from each cell type
    # mutate(markers = map(markers, ~ select_markers_for_each_cell_type(.x, sig_size))) %>% 
    
    # Add original data info to markers
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>%
    select(markers) %>%
    unnest(markers) %>%
    
    # make contrasts pretty
    mutate(contrast_pretty = str_replace(contrast, "cell_type", "") %>% str_replace("cell_type", ""))
}

LEVEL <- 2
all_contrasts <- contrast(tt, LEVEL)
all_contrasts

all_contrasts %>% group_by(cell_type) %>% select(symbol) %>% count()
all_contrasts %>% select(cell_type, symbol) %>% group_by(cell_type) %>% count()

v <- rep(c(1, 2, 3), c(4, 5, 6))
group <- 1:length(v)
df <- tibble(index = index, v = v)
df %>% group_by(v) %>% count()

group <- c("group1, group2", "group3")

group1 <- c(1, 1, 3, 3, 2)
group2 <- c(1, 2, 2, 3, 3)
group3 <- c(1, 1, 1, 2, 3)

## find difference between orig_counts.rda and new counts.rda
library(arsenal)
library(sets)
library(magrittr)
library(VennDiagram)

load("dev/orig_counts.rda")
df_old <- counts
df_old

load("data/counts.rda")
df_new <- counts
df_new

compare <- comparedf(df_old, df_new)

sym_old <- df_old %>% select(symbol) %>% as.set()
sym_new <- df_new %>% select(symbol) %>% as.set()

sets::set_complement(sym_new, sym_old)

counts %>% 
  distinct(level, sample, cell_type) %>%
  count(level, cell_type) %>%
  arrange(level, n) %>% 
  filter(n<9) %>% 
  filter(!grepl("PDGFD", cell_type)) %>% 
  arrange(n)

# Try this:
# count(sample) %>% count(n)
# Should see: nrow == 1

# The Try:
# count(sample, symbol) %>% count(n)
# Should see: # nn == 1
# NO DUPLICATION
# ALL samples have same number of genes


# Debug

all_contrasts_L2 <- readRDS("dev/all_contrasts_L2.rds")

all_contrasts_L2 %>% 
  slice(1:10) %>% 
  mutate(sil_df = map(sil_df, ~ .x %>% sil_func(LEVEL)))

x <- all_contrasts_L2 %>% 
  slice(3) %>% 
  unnest(sil_df)

y <- x %>% 
  nest(pca = - !!as.symbol(pre(LEVEL))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     reduce_dimensions(sample, symbol, count_scaled,
                                       method = "PCA",
                                       action = "add",
                                       transform = log1p)
  ))

y %>% 
  unnest(pca) %>% 
  select(contains("PC")) %>% 
  # str()
  dist()

z %>% 
  unnest(pca) %>% 
  select(contains("PC")) %>% 
  # str()
  dist()

z <- all_contrasts_L2 %>% 
  slice(2) %>% 
  unnest(sil_df) %>% 
  nest(pca = - !!as.symbol(pre(LEVEL))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     reduce_dimensions(sample, symbol, count_scaled,
                                       method = "PCA",
                                       action = "add",
                                       transform = log1p)
  ))

f <- all_contrasts_L2 %>% 
  slice(4) %>% 
  unnest(sil_df) %>% 
  nest(pca = - !!as.symbol(pre(LEVEL))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     distinct(sample, symbol, count_scaled, !!as.symbol(LEVEL)))) %>%
  mutate(pca = map(pca, ~ .x %>%
                     reduce_dimensions(sample, symbol, count_scaled,
                                       method = "PCA",
                                       action = "add",
                                       transform = log1p) ))

a <- z %>% 
  mutate(distance = map(pca, ~ .x %>%
                          select(contains("PC")) %>%
                          dist()
  ))

b <- f %>% mutate(distance = map(pca, ~ .x %>%
                                   select(contains("PC")) %>%
                                   dist() ))


# Plotly =================================================================

df <- tt_L2 %>% 
  contrast_MC("level_2")%>% 
  sig_select("level_2", 41) %>% 
  distinct(sample, symbol, count_scaled, level_2) %>% 
  reduce_dimensions(sample, symbol, count_scaled,
                    method = "PCA", 
                    action = "get",
                    transform = log1p)

df %>% 
  plot_ly(x = ~ PC1, y = ~ PC2, 
          color = ~ level_2, 
          type = "scatter", 
          mode = "markers",
          hoverinfo = "text",
          text = paste("cell_type: ", df$level_2,
                       "<br>",
                       "sample: ", df$sample)
  )

# Shiny ================================
cell_sig_select <- function(.markers) {
  .markers %>% 
    # obtain cell types in a node and nest by it to extract signatures for all cell types
    mutate(cell_type = str_extract(contrast_pretty, "([a-z]|\\_)+(?=\\s)")) %>% 
    nest(signature = - cell_type) %>% 
    mutate(signature = map(signature, ~ .x %>% 
                             pull(symbol) %>% 
                             unique()
    ))
}

x <- sig_collect %>% 
  slice(1) %>% 
  # select individual cell markers for the cell types at each node
  mutate(cell_markers_PW = map(markers_PW, ~ .x %>% cell_sig_select() ))

z <- sig_collect %>% 
  select(level, ancestor_type, markers_PW) %>% 
  mutate(cell_sig_PW = map(markers_PW, ~ cell_sig_select(.x))) 

x <- contrast_all %>% 
  mutate(contrast_PW = map(contrast_PW, ~ .x %>% 
                             mutate(node = 
                                      select(., contains("level_")) %>% 
                                      as_vector()) %>% 
                             select(-contains("level_"))
  )) %>% 
  mutate(contrast_MC = map(contrast_MC, ~ .x %>% 
                             mutate(node = 
                                      select(., contains("level_")) %>% 
                                      as_vector()) %>% 
                             select(-contains("level_"))
  )) %>% 
  unnest(contrast_PW, contrast_MC) %>% 
  select(level, ancestor_type, tt_data = data, 
         contrast_PW = markers, contrast_MC = markers1) %>% 
  mutate(opPCA_sig_PW = opPCA_sig_PW) %>% 
  mutate(opPCA_sig_MC = opPCA_sig_MC) %>% 
  mutate(markers_PW = pmap(list(contrast_PW, level, opPCA_sig_PW), ~ sig_select(..1, ..2, ..3))) %>% 
  mutate(markers_MC = pmap(list(contrast_MC, level, opPCA_sig_MC), ~ sig_select(..1, ..2, ..3)))


x <- contrast_all %>% 
  unnest(tt) %>% 
  mutate(ancestor_type = select(., contains("level_")) %>% 
           pivot_longer(contains("level_"), values_to="cell_type") %>% 
           drop_na() %>% 
           pull(cell_type) ) %>% 
  select(-contains("level_")) %>% 
  
  mutate(contrast_PW = map2(contrast_PW, ancestor_type, ~ .x %>% 
                              mutate(node = 
                                       select(., contains("level_")) %>% 
                                       as_vector()) %>% 
                              filter(node == .y)
  ))

# Plotly

signature_cell <- sig_collect %>% 
  pluck("cell_markers_MC", 1) 

sig_collect %>% 
  pluck("sil_MC_tSNE", 1) %>% 
  pluck("rdim", 1) %>%
  left_join(signature_cell, by = c("level_1" = "cell_type")) %>% 
  plot_ly(x = ~ tSNE1, y = ~ tSNE2) %>% 
  add_markers(
    color = ~ level_1,
    colors = "Set1",
    hoverinfo = "text",
    text = ~ paste("</br>Sample: ", sample,
                   "</br>Signature: ", signature
    )
  ) %>% 
  layout(
    title = sig_collect$ancestor_type[1],
    xaxis = list(zeroline = FALSE),
    yaxis = list(zeroline = FALSE)
  )

# GGplotly

signature_cell <- sig_collect %>% 
  pluck("cell_markers_MC", 1) 

p <- sig_collect %>% 
  pluck("sil_MC_tSNE", 1) %>% 
  pluck("rdim", 1) %>%
  left_join(signature_cell, by = c("level_1" = "cell_type")) %>% 
  ggplot(aes(x = tSNE1, y = tSNE2, colour = level_1, label = signature)) +
  geom_point() +
  ggtitle("cell") +
  theme(plot.title = element_text(hjust = 0.5)  )

ggplotly(p, tooltip = c("level_1", "label"))

plot_ly(x = ~ tSNE1, y = ~ tSNE2) %>% 
  add_markers(
    color = ~ level_1,
    colors = "Set1",
    hoverinfo = "text",
    text = ~ paste("</br>Sample: ", sample,
                   "</br>Signature: ", signature
    )
  ) %>% 
  layout(
    title = sig_collect$ancestor_type,
    xaxis = list(zeroline = FALSE),
    yaxis = list(zeroline = FALSE)
  )


x <- sig_collect %>% pluck("sil_MC_tSNE", 1) %>% pluck("rdim", 1)
y <- sig_collect %>% pluck("cell_markers_MC", 1) %>% pluck("signature", 1)



x %>% nest(info = - level_1) %>% 
  left_join(y, by = c("level_1"="cell_type"))

x %>% left_join(y, by = c("level_1"="cell_type"))

y <- signature_cell %>% 
  mutate(signature = map_chr(signature, ~ .x %>% str_c(collapse = "\n")  ))

x <- signature_cell %>% pluck("signature", 1)


str_c(y, sep = "\n") %>% cat()
x <- str_c(y, sep = "\n")
x <-  c("a", "b", "c", "d", "e")

sig_size <- 2

PCA2_immune_sil <- sil_tb %>% 
  pluck("sil_df", 19) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  ggtitle("sig_size=19; real_size=290; sil=0.493") +
  theme_bw()

immune_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point() +
  annotate("text", x = c(46, 176, 201, 261, 290),
           y = c(0.471, 0.473, 0.472, 0.493, 0.493),
           label = c("46", "176", "201", "261", "290"),
           vjust = -1) +
  annotate(geom = "point", x = c(46, 176, 201, 261, 290), 
           y = c(0.471, 0.473, 0.472, 0.493, 0.493),
           colour = "red", size = 3, alpha=0.5) +
  ggtitle("level_2, PW+H")


tCD4_sil <- sil_data %>%
  pluck("plot_data", 1) %>% 
  ggplot(aes(x=real_size, y=sil)) +
  geom_line() +
  geom_point() +
  annotate("text", x = c(2, 10, 18),
           y = c(0.95, 0.924, 0.891),
           label = c("2", "10", "18"),
           hjust = -1) +
  annotate(geom = "point", x = c(2, 10, 18), 
           y = c(0.95, 0.924, 0.891),
           colour = "red", size = 3, alpha=0.5) +
  ggtitle("level_4, PW+H, tCD4")

PCA4_tCD4_sil <- sil_tb %>% 
  pluck("sil_df", 9) %>% 
  pluck("rdim", 1) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("sig_size=9, real_size=18, sil=0.891")

macro_sil <- sil_data %>%
  pluck("plot_data", 2) %>% 
  ggplot(aes(real_size, sil)) +
  geom_line() +
  geom_point() +
  annotate("text", x = c(2, 16, 28, 30),
           y = c(0.652, 0.582, 0.539, 0.57),
           label = c("2", "16", "28", "30"),
           hjust = -1) +
  annotate(geom = "point", x = c(2, 16, 28, 30), 
           y = c(0.652, 0.582, 0.539, 0.57),
           colour = "red", size = 3, alpha=0.5) +
  ggtitle("level_4, PW+H, macrophage")

PCA4_macro_sil <- sil_tb %>% 
  pluck("sil_df", 15) %>% 
  pluck("rdim", 2) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = !!as.symbol(LEVEL), label = sample)) + 
  geom_point() +
  stat_ellipse(type = 't')+
  theme_bw() +
  ggtitle("sig_size=15, real_size=30, sil=0.57")

final_tSNE %>% 
  ggplot(aes(sig_size, sil, 
             color=analysis,
             shape = analysis) ) +
  geom_line(position = position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width=0.5)) +
  facet_wrap(~ ancestor_type) +
  ggtitle("Compare four marker selection strategies by silhouette score over 20 signature sizes") +
  theme(plot.title = element_text(hjust = 0.5))

a <- tibble()
a %>% add_row(y %>% slice(1))

# unused & old functions

#  1 base each contrast + 1 new marker at a time
single_marker_pw_select2 <- function(.contrast, .target_size, .level, .method) {
  
  # initialize variables
  contrast_copy <- .contrast %>% 
    mutate(ancestor = !!as.symbol(pre(.level)))
  
  signature <- tibble(ancestor = contrast_copy %>% 
                        pull(ancestor)) %>% 
    mutate(markers_cumu = map(ancestor, ~ vector()))
  
  contrast_summary_tb <- tibble()
  
  # set the base markers
  contrast_pair_tb0 <- 
    
    # contrast_copy contains all the statistics of all cell_type contrasts for each gene
    contrast_copy %>% 
    
    # select top 1 markers from each contrast, output is an unnested tibble
    sig_select(.level, 1) %>% 
    
    mutate(ancestor = !!as.symbol(pre(.level))) %>% 
    
    nest(data = - ancestor) %>% 
    
    # base markers
    mutate(base_markers = map(data, ~ .x %>% pull(symbol) %>% unique() ))
  
  # remove base markers from contrast_copy input before further selection
  contrast_copy <- contrast_copy %>%
    left_join(contrast_pair_tb0, by= "ancestor", suffix=c("", ".y")) %>% 
    select(-ends_with(".y")) %>% 
    mutate(markers = map2(markers, base_markers, ~.x %>% 
                            filter(!symbol %in% .y) ))
  
  while (any(map_int(signature$markers_cumu, ~ length(.x)) < .target_size)) {
    
    contrast_pair_tb <- 
      
      # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(.level, 1) %>% 
      
      # combine markers at the base level before iteration
      bind_rows(contrast_pair_tb0 %>% 
                  select(- base_markers) %>% 
                  unnest(data) ) %>% 
      
      mutate(ancestor = !!as.symbol(pre(.level))) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(data = map(data, ~ .x %>% 
                          
                          # nest by pairwise contrast instead of ancestor level before
                          nest(markers = - contrast_pretty) %>% 
                          
                          # calculate silhouette score for PCA of cell types resolved by the selected 1 markers
                          mutate(sil_pair = map(markers, ~ sil_func(.x, .level, .method))) %>% 
                          
                          # make the silhouette score explicit
                          mutate(sil_pair = map_dbl(sil_pair, ~ .x$sil)) %>% 
                          
                          # rank silhouette score in a descending manner
                          arrange(desc(sil_pair)) %>% 
                          
                          # extract all the new markers from each contrast
                          mutate(markers_new = map(markers, ~ .x %>% pull(symbol) %>% unique())) %>%
                          
                          # create a column of empty vectors to store cumulative markers below
                          mutate(markers_cumu = map(contrast_pretty, ~ vector()) )
      )) %>% 
      mutate(data = map(data, ~ .x %>% slice(1))) %>% 
      unnest(data)
    
    # append the base + 1 markers that result in highest silhouette score
    signature <- signature %>% 
      mutate(markers_cumu = map2(markers_cumu, ancestor, ~ .x %>% 
                                   append(
                                     contrast_pair_tb$markers_new[[which(contrast_pair_tb$ancestor==.y)]]
                                   ) %>% 
                                   unique()
      ))
    
    # store the selected contrast pair, markers and cumulative markers as a tibble
    contrast_summary_tb <- contrast_summary_tb %>% 
      bind_rows(
        contrast_pair_tb %>% 
          mutate(markers_cumu = map2(markers_cumu, ancestor, 
                                     ~ .x %>% 
                                       append(signature$markers_cumu[[which(signature$ancestor==.y)]]) 
          ))
      )
    
    # remove the selected markers from contrast_PW_L1 so that other markers can be selected in later iterations
    contrast_copy <- contrast_copy %>%
      mutate(markers = map2(markers, ancestor, ~.x %>% 
                              filter(!symbol %in% signature$markers_cumu[[which(signature$ancestor==.y)]]) ))
  }
  
  return(contrast_summary_tb %>% nest(sig_df = - ancestor))
}

# 2 markers at a time
pair_marker_pw_select <- function(.contrast, .target_size, .level, .method) {
  
  # initialize variables
  contrast_copy <- .contrast %>% 
    mutate(ancestor = !!as.symbol(pre(.level)))
  
  signature <- tibble(ancestor = contrast_copy %>% 
                        pull(ancestor)) %>% 
    mutate(markers_cumu = map(ancestor, ~ vector()))
  
  contrast_summary_tb <- tibble()
  
  while (any(map_int(signature$markers_cumu, ~ length(.x)) < .target_size)) {
    
    contrast_pair_tb <- 
      
      # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(.level, 2) %>% 
      
      mutate(ancestor = !!as.symbol(pre(.level))) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(data = map(data, ~ .x %>% 
                          
                          # nest by pairwise contrast instead of ancestor level before
                          nest(markers = - contrast_pretty) %>% 
                          
                          # calculate silhouette score for PCA of cell types resolved by the selected 1 markers
                          mutate(sil_pair = map(markers, ~ sil_func(.x, .level, .method))) %>% 
                          
                          # make the silhouette score explicit
                          mutate(sil_pair = map_dbl(sil_pair, ~ .x$sil)) %>% 
                          
                          # rank silhouette score in a descending manner
                          arrange(desc(sil_pair)) %>% 
                          
                          # extract all the new markers from each contrast
                          mutate(markers_new = map(markers, ~ .x %>% pull(symbol) %>% unique())) %>%
                          
                          # create a column of empty vectors to store cumulative markers below
                          mutate(markers_cumu = map(contrast_pretty, ~ vector()) )
      )) %>% 
      mutate(data = map(data, ~ .x %>% slice(1))) %>% 
      unnest(data)
    
    # append the base + 1 markers that result in highest silhouette score
    signature <- signature %>% 
      mutate(markers_cumu = map2(markers_cumu, ancestor, ~ .x %>% 
                                   append(
                                     contrast_pair_tb$markers_new[[which(contrast_pair_tb$ancestor==.y)]]
                                   ) %>% 
                                   unique()
      ))
    
    # store the selected contrast pair, markers and cumulative markers as a tibble
    contrast_summary_tb <- contrast_summary_tb %>% 
      bind_rows(
        contrast_pair_tb %>% 
          mutate(markers_cumu = map2(markers_cumu, ancestor, 
                                     ~ .x %>% 
                                       append(signature$markers_cumu[[which(signature$ancestor==.y)]]) 
          ))
      )
    
    # remove the selected markers from contrast_PW_L1 so that other markers can be selected in later iterations
    contrast_copy <- contrast_copy %>%
      mutate(markers = map2(markers, ancestor, ~.x %>% 
                              filter(!symbol %in% signature$markers_cumu[[which(signature$ancestor==.y)]]) ))
  }
  
  return(contrast_summary_tb %>% nest(sig_df = - ancestor))
}

# all bases + 1 new marker at a time
single_marker_pw_select <- function(.contrast, .target_size, .level, .method) {
  
  # initialize variables
  contrast_copy <- .contrast %>% 
    mutate(ancestor = !!as.symbol(pre(.level)))
  
  signature <- tibble(ancestor = contrast_copy %>% 
                        pull(ancestor)) %>% 
    mutate(markers_cumu = map(ancestor, ~ vector()))
  
  contrast_summary_tb <- tibble()
  
  # set the base markers
  contrast_pair_tb0 <- 
    
    # contrast_copy contains all the statistics of all cell_type contrasts for each gene
    contrast_copy %>% 
    
    # select top 1 markers from each contrast, output is an unnested tibble
    sig_select(.level, 1) %>% 
    
    mutate(ancestor = !!as.symbol(pre(.level))) %>% 
    
    nest(data = - ancestor) %>% 
    
    # enables markers from each contrast pair to be entitled to every other contrast so that 
    # all the base genes rather than one base gene per contrast can be incorporated in each iteration
    mutate(data = map(data, ~ .x %>% 
                        nest(markers = - contrast_pretty) %>% 
                        expand(contrast_pretty, markers) %>% 
                        unnest(markers)
    )) %>% 
    # base markers
    mutate(base_markers = map(data, ~ .x %>% pull(symbol) %>% unique() ))
  
  # remove base markers from contrast_copy input before further selection
  contrast_copy <- contrast_copy %>%
    left_join(contrast_pair_tb0, by= "ancestor", suffix=c("", ".y")) %>% 
    select(-ends_with(".y")) %>% 
    mutate(markers = map2(markers, base_markers, ~.x %>% 
                            filter(!symbol %in% .y) ))
  
  while (any(map_int(signature$markers_cumu, ~ length(.x)) < .target_size)) {
    
    contrast_pair_tb <- 
      
      # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(.level, 1) %>% 
      
      # combine markers at the base level before iteration
      bind_rows(contrast_pair_tb0 %>% 
                  select(- base_markers) %>% 
                  unnest(data) ) %>% 
      
      mutate(ancestor = !!as.symbol(pre(.level))) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(data = map(data, ~ .x %>% 
                          
                          # nest by pairwise contrast instead of ancestor level before
                          nest(markers = - contrast_pretty) %>% 
                          
                          # calculate silhouette score for PCA of cell types resolved by the selected 1 markers
                          mutate(sil_pair = map(markers, ~ sil_func(.x, .level, .method))) %>% 
                          
                          # make the silhouette score explicit
                          mutate(sil_pair = map_dbl(sil_pair, ~ .x$sil)) %>% 
                          
                          # rank silhouette score in a descending manner
                          arrange(desc(sil_pair)) %>% 
                          
                          # extract all the new markers from each contrast
                          mutate(markers_new = map(markers, ~ .x %>% pull(symbol) %>% unique())) %>%
                          
                          # create a column of empty vectors to store cumulative markers below
                          mutate(markers_cumu = map(contrast_pretty, ~ vector()) )
      )) %>% 
      mutate(data = map(data, ~ .x %>% slice(1))) %>% 
      unnest(data)
    
    # append the base + 1 markers that result in highest silhouette score
    signature <- signature %>% 
      mutate(markers_cumu = map2(markers_cumu, ancestor, ~ .x %>% 
                                   append(
                                     contrast_pair_tb$markers_new[[which(contrast_pair_tb$ancestor==.y)]]
                                   ) %>% 
                                   unique()
      ))
    
    # store the selected contrast pair, markers and cumulative markers as a tibble
    contrast_summary_tb <- contrast_summary_tb %>% 
      bind_rows(
        contrast_pair_tb %>% 
          mutate(markers_cumu = map2(markers_cumu, ancestor, 
                                     ~ .x %>% 
                                       append(signature$markers_cumu[[which(signature$ancestor==.y)]]) 
          ))
      )
    
    # remove the selected markers from contrast_PW_L1 so that other markers can be selected in later iterations
    contrast_copy <- contrast_copy %>%
      mutate(markers = map2(markers, ancestor, ~.x %>% 
                              filter(!symbol %in% signature$markers_cumu[[which(signature$ancestor==.y)]]) ))
  }
  
  return(contrast_summary_tb %>% nest(sig_df = - ancestor))
}

# old sil_score_for_markers function
sil_score_for_markers <-function(.sig_df, .contrast, .level, .method) {
  
  .contrast %>%
    
    mutate(ancestor = !!as.symbol(pre(.level))) %>% 
    
    left_join(.sig_df, by="ancestor") %>% 
    
    unnest(sig_df, names_repair = "universal") %>% 
    
    # filter markers that are in the signature
    mutate(markers...3 = map2(markers...3, markers_cumu, 
                              ~.x %>% 
                                filter(symbol %in% .y))) %>% 
    
    # format statistics from pairwise contrast
    mutate(markers...3  = map(markers...3, 
                              ~ .x %>% 
                                # Group by contrast. Comparisons both ways.
                                pivot_longer(
                                  cols = contains("___"),
                                  names_to = c("stats", "contrast"), 
                                  values_to = ".value", 
                                  names_sep="___"
                                ) %>% 
                                
                                # Markers selection within each pair of contrast
                                nest(stat_df = -contrast) %>%
                                
                                # Reshape inside each contrast
                                mutate(stat_df = map(stat_df, ~.x %>% 
                                                       pivot_wider(names_from = stats, 
                                                                   values_from = .value))) %>% 
                                
                                unnest(stat_df) )) %>% 
    
    # Add original data data to the markers selected
    mutate(markers = map2(markers...3, data, ~ left_join(.x, .y))) %>% 
    
    # select only columns needed
    select(-c(data, markers...3, markers...6)) %>% 
    
    nest(sil_df = c(!!as.symbol(pre(.level)), markers)) %>% 
    
    mutate(sil_df = map(sil_df, ~ .x %>% unnest(markers))) %>% 
    
    # Calculate silhouette score for PCA plot resulted from the markers selected
    mutate(sil_df = map(sil_df, ~ sil_func(.x, .level, .method))) %>% 
    
    mutate(sil = map_dbl(sil_df, ~ .x$sil)) %>% 
    mutate(real_size = map_int(sil_df, ~ .x$real_size)) %>% 
    
    nest(data = - ancestor)
}



# to be modified
single_marker_pw_select <- function(.contrast, .target_size, .level, .method) {
  
  # initialize variables
  
  contrast_copy <- contrast_PW_L2 %>% 
    mutate(ancestor = !!as.symbol(pre(LEVEL)))
  
  # initialise a signature tibble to store signature markers for each cell type in each iteration
  signature <- tibble(
    ancestor = contrast_copy %>% 
      pull(ancestor)
  ) %>% 
    mutate(markers_cumu = map(ancestor, ~ vector()))
  
  # initialise an output tibble containing all results of interest
  contrast_summary_tb <- vector()
  
  # set the base markers
  contrast_pair_tb0 <- 
    
    # contrast_copy contains all the statistics of all cell_type contrasts for each gene
    contrast_copy %>% 
    
    # select top 1 markers from each contrast, output is an unnested tibble
    sig_select(LEVEL, 1) %>% 
    
    mutate(ancestor = !!as.symbol(pre(LEVEL))) %>% 
    
    nest(data = - ancestor) %>% 
    
    # enables markers from each contrast pair to be entitled to every other contrast so that 
    # all the base genes rather than one base gene per contrast can be incorporated in each iteration
    mutate(data = map(data, ~ .x %>% 
                        nest(markers = - contrast_pretty) %>% 
                        expand(contrast_pretty, markers) %>% 
                        unnest(markers)
    )) %>% 
    # base markers
    mutate(base_markers = map(data, ~ .x %>% pull(symbol) %>% unique() ))
  
  # remove base markers from contrast_copy input before further selection
  contrast_copy <- contrast_copy %>%
    left_join(contrast_pair_tb0, by= "ancestor", suffix=c("", ".y")) %>% 
    select(-ends_with(".y")) %>% 
    mutate(markers = map2(markers, base_markers, ~.x %>% 
                            filter(!symbol %in% .y) ))
  
  sil_score <- tibble(ancestor = contrast_copy %>% 
                        pull(ancestor)) %>% 
    mutate(sil_pre = map_dbl(ancestor, ~ 0))
  
  # i <- rep(0, nrow(contrast_copy))
  
  j <- map_int(contrast_pair_tb0$base_markers, ~ length(.x))
  
  while (any(map_int(signature$markers_cumu, ~ length(.x)) < SIZE) & 
         any(j < 500) & (!identical(map_int(contrast_copy$markers, ~ nrow(.x)), integer(0))) ) {
    
    contrast_pair_tb <- 
      
      # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>% 
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(LEVEL, 1) %>% 
      
      mutate(ancestor = !!as.symbol(pre(LEVEL))) %>% 
      
      # combine markers at the base level before iteration
      bind_rows(contrast_pair_tb0 %>% 
                  select(- base_markers) %>% 
                  unnest(data) ) %>% 
      
      nest(data = - ancestor) %>% 
      
      mutate(data = map(data, ~ .x %>% sil_for_each_contrast(LEVEL, METHOD) )) %>% 
      
      mutate(is_sil_greater = 
               map2_lgl(data, ancestor, 
                        ~ if (.x$sil_pair[1]> with(sil_score, sil_pre[ancestor==.y])) 
                        {TRUE} else {FALSE} )) %>% 
      
      mutate(sil_current = map_dbl(data, ~ .x[[1, "sil_pair"]] )) %>% 
      
      mutate(markers_to_filter = map2(data, is_sil_greater, ~
                                        if(.y == TRUE) {
                                          .x %>% 
                                            # select the first marker (with highest silhouette score)
                                            slice(1) %>% 
                                            pull(markers_new) %>% 
                                            unlist()
                                        } else {
                                          .x %>% 
                                            # select all the unique markers from all contrasts for removal
                                            pull(markers_new) %>% 
                                            unlist() %>% 
                                            unique()
                                        } )) %>% 
      # if we want the marker, keep contrast information, else, we don't care if the contrast is correct
      mutate(contrast = map_chr(data, ~ .x$contrast_pretty[[1]])) %>% 
      
      mutate(markers_cumu = map(ancestor, ~ vector())) %>% 
      
      select(-data)
    
    sil_score <- sil_score %>% 
      mutate(sil_pre = map2_dbl(sil_pre, ancestor, ~ 
                                  if (with(contrast_pair_tb, is_sil_greater[ancestor == .y]) )
                                  {with(contrast_pair_tb, sil_current[ancestor == .y])} else{
                                    .x} ))
    
    # append the base + 1 markers that result in highest silhouette score
    signature <- signature %>% 
      mutate(markers_cumu = map2(markers_cumu, ancestor, ~ 
                                   if (with(contrast_pair_tb, is_sil_greater[ancestor == .y]) ) 
                                   {.x %>% 
                                       append(
                                         with(contrast_pair_tb, markers_to_filter[ancestor == .y][[1]])
                                       ) %>% 
                                       unique()} else {.x}
      ))
    
    contrast_summary_tb <- contrast_summary_tb %>% 
      append(
        contrast_pair_tb %>% 
          mutate(markers_cumu = map2(markers_cumu, ancestor,
                                     ~ .x %>% 
                                       append(
                                         with(signature, markers_cumu[ancestor==.y][[1]])
                                       ))) %>% 
          nest(data = - is_sil_greater) %>% 
          with(data[is_sil_greater==TRUE])
      )
    
    contrast_copy <- contrast_copy %>% 
      mutate(markers = map2(markers, ancestor,
                            ~ if(with(contrast_pair_tb, is_sil_greater[ancestor==.y])) {
                              .x %>% 
                                filter(!symbol %in% with(signature, markers_cumu[ancestor==.y][[1]]))
                            } else {
                              .x %>% 
                                filter(!symbol %in% with(contrast_pair_tb, markers_to_filter[ancestor==.y][[1]]))
                            } 
      ))
    
    # number of genes selected for each node (THIS IS WRONG, such cumulative counter contains duplicates of base markers)
    # i <- i + map_int(contrast_pair_tb$markers_to_filter, ~ length(.x)) * contrast_pair_tb$is_sil_greater
    
    # number of genes discarded for each node
    j <- j + 
      (map_int(contrast_pair_tb$markers_to_filter, ~ length(.x)) -
         map_int(contrast_pair_tb0$base_markers, ~ length(.x)) ) * 
      (!contrast_pair_tb$is_sil_greater)
    
    cat("genes discarded for each node: ", j, "\n")
    cat("genes selected for each node: ", map_int(signature$markers_cumu, ~ length(.x)),  "\n")
  }
  
  return(contrast_summary_tb %>% 
           bind_rows() %>% 
           nest(sig_df = - ancestor))
}

sil_score_for_markers <-function(.sig_df, .contrast, .level, .method) {
  
  contrast_PW_L3 %>%
    
    mutate(ancestor = !!as.symbol(pre(LEVEL))) %>% 
    
    left_join(x, by="ancestor") %>% 
    
    unnest(sig_df) %>% 
    
    # filter markers that are in the signature (markers_cumu here)
    mutate(markers = map2(markers, markers_cumu, 
                          ~.x %>% 
                            filter(symbol %in% .y))) %>% 
    
    # format statistics from pairwise contrast
    mutate(markers  = map(markers, 
                          ~ .x %>% 
                            # Group by contrast. Comparisons both ways.
                            pivot_longer(
                              cols = contains("___"),
                              names_to = c("stats", "contrast"), 
                              values_to = ".value", 
                              names_sep="___"
                            ) %>% 
                            
                            # Markers selection within each pair of contrast
                            nest(stat_df = -contrast) %>%
                            
                            # Reshape inside each contrast
                            mutate(stat_df = map(stat_df, ~.x %>% 
                                                   pivot_wider(names_from = stats, 
                                                               values_from = .value))) %>% 
                            
                            unnest(stat_df) )) %>% 
    
    # Add original data data to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>% 
    
    # select only columns needed
    select(-data) %>% 
    
    nest(sil_df = c(!!as.symbol(pre(LEVEL)), markers)) %>% 
    
    mutate(sil_df = map(sil_df, ~ .x %>% unnest(markers))) %>% 
    
    # Calculate silhouette score for PCA plot resulted from the markers selected
    mutate(sil_df = map(sil_df, ~ sil_func(.x, LEVEL, METHOD))) %>% 
    
    mutate(sil = map_dbl(sil_df, ~ .x$sil)) %>% 
    mutate(real_size = map_int(sil_df, ~ .x$real_size)) %>% 
    
    nest(data = - ancestor)
}

sig_select(contrast_PW_L4, LEVEL, 1)

contrast_PW_L4 %>% 
  pluck("markers", 4) %>% 
  select_markers_for_each_contrast(1)
# new method

# pairwise selection 1 marker at a time
single_marker_pw_select <- 
  
  function(.contrast, .level, .target_size=NULL, .discard_num=NULL, .method="PCA") {
    
    .target_size = 100
    .discard_num = NULL
    
    # initialize variables
    
    contrast_copy <- contrast_PW_L3 %>% 
      mutate(ancestor = !!as.symbol(pre(LEVEL)))
    
    # initialise a signature tibble to store signature markers for each cell type in each iteration
    signature <- tibble(
      ancestor = contrast_copy %>% 
        pull(ancestor)
    ) %>% 
      mutate(cumulative_signature = map(ancestor, ~ vector())) %>% 
      mutate(last_silhouette = 0)
    
    # initialise an output tibble containing all results of interest
    summary_tb <- tibble(
      ancestor = character(),
      data = list(),
      new_challengers = list(),
      winner = list(),
      cumulative_signature = list()
    )
    
    # set the base markers
    contrast_pair_tb0 <-
      
      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>%
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(LEVEL, 1) %>%
      
      mutate(ancestor = !!as.symbol(pre(LEVEL))) %>%
      
      nest(data = - ancestor) %>%
      
      mutate(data = map(data, ~ .x %>% 
                          nest(markers = - contrast_pretty) %>% 
                          
                          mutate(new_challengers = map_chr(
                            markers, ~ .x %>% pull(symbol) %>% unique()))
      )) %>% 
      
      mutate(new_challengers = map(data, ~.x %>% pull(new_challengers))) %>% 
      
      mutate(data = map2(data, new_challengers,
                         ~ .x %>% 
                           mutate(challengers_for_silhouette = list(.y))
      )) %>% 
      
      mutate(data = map2(
        data, ancestor,
        ~ .x %>% 
          mutate(silhouette = map2(
            challengers_for_silhouette, .y,
            ~ silhouette_for_markers(.x, .y, contrast_PW_L3, LEVEL, METHOD)
          )) %>% 
          unnest(silhouette)
      )) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          mutate(is_greater = map2_lgl(
            sil, ancestor,
            ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
          ))
      )) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% select(-ancestor)
      )) %>% 
      
      mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
      
      mutate(cumulative_signature = winner)
    
    
    signature <- signature %>%
      mutate(cumulative_signature = map2(
        cumulative_signature, ancestor,
        ~ .x %>% 
          append(with(contrast_pair_tb0, cumulative_signature[ancestor==.y][[1]]))
      )) %>% 
      mutate(last_silhouette = map_dbl(
        ancestor,
        ~ with(contrast_pair_tb0$data[ancestor==.x][[1]], unique(sil))
      ))
    
    summary_tb <- summary_tb %>% 
      bind_rows(contrast_pair_tb0)
    
    # remove base markers from contrast_copy input before further selection
    contrast_copy <- contrast_copy %>%
      mutate(markers = map2(
        markers, ancestor, 
        ~ .x %>%
          filter(!symbol %in% with(signature, cumulative_signature[ancestor==.y][[1]]))
      ))
    
    # counter for number of genes discarded
    j <- map_int(signature$cumulative_signature, ~ length(.x))
    
    condition <- when(
      !is.null(.target_size) & is.null(.discard_num) ~ any(map_int(signature$cumulative_signature, ~ length(.x)) < .target_size),
      !is.null(.discard_num) & is.null(.target_size) ~ any(j < .discard_num),
      ~ stop("please supply either target signature size or the numer of genes to be discarded")
    )
    
    while (condition & all(map_int(contrast_copy$markers, ~ nrow(.x)) > 1)) {
      
      contrast_pair_tb <- 
        
        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        my.list$contrast_copy %>% 
        
        # select top 1 markers from each contrast, output is an unnested tibble
        sig_select(LEVEL, 1) %>% 
        
        mutate(ancestor = !!as.symbol(pre(LEVEL))) %>% 
        
        nest(data = - ancestor) %>% 
        
        mutate(data = map(data, ~ .x %>% 
                            nest(markers = - contrast_pretty) %>% 
                            mutate(new_challengers = map_chr(
                              markers, 
                              ~.x %>% pull(symbol) %>% unique()
                            ))
        )) %>% 
        
        unnest(data) %>% 
        
        mutate(challengers_for_silhouette = map2(
          new_challengers, ancestor, 
          ~ with(signature, cumulative_signature[ancestor==.y][[1]]) %>% 
            append(.x)
        )) %>% 
        
        mutate(silhouette = map2(
          challengers_for_silhouette, ancestor, 
          ~ silhouette_for_markers(.x, .y, contrast_PW_L3, LEVEL, METHOD) %>% 
            select(-ancestor)
        )) %>% 
        
        unnest(silhouette) %>% 
        
        group_by(ancestor) %>% 
        
        arrange(desc(sil), .by_group = TRUE) %>% 
        
        ungroup() %>% 
        
        mutate(is_greater = map2_lgl(
          sil, ancestor, 
          ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
        )) %>% 
        
        nest(data = - ancestor) %>% 
        
        mutate(new_challengers = map(
          data,
          ~ .x %>% 
            pull(new_challengers)
        )) %>% 
        
        mutate(winner = map(
          data,
          ~ if(.x[1, ]$is_greater){
            .x[1, ]$new_challengers
          } else {NA}
        )) %>% 
        
        mutate(cumulative_signature = pmap(
          list(data, winner, ancestor),
          ~ if(!is.na(..2)) {
            with(..1[1, ], challengers_for_silhouette[[1]])
          } else {
            with(signature, cumulative_signature[ancestor==..3][[1]])
          }
        ))
      
      # append the base + 1 markers that result in highest silhouette score
      signature <- signature %>% 
        
        mutate(cumulative_signature = map(
          ancestor,
          ~ with(contrast_pair_tb, cumulative_signature[ancestor==.x][[1]])
        )) %>% 
        
        mutate(last_silhouette = map2_dbl(
          ancestor, last_silhouette,
          ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x]))) {
            with(contrast_pair_tb, data[ancestor==.x][[1]][[1, "sil"]])
          } else {.y}
        ))
      
      summary_tb <- summary_tb %>% 
        bind_rows(
          contrast_pair_tb %>% 
            filter(!is.na(winner))
        )
      
      contrast_copy <- contrast_copy %>% 
        mutate(markers = map2(
          markers, ancestor, 
          ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y]))){
            .x %>% 
              filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))
          } else {
            .x %>% 
              filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]]))
          }
        ))
      
      # number of genes discarded for each node
      j <- j + 
        
        # unsuccessful candidates
        map_int(contrast_pair_tb$new_challengers, ~length(.x)) *
        is.na(contrast_pair_tb$winner) +
        
        # winning candidates
        map_int(contrast_pair_tb$winner, ~length(.x)) *
        !is.na(contrast_pair_tb$winner)
      
      cat("genes discarded for each node: ", j, "\n")
      cat("genes selected for each node: ", map_int(signature$cumulative_signature, ~ length(.x)),  "\n")
    }
    
    return(summary_tb %>% 
             nest(sig_df = - ancestor))
  }

# calculates silhouette score for each set of signature (cumulative markers) at a signature size
silhouette_for_markers <-function(.signature, .ancestor, .contrast, .level, .method) {
  
  .contrast %>%
    
    mutate(ancestor = !!as.symbol(pre(.level))) %>% 
    
    filter(ancestor == .ancestor) %>% 
    
    # filter markers that are in the signature
    mutate(markers = map2(markers, ancestor, ~.x %>% 
                            filter(symbol %in% .signature))) %>% 
    
    # format statistics from pairwise contrast
    mutate(markers  = map(markers, 
                          ~ .x %>% 
                            # Group by contrast. Comparisons both ways.
                            pivot_longer(
                              cols = contains("___"),
                              names_to = c("stats", "contrast"), 
                              values_to = ".value", 
                              names_sep="___"
                            ) %>% 
                            
                            # Markers selection within each pair of contrast
                            nest(stat_df = -contrast) %>%
                            
                            # Reshape inside each contrast
                            mutate(stat_df = map(stat_df, ~.x %>% 
                                                   pivot_wider(names_from = stats, 
                                                               values_from = .value))) %>% 
                            
                            unnest(stat_df) )) %>% 
    
    # Add original data data to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>% 
    
    # select only columns needed
    select(-data) %>% 
    
    nest(sil_df = c(!!as.symbol(pre(.level)), markers)) %>% 
    
    mutate(sil_df = map(sil_df, ~ .x %>% unnest(markers))) %>% 
    
    # Calculate silhouette score for PCA plot resulted from the markers selected
    mutate(sil_df = map(sil_df, ~ sil_func(.x, .level, .method))) %>% 
    
    mutate(sil = map_dbl(sil_df, ~ .x$sil)) %>% 
    mutate(real_size = map_int(sil_df, ~ .x$real_size))
  
  # nest(data = - ancestor)
}

# below are real functions ===================================

# pairwise selection 1 marker at a time
single_marker_pw_select <- 
  
  function(.contrast, .level, .discard_num, .method="PCA") {
    
    # initialize variables
    
    contrast_copy <- .contrast %>% 
      mutate(ancestor = !!as.symbol(pre(.level)))
    
    # initialise a signature tibble to store signature markers for each cell type in each iteration
    signature <- tibble(
      ancestor = contrast_copy %>% 
        pull(ancestor)
    ) %>% 
      mutate(cumulative_signature = map(ancestor, ~ vector())) %>% 
      mutate(last_silhouette = 0)
    
    # initialise an output tibble containing all results of interest
    summary_tb <- tibble(
      ancestor = character(),
      data = list(),
      new_challengers = list(),
      winner = list(),
      cumulative_signature = list()
    )
    
    # set the base markers
    contrast_pair_tb0 <-
      
      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      contrast_copy %>%
      
      # select top 1 markers from each contrast, output is an unnested tibble
      sig_select(.level, 1) %>%
      
      mutate(ancestor = !!as.symbol(pre(.level))) %>%
      
      nest(data = - ancestor) %>%
      
      mutate(data = map(data, ~ .x %>% 
                          nest(markers = - contrast_pretty) %>% 
                          
                          mutate(new_challengers = map_chr(
                            markers, ~ .x %>% pull(symbol) %>% unique()))
      )) %>% 
      
      mutate(new_challengers = map(data, ~.x %>% pull(new_challengers))) %>% 
      
      mutate(data = map2(data, new_challengers,
                         ~ .x %>% 
                           mutate(challengers_for_silhouette = list(.y))
      )) %>% 
      
      mutate(data = map2(
        data, ancestor,
        ~ .x %>% 
          mutate(silhouette = map2(
            challengers_for_silhouette, .y,
            ~ silhouette_for_markers(.x, .y, .contrast, .level, .method)
          )) %>% 
          unnest(silhouette)
      )) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% 
          mutate(is_greater = map2_lgl(
            sil, ancestor,
            ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
          ))
      )) %>% 
      
      mutate(data = map(
        data,
        ~ .x %>% select(-ancestor)
      )) %>% 
      
      mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
      
      mutate(cumulative_signature = winner)
    
    
    signature <- signature %>%
      mutate(cumulative_signature = map2(
        cumulative_signature, ancestor,
        ~ .x %>% 
          append(with(contrast_pair_tb0, cumulative_signature[ancestor==.y][[1]]))
      )) %>% 
      mutate(last_silhouette = map_dbl(
        ancestor,
        ~ with(contrast_pair_tb0$data[ancestor==.x][[1]], unique(sil))
      ))
    
    summary_tb <- summary_tb %>% 
      bind_rows(contrast_pair_tb0)
    
    # remove base markers from contrast_copy input before further selection
    contrast_copy <- contrast_copy %>%
      mutate(markers = map2(
        markers, ancestor, 
        ~ .x %>%
          filter(!symbol %in% with(signature, cumulative_signature[ancestor==.y][[1]]))
      ))
    # counter for number of genes discarded
    j <- map_int(signature$cumulative_signature, ~ length(.x))
    
    while (any(j < .discard_num) & all(map_int(contrast_copy$markers, ~ nrow(.x)) > 1)) {
      
      contrast_pair_tb <- 
        
        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        contrast_copy %>% 
        
        # select top 1 markers from each contrast, output is an unnested tibble
        sig_select(.level, 1) %>% 
        
        mutate(ancestor = !!as.symbol(pre(.level))) %>% 
        
        nest(data = - ancestor) %>% 
        
        mutate(data = map(data, ~ .x %>% 
                            nest(markers = - contrast_pretty) %>% 
                            mutate(new_challengers = map_chr(
                              markers, 
                              ~.x %>% pull(symbol) %>% unique()
                            ))
        )) %>% 
        
        unnest(data) %>% 
        
        mutate(challengers_for_silhouette = map2(
          new_challengers, ancestor, 
          ~ with(signature, cumulative_signature[ancestor==.y][[1]]) %>% 
            append(.x)
        )) %>% 
        
        mutate(silhouette = map2(
          challengers_for_silhouette, ancestor, 
          ~ silhouette_for_markers(.x, .y, .contrast, .level, .method) %>% 
            select(-ancestor)
        )) %>% 
        
        unnest(silhouette) %>% 
        
        group_by(ancestor) %>% 
        
        arrange(desc(sil), .by_group = TRUE) %>% 
        
        ungroup() %>% 
        
        mutate(is_greater = map2_lgl(
          sil, ancestor, 
          ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
        )) %>% 
        
        nest(data = - ancestor) %>% 
        
        mutate(new_challengers = map(
          data,
          ~ .x %>% 
            pull(new_challengers)
        )) %>% 
        
        mutate(winner = map(
          data,
          ~ if(.x[1, ]$is_greater){
            .x[1, ]$new_challengers
          } else {NA}
        )) %>% 
        
        mutate(cumulative_signature = pmap(
          list(data, winner, ancestor),
          ~ if(!is.na(..2)) {
            with(..1[1, ], challengers_for_silhouette[[1]])
          } else {
            with(signature, cumulative_signature[ancestor==..3][[1]])
          }
        ))
      
      # append the base + 1 markers that result in highest silhouette score
      signature <- signature %>% 
        
        mutate(cumulative_signature = map(
          ancestor,
          ~ with(contrast_pair_tb, cumulative_signature[ancestor==.x][[1]])
        )) %>% 
        
        mutate(last_silhouette = map2_dbl(
          ancestor, last_silhouette,
          ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x]))) {
            with(contrast_pair_tb, data[ancestor==.x][[1]][[1, "sil"]])
          } else {.y}
        ))
      
      summary_tb <- summary_tb %>% 
        bind_rows(
          contrast_pair_tb %>% 
            filter(!is.na(winner))
        )
      
      contrast_copy <- contrast_copy %>% 
        mutate(markers = map2(
          markers, ancestor, 
          ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y]))){
            .x %>% 
              filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))
          } else {
            .x %>% 
              filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]]))
          }
        ))
      
      # number of genes discarded for each node
      j <- j + 
        
        # unsuccessful candidates
        map_int(contrast_pair_tb$new_challengers, ~length(.x)) *
        is.na(contrast_pair_tb$winner) +
        
        # winning candidates
        map_int(contrast_pair_tb$winner, ~length(.x)) *
        !is.na(contrast_pair_tb$winner)
      
      cat("genes discarded for each node: ", j, "\n")
      cat("genes selected for each node: ", map_int(signature$cumulative_signature, ~ length(.x)),  "\n")
      
    }
    
    my.list <- list(
      "contrast_pair_tb" = contrast_pair_tb, 
      "signature" = signature, 
      "summary_tb" = summary_tb, 
      "contrast_copy" = contrast_copy, 
      "j" = j
    )
    
    return(my.list)
  }

# calculates silhouette score for each set of signature (cumulative markers) at a signature size
silhouette_for_markers <-function(.signature, .ancestor, .contrast, .level, .method) {
  
  .contrast %>%
    
    mutate(ancestor = !!as.symbol(pre(.level))) %>% 
    
    filter(ancestor == .ancestor) %>% 
    
    # filter markers that are in the signature
    mutate(markers = map2(markers, ancestor, ~.x %>% 
                            filter(symbol %in% .signature))) %>% 
    
    # format statistics from pairwise contrast
    mutate(markers  = map(markers, 
                          ~ .x %>% 
                            # Group by contrast. Comparisons both ways.
                            pivot_longer(
                              cols = contains("___"),
                              names_to = c("stats", "contrast"), 
                              values_to = ".value", 
                              names_sep="___"
                            ) %>% 
                            
                            # Markers selection within each pair of contrast
                            nest(stat_df = -contrast) %>%
                            
                            # Reshape inside each contrast
                            mutate(stat_df = map(stat_df, ~.x %>% 
                                                   pivot_wider(names_from = stats, 
                                                               values_from = .value))) %>% 
                            
                            unnest(stat_df) )) %>% 
    
    # Add original data data to the markers selected
    mutate(markers = map2(markers, data, ~ left_join(.x, .y))) %>% 
    
    # select only columns needed
    select(-data) %>% 
    
    nest(sil_df = c(!!as.symbol(pre(.level)), markers)) %>% 
    
    mutate(sil_df = map(sil_df, ~ .x %>% unnest(markers))) %>% 
    
    # Calculate silhouette score for PCA plot resulted from the markers selected
    mutate(sil_df = map(sil_df, ~ sil_func(.x, .level, .method))) %>% 
    
    mutate(sil = map_dbl(sil_df, ~ .x$sil)) %>% 
    
    mutate(real_size = map_int(sil_df, ~ .x$real_size))
  
}

# Testing =====

my.list <- contrast_PW_L3 %>% 
  single_marker_pw_select(LEVEL, .discard_num = 376)

saveRDS(pw_markers_L3, "pw_markers_L3.rds")


full_df %>% 
  filter(str_detect(method, 'silhouette')) %>% 
  unnest(data) %>% 
  ggplot(aes(real_size, silhouette, colour = method)) +
  geom_point(size = 0.1) +
  geom_line() +
  xlim(0, 50)+
  facet_wrap(~ ancestor) +
  guides(color = guide_legend(
    title.position = "left",
    nrow = 2,
    byrow = TRUE))+
  theme(legend.position = "bottom",
        legend.title = element_text(size=10),
        legend.title.align = 0.5,
        plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("silhouette selecttion")


tt_non_hierarchy %>% 
  unnest(tt) %>% 
  pluck("data", 1) %>% 
  distinct(cell_type, level_1, level_2, level_3, level_4, level_5) %>% 
  as.data.frame() %>% 
  FromDataFrameTable()

FromDataFrameTable(x, 
                   pathName = paste('root', 
                                    x$level_1, 
                                    x$level_2, 
                                    x$level_3, 
                                    x$level_4, 
                                    x$level_5, 
                                    x$cell_type),
                   
)

x <- plot(cellsig::tree)
saveWidget(x, "tree.html")
savePlot("tree", "png")
webshot("tree.html", "tree.png", expand = c(10, 50, 0, 50), zoom = 5)

single_marker_pw_selection_using_silhouette <- 
  function(.ranked, .discard_number=NULL, .reduction_method="PCA") {
    
    # initialize variables
    
    # ranked_copy is created as a pool of markers for selection, 
    # which continuously decrease with each iterative selection,
    # input .rank is used for calculating silhouette score for the selected markers
    ranked_copy <- ranked_PW_L4
    
    # initialise a signature tibble to store signature markers for each cell type in each iteration
    signature <- ranked_PW_L4 %>% 
      select(level, ancestor) %>% 
      mutate(signature = map(ancestor, ~ vector())) %>% 
      mutate(last_silhouette = 0)
    
    # initialise an output tibble containing all results of interest
    summary_tb <- tibble(
      level = character(),
      ancestor = character(),
      new_challengers = list(),
      winner = list(),
      children = list(),
      signature = list(),
      # reduced_dimensions = list(),
      silhouette = double()
    )
    
    # set the base markers
    contrast_pair_tb0 <- 
      
      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      ranked_PW_L4 %>%
      
      # select top 1 markers from each contrast
      naive_selection(1) %>%
      
      dplyr::rename("new_challengers" = "signature") %>% 
      
      # there might be a bug (although unlikely) when there are only two contrasts and the top ranked gene are in both contrasts
      mutate(winner = map(new_challengers, ~ unique(.x))) %>% 
      
      # mutate(winning_contrast = map(
      #   markers,
      #   ~ .x %>% pull(contrast) %>% unique()
      #   # distinct(contrast, symbol) %>% 
      #   # mutate(contrast = contrast %>% str_extract(".*(?=\\s\\-)")) %>% 
      #   # mutate(contrast_symbol = map2_chr(contrast, symbol, ~ paste(.x, .y, sep = "."))) %>% 
      #   # pull(contrast_symbol)
      # )) %>%
      
      mutate(signature = winner) %>% 
      
      silhouette_function("PCA") %>% 
      
      select(-c(reduced_dimensions, real_size))
    
    
    signature <- signature %>%
      
      # append cumulative markers
      mutate(signature = map2(
        signature, ancestor,
        ~ .x %>% 
          append(with(contrast_pair_tb0, signature[ancestor==.y][[1]]))
      )) %>% 
      
      # append silhouette scores for these markers
      mutate(last_silhouette = map_dbl(
        ancestor,
        ~ with(contrast_pair_tb0, silhouette[ancestor==.x])
      ))
    
    summary_tb <- summary_tb %>% 
      bind_rows(contrast_pair_tb0)
    
    # remove base markers from contrast_copy input before further selection
    ranked_copy <- ranked_copy %>%
      mutate(markers = map2(
        markers, ancestor, 
        ~ .x %>%
          unnest(stat_df) %>% 
          filter(!symbol %in% with(signature, signature[ancestor==.y][[1]])) %>% 
          nest(stat_df = - contrast)
      ))
    
    # counter for number of genes discarded
    j <- map_int(signature$signature, ~ length(.x))
    
    # count the number of iterations
    i <- 0
    while (any(j < .discard_number) &
           # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
           all(map_int(ranked_copy$markers, 
                       # hence the boundary should be the number of satisfactory genes selected
                       ~ .x %>% unnest(stat_df) %>% nrow()) > 0)) {
      
      contrast_pair_tb <- 
        
        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        ranked_copy %>% 
        
        # select top 1 markers from each contrast, ignore the signature output
        naive_selection(1) %>% 
        
        # # pick the one new challenger from each contrast
        # mutate(markers = map(
        #   markers,
        #   ~ .x %>% 
        #     nest(new_challenger = - contrast) %>% 
        #     mutate(new_challenger = map_chr(new_challenger, ~.x %>% distinct(symbol) %>% pull()))
        # )) %>% 
        # unnest(markers) %>% 
        # select(-c(signature, real_size)) %>% 
        select(-c(markers, signature, real_size)) %>% 
        unnest(children) %>% 
        unnest(enriched) %>% 
        dplyr::rename(new_challenger = symbol) %>% 
        
        # append the new challenger from each contrast to the base markers for that ancestor node
        mutate(challengers_for_silhouette = map2(
          new_challenger, ancestor, 
          ~ with(signature, signature[ancestor==.y][[1]]) %>% 
            append(.x)
        )) %>% 
        
        # calculate silhouette score for the challengers from each contrast
        mutate(silhouette = map2_dbl(
          challengers_for_silhouette, ancestor, 
          ~ silhouette_for_markers(ranked_PW_L4, .x, .y, "PCA") %>% 
            pull(silhouette)
        )) %>% 
        
        # arrange silhouette score in a descending manner within each ancestor node
        group_by(ancestor) %>% 
        arrange(desc(silhouette), .by_group = TRUE) %>% 
        ungroup() %>% 
        
        # check if the silhouette score for the challengers is greater than previous silhouette score
        mutate(is_greater = map2_lgl(
          silhouette, ancestor, 
          ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
        )) %>% 
        
        # nest under ancestor node to select markers that is TRUE for is_greater
        nest(data = - c(level, ancestor)) %>% 
        
        # record new_challengers
        mutate(new_challengers = map(data, ~ .x %>% pull(new_challenger))) %>% 
        
        # check if the biggest silhouette score is greater than previous score, if true we have a winner, else no winner
        mutate(winner = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ]$new_challenger
        } else {NA}
        )) %>% 
        
        # record which contrast the winner comes from
        mutate(children = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ] %>% 
            select(contrast, symbol=new_challenger, rank) %>% 
            nest(enriched = -contrast)
        } else {NA}
        )) %>%

        
        # cummulative signature: winner + previously selected
        mutate(signature = pmap(
          list(data, winner, ancestor),
          ~ if(!is.na(..2)) {
            with(..1[1, ], challengers_for_silhouette[[1]])
          } else {
            with(signature, signature[ancestor==..3][[1]])
          }
        )) %>% 
        
        # silhouette score
        mutate(silhouette = map_dbl(data, ~ .x[[1, "silhouette"]]))
      
      
      # append the base + 1 markers that result in highest silhouette score
      signature <- signature %>% 
        
        mutate(signature = map(
          ancestor,
          ~ with(contrast_pair_tb, signature[ancestor==.x][[1]])
        )) %>% 
        
        mutate(last_silhouette = map2_dbl(
          ancestor, last_silhouette,
          ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x]))) {
            with(contrast_pair_tb, silhouette[ancestor==.x])
          } else {.y}
        )) 
      
      # append the winning signatures into the output summary table
      summary_tb <- summary_tb %>% 
        bind_rows(
          contrast_pair_tb %>% 
            filter(!is.na(winner)) %>% 
            # mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>% 
            select(-data)
        )
      
      # remove the signatures and unsuccessful genes from the selection list(ranked_copy)
      ranked_copy <- ranked_copy %>% 
        mutate(markers = map2(
          markers, ancestor, 
          ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y]))){
            .x %>% 
              unnest(stat_df) %>% 
              filter(!symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]])) %>% 
              nest(stat_df = -contrast)
          } else {
            .x %>% 
              unnest(stat_df) %>% 
              filter(symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]])) %>% 
              nest(stat_df = -contrast)
          }
        ))
      
      # number of genes discarded for each node
      j <- j + 
        
        # unsuccessful candidates
        map_int(contrast_pair_tb$new_challengers, ~length(.x)) *
        is.na(contrast_pair_tb$winner) +
        
        # winning candidates
        map_int(contrast_pair_tb$winner, ~length(.x)) *
        !is.na(contrast_pair_tb$winner)
      
      cat("genes discarded for each node: ", j, "\n")
      cat("genes selected for each node: ", map_int(signature$signature, ~ length(.x)),  "\n")
      
      i <- i + 1
      cat("iteration: ", i, "\n")
      
    }
    
    # format output for optimisation
    output <- summary_tb %>% 
      mutate(real_size = map_int(signature, ~ length(.x))) %>% 
      nest(data = - c(level, ancestor))
    
    return(output)
  }


signature.optimised <- signature.unoptimised %>% do_optimisation("curvature")

# for naive optimisation
optimised1 <- xx %>% 
  
  mutate(optimal_size = map_int(data, ~ penalised_silhouette(.x))) %>% 
  
  unnest(data) %>% 
  
  filter(real_size == optimal_size) %>% 
  
  select(level, ancestor, children, signature, silhouette)


optimised <- xx %>% 
  
  curvature_of_kernel_smoothed_trend() %>% 
  
  # mutate(optimal_size = map_int(data, ~ penalised_silhouette(.x))) %>%
  
  unnest(data) %>% 
  
  filter(real_size <= optimal_size) %>% 
  
  nest(children = -c(level, ancestor)) %>% 
  
  mutate(signature = map(children, ~ .x %>% tail(1) %>% pull(signature) %>% unlist() )) %>% 
  
  mutate(silhouette = map_dbl(children, ~ .x %>% tail(1) %>% pull(silhouette) %>% unlist() )) %>% 
  
  mutate(children = map(
    children,
    ~ .x %>% 
      unnest(children) %>% 
      unnest(enriched) %>% 
      distinct(contrast, symbol, rank) %>% 
      nest(enriched = -contrast)
  ))


x <- signature.unoptimised %>% 
  mutate(data = map(data, 
                    ~ .x %>% 
                      mutate(winning_contrast = map(
                        winning_contrast, ~ 
                          if(str_detect(.x, "\\.")){
                            .x %>% 
                              str_extract(".*(?=\\.)") %>% 
                              unlist()
                          }else{.x}
                        ))
                      ))


silhouette_score <- function(.reduced_dimensions, .distance, .level){
  
  .reduced_dimensions %>% 
    
    pull(!!as.symbol(.level)) %>% 
    
    as.factor() %>% 
    
    as.numeric() %>% 
    
    silhouette(.distance) %>% 
    
    summary()
  
}

silhouette_function <- function(.selected, .reduction_method){
  
  .selected %>% 
    
    # reduce dimensions
    mutate(reduced_dimensions = map2(
      markers, level, 
      ~ dimension_reduction(.x, .y, .reduction_method)
    )) %>% 
    
    # calculate distance matrix using PC1 & PC2
    mutate(distance = map(
      reduced_dimensions,
      ~ distance_matrix(.x, .reduction_method)
    )) %>% 
    
    # calculate silhouette score
    mutate(silhouette = pmap(
      list(reduced_dimensions, distance, level),
      ~ silhouette_score(..1, ..2, ..3)
    )) %>% 
    
    # remove unnecessary columns
    select(-c(markers, distance))
  
}

yy <- full_df %>% 
  nest(signature = -method) %>% 
  mutate(signature = map(signature, ~.x %>% pull(signature) %>% unlist() %>% unique())) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score for all signatures combined in each method
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

cibersortx <- readRDS("dev/topInf_scaleFALSE/cibersortx.new.rds")
cibersort_signature <- cibersortx$signature[[1]]

cibersortx <- tibble(method = "cibersortx") %>% 
  mutate(signature = list(cibersort_signature)) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

yy <- yy %>% bind_rows(cibersortx)

yy %>% 
  mutate(cluster.silhouette = map(silhouette, ~ .x$clus.avg.widths)) %>% 
  mutate(avg.silhouette = map_dbl(silhouette, ~ .x$avg.width)) %>% 
  select(-c(reduced_dimensions, silhouette)) %>% 
  unnest(cluster.silhouette) %>% 
  ggplot(aes(x=reorder(method, avg.silhouette), y=cluster.silhouette, colour=method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  theme(axis.text.x = element_blank())



all_methods_silhouette <- full_df %>%
  mutate(data = map(data, ~ if ("cumulative_signature" %in% names(.x))
  {rename(.x, signature=cumulative_signature)}else(.x))) %>% 
  nest(signature = -method) %>% 
  mutate(signature = map(signature, ~ .x %>% do_optimisation("penalised", 0.4))) %>% 
  mutate(signature = map(signature, ~.x %>% pull(signature) %>% unlist() %>% unique())) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score for all signatures combined in each method
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)

cibersortx <- readRDS("topInf_scaleFALSE/cibersortx.new.rds")

all_methods_comparison <- all_methods_silhouette %>% 
  bind_rows(cibersortx) %>% 
  arrange(desc(silhouette))

all_methods_comparison

all_methods_comparison %>% 
  ggplot(aes(reorder(method, silhouette), silhouette, fill = method)) +
  geom_col() +
  geom_text(aes(label = round(silhouette, 3)), vjust = 1.5) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("All methods comparison using silhouette score, penalty_rate=1.2")

signature_hierarchial_methods %>% 
  mutate(signature = map(benchmark, ~ .x$signature %>% unlist() %>% unique()))

evaluation <- function(.signature, .reduction_method, .preprocessed_non_hierarchy){
  
  # modified silhouette_function and silhouette_score for evaluation
  
  silhouette_score <- function(.reduced_dimensions, .distance, .level){
    
    .reduced_dimensions %>% 
      
      pull(!!as.symbol(.level)) %>% 
      
      as.factor() %>% 
      
      as.numeric() %>% 
      
      silhouette(.distance) %>% 
      
      summary()
    
  }
  
  silhouette_function <- function(.selected, .reduction_method){
    
    .selected %>% 
      
      # reduce dimensions
      mutate(reduced_dimensions = map2(
        markers, level, 
        ~ dimension_reduction(.x, .y, .reduction_method)
      )) %>% 
      
      # calculate distance matrix using PC1 & PC2
      mutate(distance = map(
        reduced_dimensions,
        ~ distance_matrix(.x, .reduction_method)
      )) %>% 
      
      # calculate silhouette score
      mutate(silhouette = pmap(
        list(reduced_dimensions, distance, level),
        ~ silhouette_score(..1, ..2, ..3)
      )) %>% 
      
      # remove unnecessary columns
      select(-c(markers, distance))
    
  }
  
  .preprocessed_non_hierarchy %>% 
    unnest(tt) %>% 
    unnest(data) %>% 
    filter(symbol %in% .signature) %>% 
    nest(markers = -c(level, ancestor)) %>% 
    # calculate silhouette score for all signatures combined in each method
    silhouette_function(.reduction_method = .reduction_method) %>% 
    select(reduced_dimensions, silhouette)
}

x <- signature_hierarchial_methods %>% 
  mutate(stream = "mean_contrast_edgR_logFC_naive_penalty") %>%
  select(benchmark, stream) %>% 
  mutate(signature = map(benchmark, ~ .x$signature %>% unlist() %>% unique())) %>% 
  mutate(silhouette = map(signature, ~ evaluation(.x, "PCA", counts_non_hierarchy)))

  unnest(silhouette) %>% 
  select(-reduced_dimensions) %>% 
  mutate(cluster_silhouette = map(silhouette, ~ .x$clus.avg.widths)) %>% 
  mutate(avg_silhouette = map_dbl(silhouette, ~ .x$avg.width)) %>% 
  unnest(cluster_silhouette) %>% 
  
  ggplot(aes(x=reorder(stream, avg_silhouette), y=cluster_silhouette, colour=stream)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.5) +
  theme(axis.text.x = element_blank())

  
cibersortx <- tibble(method = "cibersortx") %>% 
  mutate(signature = list(cibersort_signature)) %>% 
  mutate(silhouette = map(
    signature, 
    ~ tt_non_hierarchy %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      filter(symbol %in% .x) %>% 
      nest(markers = -c(level, ancestor)) %>% 
      # calculate silhouette score
      silhouette_function(METHOD) %>% 
      select(reduced_dimensions, silhouette)
  )) %>% 
  unnest(silhouette)
  
x <- 
tibble(
  contrast = c(mean_contrast, pairwise_contrast, NA),
  contrast_name = c("mean_contrast", "pairwise_contrast", NA), 
  rank = c(rank_edgR_quasi_likelihood, rank_edgR_robust_likelihood_ratio, rank_bayes),
  rank_name = c("edgR", "edgR_robust", "bayes"),
  rank_stat = c("logFC", "PValue", NA),
  selection = c("silhouette", "naive", NA),
  optimisation = c("penalty", "curvature", NA)
  ) %>% 
  tidyr::expand(nesting(contrast, contrast_name), nesting(rank, rank_name), rank_stat, selection, optimisation) %>%
  
  # Drop arguments for some methods
  filter(!( rank_name == "edgR_robust" & rank_stat == "logFC")) %>%
  filter(!(selection == "naive" & optimisation == "curvature")) %>% 
  mutate(rank_stat = map2(
    rank_stat, rank_name,
    ~ if (.y == "bayes"){.x = NULL} else {.x}
  )) %>%
  
  filter(!(is.na(contrast) | is.na(rank_stat) | is.na(selection) | is.na(optimisation)) ) %>% 
  distinct()

# How to call this file
# From bash
# Rscript database_plus_tree_as_input_benchmark_plot_as_output.R
#
# Output -> a PDF file
# Input nothing
#
# What does it do
# Produce a pdf running the whole pipeline, including all methods, and testing with sihuette and deconvolution

naive_selection <- function(.ranked, .k) {
  
  # Args:
  # .ranked: output from do_ranking()
  # .k: the number of genes selected from each cell_type contrast
  x <- 
  ranked_PW_L4 %>% 
    
    # selection markers from each contrast
    mutate(markers = map(
      markers,
      ~ .x %>% 
        mutate(stat_df = map(stat_df, ~ .x %>% dplyr::slice(1: 7)
                             )) %>% 
        unnest(stat_df)
    )) %>% 
    
    # collect from which contrasts signature genes are extracted
    mutate(children = map(markers, ~ .x %>% 
                            select(contrast, symbol, rank) %>% 
                            nest(enriched = - contrast)
                            )) %>% 
    
    # Add original data info to the markers selected, 
    # use inner_join to ensure symbols are present in both markers and data
    mutate(markers = map2(markers, data, ~ inner_join(.x, .y, by="symbol"))) %>%
    
    # remove unnecessary column
    select(-data) %>% 
    
    # collect signature genes selected
    mutate(signature = map(markers, ~ .x$symbol %>% unique())) %>% 
    
    # number of the signature genes
    mutate(real_size = map_int(signature, ~ length(.x)))
  
}

do_naive_selection <- function(.ranked, .kmax, .reduction_method) {
  
  # Args:
  # .ranked: output from do_ranking
  # .kmax: maximum number of markers selected from each cell type contrast
  # .reduction_method: method used to reduce dimensions such as "PCA", "tSNE", "MSA"
  
  tibble(number_of_markers_from_each_contrast = 1: .kmax) %>% 
    
    # select signature and calculate silhouette score 
    mutate(data = map(
      number_of_markers_from_each_contrast,
      ~ naive_selection(.ranked=.ranked, .x) %>% 
        silhouette_function(.reduction_method=.reduction_method)
    )) %>% 
    
    # nest by ancestor nodes/cell types
    unnest(data) %>%
    nest(data = - c(level, ancestor))
  
}



signature_stefano <- 
signature_stefano_unoptimised %>% 
  mutate(method = "unoptimised", .before = 1) %>% 
  mutate(data = map(data, ~ .x %>% 
                      unnest(c(winner, winning_contrast)) %>% 
                      select(winning_contrast, winner, rank) %>% 
                      nest(enriched = -winning_contrast)
                      )) %>% 
  bind_rows(
    signature_stefano_optimised_by_curvature %>% 
      mutate(method = "curvature", .before = 1) %>% 
      select(-c(signature, silhouette))
  ) %>% 
  
  bind_rows(
    signature_stefano_optimised_by_penalty %>% 
      mutate(method = "penalty", .before = 1) %>% 
      select(-c(signature, silhouette))
  )
 
signature_stefano %>% 
  unnest(data) %>% 
  unnest(enriched) %>% 
  pivot_wider(names_from = method, values_from = rank) %>% 
  saveRDS("./dev/signature_stefano.rds", compress = "xz")
  # nest(data = -c(level, ancestor, winning_contrast)) %>% 
  # pluck("data", 2) %>% print(n=30)
  filter((level=="level_1" & unoptimised < 100) | (level=="level_2" & !is.na(curvature)))
  
 f <-  function(x){
    when(x,
         (.) < 0 ~ x+1,
         (.) >= 0 ~ -x+1
    )
  }

 f(3)
 
readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_non_hierarchy.rds") %>% 
 main(.is_hierarchy=FALSE,
      .contrast_method = mean_contrast, .ranking_method = rank_bayes, .bayes = bayes,
      .selection_method = "silhouette", .kmax = 60, .discard_number = 100, .reduction_method = "PCA",
      .optimisation_method= "curvature") %>% 
 saveRDS(file = glue("{output_file}", compress = "xz"))

ranked_NH_bayes <- readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_non_hierarchy.rds") %>% 
  do_ranking(.ranking_method=rank_bayes, .contrast_method=mean_contrast, .rank_stat=NULL, .bayes=bayes)

selected_NH_bayes <- ranked_NH_bayes %>% 
  do_selection(.selection_method = "naive", "PCA", .kmax = 60, .discard_number = NULL)

selected_NH_bayes %>% 
  do_optimisation("curvature")

x <- dir(input_directory) %>%
  `names<-`(dir(input_directory)) %>%
  .[1:5] %>% 
  map_dfr(~ readRDS(glue("{input_directory}{.x}")), .id = "stream") %>% 
  mutate(stream = str_remove(stream, "\\.rds")) %>% 
  nest(signature = - stream) %>% 
  mutate(signature = map(signature, ~ .x$signature %>% unlist() %>% unique())) %>% 
  
  # silhouette evaluation
  mutate(silhouette = map(
    signature, 
    ~ silhouette_evaluation(
      .signature = .x,
      .reduction_method = "PCA",
      .preprocessed_non_hierarchy = counts_non_hierarchy)
  ))

yy <- x %>% 
  mutate(silhouette = map(signature, ~ silhouette_evaluation(.x, "PCA", counts_non_hierarchy)))

silhouette_score <- function(.reduced_dimensions, .distance, .level){
  
  .reduced_dimensions %>% 
    
    pull(!!as.symbol(.level)) %>% 
    
    as.factor() %>% 
    
    as.numeric() %>% 
    
    silhouette(.distance) %>% 
    
    summary()
  
}

silhouette_function <- function(.selected, .reduction_method){
  
  .selected %>% 
    
    # reduce dimensions
    mutate(reduced_dimensions = map2(
      markers, level, 
      ~ dimension_reduction(.x, .y, .reduction_method)
    )) %>% 
    
    # calculate distance matrix using PC1 & PC2
    mutate(distance = map(
      reduced_dimensions,
      ~ distance_matrix(.x, .reduction_method)
    )) %>% 
    
    # calculate silhouette score
    mutate(silhouette = pmap(
      list(reduced_dimensions, distance, level),
      ~ silhouette_score(..1, ..2, ..3)
    )) %>% 
    
    # remove unnecessary columns
    select(-c(markers, distance))
  
}

non_hierarchical_mean_contrast_edgR_PValue_naive_penalty %>% 
  expand_grid(mix100 %>% dplyr::slice(1), .) %>% 
  mutate(deconvolution = map2(
    signature, mix,
    ~ deconvolution_evaluation(
      .signature = .x, 
      .mix= .y, 
      .preprocessed_non_hierarchy = counts_non_hierarchy)
  ))



  
  # Produce the plot from the results
c(sprintf(
  "dev/benchmark_results/benchmark_plot.pdf:\n\tRscript dev/benchmark_code/produce_plot_from_results.R %s %s", 
  result_directory, 
  result_directory)) %>%
  
  # Add SLURM requirements
  purrr::prepend("CATEGORY=yes_no_hierarchy\nMEMORY=80000\nCORES=2\nWALL_TIME=86400") %>% 
  
  write_lines("./dev/benchmark_code/benchmark_plot.makeflow")

source("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/function_jian.R")

non_hierarchical_mean_contrast_edgR_PValue_naive_penalty <- 
  
  counts_non_hierarchy %>% 
  
  do_ranking(.ranking_method=rank_edgR_quasi_likelihood, .contrast_method=mean_contrast, .rank_stat="PValue") %>%  
  
  # Input: data.frame columns_1 <int> | ...
  # Output: 
  do_selection(.selection_method="naive", .reduction_method="PCA", .kmax=60, .discard_number = NULL) %>% 
  
  format_output(.is_complete = TRUE)

non_hierarchical_mean_contrast_edgR_PValue_naive_penalty %>% 
  
  # saveRDS("dev/intermediate_data/non_hierarchical_mean_contrast_edgR_PValue_naive_unoptimised.rds", compress = "xz")
  
  unnest(data) %>% 
  
  ggplot(aes(real_size, silhouette)) +
  
  geom_point() +
  
  geom_line()


x <- non_hierarchical_mean_contrast_edgR_PValue_naive_penalty %>% 
  
  do_optimisation("penalty")

table_of_commands %>%
  
  filter((is_hierarchy == "non_hierarchical" & contrast_name == "pairwise_contrast" & rank_name=="edgR_robust")| 
           (is_hierarchy == "non_hierarchical" & selection == "silhouette")) %>% 
  
  pull(makeflow_command) %>%
  
  # Add SLURM requirements
  purrr::prepend("CATEGORY=yes_no_hierarchy\nMEMORY=120000\nCORES=2\nWALL_TIME=172800") %>% 
  # 
  # mutate(SLURM_command = glue::glue("sbatch Rscript {R_command}")) %>% 
  # pull(SLURM_command) %>%
  write_lines("./dev/benchmark_code/benchmark_non_hierarchical.makeflow")

x <- 
counts_bayes_imputed_hierarchy %>% 
  
  unnest(tt) %>% 
  
  mutate(contrast = map2(data, level, ~ mean_contrast(.x, .y))) %>% 
  
  unnest(contrast) %>% 
  mutate(contrast = str_remove_all(contrast, glue("{level}"))) %>% 
  
  mutate(lower_quantile = map2(
    data, contrast,
    ~ .x %>% 
      filter(cell_type == str_extract(.y, ".*(?=\\s\\-)")) %>% 
      # filter out genes with imputation ratio greater than 0.2
      filter(ratio_imputed_samples < 0.2) %>% 
      select(symbol, lower_quantile='25%')
    # arrange(symbol)
  )) %>% 
  
  mutate(mean_upper_quantile = map2(
    data, contrast,
    ~ {
      # obtain background cell type(s)
      background <- (.y) %>% 
        str_extract("(?<=\\-\\s).*") %>% 
        str_split("\\+") %>% 
        unlist() %>% 
        str_remove_all("(?<=\\/).*|\\W")
      
      (.x) %>% 
        # calculate the mean 75% quantile of each gene over all background cell types
        filter(cell_type %in% background) %>% 
        group_by(symbol) %>% 
        summarise(symbol, mean_upper_quantile = mean(`75%`)) %>% 
        distinct() %>% 
        ungroup()
    }
  )) %>% 
  
  mutate(stat_df = map2(
    lower_quantile, mean_upper_quantile,
    ~ left_join(.x, .y, by= "symbol")
  )) %>% 
  select(-c(lower_quantile, mean_upper_quantile)) %>% 
  
  mutate(stat_df = map(
    stat_df,
    ~ .x %>% 
      mutate(difference = lower_quantile - mean_upper_quantile) %>% 
      arrange(desc(difference))
  ))

x %>% 
  
  nest(markers = -c(level, ancestor, data))

counts_imputed %>% 
  with_groups(c(cell_type, .feature), ~ .x %>% summarise(`50%` = median(count)))


counts_bayes_imputed_hierarchy %>% 
  slice(5) %>% 
  unnest(tt) %>% 
  unnest(data) %>% 
  distinct(level_5, symbol, ratio_imputed_samples) %>% 
  rename("ratio_bayes" = "ratio_imputed_samples") %>% 
  
  left_join(
    counts_imputed_hierarchy %>% 
      slice(5) %>% 
      unnest(tt) %>% 
      unnest(data) %>% 
      distinct(level_5, symbol, ratio_imputed_samples) %>% 
      rename("ratio_non_bayes" = "ratio_imputed_samples"),
    
    by = c("level_5", "symbol")
  ) %>% 
  
  ggplot(aes(ratio_bayes, ratio_non_bayes)) +
  geom_point()

z <- tree %>%
  ToDataFrameTypeColFull(fill=NA) %>% 
  unite("ancestors", contains("level"), sep=".", na.rm = TRUE) %>% 
  mutate(ancestors = map2(
    ancestors, cell_type,
    ~.x %>% 
      paste("Tissue", ., sep=".") %>% 
      str_split("\\.") %>% 
      unlist() %>% 
      .[. != .y] %>% 
      rev()
  )) %>% 
  mutate(cell_type_allowed = t_helper_tree$Get("level") %>% names() %>% list()) %>% 
  mutate(cell_type_to_be = pmap_chr(
    list(ancestors, cell_type, cell_type_allowed),
    ~ if (! ..2 %in% ..3){
      ..1[which(..1 %in% ..3)[1]]
    } else {..2}
  )) %>% 
  select(cell_type, cell_type_to_be)

counts_t_helper_tree %>% 
  left_join(z, by = "cell_type") %>% 
  select(-cell_type) %>% 
  rename("cell_type" = "cell_type_to_be") %>% 
  tree_and_signatures_to_database(t_helper_tree, ., sample, cell_type, symbol, count)

counts_t_helper_tree %>% 
  # slice(1:100) %>% 
  mutate(cell_type = map_chr(cell_type, ~ with(z, cell_type_to_be[which(cell_type == .x)]) ))

# without enquo()
preprocess <- function(.transcript, .level) {
  # this preproces function ranged data in hierarchy(or non_hierarchy) and
  # calculates the imputation ratio for genes in each hierarchy
  
  # load data
  .transcript %>%
    
    dplyr::rename("symbol" = "feature") %>% 
    
    # tidybulk(sample, symbol, count_scaled) %>% for imputed counts data
    tidybulk(.sample = sample, .transcript = symbol, .abundance = count_scaled) %>%
    
    # filter for cells at the level of interest. .level == level_1
    filter(!is.na(!!as.symbol(.level))) %>%
    
    # calculate the ratio of imputation for genes in a cell type
    nest(data = -c(symbol, !!as.symbol(.level))) %>%
    
    # for a cell type some samples may miss genes in other samples: so for the same cell type genes may have different number of samples
    mutate(n_samples_per_gene = map_int(
      data,
      ~ .x %>%
        distinct(sample) %>%
        nrow)) %>%
    
    unnest(data) %>%
    nest(data = -!!as.symbol(.level)) %>%
    
    mutate(n_samples= map_int(
      data,
      ~ .x %>%
        distinct(sample) %>%
        nrow
    )) %>%
    
    unnest(data) %>% 
    
    mutate(ratio_imputed_samples = 1 - n_samples_per_gene / n_samples) %>% 
    
    # nest by ancestor
    nest(data = - !!as.symbol(pre(.level)))
  
}

counts_imputed_t_helper_tree %>% 
  dplyr::rename(symbol = feature) %>% 
  main(.is_hierarchy=TRUE, 
       .sample = sample, 
       .symbol = symbol,
       .contrast_method=mean_contrast, 
       .ranking_method=rank_bayes, 
       .rank_stat=NULL,
       .bayes=counts_bayes_imputed, 
       .selection_method="silhouette", .kmax=60, .discard_number=2000, .reduction_method = "PCA", .dims=2,
       .optimisation_method="curvature", .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
       .is_complete = TRUE) 

counts_imputed_t_helper_tree %>% 
  dplyr::rename(symbol = feature) %>%
  do_hierarchy()

x %>% 
  # slice(9) %>% 
  mutate(level.copy = level) %>% 
  nest(data = -level.copy) %>% 
  
  mutate(data = map(
    data,
    ~ .x %>% 
      
      do_optimisation(.optimisation_method = "curvature", .symbol=symbol)
  )) %>% 
  
  unnest(data) %>% 
  select(-level.copy)

hierarchical pairwise_contrast bayes _ naive penalty 10 1

x <- 
counts_imputed_hierarchy %>% 
  do_ranking(.sample=sample, .symbol=symbol, .ranking_method = rank_bayes, 
             .contrast_method = pairwise_contrast, .bayes = counts_bayes_imputed_hierarchy)

y <- x %>% 
  mutate(level.copy = level) %>% 
  nest(data = -level.copy) %>% 
  
  mutate(data = map(
    data,
    ~ .x %>% 
      
      do_selection(.sample = sample, 
                   .symbol = symbol,
                   .selection_method="naive", 
                   .reduction_method="PCA", 
                   .kmax=10,
                   .dims=10)
    
  )) %>% 
  
  mutate(data = map(
    data,
    ~ .x %>% 
      
      do_optimisation(.optimisation_method = .optimisation_method, .symbol=!!.symbol) %>%
      
      format_output(.is_complete = .is_complete)
    
  )) %>%
  
  unnest(data) %>% 
  select(-level.copy)

do_naive_selection <- function(.ranked, .sample, .symbol, .kmax, .reduction_method, .dims=2) {
  
  # Args:
  # .ranked: output from do_ranking
  # .kmax: maximum number of markers selected from each cell type contrast
  # .reduction_method: method used to reduce dimensions such as "PCA", "tSNE", "MSA", "UMAP"
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  
  .ranked %>% 
    
    # calculate the minimum number of genes need to be selected for feasible dimension reduction
    mutate(k0 = map_int(markers, ~ min_markers_per_contrast(.x, .dims=.dims))) %>% 
    
    # expand the number of markers selected to .kmax
    mutate(n_markers_from_each_contrast = map(k0, ~ .x: .kmax)) %>% 
    
    # clean unnecessary data
    select(-k0) %>%
    
    # expand the nuber of markers need to be selected
    unnest(n_markers_from_each_contrast) %>% 
    nest(data = -n_markers_from_each_contrast) %>% 
    
    # do naive selection according to the number of markers from each contrast and calculate silhouette score for the signature
    mutate(data = map2(
      data, n_markers_from_each_contrast,
      
      ~ naive_selection(.ranked=.x, .k=.y, .symbol=!!.symbol) %>%
        
        silhouette_function(.sample=!!.sample,
                            .symbol=!!.symbol,
                            .reduction_method=.reduction_method,
                            .dims=.dims) %>%
        
        # optional: remove reduced_dimensions matrix to save memory, if kept it can be used to plot PCA
        select(-reduced_dimensions)
    )) %>%
    
    # nest by ancestor nodes/cell types
    unnest(data) %>%
    nest(data = - c(level, ancestor))
  
}

min_markers_per_contrast <- function(.markers, .dims){
  k = 1L
  n_unique_markers = .markers %>% 
    mutate(top_k = map(stat_df, ~ .x %>% slice(1:k) %>% pull(symbol))) %>% 
    pull(top_k) %>% 
    unlist %>% 
    n_distinct
  
  while (n_unique_markers < .dims) {
    k = k + 1L
    n_unique_markers = .markers %>% 
      mutate(top_k = map(stat_df, ~ .x %>% slice(1:k) %>% pull(symbol))) %>% 
      pull(top_k) %>% 
      unlist %>% 
      n_distinct
  }
  
  return(k)
}

hierarchical pairwise_contrast edgR_robust PValue silhouette curvature 10 0

hierarchical_mean_contrast_edgR_PValue_silhouette_curvature_10.rds

non_hierarchical_mean_contrast_bayes___silhouette_curvature_2.rds

hierarchical_mean_contrast_edgR_PValue_silhouette_curvature_2.rds

non_hierarchical_pairwise_contrast_edgR_logFC_naive_penalty_4.rds

xx <- counts_imputed_non_hierarchy %>% 
  # counts_imputed_non_hierarchy %>% 
  do_ranking(.sample=sample, .symbol=symbol, .ranking_method = rank_edgR_quasi_likelihood, 
             .contrast_method = pairwise_contrast, .rank_stat = "logFC")

y <- xx %>% 
  do_selection(.sample=sample, 
               .symbol=symbol,
               .selection_method="naive", 
               .reduction_method="PCA", 
               .kmax=60,
               .dims=4)

y <- xx %>% 
  mutate(level.copy = level) %>% 
  nest(data = -level.copy) %>% 
  
  mutate(data = map(
    data,
    ~ .x %>% 
      
      do_selection(.sample = sample, 
                   .symbol = symbol,
                   .selection_method="silhouette", 
                   .reduction_method="PCA", 
                   .discard_number = 2000,
                   .dims=2)
    
  ))

z <- y %>% 
  
  # slice(3) %>% 
  
  mutate(data = map(
    data,
    ~ .x %>% 
      
      do_optimisation(.optimisation_method = "curvature", .symbol=symbol) %>%
      
      format_output(.is_complete = TRUE)
    
  )) %>%
  
  unnest(data) %>% 
  select(-level.copy)


do_scaling <- function(.data_tree, .sample, .symbol, .count, .cell_type) {
  
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .count = enquo(.count)
  .cell_type = enquo(.cell_type)
  
  .data_tree %>% 
    
    nest(data = -c(level_1, !!.symbol)) %>%
    add_count(!!.symbol) %>%
    filter(n==n_distinct(.data_tree$level_1)) %>%
    select(-n) %>%
    unnest(data) %>%
    
    # Convert to SE
    as_SummarizedExperiment(!!.sample, !!.symbol, !!.count)  %>%
    
    # Scale with first degree imputation. 
    # This because there are no common genes to all samples
    impute_missing_abundance(~ !!.cell_type, .sample=sample, .transcript = symbol, .abundance = count) %>%
    identify_abundant() %>%
    scale_abundance() %>%
    filter(!.imputed) %>% 
    
    select(-.imputed)  %>%
    
    # Just needed for the old version
    # select(-one_of("exposure_rate")) %>%
    
    # Calculate exposure for Bayes model
    mutate(exposure_rate = -log(multiplier))
}

do_imputation <- function(.scaled_counts, .sample, .symbol, .cell_type){
  .sample = enquo(.sample)
  .symbol = enquo(.symbol)
  .cell_type = enquo(.cell_type)
  
  .scaled_counts %>% 
    # Convert to SE
    # as_SummarizedExperiment(.sample, .feature, count) %>%
    as_SummarizedExperiment(!!.sample, !!.symbol, count_scaled) %>%
    
    # Hierarchical imputation. Suffix = "" equated to overwrite counts
    impute_missing_abundance(~ !!.cell_type, suffix="") %>%
    
    {
      for(level in (.) %>% colnames %>% str_extract("level\\_\\d") %>% .[!is.na(.)]) {
        (.) <- (.) %>% 
          impute_missing_abundance(~ !!as.symbol(level), suffix="")
      }
      
    } %>% 
    
    # Convert back to tibble
    as_tibble() %>%
    
    mutate(.imputed = if_any(contains("imputed"), ~ .x != 0)) %>% 
    
    select(-matches("imputed\\.\\d"))
}

tibble(command = sprintf("error_file:\n\tRscript dev/error_files/SLURM_cellsig_error.R > dev/AAA_err.stderr  2>&1")) %>% 
  
  pull(command) %>% 
  purrr::prepend("CATEGORY=yes_no_hierarchy\nMEMORY=10000\nCORES=2\nWALL_TIME=86400") %>% 
  
  write_lines("dev/error_files/slurm_cellsig_error.makeflow")

tree = cellsig::tree

counts_scaled_old_tree <- counts %>% 
  
  adapt_tree(tree) %>%
  
  # Parse into hierarchical dataset
  tree_and_signatures_to_database(tree, ., sample, cell_type, symbol, count)  %>%
  
  # Remove redundant samples
  remove_redundancy(sample, symbol, count, correlation_threshold = 0.999, top = 500, method = "correlation") %>%
  droplevels() %>% 
  
  # Eliminate suspicious samples
  filter(!grepl("GSM3722278|GSM3722276|GSM3722277", sample))

counts_scaled_old_tree <- counts_scaled_old_tree %>% 
  
  nest(data = -c(level_1, symbol)) %>%
  add_count(symbol) %>%
  filter(n==n_distinct(.$level_1)) %>%
  select(-n) %>%
  unnest(data)

counts_scaled_old_tree %>% 
  
  # impute_missing_abundance(~ cell_type, .sample=sample, .transcript = symbol, .abundance = count)
  
  # Convert to SE
  # as_SummarizedExperiment(sample, symbol, count) %>% 
  
  do_scaling(.sample = sample, .symbol=symbol, .count = count, .cell_type = cell_type)

counts_imputed_old_tree <- counts_scaled_old_tree %>% 
  
  do_imputation(.sample = sample, .symbol=symbol, .cell_type = cell_type)

x <- counts_scaled_cellsig_tree %>% 
  
  # Convert to SE
  as_SummarizedExperiment(sample, feature, c(count, count_scaled))

y0 <- x %>% 
  
  # Hierarchical imputation. Suffix = "" equated to overwrite counts
  impute_missing_abundance(~ cell_type, .abundance = c(count, count_scaled))
y1 <- y0 %>% 
  impute_missing_abundance(~ level_5, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_4, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_3, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_2, .abundance = c(count, count_scaled)) %>%
  impute_missing_abundance(~ level_1, .abundance = c(count, count_scaled)) %>% 
  
levelwise_imputation <- function(.scaled_data, level){
  if (level == "level_6") {
    return(.scaled_data %>% 
             impute_missing_abundance(.formula = as.formula(sprintf("~ %s", pre(level))),
                                      .abundance = c(count, count_scaled))
    )
  }else{
    
    level <- level %>% 
      str_split("_") %>%
      {as.numeric(.[[1]][2]) + 1} %>%
      paste("level", ., sep = "_")
      
    return(
        levelwise_imputation(.scaled_data, level)
      )
    }
}
  
  
  # Convert back to tibble
  as_tibble() %>%
  
  mutate(.imputed = if_any(contains("imputed"), ~ .x != 0)) %>% 
  
  select(-matches("imputed\\.\\d")) %>% 
  
  # Merge the imputed column
  mutate(.imputed = .imputed | .imputed.1 | .imputed.2 | .imputed.3 |.imputed.4 |.imputed.5  ) %>%
  select(-c( .imputed.1 , .imputed.2 , .imputed.3 ,.imputed.4 ,.imputed.5  )) %>%
  
  # Save
  saveRDS("dev/intermediate_data/counts_imputed.rds", compress = "xz")

## create reference file


debugonce(produce_cibersortx_bulk_rnaseq_input)

counts_imputed %>% 
  rename(symbol = feature) %>% 
  produce_cibersortx_bulk_rnaseq_input(.transcript=symbol, .sample=sample, .cell_type=cell_type, .count=count_scaled, 
                                         .dir="dev/jian_R_files/cibersortx", .suffix="_zijie")


toy_data <- counts_imputed %>%  
  select(feature, sample, count_scaled, cell_type) %>% 
  filter(cell_type %in% c("epithelial", "b_cell", "t_cell")) %>% 
  distinct(cell_type, sample) %>% 
  nest(sample = -cell_type) %>% 
  mutate(sample = map(sample, 
                      ~ .x %>% slice(
                        sample(1:nrow(.x), 10, replace = FALSE)
                      ))) %>% 
  unnest(sample) %>% 
  left_join(counts_imputed %>% 
              select(feature, sample, count_scaled, cell_type), 
            by = c("cell_type", "sample"))

saveRDS(toy_data, "dev/intermediate_data/toy_data.rds", compress = "xz")

tree <- read_yaml("dev/tree.yaml") %>% as.Node
Prune(tree$epithelial, function(node) node$name == "epithelial")

Prune(tree, function(node) node$level <= 3)
treeClone = Clone(tree, pruneFun = function(node) node$level <= 3)
LEVEL <- "level_1"
L = LEVEL %>% str_split("_") %>% {as.numeric(.[[1]][2])+1}  

treeClone = Clone(tree, pruneFun = function(node) node$level <= L)

get_leaf_nodes_at_a_level <- function(.tree, .level){
  
  L  = .level %>% str_split("_") %>% {as.numeric(.[[1]][2])+1}
  
  Clone(.tree, pruneFun = function(node) node$level <= L) %>% 
    as.phylo %>% 
    .$tip.label
  
}

tree %>% 
  get_leaf_nodes_at_a_level("level_5")


# deconvolve cellularity problem debug
x <- deconvolve_cellularity(
  .data = mix100 %>% pluck("mix", 1),
  .sample = replicate,
  .transcript = symbol, # the column in the mixture is called symbol not provided by the pipeline
  .abundance = count_mix,
  reference = reference,
  method = "llsr",
  prefix = "llsr_",
  action = "get",
  intercept=FALSE)

plot_data = 
tibble(levels = if(TRUE){"root"}else{sprintf("level_%s", 1:5)}
       
       ) %>%
  mutate(data = map(levels, ~ non_hierarchical_pairwise_contrast_bayes___naive_penalty_10 %>%
                      nest(data = -level) %>% 
                      mutate(markers = map(data, ~ .x$signature %>% unlist %>% unique)) %>% 
                      mutate(level = as.ordered(level)) %>% 
                      filter(level <= .x)
  )) %>% 
  mutate(markers = map(data, ~ .x$markers %>% unlist %>% unique)) %>% 
  mutate(leaf_nodes = map(levels, ~ get_leaf_nodes_at_a_level(new_tree, .x))) %>% 
  mutate(reduced_dimensions = map2(
    markers, leaf_nodes, 
    ~ expression %>% 
      filter(cell_type2 %in% .y) %>% 
      filter(symbol %in% .x) %>% 
      pivot_wider(names_from = "level", values_from = "cell_type2") %>% 
      reduce_dimensions(.element = sample,
                        .feature = symbol,
                        .abundance = count_scaled,
                        action = "get",
                        method = "tSNE",
                        .dims = 2,
                        log_transform = TRUE,
                        top = Inf,
                        scale = FALSE,
                        check_duplicates = FALSE) %>% 
      unite(cell_type2, contains("level"), na.rm = TRUE) %>% 
      # because action is "get" no symbols all samples should be distinct
      select(sample, cell_type2, contains(str_sub("tSNE", end = -2L)))
  )) %>% 
  mutate(marker_text = map(data, ~ .x %>% get_target_cell_markers())) %>% 
  mutate(reduced_dimensions = map2(
    reduced_dimensions, marker_text,
    ~ left_join(.x, .y, by = c("cell_type2" = "target"))
  )) %>% 
  mutate(data = map2(
    data, leaf_nodes,
    ~ .x %>% 
      unnest(data) %>% 
      unnest(children) %>% 
      mutate(contrast = str_extract(contrast, ".*(?=\\s\\-)")) %>% 
      rename(target = contrast) %>% 
      select(target, enriched) %>% 
      filter(target %in% .y) %>% 
      unnest(enriched) %>% 
      distinct(target, symbol, rank) %>% 
      group_by(target) %>% 
      arrange(rank, .by_group = TRUE) %>% 
      ungroup
  )) %>% 
  select(levels, data, leaf_nodes, reduced_dimensions)



c("\tRscript dev/jian_R_files/for_zijie.R dev/intermediate_data/zijie_bulk_signature.rds") %>% 
  
  purrr::prepend("CATEGORY=yes_no_hierarchy\nMEMORY=30000\nCORES=2\nWALL_TIME=86400") %>%
  
  write_lines("dev/benchmark_code/zijie_bulk_signaure.makeflow")

counts %>% 
  
  adapt_tree(.tree = new_tree) %>%

  tree_and_signatures_to_database(tree=new_tree, ., .sample=sample, .cell_type=cell_type,
                                 .symbol=symbol, .count=count) %>%

  # Remove redundant samples
  remove_redundancy(.element=sample, .feature=symbol, .abundance=count, correlation_threshold = 0.999, top = 500, method = "correlation") %>%
  droplevels() %>%
  
  # Eliminate suspicious samples
  filter(!grepl("GSM3722278|GSM3722276|GSM3722277", sample)) %>%

  do_scaling(.sample = sample, .symbol= symbol , .count= count, .cell_type= cell_type) %>%
  
  saveRDS("dev/intermediate_data/counts_scaled.rds", compress = "xz")

counts_imputed_hierarchy %>% 
  
  main(.sample=sample, .symbol=symbol, .count = NULL, .cell_type = cell_type,
       .is_hierarchy=TRUE, 
       .contrast_method=pairwise_contrast, 
       .ranking_method=rank_bayes, 
       .rank_stat=NULL, 
       .bayes=counts_bayes_imputed_hierarchy, 
       .tree = NULL,
       .selection_method="silhouette", .kmax=60, .discard_number=2000, .reduction_method = "tSNE",
       .dims=2,
       .optimisation_method = "penalty", .penalty_rate = 0.2, .kernel = "normal", .bandwidth = 0.05, .gridsize = 100,
       .is_complete = TRUE)


do_silhouette_selection <-
  function(.ranked, .sample, .symbol, .discard_number, .reduction_method, .dims=2) {
    
    .sample = enquo(.sample)
    .symbol = enquo(.symbol)
    
    # initialize variables
    
    # ranked_copy is created as a pool of markers for selection,
    # which continuously decrease with each iterative selection,
    # input .rank is used for calculating silhouette score for the selected markers
    ranked_copy <- .ranked
    
    # initialise a signature tibble to store signature markers for each cell type in each iteration
    signature <- .ranked %>%
      select(level, ancestor) %>%
      mutate(signature = map(ancestor, ~ vector())) %>%
      mutate(last_silhouette = 0)
    
    # initialise an output tibble containing all results of interest
    summary_tb <- tibble(
      level = character(),
      ancestor = character(),
      new_challengers = list(),
      winner = list(),
      children = list(),
      signature = list(),
      # reduced_dimensions = list(),
      silhouette = double()
    )
    
    # set the base markers
    contrast_pair_tb0 <-
      
      # contrast_copy contains all the statistics of all cell_type contrasts for each gene
      .ranked %>%
      
      # calculate the minimum number of genes need to be selected for feasible dimension reduction
      mutate(k0 = map_int(markers, ~ min_markers_per_contrast(.x, .dims=.dims, .symbol = !!.symbol))) %>% 
      
      # nest by k0 because naive selection takes the whole ranked data frame
      nest(data = -k0) %>% 
      
      # select the minimum number of base markers for each ancestor node
      mutate(data = map2(data, k0, ~ naive_selection(.ranked=.x, .k=.y, .symbol=!!.symbol))) %>% 
      
      unnest(data) %>% 
      
      # clean up data
      select(-k0) %>% 
      
      dplyr::rename(new_challengers = signature) %>%
      
      mutate(winner = map(new_challengers, ~ unique(.x))) %>%
      
      mutate(signature = winner) %>%
      
      silhouette_function(.sample=!!.sample, .symbol=!!.symbol, .reduction_method=.reduction_method, .dims=.dims) %>%
      
      select(-c(reduced_dimensions, real_size))
    
    
    signature <- signature %>%
      
      # append cumulative markers
      mutate(signature = map2(
        signature, ancestor,
        ~ .x %>%
          append(with(contrast_pair_tb0, signature[ancestor==.y][[1]]))
      )) %>%
      
      # append silhouette scores for these markers
      mutate(last_silhouette = map_dbl(
        ancestor,
        ~ with(contrast_pair_tb0, silhouette[ancestor==.x])
      ))
    
    summary_tb <- summary_tb %>%
      bind_rows(contrast_pair_tb0)
    
    # remove base markers from contrast_copy input before further selection
    ranked_copy <- ranked_copy %>%
      mutate(markers = map2(
        markers, ancestor,
        ~ .x %>%
          unnest(stat_df) %>%
          filter(!(!!.symbol %in% with(signature, signature[ancestor==.y][[1]]))) %>%
          nest(stat_df = - contrast)
      ))
    
    # counter for number of genes discarded
    j <- map_int(signature$signature, ~ length(.x))
    
    # count the number of iterations
    i <- 0L
    while (any(j < .discard_number) &
           # markers contains genes including many that do not satisfy logFC > 2 & FDR < 0.05 & logCPM > mean(logCPM)
           all(map_int(ranked_copy$markers,
                       # hence the boundary should be the number of satisfactory genes selected
                       ~ .x %>% unnest(stat_df) %>% nrow()) > 0)) {
      
      cat("step_a: ", ranked_copy$level, "\n")
      
      contrast_pair_tb <-
        
        # contrast_PW_L1 contains all the statistics of all cell_type contrasts for each gene
        ranked_copy %>%
        
        # select top 1 markers from each contrast, ignore the signature output
        naive_selection(.k=1, .symbol=!!.symbol) %>%
        
        # pick the one new challenger from each contrast
        select(-c(markers, signature, real_size)) %>%
        unnest(children) %>%
        unnest(enriched) %>%
        dplyr::rename(new_challenger = !!.symbol) %>%
        
        # append the new challenger from each contrast to the base markers for that ancestor node
        mutate(challengers_for_silhouette = map2(
          new_challenger, ancestor,
          ~ with(signature, signature[ancestor==.y][[1]]) %>%
            append(.x)
        )) %>%
        
        # calculate silhouette score for the challengers from each contrast
        mutate(silhouette = map2_dbl(
          challengers_for_silhouette, ancestor,
          ~ silhouette_for_markers(.ranked=.ranked, .sample=!!.sample, .symbol=!!.symbol,
                                   .signature=.x, .ancestor=.y, 
                                   .reduction_method=.reduction_method, .dims=.dims) %>%
            pull(silhouette)
        )) %>%
        
        # arrange silhouette score in a descending manner within each ancestor node
        group_by(ancestor) %>%
        arrange(desc(silhouette), .by_group = TRUE) %>%
        ungroup() %>%
        
        # check if the silhouette score for the challengers is greater than previous silhouette score
        mutate(is_greater = map2_lgl(
          silhouette, ancestor,
          ~ if(.x > with(signature, last_silhouette[ancestor==.y])){TRUE}else{FALSE}
        )) %>%
        
        # nest under ancestor node to select markers that is TRUE for is_greater
        nest(data = - c(level, ancestor)) %>%
        
        # record new_challengers
        mutate(new_challengers = map(data, ~ .x %>% pull(new_challenger))) %>%
        
        # check if the biggest silhouette score is greater than previous score, if true we have a winner, else no winner
        mutate(winner = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ]$new_challenger
        } else {NA}
        )) %>%
        
        # record which contrast the winner comes from
        mutate(children = map(data, ~ if(.x[1, ]$is_greater){
          .x[1, ] %>%
            select(contrast, !!.symbol := new_challenger, rank) %>%
            nest(enriched = -contrast)
        } else {NA}
        )) %>%
        
        
        # cummulative signature: winner + previously selected
        mutate(signature = pmap(
          list(data, winner, ancestor),
          ~ if(!is.na(..2)) {
            with(..1[1, ], challengers_for_silhouette[[1]])
          } else {
            with(signature, signature[ancestor==..3][[1]])
          }
        )) %>%
        
        # silhouette score
        mutate(silhouette = map_dbl(data, ~ .x[[1, "silhouette"]]))
      
      cat("step_b: ", ranked_copy$level, "\n")
      
      
      # append the base + 1 markers that result in highest silhouette score
      signature <- signature %>%
        
        mutate(signature = map(
          ancestor,
          ~ with(contrast_pair_tb, signature[ancestor==.x][[1]])
        )) %>%
        
        mutate(last_silhouette = map2_dbl(
          ancestor, last_silhouette,
          ~ if(!is.na(with(contrast_pair_tb, winner[ancestor==.x]))) {
            with(contrast_pair_tb, silhouette[ancestor==.x])
          } else {.y}
        ))
      
      cat("step_c: ", ranked_copy$level, "\n")
      
      # append the winning signatures into the output summary table
      summary_tb <- summary_tb %>%
        bind_rows(
          contrast_pair_tb %>%
            filter(!is.na(winner)) %>%
            # mutate(reduced_dimensions = map(data, ~ .x$reduced_dimensions[[1]])) %>%
            select(-data)
        )
      
      cat("step_d: ", ranked_copy$level, "\n")
      
      # remove the signatures and unsuccessful genes from the selection list(ranked_copy)
      ranked_copy <- ranked_copy %>%
        mutate(markers = map2(
          markers, ancestor,
          ~ if(is.na(with(contrast_pair_tb, winner[ancestor==.y]))){
            .x %>%
              unnest(stat_df) %>%
              filter(!(!!.symbol %in% with(contrast_pair_tb, new_challengers[ancestor==.y][[1]]))) %>%
              nest(stat_df = -contrast)
          } else {
            .x %>%
              unnest(stat_df) %>%
              filter(!!.symbol != with(contrast_pair_tb, winner[ancestor==.y][[1]])) %>%
              nest(stat_df = -contrast)
          }
        ))
      
      cat("step_e: ", ranked_copy$level, "\n")
      
      # number of genes discarded for each node
      j <- j +
        
        # unsuccessful candidates
        map_int(contrast_pair_tb$new_challengers, ~length(.x)) *
        is.na(contrast_pair_tb$winner) +
        
        # winning candidates
        map_int(contrast_pair_tb$winner, ~length(.x)) *
        !is.na(contrast_pair_tb$winner)
      
      cat("genes discarded for each node: ", j, "\n")
      cat("genes selected for each node: ", map_int(signature$signature, ~ length(.x)),  "\n")
      
      cat("step_f: ", ranked_copy$level, "\n")
      
      i <- i + 1L
      cat("iteration: ", i, "\n")
      
      cat("step_g: ", ranked_copy$level, "\n")
      
    }
    
    # format output for optimisation
    # output <- summary_tb %>%
    #   mutate(real_size = map_int(signature, ~ length(.x))) %>%
    #   nest(data = - c(level, ancestor))
    
    return(ranked_copy)
  }

packages <- c("library(yaml)", "library(tidytext)",
"library(data.tree)",
"library(tidytree)",
"library(ape)",
"library(glue)",
"library(rlang)",
"library(factoextra)",
"library(stringr)",
"library(scales)",
"library(KernSmooth)",
"library(splus2R)",
"library(data.tree)",
"library(cluster)",
"library(tidyverse)",
"library(tidybulk)",
"library(cellsig)",
"library(patchwork)",
"library(tidySummarizedExperiment)",
"library(networkD3)",
"library(htmlwidgets)",
"library(webshot)",
"library(kableExtra)",
"library(shiny)",
"library(shinyjs)",
"library(DT)",
"library(plotly)")

packages %>% 
  str_extract("(?<=\\().*(?=\\))") %>% 
  unique() %>% 
  str_sort()

org <- data.frame(
  Manager = c(
    NA, "Ana", "Ana", "Bill", "Bill", "Bill", "Claudette", "Claudette", "Danny",
    "Fred", "Fred", "Grace", "Larry", "Larry", "Nicholas", "Nicholas"
  ),
  Employee = c(
    "Ana", "Bill", "Larry", "Claudette", "Danny", "Erika", "Fred", "Grace",
    "Henri", "Ida", "Joaquin", "Kate", "Mindy", "Nicholas", "Odette", "Peter"
  ),
  Title = c(
    "President", "VP Operations", "VP Finance", "Director", "Director", "Scientist",
    "Manager", "Manager", "Jr Scientist", "Operator", "Operator", "Associate",
    "Analyst", "Director", "Accountant", "Accountant"
  )
)

cellsig::tree


new_tree <- read_yaml("dev/jian_R_files/new_tree.yaml") %>% as.Node

tally <- expression %>% 
  mutate(across(where(is.factor), as.character)) %>% 
  group_by(level, cell_type2) %>% 
  summarise(n_sample = n_distinct(sample), n_study = n_distinct(database), .groups = "drop") %>% 
  drop_na() %>% 
  distinct()

new_tree %>% 
  ToDataFrameNetwork() %>%
  left_join(tally, by = c("to" = "cell_type2")) %>%
  # create the root node
  rbind(c(from=NA, to="Tissue", level="level_0", n_sample=1053L, n_study=41L), .) %>%
  # rename level because it's a reserved name for collapseTreeNetwork
  rename(levels = level) %>% 
  mutate(across(starts_with("n_"), as.integer)) %>% 
  mutate(color = colorspace::rainbow_hcl(41)) %>% 
  mutate(tooltip = paste0(to, "<br>",
                          levels,
                          "<br>n_sample: ", n_sample,
                          "<br>n_study: ", n_study)
         ) %>% 
  collapsibleTreeNetwork(attribute = "n_sample", 
                         aggFun = identity,
                         nodeSize = "n_sample",
                         fill = "color",
                         tooltipHtml = "tooltip",
                         collapsed=TRUE)


job::job({saveRDS(pca_df, "dev/intermediate_data/pca_df.rds", compress = "xz")})

# check object size
object.size(pca_with_counts )

object.size(pca_with_counts ) %>% `/` (1e6)

object.size(pca_with_counts_se)

object.size(pca_with_counts_se) %>% `/` (1e6)

expression <-
  # counts_imputed %>%
  readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_imputed.rds") %>%
  rename(symbol = feature) %>%
  pivot_longer(contains("level"), names_to = "level", values_to="cell_type2")

x = expression %>%
  select(symbol, sample, count_scaled, cell_type, level, cell_type2) %>%
  mutate(count_scaled = as.integer(count_scaled)) %>%
  # tidybulk(sample, symbol, count_scaled) %>%
  tidybulk::as_SummarizedExperiment(sample, symbol, count_scaled)

expression %>% 
  assay()

input_address = "/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig"

sprintf("%s/output.txt:\n\tRscript %s/test.R %s/output.R", input_address, input_address, input_address) %>% cat()
  prepend("CATEGORY=yes_no_hierarchy\nMEMORY=80000\nCORES=2\nWALL_TIME=172800") %>% 
  write_lines(glue("{input_address}makeflow.makeflow"))
  
  
path = "/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/jian_R_files/error_files/"
sprintf("%s:\n\tRscript %screate_input.R %s", path, path, path) %>% 
  prepend("CATEGORY=yes_no_hierarchy\nMEMORY=1000\nCORES=1\nWALL_TIME=200") %>% 
  write_lines(glue("{path}test.makeflow"))

raw_counts_kam %>%
  
  adapt_tree(.tree = kamran_tree, .node = NULL) %>%
  
  tree_and_signatures_to_database(tree=kamran_tree, ., .sample=sample, .cell_type=cell_type,
                                  .symbol=symbol, .count=count) %>% 
  
  # Remove redundant samples
  remove_redundancy(.element=sample, .feature=symbol, .abundance=count, 
                    correlation_threshold = 0.999, top = 500, method = "correlation") %>%
  droplevels() %>% 
  
  # Eliminate suspicious samples
  filter(!grepl("GSM3722278|GSM3722276|GSM3722277", sample)) %>%
  
  do_scaling(.sample = sample, .symbol= symbol , .count= count, .cell_type=cell_type) %>% 
  
  saveRDS("dev/intermediate_data/counts_scaled_kamran.rds", compress = "xz")

counts_scaled_kamran %>% 
  
  do_imputation(.sample = .sample, .symbol = .feature, .count = count, .cell_type = cell_type) %>% 
  
  saveRDS("dev/intermediate_data/counts_imputed_kamran.rds", compress = "xz")

counts_hierarchy_kamran <- counts_imputed_kamran %>% 
  dplyr::rename(symbol = .feature, sample = .sample) %>% 
  do_hierarchy(.sample=sample,
               .symbol=symbol,
               .cell_type = cell_type,
               .tree = kamran_tree,
               .is_hierarchy=TRUE)

counts_ranked_kamran <- counts_hierarchy_kamran %>% 
  do_ranking(.sample=sample, 
             .symbol=symbol,
             .cell_type = cell_type,
             .ranking_method=rank_edgR_quasi_likelihood, 
             .contrast_method=pairwise_contrast, 
             .rank_stat="PValue", 
             .tree = kamran_tree)


de_kamran <- counts_hierarchy_kamran %>%
  unnest(tt) %>%
  
  # Differential transcription: generate contrast
  mutate(markers = map2(
    data, level,
    ~ .x %>%
      test_differential_abundance(
        .formula = as.formula(sprintf("~ 0 + %s", .y)),
        .sample = sample,
        .transcript = symbol,
        .abundance = count_scaled,
        .contrasts = pairwise_contrast(.x, .y),
        method = "edger_robust_likelihood_ratio",
        test_above_log2_fold_change = 1,
        action="only")
  ))
 

de_kamran %>% 
  slice_tail(n=2) %>% 
  saveRDS("dev/intermediate_data/kamran_memory_de.rds", compress = "xz")

x <- kamran_memory_de %>% 
  mutate(markers = map(
    markers,
    ~ rank_by_stat(.x, "PValue")
  ))
  
de_jian %>% 
  mutate(markers = map(
    markers,
    ~ rank_by_stat(.x, "PValue")
  ))
  
counts_imputed_kamran %>% 
  as_tibble() %>% 
  group_by(cell_type) %>% 
  summarise(n_sample = n_distinct(.sample), 
            n_gene = n_distinct(.feature),
            .groups = "drop") %>% 
  saveRDS("dev/intermediate_data/n_sample_after_imputation_kamran.rds", compress = "xz")

counts_imputed %>% 
  group_by(cell_type) %>% 
  summarise(n_sample = n_distinct(sample), 
            n_gene = n_distinct(feature),
            .groups = "drop") %>% 
  saveRDS("dev/intermediate_data/n_sample_after_imputation_jian.rds", compress = "xz")

n_sample_after_imputation_kamran$cell_type %>% 
  .[!n_sample_after_imputation_kamran$cell_type %in% n_sample_after_imputation_jian$cell_type]

n_sample_after_imputation_jian %>% 
  left_join(n_sample_after_imputation_kamran, by = "cell_type", suffix = c(".jian", ".kamran")) %>% 
  drop_na() %>% 
  print(n=37)


de_jian <- counts_imputed_hierarchy %>%
  unnest(tt) %>%
  slice(13, 15) %>% 
  # Differential transcription: generate contrast
  mutate(markers = map2(
    data, level,
    ~ .x %>%
      test_differential_abundance(
        .formula = as.formula(sprintf("~ 0 + %s", .y)),
        .sample = sample,
        .transcript = symbol,
        .abundance = count_scaled,
        .contrasts = pairwise_contrast(.x, .y),
        method = "edger_robust_likelihood_ratio",
        test_above_log2_fold_change = 1,
        action="only")
  ))


de_jian %>% 
  pluck("markers", 1) %>% 
  pivot_longer(contains("___"), 
               names_to = c("stats", "contrast"), 
               values_to = ".value", 
               names_sep = "___") %>% 
  nest(stat_df = - contrast) %>% 
  mutate(stat_df = map(stat_df, ~.x %>%
                         pivot_wider(names_from = stats, values_from = ".value" )
                       ))

de_jian %>% 
  pluck("markers", 2) %>% 
  pivot_longer(contains("___"), 
               names_to = c("stats", "contrast"), 
               values_to = ".value", 
               names_sep = "___") %>% 
  nest(stat_df = - contrast) %>% 
  mutate(stat_df = map(stat_df, ~.x %>%
                         pivot_wider(names_from = stats, values_from = ".value" )
  )) %>% 
  pluck("stat_df", 1) %>% 
  arrange(FDR)


x = kamran_memory_de %>% 
  # pluck("markers", 1) %>% 
  mutate(markers = map(
    markers,
    ~ .x %>% 
      pivot_longer(
        cols = contains("___"),
        names_to = c("stats", "contrast"),
        values_to = ".value",
        names_sep="___"
      ) %>%
      # Markers selection within each pair of contrast
      nest(stat_df = -contrast) %>%
      # Reshape inside each contrast
      mutate(stat_df = map(stat_df, ~.x %>% pivot_wider(names_from = stats, values_from = .value)))
  ))


x %>% 
  pluck("markers", 2) %>% 
  pluck("stat_df", 1) %>% 
  filter(FDR < 0.05 & logFC > 2) %>%
  filter(logCPM > mean(logCPM))
