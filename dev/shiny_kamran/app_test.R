# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# specify the repository for packages from Bioconductor
options(repos = BiocManager::repositories())

library(shiny)
library(ggrepel)
library(ggnewscale)
library(rsconnect)
library(yaml)
library(tidytext)
library(data.tree)
library(tidytree)
library(ape)
library(scales)
library(data.tree)
library(tidyverse)
library(tidybulk)
library(tidySummarizedExperiment)
# library(profvis)

# profvis({

# load data locally
# new_tree <- read_yaml("dev/shiny_pub/new_tree.yaml") %>%
#   as.Node
# load("dev/shiny_pub/geneNames.RData")
# load("dev/shiny_pub/expression.RData")
# load("dev/shiny_pub/pca_df.RData")

# setwd("/stornext/Home/data/allstaff/k/khan.k/cellsig/dev/shiny_kamran")
#  load data for webpage deployment
new_tree <- read_yaml("tree.yaml") %>% as.Node
load("geneNames.RData")
load("expression.RData")
load("pca_df.RData")
load("data_for_download.RData")
load("celltypes.RData")
load("bayes.RData")


# expression = expression@assays@data$count_scaled %>% Matrix(sparse = TRUE)


# Define UI for application that draws a histogram
ui <- navbarPage(
  
  title = "HBCC",
  
  tags$style(type='text/css', 
             ".selectize-input {font-size: 20px; line-height: 22px;} 
             .selectize-dropdown {font-size: 20px; line-height: 22px;}"),
  
  
  tabPanel(
    "Gene Expression Comparison",
    
    sidebarLayout(
      
      sidebarPanel(
        selectInput("cell",
                    h3(strong("Cellular Compartment")),
                    choices = list("All", "Immune", "Non-immune"),
                    selected = "All"
        ),
        selectizeInput("gene",
                       h3(strong("Gene")),
                       choices = NULL,
                       selected = 1
        ),
        
        width = 3
      ),
      
      mainPanel(
        
        fluidRow(
          
          h3("PCA"),
          
          column(width = 6, 
                 span(textOutput("name"), style="font-size: 16px"),
                 plotOutput("geneExp_pca")
          ),
          column(width = 6,
                 h4("Cell type"),
                 plotOutput("cellType_pca")
          )
        ),
        
        fluidRow(
          h3("Gene Expression"),
          plotOutput("boxplot")
        ),
        
        width = 9
      ) # mainPanel
    ) # sidebarLayout
    
  ) , # tabPanel
  
  tabPanel(
    "Download Dataset",
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Choose dataset ----
        selectInput("dataset", 
                    h3(strong("Download Cell-type wise datasets:")),
                    choices = celltypes),
        
        # Button
        downloadButton("downloadData", "Download"),
        
        
        # adding the new div tag to the sidebar            
        tags$div(class="header", checked=NA,
                 tags$br(),
                 tags$a(href="https://doi.org/10.5281/zenodo.7582421", "Alternatively retrieve the entire database from here", target="_blank")
        )
        
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        
        tableOutput("table")
        
      )
      
    )
    
  )
  
) # navbarPage


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  updateSelectizeInput(session = session, inputId = "gene", choices = geneNames, server = TRUE)
  
  output$name <- renderText({
    req(input$gene)
    input$gene
  })
  
  geneExpPCA <- reactive({
    req(input$gene)
    
    pca_df %>% 
      left_join(expression[input$gene, ] %>%
                  as_tibble(),
                by = c(".sample", "cell_type")
      ) %>%
      
      {
        if (input$cell == "All") {
          (.)
        } else if (input$cell == "Non-immune") {
          (.) %>% 
            filter(is.na(level_2))
        } else { # input$cell == "Immune"
          (.) %>% 
            filter(!is.na(level_2))
        }
      } 
  })
  
  
  output$geneExp_pca <- renderPlot({
    
    geneExpPCA() %>% 
      ggplot(aes(PC1, PC2, colour = count_scaled)) + 
      geom_point()+
      scale_colour_viridis_c(trans = "log10") +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
      )
    
  })
  
  
  
  output$cellType_pca <- renderPlot({
    
    if (input$cell == "All") {
      
      pca_df %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(is.na(level_2)),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Greys") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "t_cell"),
                   aes(color = cell_type)
        ) +
        scale_color_manual(values = gradient_n_pal(brewer_pal(palette="Reds")(9))(seq(0, 1, length.out=10))) +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "b_cell"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Greens") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "granulocyte"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Oranges") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "mono_derived"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Blues") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "natural_killer"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Purples") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "mast_cell"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "YlOrBr") +
        new_scale_color() +
        
        geom_text_repel(data = pca_df %>% 
                          group_by(cell_type) %>% 
                          summarise(PC1 = median(PC1), PC2 = median(PC2), .groups = "drop"),
                        aes(label = cell_type)
        ) +
        theme_bw() +
        theme(plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.position = "none")
      
    }else if(input$cell == "Non-immune") {
      
      pca_df %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(is.na(level_2)),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Greys") +
        geom_text_repel(data = pca_df %>% 
                          # mutate(across(where(is.factor), as.character)) %>% 
                          filter(is.na(level_2)) %>% 
                          group_by(cell_type) %>% 
                          summarise(PC1 = median(PC1), PC2 = median(PC2), .groups = "drop"),
                        aes(label = cell_type)
        ) +
        theme_bw() +
        theme(plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.position = "none")
      
    } else{ # cell == "Immune"
      
      pca_df %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "t_cell"),
                   aes(color = cell_type)
        ) +
        scale_color_manual(values = gradient_n_pal(brewer_pal(palette="Reds")(9))(seq(0, 1, length.out=10))) +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "b_cell"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Greens") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "granulocyte"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Oranges") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "mono_derived"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Blues") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "natural_killer"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "Purples") +
        new_scale_color() +
        geom_point(data = pca_df %>% 
                     # mutate(across(where(is.factor), as.character)) %>% 
                     filter(level_2 == "mast_cell"),
                   aes(color = cell_type)
        ) +
        scale_color_brewer(palette = "YlOrBr") +
        new_scale_color() +
        
        geom_text_repel(data = pca_df %>% 
                          # mutate(across(where(is.factor), as.character)) %>% 
                          filter(!is.na(level_2)) %>%
                          group_by(cell_type) %>% 
                          summarise(PC1 = median(PC1), PC2 = median(PC2), .groups = "drop"),
                        aes(label = cell_type)
        ) +
        theme_bw() +
        theme(plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.position = "none")
      
    }
    
  })
  
  
  
  output$boxplot <- renderPlot({
    req(input$gene)
    
    expression[input$gene, ] %>% 
      as_tibble() %>% 
      # left_join(pca_df %>% select(.sample, cell_type), by = ".sample") %>% 
      # filter(!is.na(cell_type)) %>% 
      filter(cell_type %in% (new_tree %>% as.phylo %>% .$tip.label)) %>%
      left_join(bayes) %>%
      # filter(.feature == gene()) %>% 
      ggplot(aes(cell_type, count_scaled + 1)) +
      scale_y_log10() +
      labs(y="expression") +
      geom_violin() +
      geom_jitter(alpha=0.3) +
      
      ### Add a segment based on 10% and 90% values from bayes data
      geom_errorbar(aes(x = cell_type, ymin = tenth, ymax = ninetieth), width=.2, color = "red") +
      
      
      theme_bw() +
      theme(text = element_text(size=20),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle=45, vjust=1, hjust=1),
            axis.title.x=element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
      )
    
  })
  
  
  datasetInput <- reactive({
    req(input$dataset)
    
    data_for_download %>%
      filter(cell_type == input$dataset)
    
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if (input$dataset == "ALL"){
        paste(input$dataset, ".rds", sep = "")
      } else {
        paste(input$dataset, ".csv", sep = "")
      }
    },
    content = function(file) {
      if (input$dataset == "ALL"){
        saveRDS(data_for_download, file, compress = "xz")
      } else {
        write.csv(datasetInput(), file, row.names = FALSE)
      }
    }
  )
}

# Run the application 

shinyApp(ui = ui, server = server)

# })


# Launnch the shiny app
# rsconnect::deployApp("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/shiny_pub")


# test code

# show gene expression level with PCA
# pca_with_counts_se %>% 
#   # mutate(expression = log1p(count_scaled)) %>% 
#   filter(is.na(level_2)) %>% 
#   filter(.feature == "FOXP3") %>% 
#   ggplot(aes(PC1, PC2, colour = expression, label=sample)) + 
#   geom_point()+
#   scale_colour_viridis_c() +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank()
#   )
# 
# # show annotated cell clusters on PCA
# pca_with_counts %>% 
#   mutate(expression = log1p(count_scaled)) %>% 
#   ggplot(aes(PC1, PC2)) +
#   geom_point(data = pca_with_counts %>% 
#                mutate(across(where(is.factor), as.character)) %>% 
#                filter(is.na(level_2)),
#              aes(color = cell_type)
#   ) +
#   scale_color_brewer(palette = "Greys") +
#   new_scale_color() +
#   geom_point(data = pca_with_counts %>% 
#                mutate(across(where(is.factor), as.character)) %>% 
#                filter(level_2 == "t_cell"),
#              aes(color = cell_type)
#   ) +
#   scale_color_manual(values = gradient_n_pal(brewer_pal(palette="Reds")(9))(seq(0, 1, length.out=10))) +
#   new_scale_color() +
#   geom_point(data = pca_with_counts %>% 
#                mutate(across(where(is.factor), as.character)) %>% 
#                filter(level_2 == "b_cell"),
#              aes(color = cell_type)
#   ) +
#   scale_color_brewer(palette = "Greens") +
#   new_scale_color() +
#   geom_point(data = pca_with_counts %>% 
#                mutate(across(where(is.factor), as.character)) %>% 
#                filter(level_2 == "granulocyte"),
#              aes(color = cell_type)
#   ) +
#   scale_color_brewer(palette = "Oranges") +
#   new_scale_color() +
#   geom_point(data = pca_with_counts %>% 
#                mutate(across(where(is.factor), as.character)) %>% 
#                filter(level_2 == "mono_derived"),
#              aes(color = cell_type)
#   ) +
#   scale_color_brewer(palette = "Blues") +
#   new_scale_color() +
#   geom_point(data = pca_with_counts %>% 
#                mutate(across(where(is.factor), as.character)) %>% 
#                filter(level_2 == "natural_killer"),
#              aes(color = cell_type)
#   ) +
#   scale_color_brewer(palette = "Purples") +
#   new_scale_color() +
#   geom_point(data = pca_with_counts %>% 
#                mutate(across(where(is.factor), as.character)) %>% 
#                filter(level_2 == "mast_cell"),
#              aes(color = cell_type)
#   ) +
#   scale_color_brewer(palette = "YlOrBr") +
#   new_scale_color() +
#   
#   geom_text_repel(data = pca_with_counts %>% 
#                     group_by(cell_type) %>% 
#                     summarise(PC1 = median(PC1), PC2 = median(PC2), .groups = "drop"),
#                   aes(label = cell_type)
#   ) +
#   theme_bw() +
#   theme(plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "none")
# 
# 
# # show gene expression boxplot
# expression %>% 
#   filter(cell_type2 %in% (new_tree %>% as.phylo %>% .$tip.label)) %>% 
#   filter(symbol == "FOXP3") %>% 
#   ggplot(aes(cell_type2, log1p(count_scaled))) +
#   labs(y="expression") +
#   geom_violin() +
#   geom_jitter(alpha=0.3) +
#   theme_bw() +
#   theme(text = element_text(size=20),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle=45, vjust=1, hjust=1),
#         axis.title.x=element_blank(),
#         plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#   )

# non_hierarchical_mean_contrast_edgR_PValue_silhouette_penalty_2 <- 
#   readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/benchmark_results_multiPC_NH/non_hierarchical_mean_contrast_edgR_PValue_silhouette_penalty_2.rds") %>% 
#   pull(signature) %>% 
#   unlist
# 
# geneNames <- 
#   readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_imputed.rds") %>% 
#   rename(symbol = feature) %>%
#   distinct(symbol) %>% 
#   pull
# 
# expression <- 
#   # counts_imputed %>% 
#   readRDS("/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/intermediate_data/counts_imputed.rds") %>%
#   rename(symbol = feature) %>%
#   pivot_longer(contains("level"), names_to = "level", values_to="cell_type2")
# 
# expression = expression %>% 
#   select(symbol, sample, count_scaled, cell_type, level, cell_type2) %>% 
#   mutate(count_scaled = as.integer(count_scaled)) %>% 
#   tidybulk(sample, symbol, count_scaled) %>% 
#   tidybulk::as_SummarizedExperiment(sample, symbol, count_scaled)
# 
# signature_directory = "/stornext/Home/data/allstaff/w/wu.j/Master_Project/cellsig/dev/benchmark_results_multiPC_NH/"
# 
# pca_with_counts <- expression %>% 
#   filter(cell_type2 %in% (new_tree %>% as.phylo %>% .$tip.label)) %>% 
#   filter(symbol %in% non_hierarchical_mean_contrast_edgR_PValue_silhouette_penalty_2) %>% 
#   pivot_wider(names_from = "level", values_from = "cell_type2") %>% 
#   reduce_dimensions(.element = sample,
#                     .feature = symbol,
#                     .abundance = count_scaled,
#                     action = "only",
#                     method = "PCA",
#                     .dims = 2,
#                     log_transform = TRUE,
#                     top = Inf,
#                     scale = FALSE,
#                     check_duplicates = FALSE) %>% 
#   left_join(
#     expression %>% 
#       pivot_wider(names_from = "level", values_from = "cell_type2") %>% 
#       filter(cell_type %in% (new_tree %>% as.phylo %>% .$tip.label)),
#     by="sample"
#   ) %>% 
#   select(sample, contains("PC"), symbol, count_scaled, cell_type, level_2) %>% 
#   # mutate(across(where(is.factor), as.character)) %>% 
#   mutate(count_scaled = count_scaled %>% as.integer) %>% 
#   tidybulk(sample, symbol, count_scaled) %>% 
#   tidybulk::as_SummarizedExperiment(sample, symbol, count_scaled)
# 
# pca_with_counts %>% saveRDS("dev/intermediate_data/pca_with_counts.rds", compress = "xz")
# 
# 
# # OPTIMISATION
# 
# pca_df = 
#   pca_with_counts %>% 
#   select(.sample, PC1, PC2, cell_type,  level_2) %>%
#   distinct() 
# 
# 
# 
# user_choose_a_gene = "CD3G"
# 
# # PCA plotting
# pca_df %>%
#   left_join(
#     expression[user_choose_a_gene,] %>%
#       as_tibble(),
#     by = c(".sample", "cell_type")
#   )
# 
# 
# # Box plot
# expression[user_choose_a_gene,] %>%
#   as_tibble()
# 
# 
# 
# # pca_with_counts_se = pca_with_counts %>% 
# #   mutate(count_scaled = count_scaled + 1L) %>% 
# #   tidybulk(sample, symbol, count_scaled) %>% 
# #   tidybulk::as_SummarizedExperiment(sample, symbol, count_scaled) 
# 