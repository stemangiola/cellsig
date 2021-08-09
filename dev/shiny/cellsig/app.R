#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(tidyverse)
library(tidybulk)
library(ggplot2)
library(plotly)

load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")

# Functions

format_name <- function(.method) {
    paste("markers", .method, sep = "_")
}

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


# Obsolete
# processed <- counts %>% 
#     select(-cell_type) %>% 
#     pivot_longer(contains("level_"), names_to = "level", values_to="cell_type") %>% 
#     filter(is.na(cell_type)==F & is.na(symbol)==F) %>%
#     tidybulk(sample, symbol, count) %>%
#     nest(data = -level) %>%
#     mutate(data = map(data, ~ droplevels(.x))) %>%
#     # mutate(data = map(data, ~aggregate_duplicates(.x))) %>%
#     # mutate(data = map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
#     mutate(data = map(data, ~ .x %>%
#              identify_abundant(factor_of_interest = cell_type) %>%
#              scale_abundance(sample, symbol, count)
#              ))

processed <- counts %>% 
    
    # filling all the NAs in level_*
    pivot_longer(contains("level_"), names_to = "level", values_to="cell_type2") %>% 
    mutate(cell_type2 = ifelse(is.na(cell_type2), cell_type, cell_type2)) %>% 
    
    # process to scale abundance
    filter(is.na(cell_type2)==F & is.na(symbol)==F) %>%
    tidybulk(sample, symbol, count) %>%
    nest(data = -level) %>%
    mutate(data = map(data, ~ droplevels(.x))) %>%
    # mutate(data = map(data, ~aggregate_duplicates(.x))) %>%
    # mutate(data = map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
    mutate(data = map(data, ~ .x %>%
                          identify_abundant(factor_of_interest = cell_type2) %>%
                          scale_abundance(sample, symbol, count)
    ))


geneNames <- counts %>%
    select(symbol) %>%
    filter(is.na(symbol)==F) %>%
    distinct() %>%
    pull()

tt_all <- tibble(level = 1:5) %>%
    # slice(1) %>%
    mutate(level = paste("level", level, sep = "_")) %>%
    
    # preprocess data
    mutate(tt = map(level, ~ counts %>%
                        mutate(level_0 = "cell") %>%
                        preprocess(.x)))

contrast_all <- tt_all %>% 
    
    # generate contrast by pairwise comparison
    mutate(contrast_PW = map2(tt, level, ~ contrast_PW(.x, .y) )) %>% 
    
    # generate contrast by mean contrast method
    mutate(contrast_MC = map2(tt, level, ~ contrast_MC(.x, .y) )) 


## Identify optimal signature sizes for each cell type (of the ancestor nodes)
opPCA_sig_PW <- c(44, 17, 60, 5, 1, 4, 15, 1, 1, 6, 0, 3, 4, 1, 3)#[1]
opPCA_sig_MC <- c(60, 41, 29, 31, 1, 4, 15, 1, 1, 6, 0, 3, 4, 1, 2)#[1]

## create a collection of all signature genes from each node at each level
sig_collect <- contrast_all %>% 
    
    # reveal all nodes as the signature selection has a specific size for each node
    unnest(tt) %>% 
    mutate(ancestor_type = select(., contains("level_")) %>% 
               pivot_longer(contains("level_"), values_to="cell_type") %>% 
               drop_na() %>% 
               pull(cell_type) ) %>% 
    select(-contains("level_")) %>% 
    
    # select the actual contrast data that corresponds to the ancestor_type 
    # while keeping the structure for sig_select function
    mutate(contrast_PW = map2(contrast_PW, ancestor_type, 
                              ~ .x %>% 
                                  
                                  # create a uniform naming variable for filteration
                                  mutate(node = select(., contains("level_")) %>% 
                                             as_vector()) %>% 
                                  
                                  # select the actual data that corresponds to the ancestor_type outside
                                  filter(node == .y)
    )) %>% 
    mutate(contrast_MC = map2(contrast_MC, ancestor_type, 
                              ~ .x %>% 
                                  
                                  # create uniform naming variable for filteration
                                  mutate(node = select(., contains("level_")) %>% 
                                            as_vector()) %>% 
                                 
                                 # select the actual data that corresponds to the ancestor_type outside
                                 filter(node == .y)
    )) %>% 
    
    # Enter the optimal signature size for each ancestor cell type (or node)
    mutate(opPCA_sig_PW = opPCA_sig_PW) %>% 
    mutate(opPCA_sig_MC = opPCA_sig_MC) %>% 
    
    # select markers of the optimal size for that ancestor cell type
    mutate(markers_PW = pmap(list(contrast_PW, level, opPCA_sig_PW), ~ sig_select(..1, ..2, ..3))) %>% 
    mutate(markers_MC = pmap(list(contrast_MC, level, opPCA_sig_MC), ~ sig_select(..1, ..2, ..3))) %>% 
    
    # select individual cell markers for the cell types at each node
    mutate(cell_markers_PW = map(markers_PW, ~ cell_sig_select(.x) )) %>% 
                                 
    mutate(cell_markers_MC = map(markers_MC, ~ cell_sig_select(.x) )) %>% 
    
    # obtain plot data for either PCA or tSNE method
    mutate(sil_PW_PCA = pmap(list(markers_PW, level, "PCA"), ~ sil_func(..1, ..2, ..3))) %>% 
    mutate(sil_MC_PCA = pmap(list(markers_MC, level, "PCA"), ~ sil_func(..1, ..2, ..3))) %>% 
    mutate(sil_PW_tSNE = pmap(list(markers_PW, level, "tSNE"), ~ sil_func(..1, ..2, ..3))) %>% 
    mutate(sil_MC_tSNE = pmap(list(markers_MC, level, "tSNE"), ~ sil_func(..1, ..2, ..3))) %>% 

    select(-c(data, contrast_PW, contrast_MC, sil_MC_PCA, sil_PW_tSNE))



# Define UI for application that draws a histogram
ui <- fluidPage(
    
    tabsetPanel(type = "tabs",
                
                # First Tab - signature expression =========================

                tabPanel("Gene Expression Visualization",

                         # Application title
                         titlePanel("Gene Expression Visualization"),

                         # Sidebar with a slider input for number of bins
                         sidebarLayout(

                             sidebarPanel(

                                 selectInput("gene",
                                             h3("Enter a gene: "),
                                             choices = geneNames
                                 ),

                                 actionButton("submit",
                                              "Submit"),

                                 width = 3

                             ),

                             # Show a plot of the generated distribution
                             mainPanel(
                                 fluidRow(
                                     h4("Level 1"),
                                     plotOutput("geneExprL1")),
                                 fluidRow(
                                     h4("Level 2"),
                                     plotOutput("geneExprL2")),
                                 fluidRow(
                                     h4("Level 3"),
                                     plotOutput("geneExprL3")),
                                 fluidRow(
                                     h4("Level 4"),
                                     plotOutput("geneExprL4")),
                                 fluidRow(
                                     h4("Level 5"),
                                     plotOutput("geneExprL5")),
                                 width = 9
                             )
                         )
                 ),

                # Second Tab Plots & Signature =================================

                tabPanel("Plots & Signatures",
                         
                         # Application title
                         titlePanel("Plots & Signatures"),
                         
                         fluidRow(
                             
                             wellPanel(
                                 
                                 h3("Settings"),
                                 
                                 fluidRow(
                                     
                                     column(6, 
                                            selectInput("rdim_method",
                                                        h4("Select method of dimensionality reduction: "),
                                                        choices = list("PCA", "tSNE"),
                                                        selected = "PCA")
                                            ),
                                     column(6, 
                                            selectInput("sig_select_method",
                                                        h4("Choose method of signature selection: "),
                                                        choices = list("Pairwise Comparison" = "PW", 
                                                                       "Mean Contrast" = "MC"),
                                                        selected = "PW")
                                            )
                                     
                                        ),
                                 
                                 actionButton("run", "Run")
                                 
                                        )
                             
                                ),
                         
                         mainPanel(
                             
                             fluidRow(
                                 
                                 h4("Level 1"),
                                 
                                 column(9,
                                        plotlyOutput("plot_cell")  
                                        ),
                                 
                                 column(3, 
                                        uiOutput("markers_count"),
                                        DT::dataTableOutput("markers_cell")
                                        )
                                 
                                    )

                                 )
                             
                        )
                         
              )
                
)

    
 


# Define server logic required to draw a histogram============

server <- function(input, output) {
    
    # Tab 1 output =============================
    gene <- reactive({input$gene})

    plotData1 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_1") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })

    output$geneExprL1 <- renderPlot({

        plotData1() %>%
            ggplot(aes(cell_type2, count_scaled+1, color=cell_type2)) +
            labs(title=isolate({gene()}), y="Scaled count", color="Cell type") +
            scale_y_log10() +
            geom_violin() +
            geom_jitter(alpha=0.3) +
            theme_bw() +
            theme(text = element_text(size=20),
                  plot.title = element_text(hjust = 0.5),
                  axis.title.x=element_blank(),
                  axis.text.x = element_text(angle=45, vjust=1, hjust=1))

    })

    plotData2 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_2") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })

    output$geneExprL2 <- renderPlot({

        plotData2() %>%
            ggplot(aes(cell_type2, count_scaled+1, color=cell_type2)) +
            labs(title=isolate({gene()}), y="Scaled count", color="Cell type") +
            scale_y_log10() +
            geom_violin() +
            geom_jitter(alpha=0.3) +
            theme_bw() +
            theme(text = element_text(size=20),
                  plot.title = element_text(hjust = 0.5),
                  axis.title.x=element_blank(),
                  axis.text.x = element_text(angle=45, vjust=1, hjust=1))

    })

    plotData3 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_3") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })

    output$geneExprL3 <- renderPlot({

        plotData3() %>%
            ggplot(aes(cell_type2, count_scaled+1, color=cell_type2)) +
            labs(title=isolate({gene()}), y="Scaled count", color="Cell type") +
            scale_y_log10() +
            geom_violin() +
            geom_jitter(alpha=0.3) +
            theme_bw() +
            theme(text = element_text(size=20),
                  plot.title = element_text(hjust = 0.5),
                  axis.title.x=element_blank(),
                  axis.text.x = element_text(angle=45, vjust=1, hjust=1))

    })

    plotData4 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_4") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })

    output$geneExprL4 <- renderPlot({

        plotData4() %>%
            ggplot(aes(cell_type2, count_scaled+1, color=cell_type2)) +
            labs(title=isolate({gene()}), y="Scaled count", color="Cell type") +
            scale_y_log10() +
            geom_violin() +
            geom_jitter(alpha=0.3) +
            theme_bw() +
            theme(text = element_text(size=20),
                  plot.title = element_text(hjust = 0.5),
                  axis.title.x=element_blank(),
                  axis.text.x = element_text(angle=45, vjust=1, hjust=1))

    })

    plotData5 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_5") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })

    output$geneExprL5 <- renderPlot({

        plotData5() %>%
            ggplot(aes(cell_type2, count_scaled+1, color=cell_type2)) +
            labs(title=isolate({gene()}), y="Scaled count", color="Cell type") +
            scale_y_log10() +
            geom_violin() +
            geom_jitter(alpha=0.3) +
            theme_bw() +
            theme(text = element_text(size=20),
                  plot.title = element_text(hjust = 0.5),
                  axis.title.x=element_blank(),
                  axis.text.x = element_text(angle=45, vjust=1, hjust=1))

    })


    # Tab 2 output ===========================
    
    rdim_method <- reactive({input$rdim_method})
    
    sig_select_method <- reactive({input$sig_select_method})
    
    # Plot output
    plotCell <- eventReactive(input$run, {
        if(isolate({rdim_method()}) =="PCA") {
            
            which_sil <- paste("sil", sig_select_method(), rdim_method(), sep = "_")
            which_cell_markers <- paste("cell_markers", sig_select_method(), sep = "_")
            
            signature_cell <- sig_collect %>% 
                pluck(which_cell_markers, 1) 
            
            p_cell <- sig_collect %>% 
                pluck(which_sil, 1) %>% 
                pluck("rdim", 1) %>%
                left_join(signature_cell, by = c("level_1" = "cell_type")) %>% 
                plot_ly(x = ~ PC1, y = ~ PC2) %>% 
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
        } else if (isolate({rdim_method()}) =="tSNE") {
            
            which_sil <- paste("sil", sig_select_method(), rdim_method(), sep = "_")
            which_cell_markers <- paste("cell_markers", sig_select_method(), sep = "_")
            
            signature_cell <- sig_collect %>% 
                pluck(which_cell_markers, 1) 
            
            p_cell <- sig_collect %>% 
                pluck(which_sil, 1) %>% 
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
        }
        return(p_cell)
    })
    
    output$plot_cell <- renderPlotly({ plotCell() })
    
    # Markers output
    markersCell <- eventReactive(input$run, {
        if(isolate({sig_select_method()}) =="PW") {
            marker <- sig_collect %>% 
                pluck("markers_PW", 1) %>% 
                select(Gene = symbol) %>% 
                distinct()
            
        } else if (isolate({sig_select_method()}) =="MC") {
            marker <- sig_collect %>% 
                pluck("markers_MC", 1) %>% 
                select(Gene = symbol) %>% 
                distinct()
        }
        return(marker)
    })
    
    
    output$markers_count <- renderUI({  
        tags$b(
            paste(nrow(markersCell()), 
                  "markers are selected for", 
                  sig_collect$ancestor_type[1], 
                  ": ")
        )
    })
    
    output$markers_cell <- DT::renderDataTable({ markersCell() })
    
 }

# Run the application 
shinyApp(ui = ui, server = server)
