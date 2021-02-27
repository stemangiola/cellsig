#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(tidybulk)
library(ggplot2)
library(plotly)

load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")

format_name <- function(.method) {
    paste("markers", .method, sep = "_")
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
    pivot_longer(contains("level_"), names_to = "level", values_to="cell_type2") %>% 
    mutate(cell_type2 = ifelse(is.na(cell_type2), cell_type, cell_type2)) %>% 
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
    mutate(level = paste("level", level, sep = "_")) %>% 
    
    # preprocess data
    mutate(tt = map(level, ~ counts %>% 
                        mutate(level_0 = "cell") %>% 
                        preprocess(.x))) %>% 
    
    # generate contrast by pairwise comparison
    mutate(contrast_PW = map2(tt, level, ~ contrast_PW(.x, .y) )) %>% 
    
    # generate contrast by mean contrast method
    mutate(contrast_MC = map2(tt, level, ~ contrast_MC(.x, .y) ))


## Identify optimal signature sizes for each cell type (of the ancestor nodes)
optim_sig_PW <- c(44, 17, 60, 5, 1, 4, 15, 1, 1, 6, 0, 3, 4, 1, 3)
optim_sig_MC <- c(60, 41, 29, 31, 1, 4, 15, 1, 1, 6, 0, 3, 4, 1, 2)

## create a collection of all signature genes from each node at each level
markers_collect <- tt_all %>% 
    unnest(tt) %>% 
    mutate(optim_sig_PW = optim_sig_PW) %>% 
    mutate(optim_sig_MC = optim_sig_MC) %>% 
    # nest(tt = -c(level, contrast_PW, contrast_MC, optim_sig_PW, optim_sig_MC)) %>% 
    mutate(markers_PW = pmap(list(contrast_PW, level, optim_sig_PW), ~ sig_select(..1, ..2, ..3))) %>% 
    mutate(markers_MC = pmap(list(contrast_MC, level, optim_sig_MC), ~ sig_select(..1, ..2, ..3))) %>% 
    select(-c(level_0, level_1, level_2, level_3, level_4, data))


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
                                 column(6,
                                        plotOutput("plot_cell")  
                                        ),
                                 textOutput("sig_cell")),
                             
                             fluidRow(
                                 h4("Level 2"),
                                 column(6,
                                        plotOutput("plot_immune")
                                        ),
                                 textOutput("sig_immune")),
                             
                             fluidRow(
                                 h4("Level 3"),
                                 fluidRow(
                                     column(6, 
                                            plotOutput("plot_mono_derived"),
                                            textOutput("sig_mono_derived")
                                            ),
                                     column(6, 
                                            plotOutput("plot_t_cell"),
                                            textOutput("sig_t_cell"))
                                     ),
                                 
                                 fluidRow(
                                     column(6, 
                                            plotOutput("plot_granulocyte"),
                                            textOutput("sig_granulocyte")
                                            ),
                                     column(6,
                                            plotOutput("plot_b_cell"),
                                            textOutput("sig_b_cell")
                                            )
                                     ),
                                 
                                 fluidRow(
                                     column(6, 
                                            plotOutput("plot_NK"),
                                            textOutput("sig_NK")
                                            )
                                        )
                                    ),
                             
                             fluidRow(
                                 h4("Level 4"),
                                 fluidRow(
                                     column(6,
                                            plotOutput("plot_t_CD4"),
                                            textOutput("sig_t_CD4")
                                            ),
                                     column(6,
                                            plotOutput("plot_macrophage"),
                                            textOutput("sig_macrophage")
                                            )
                                 ),
                                 
                                 fluidRow(
                                     column(6,
                                            plotOutput("plot_t_CD8"),
                                            textOutput("sig_t_CD8")
                                            ),
                                     column(6,
                                            plotOutput("plot_DC_myeloid"),
                                            textOutput("sig_DC_myeloid")
                                     )
                                 ),
                                 
                                 fluidRow(
                                     column(6,
                                            plotOutput("plot_NK_primed"),
                                            textOutput("sig_NK_primed")
                                            )
                                    )
                                 ),
                             
                             fluidRow(
                                 h4("Level 5"),
                                 fluidRow(
                                     column(6,
                                            plotOutput("plot_t_CD4_memory"),
                                            textOutput("sig_t_CD4_memory")
                                            ),
                                     column(6,
                                            plotOutput("plot_t_CD8_memory"),
                                            textOutput("sig_t_CD8_memory")
                                            )
                                 ),
                                 
                                 fluidRow(
                                     column(6,
                                            plotOutput("plot_t_helper"),
                                            textOutput("sig_t_helper")
                                            )
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

    sig_data <- eventReactive(input$run, {
        markers_collect %>% 
            mutate(sil = 
                       pmap(list(!!as.symbol(format_name(sig_select_method())), level, rdim_method()),
                            ~ sil_func(..1, ..2, ..3))) %>% 
            mutate(sig =
                       map(!!as.symbol(format_name(sig_select_method())), 
                           ~.x %>% 
                               pull(symbol) %>% 
                               unique()))
    })
    
    
    output$plot_cell <- renderPlot({
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 1) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_1, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 1) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_1, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_cell <- renderText({
        sig_data() %>% 
            pluck("sig", 1) %>% 
            str_pad(10, "both")
    })
    
    output$plot_immune <- renderPlot({
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 2) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_2, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw() 
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 2) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_2, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_immune <- renderText({
        sig_data() %>% 
            pluck("sig", 2) %>% 
            str_pad(10, "both")
    })
    
    output$plot_mono_derived <- renderPlot({
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 3) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 3) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_mono_derived <- renderText({
        sig_data() %>% 
            pluck("sig", 3) %>% 
            str_pad(10, "both")
    })
    
    output$plot_t_cell <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 4) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 4) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_t_cell <- renderText({
        sig_data() %>% 
            pluck("sig", 4) %>% 
            str_pad(10, "both")
    })
    
    output$plot_granulocyte <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 5) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 5) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_granulocyte <- renderText({
        sig_data() %>% 
            pluck("sig", 5) %>% 
            str_pad(10, "both")
    })
    
    output$plot_b_cell <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 6) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 6) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_b_cell <- renderText({
        sig_data() %>% 
            pluck("sig", 6) %>% 
            str_pad(10, "both")
    })
    
    output$plot_NK <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 7) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 7) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_3, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_NK <- renderText({
        sig_data() %>% 
            pluck("sig", 7) %>% 
            str_pad(10, "both")
    })
    
    output$plot_t_CD4 <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 8) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 8) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_t_CD4 <- renderText({
        sig_data() %>% 
            pluck("sig", 8) %>% 
            str_pad(10, "both")
    })
    
    output$plot_macrophage <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 9) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 9) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_macrophage <- renderText({
        sig_data() %>% 
            pluck("sig", 9) %>% 
            str_pad(10, "both")
    })
    
    output$plot_t_CD8 <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 10) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 10) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_t_CD8 <- renderText({
        sig_data() %>% 
            pluck("sig", 10) %>% 
            str_pad(10, "both")
    })
    
    output$plot_DC_myeloid <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 11) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 11) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_DC_myeloid <- renderText({
        sig_data() %>% 
            pluck("sig", 11) %>% 
            str_pad(10, "both")
    })
    
    output$plot_NK_primed <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 12) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 12) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_4, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_NK_primed <- renderText({
        sig_data() %>% 
            pluck("sig", 12) %>% 
            str_pad(10, "both")
    })
    
    output$plot_t_CD4_memory <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 13) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_5, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 13) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_5, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_t_CD4_memory <- renderText({
        sig_data() %>% 
            pluck("sig", 13) %>% 
            str_pad(10, "both")
    })
    
    output$plot_t_CD8_memory <- renderPlot({
        
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 14) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_5, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 14) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_5, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_t_CD8_memory <- renderText({
        sig_data() %>% 
            pluck("sig", 14) %>% 
            str_pad(10, "both")
    })
    
    output$plot_t_helper <- renderPlot({
        if(isolate({rdim_method()}) =="PCA") {
            sig_data() %>% 
                pluck("sil", 15) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = PC1, y = PC2, colour = level_5, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        } else if (isolate({rdim_method()}) =="tSNE") {
            sig_data() %>% 
                pluck("sil", 15) %>% 
                pluck("Rdim", 1) %>% 
                ggplot(aes(x = tSNE1, y = tSNE2, colour = level_5, label = sample)) + 
                geom_point() +
                stat_ellipse(type = 't')+
                theme_bw()
        }
    })
    
    output$sig_t_helper <- renderText({
        sig_data() %>% 
            pluck("sig", 15) %>% 
            str_pad(10, "both")
    })
 }

# Run the application 
shinyApp(ui = ui, server = server)
