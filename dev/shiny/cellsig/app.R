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

load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/data/counts.rda")

processed <- counts %>% 
    tidybulk(sample, symbol, count) %>%
    identify_abundant(factor_of_interest = cell_type) %>%
    scale_abundance()
    
    select(-cell_type) %>% 
    pivot_longer(contains("level"), names_to="level", values_to="cell_type") %>% 
    filter(is.na(cell_type)==F & is.na(symbol)==F) %>% 
    tidybulk(sample, symbol, count) %>%
    nest(data = -level) %>%
    mutate(data = map(data, ~ droplevels(.x))) %>%
    # mutate(data = map(data, ~aggregate_duplicates(.x))) %>%
    # mutate(data = map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
    mutate(data = map(
        data, ~ .x %>% 
            identify_abundant(factor_of_interest = cell_type) %>%
            scale_abundance(sample, symbol, count) 
    ))


# mutate(data = map(
#     data, ~ .x %>% 
#         identify_abundant(factor_of_interest = cell_type) %>%
#         scale_abundance(.sample=sample, .transcript=symbol, .abundance=count) 
# ))


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Gene Expression Visualisation"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            selectInput("gene",
                        h3("Enter a gene:"),
                        choices = counts %>% select(symbol) %>% filter(is.na(symbol)==F) %>% distinct() %>% pull()
            ),
            
            actionButton("submit",
                         "Submit")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Level 1", plotOutput("geneExprL1")),
                        tabPanel("Level 2", plotOutput("geneExprL2")),
                        tabPanel("Level 3", plotOutput("geneExprL3")),
                        tabPanel("Level 4", plotOutput("geneExprL4")),
                        tabPanel("Level 5", plotOutput("geneExprL5"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    gene <- reactive({input$gene})
    
    output$geneExprL1 <- renderPlot({
        
        observeEvent(input$submit, {
            tt %>%
                filter(level=="level_1") %>%
                unnest(data) %>%
                filter(symbol==gene) %>%
                ggplot(aes(cell_type, counts_scaled+1, color=cell_type)) +
                scale_y_log10()+
                geom_point() +
                # facet_wrap(~sample) +
                theme_bw() +
                theme(text = element_text(size=6),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        })
        
    })
    
    output$geneExprL2 <- renderPlot({
            
            observeEvent(input$submit, {
                tt %>%
                    filter(level==2) %>%
                    unnest(data) %>%
                    filter(symbol==gene) %>%
                    ggplot(aes(cell_type, counts_scaled+1, color=cell_type)) +
                    scale_y_log10()+
                    geom_point() +
                    # facet_wrap(~sample) +
                    theme_bw() +
                    theme(text = element_text(size=6),
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            })
        })
        
    output$geneExprL3 <- renderPlot({
            
        observeEvent(input$submit, {
            tt %>%
                filter(level==3) %>%
                unnest(data) %>%
                filter(symbol==gene) %>%
                ggplot(aes(cell_type, counts_scaled+1, color=cell_type)) +
                scale_y_log10()+
                geom_point() +
                # facet_wrap(~sample) +
                theme_bw() +
                theme(text = element_text(size=6),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            })
        })

    output$geneExprL4 <- renderPlot({
        
        observeEvent(input$submit, {
            tt %>%
                filter(level==4) %>%
                unnest(data) %>%
                filter(symbol==gene) %>%
                ggplot(aes(cell_type, counts_scaled+1, color=cell_type)) +
                scale_y_log10()+
                geom_point() +
                # facet_wrap(~sample) +
                theme_bw() +
                theme(text = element_text(size=6),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        })
    })
        
    output$geneExprL5 <- renderPlot({
        
        observeEvent(input$submit, {
            tt %>%
                filter(level==5) %>%
                unnest(data) %>%
                filter(symbol==gene) %>%
                ggplot(aes(cell_type, counts_scaled+1, color=cell_type)) +
                scale_y_log10()+
                geom_point() +
                # facet_wrap(~sample) +
                theme_bw() +
                theme(text = element_text(size=6),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        })
    })
   
 }

# Run the application 
shinyApp(ui = ui, server = server)
