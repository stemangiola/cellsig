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

# load("/stornext/Home/data/allstaff/w/wu.j/Master Project/cellsig/dev/counts.rda")
# 
# processed <- counts %>%
#     select(-cell_type) %>%
#     pivot_longer(contains("level"), names_to="level", values_to="cell_type") %>%
#     filter(is.na(cell_type)==F & is.na(symbol)==F) %>%
#     tidybulk(sample, symbol, count) %>%
#     nest(data = -level) %>%
#     mutate(data = map(data, ~ droplevels(.x))) %>%
#     # mutate(data = map(data, ~aggregate_duplicates(.x))) %>%
#     # mutate(data = map(data, ~ fill_missing_abundance(.x, fill_with = 0))) %>%
#     mutate(data = map(
#         data, ~ .x %>%
#             identify_abundant(factor_of_interest = cell_type) %>%
#             scale_abundance(sample, symbol, count)
#     ))
# 
# geneNames <- counts %>%
#     select(symbol) %>%
#     filter(is.na(symbol)==F) %>%
#     distinct() %>%
#     pull()


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Gene Expression Visualization"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        
        sidebarPanel(
            
            selectInput("gene",
                        h3("Enter a gene:"),
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
            
            # tabsetPanel(type = "tabs",
            #             tabPanel("Level 1", plotOutput("geneExprL1")),
            #             tabPanel("Level 2", plotOutput("geneExprL2")),
            #             tabPanel("Level 3", plotOutput("geneExprL3")),
            #             tabPanel("Level 4", plotOutput("geneExprL4")),
            #             tabPanel("Level 5", plotOutput("geneExprL5"))
            # )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    gene <- reactive({input$gene})
    
    plotData1 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_1") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })
    
    plotData2 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_2") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })
    
    plotData3 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_3") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })
    
    plotData4 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_4") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })
    
    plotData5 <- eventReactive(input$submit, {
        processed %>%
            filter(level=="level_5") %>%
            unnest(data) %>%
            filter(symbol==gene())
    })
    
    output$geneExprL1 <- renderPlot({
        
        plotData1() %>% 
            ggplot(aes(cell_type, count_scaled+1, color=cell_type)) +
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
    
    output$geneExprL2 <- renderPlot({
        
        plotData2() %>% 
            ggplot(aes(cell_type, count_scaled+1, color=cell_type)) +
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

    output$geneExprL3 <- renderPlot({
        
        plotData3() %>% 
            ggplot(aes(cell_type, count_scaled+1, color=cell_type)) +
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
    
    output$geneExprL4 <- renderPlot({
        
        plotData4() %>% 
            ggplot(aes(cell_type, count_scaled+1, color=cell_type)) +
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
    
    output$geneExprL5 <- renderPlot({
        
        plotData5() %>% 
            ggplot(aes(cell_type, count_scaled+1, color=cell_type)) +
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
 }

# Run the application 
shinyApp(ui = ui, server = server)
