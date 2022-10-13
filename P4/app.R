library(shiny)
library(ggplot2)
library(dplyr)

expressionData <- read.table("expression_data.txt", header = T, stringsAsFactors = FALSE)
expressionData$median_expression <- log(expressionData$median_expression) # Normalize the data

ui <- fluidPage(
  titlePanel("P4: Interactive and dynamic data visualization"),
  "First steps with ", code("shiny"), br(),
  strong("Alessandro Giulivo"), br(),
  em("15 October 2021"), br(),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("expressionInput", "Median expression (log)", min = -6, max = 12,
                             value = c(1, 12)),
      checkboxGroupInput("tissueInput", "Tissue",
                     choices = c("Adipose", "Brain", "Liver",
                                 "Lung", "Lymphocytes", "Muscle",
                                 "Stomach", "Testis"),
                     selected = "Brain"),
  
      radioButtons("geneAnnotInput", "Gene annotation",
               choices = c("ENSEMBL", "HAVANA"),
               selected = "HAVANA"),
  
      selectInput("geneTypeInput", "Gene type",
              choices = c("antisense", "lincRNA", "miRNA", "misc_RNA", "protein_coding",
                          "pseudogene", "rRNA", "snoRNA", "snRNA"),
              selected = "protein_coding"),
  
      checkboxGroupInput("ChrInput", "Chromosome",
                     choices = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                 "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                 "21", "22", "X", "Y"),
                     selected = c("1", "2", "3", "4"))
      ),
    
    mainPanel(plotOutput("plot"), br(), br(),
              textOutput("number"), br(), br(),
              tableOutput("results"))
    
  )
  
  
)

server <- function(input, output) {
  output$plot <- renderPlot({
    
    filtered <-
      expressionData %>%
      filter(median_expression >= input$expressionInput[1],
             median_expression <= input$expressionInput[2],
             tissue == input$tissueInput,
             gene_annotation == input$geneAnnotInput,
             gene_type == input$geneTypeInput,
             gene_chr == input$ChrInput
      )
    
    ggplot(filtered, aes(median_expression)) +
      geom_histogram() +
      scale_fill_brewer(palette = "Blues") +
      facet_wrap(.~gene_chr)
  })
  
  output$number <- renderText({
    filtered <-
      expressionData %>%
      filter(
        median_expression >= input$expressionInput[1],
        median_expression <= input$expressionInput[2],
        tissue == input$tissueInput,
        gene_annotation == input$geneAnnotInput,
        gene_type == input$geneTypeInput,
        gene_chr == input$ChrInput
      )
    paste("Number of results found =", nrow(filtered))
  })
  
  output$results <- renderTable({
    filtered <-
      expressionData %>%
      filter(
        median_expression >= input$expressionInput[1],
        median_expression <= input$expressionInput[2],
        tissue == input$tissueInput,
        gene_annotation == input$geneAnnotInput,
        gene_type == input$geneTypeInput,
        gene_chr == input$ChrInput
      )
    filtered
  })

}

shinyApp(ui = ui, server = server)
