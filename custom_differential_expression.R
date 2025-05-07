# custom_differential_expression.R

library(shiny)
library(Seurat)

custom_differential_expression_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 17: Custom Differential Expression"),
    textInput(ns("custom_genes"), "Enter Genes (comma-separated):"),
    actionButton(ns("run_custom_differential"), "Run Custom Differential Expression"),
    plotOutput(ns("custom_differential_plot")),
    actionButton(ns("finish_analysis"), "Finish Analysis")
  )
}

custom_differential_expression_server <- function(id, data_for_custom_differential) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$run_custom_differential, {
      data <- data_for_custom_differential()
      
      if (is.null(data)) return()
      
      genes <- strsplit(input$custom_genes, ",")[[1]]
      
      # Run custom differential expression analysis
      custom_results <- FindMarkers(data, features = genes)
      
      # Plot custom differential expression results
      output$custom_differential_plot <- renderPlot({
        FeaturePlot(data, features = genes[1:10])  # Example
      })
      
      shinyjs::enable("finish_analysis")
    })
  })
}
