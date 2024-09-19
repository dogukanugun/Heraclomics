# gene_coexpression.R

library(shiny)
library(Seurat)

gene_coexpression_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 10: Gene Co-expression"),
    textInput(ns("genes_coexpression"), "Enter Genes (comma-separated):"),
    actionButton(ns("plot_coexpression"), "Plot Co-expression"),
    plotOutput(ns("coexpression_plot")),
    actionButton(ns("next_step"), "Continue to Gene Regulatory Networks", disabled = TRUE)
  )
}

gene_coexpression_server <- function(id, data_for_coexpression) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$plot_coexpression, {
      data <- data_for_coexpression()
      
      if (is.null(data)) return()
      
      genes <- strsplit(input$genes_coexpression, ",")[[1]]
      
      output$coexpression_plot <- renderPlot({
        FeaturePlot(data, features = genes)
      })
      
      shinyjs::enable("next_step")
    })
  })
}
