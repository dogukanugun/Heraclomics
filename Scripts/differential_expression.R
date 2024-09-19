# differential_expression.R

library(shiny)
library(Seurat)

differential_expression_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 16: Differential Expression"),
    actionButton(ns("run_differential"), "Run Differential Expression"),
    plotOutput(ns("differential_plot")),
    actionButton(ns("next_step"), "Continue to Custom Differential Expression", disabled = TRUE)
  )
}

differential_expression_server <- function(id, data_for_differential) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$run_differential, {
      data <- data_for_differential()
      
      if (is.null(data)) return()
      
      # Run differential expression analysis
      differential_results <- FindMarkers(data)  # Example function
      
      # Plot differential expression results
      output$differential_plot <- renderPlot({
        FeaturePlot(data, features = rownames(differential_results)[1:10])  # Example
      })
      
      shinyjs::enable("next_step")
    })
  })
}
