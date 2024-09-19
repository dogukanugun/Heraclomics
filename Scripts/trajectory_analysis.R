library(shiny)
library(Seurat)
library(slingshot)

trajectory_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 13: Trajectory Analysis"),
    actionButton(ns("run_trajectory"), "Run Trajectory Analysis"),
    plotOutput(ns("trajectory_plot")),
    actionButton(ns("next_step"), "Continue to Cell Communication", disabled = TRUE)
  )
}

trajectory_analysis_server <- function(id, data_for_trajectory) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$run_trajectory, {
      data <- data_for_trajectory()
      
      if (is.null(data)) return()
      
      # Preprocess data for Slingshot
      seurat_data <- as.SingleCellExperiment(data)
      slingshot_data <- slingshot(seurat_data, clusterLabels = "seurat_clusters")
      
      # Plot trajectory results
      output$trajectory_plot <- renderPlot({
        plot(slingshot_data)  # Replace with actual plot code
      })
      
      shinyjs::enable("next_step")
    })
  })
}
