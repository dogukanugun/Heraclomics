# cell_communication.R

library(shiny)
library(Seurat)
library(CellChat)

cell_communication_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 14: Cell Communication"),
    actionButton(ns("run_cellchat"), "Run CellChat"),
    plotOutput(ns("cellchat_plot")),
    actionButton(ns("next_step"), "Continue to Multimodal Analysis", disabled = TRUE)
  )
}

cell_communication_server <- function(id, data_for_cellchat) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$run_cellchat, {
      data <- data_for_cellchat()
      
      if (is.null(data)) return()
      
      # Run CellChat analysis
      cellchat_results <- runCellChat(data)  # Placeholder function
      
      # Plot CellChat results
      output$cellchat_plot <- renderPlot({
        plot(cellchat_results)  # Replace with actual plot code
      })
      
      shinyjs::enable("next_step")
    })
  })
}
