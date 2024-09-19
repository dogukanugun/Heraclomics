# clustering.R

library(shiny)
library(Seurat)

clustering_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 6: Clustering"),
    numericInput(ns("n_neighbors"), "Number of Nearest Neighbors (k):", value = 20, min = 5, max = 100),
    numericInput(ns("resolution"), "Clustering Resolution:", value = 0.5, min = 0.1, max = 2.0),
    selectInput(ns("algorithm"), "Clustering Algorithm", choices = c("Louvain", "SLM", "Leiden")),
    actionButton(ns("run_clustering"), "Run Clustering"),
    plotOutput(ns("cluster_plot")),
    actionButton(ns("next_step"), "Continue to Gene Expression", disabled = TRUE)
  )
}

clustering_server <- function(id, data_for_clustering) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$run_clustering, {
      data <- data_for_clustering()
      
      if (is.null(data)) return()
      
      # Find Neighbors
      output$clustering_status <- renderText("Finding nearest neighbors...")
      data <- FindNeighbors(data, dims = 1:10, k.param = input$n_neighbors)
      
      # Select the clustering algorithm
      algorithm <- switch(input$algorithm, 
                          "Louvain" = 1,
                          "SLM" = 2,
                          "Leiden" = 3)
      
      # Run Clustering
      output$clustering_status <- renderText("Running clustering...")
      data <- FindClusters(data, resolution = input$resolution, algorithm = algorithm)
      
      # Plot clusters
      output$cluster_plot <- renderPlot({
        DimPlot(data, reduction = "umap", label = TRUE) + ggtitle("Clustering Result")
      })
      
      # Store clustered data for the next step
      shinyjs::enable("next_step")
    })
    
    observeEvent(input$next_step, {
      data_for_dim_reduction(data_for_clustering())  # Pass data to the next step
    })
  })
}
