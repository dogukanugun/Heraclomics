library(Seurat)
library(shiny)
library(shinyjs)
library(dplyr)
library(ggplot2)

# Function to find cluster markers
find_cluster_markers <- function(data) {
  cluster_markers <- FindAllMarkers(data, only.pos = TRUE)
  cluster_markers %>%
    group_by(cluster) %>%
    top_n(10, avg_log2FC) %>%
    pull(gene) %>%
    unique()
}

# Function to calculate module scores
calculate_module_scores <- function(data, top_markers) {
  gene_sets <- list(
    "Cluster1" = top_markers[1:10],
    "Cluster2" = top_markers[11:min(20, length(top_markers))]
  )
  
  AddModuleScore(
    object = data,
    features = gene_sets,
    name = "ModuleScore"
  )
}

# UI for cluster labeling
cluster_labeling_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 7: Cluster Labeling"),
    actionButton(ns("apply_labels"), "Apply Cluster Labels"),
    plotOutput(ns("label_plot")),
    actionButton(ns("next_step"), "Continue to Gene Expression", disabled = TRUE)
  )
}

# Server logic for cluster labeling
cluster_labeling_server <- function(id, data_for_labeling) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Create a reactive value to cache the data
    data_reactive <- reactive({
      data_for_labeling()
    })
    
    # Reactive expression for rendering the plot
    plot_reactive <- reactive({
      data <- data_reactive()
      if (!is.null(data)) {
        DimPlot(data, reduction = "umap", group.by = "ModuleScore1") + ggtitle("Cluster Labels")
      }
    })
    
    output$label_plot <- renderPlot({
      plot_reactive()
    })
    
    observeEvent(input$apply_labels, {
      tryCatch({
        data <- data_reactive()
        
        if (!is.null(data)) {
          top_markers <- find_cluster_markers(data)
          data <- calculate_module_scores(data, top_markers)
          
          # Debugging: Check the result of AddModuleScore
          print(head(data@meta.data))
          
          # Enable next step
          shinyjs::enable("next_step")
        } else {
          showNotification("No data available for cluster labeling.", type = "error")
        }
      }, error = function(e) {
        showNotification(paste("Error in cluster labeling:", e$message), type = "error")
      })
    })
    
    observeEvent(input$next_step, {
      data_for_gene_expression(data_reactive())  # Pass data to the next step
    })
  })
}
