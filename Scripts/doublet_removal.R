# doublet_removal.R

library(shiny)
library(shinyjs)
library(Seurat)
library(DoubletFinder)
library(plotly)

# UI for Doublet Removal step
doublet_removal_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 3: Doublet Removal"),
    
    # Initial plot of the data
    plotlyOutput(ns("initial_plot")),
    
    # Start and Skip buttons
    fluidRow(
      column(6, actionButton(ns("start_doublet_removal"), "Start Doublet Removal")),
      column(6, actionButton(ns("skip_doublet_removal"), "Skip Doublet Removal"))
    ),
    
    # Plot type selection
    fluidRow(
      column(6, selectInput(ns("plot_type"), "Select Plot Type", choices = c("UMAP 2D", "UMAP 3D", "t-SNE 2D", "t-SNE 3D"))),
      column(6, selectInput(ns("label_by"), "Label Cells By", choices = c("Detected Doublet", "Samples")))
    ),
    
    # Adjusters for plot settings
    fluidRow(
      column(3, sliderInput(ns("label_size"), "Label Size", min = 8, max = 20, value = 10)),
      column(3, sliderInput(ns("point_size"), "Point Size", min = 0.1, max = 10, value = 1)),
      column(3, sliderInput(ns("plot_height"), "Plot Height", min = 200, max = 1000, value = 300)),
      column(3, colourInput(ns("point_color"), "Point Color", value = "orange"))
    ),
    
    # Plot output for doublet visualization
    plotlyOutput(ns("doublet_plot")),
    
    # Simulate and Remove Doublets buttons
    fluidRow(
      column(6, actionButton(ns("simulate_doublets"), "Simulate Doublets")),
      column(6, actionButton(ns("remove_doublets"), "Remove Doublets", disabled = TRUE))
    ),
    
    # Continue to next step button
    fluidRow(
      column(12, actionButton(ns("next_step"), "Continue to Clustering", disabled = TRUE))
    )
  )
}

# Server logic for Doublet Removal step
# Server logic for Doublet Removal step
doublet_removal_server <- function(input, output, session, app_state) {
  ns <- session$ns
  
  cleaned_data <- reactiveVal(NULL)
  
  # Initial plot rendering
  output$initial_plot <- renderPlotly({
    req(app_state$data)
    data <- app_state$data
    
    # Ensure UMAP is computed
    if (!"umap" %in% names(data@reductions)) {
      data <- NormalizeData(data)
      data <- FindVariableFeatures(data)
      data <- ScaleData(data)
      data <- RunPCA(data)
      data <- RunUMAP(data, dims = 1:10)
      app_state$data <- data  # Update the app state with the modified data
    }
    
    plot <- DimPlot(data, reduction = "umap")  # You can change to t-SNE if needed
    ggplotly(plot)
  })
  
  observeEvent(input$start_doublet_removal, {
    shinyjs::enable("simulate_doublets")
    shinyjs::enable("remove_doublets")
  })
  
  observeEvent(input$skip_doublet_removal, {
    shinyjs::disable("simulate_doublets")
    shinyjs::disable("remove_doublets")
    shinyjs::enable("next_step")
  })
  
  observeEvent(input$simulate_doublets, {
    req(app_state$data)
    
    data <- app_state$data
    
    # Placeholder for determining optimal pK value for DoubletFinder
    sweep.res <- paramSweep(data, PCs = 1:10, sct = FALSE)  # Update this if needed
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
    
    # Running DoubletFinder with estimated homotypic proportion
    homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
    nExp_poi <- round(0.075 * ncol(data))  # Adjust doublet rate as needed
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    nExp_poi.adj1 <- nExp_poi.adj
    # Identifying doublets
    data <- DoubletFinder::doubletFinder(data, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj1, reuse.pANN = FALSE, sct = FALSE)
    
    # Store the data with doublets identified
    cleaned_data(data)
    
    # Enable the remove doublets button
    shinyjs::enable("remove_doublets")
  })
  
  observeEvent(input$remove_doublets, {
    req(cleaned_data())
    
    data <- cleaned_data()
    
    # Filter out doublets
    cleaned_data <- subset(data, subset = DF.classifications_0.25_0.75 == "Singlet")  # Use correct naming from DoubletFinder
    cleaned_data(cleaned_data)
    
    # Enable the next step button
    shinyjs::enable("next_step")
  })
  
  output$doublet_plot <- renderPlotly({
    req(cleaned_data())
    
    data <- cleaned_data()
    
    plot_type <- input$plot_type
    label_by <- input$label_by
    
    if (plot_type == "UMAP 2D") {
      plot <- DimPlot(data, reduction = "umap", group.by = label_by) +
        theme(text = element_text(size = input$label_size)) +
        geom_point(size = input$point_size, color = input$point_color)
      ggplotly(plot)
    } else if (plot_type == "UMAP 3D") {
      plot <- plot_ly(Embeddings(data, "umap"), x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = data@meta.data[[label_by]], type = "scatter3d", mode = "markers", marker = list(size = input$point_size, color = input$point_color))
      plot
    } else if (plot_type == "t-SNE 2D") {
      plot <- DimPlot(data, reduction = "tsne", group.by = label_by) +
        theme(text = element_text(size = input$label_size)) +
        geom_point(size = input$point_size, color = input$point_color)
      ggplotly(plot)
    } else if (plot_type == "t-SNE 3D") {
      plot <- plot_ly(Embeddings(data, "tsne"), x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, color = data@meta.data[[label_by]], type = "scatter3d", mode = "markers", marker = list(size = input$point_size, color = input$point_color))
      plot
    }
  })
  
  observeEvent(input$next_step, {
    # Navigate to the next step
    updateTabItems(session, "tabs", "clustering")  # Assuming "clustering" is the next tab
  })
  
  # Return the cleaned data for the next step
  return(cleaned_data)
}
