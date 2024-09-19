# integration_dimensionality_reduction.R

library(shiny)
library(Seurat)
library(harmony)

# UI for integration and dimensionality reduction
integration_dimensionality_reduction_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 5: Integration and Dimensionality Reduction"),
    radioButtons(ns("norm_method"), "Normalization Method", 
                 choices = c("LogNormalization", "SCTransform")),
    numericInput(ns("num_variable_features"), "Number of Variable Features", value = 2000, min = 500, max = 5000),
    selectInput(ns("integration_method"), "Integration Method", 
                choices = c("CCA", "RPCA", "Harmony")),
    numericInput(ns("k_anchors"), "k-Anchors (for CCA/RPCA)", value = 5, min = 1, max = 50),
    actionButton(ns("integrate_data"), "Integrate Datasets"),
    verbatimTextOutput(ns("integration_status")),
    actionButton(ns("run_dim_reduction"), "Run Dimensionality Reduction"),
    plotOutput(ns("dim_reduction_plot")),
    actionButton(ns("next_step"), "Continue to Gene Expression", disabled = TRUE)
  )
}

# Server logic for integration and dimensionality reduction
integration_dimensionality_reduction_server <- function(id, data_for_integration, data_for_dim_reduction) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$integrate_data, {
      data_list <- data_for_integration()
      
      if (is.null(data_list) || length(data_list) < 2) {
        output$integration_status <- renderText("Need at least two datasets to integrate.")
        return()
      }
      
      output$integration_status <- renderText("Normalizing datasets...")
      
      # Apply normalization
      if (input$norm_method == "LogNormalization") {
        data_list <- lapply(data_list, function(x) {
          NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
        })
      } else if (input$norm_method == "SCTransform") {
        data_list <- lapply(data_list, SCTransform)
      }
      
      # Select variable features
      data_list <- lapply(data_list, function(x) {
        FindVariableFeatures(x, selection.method = "vst", nfeatures = input$num_variable_features)
      })
      
      # Integration method
      output$integration_status <- renderText("Finding integration anchors...")
      if (input$integration_method %in% c("CCA", "RPCA")) {
        anchors <- FindIntegrationAnchors(object.list = data_list, reduction = tolower(input$integration_method), k.anchor = input$k_anchors)
        integrated_data <- IntegrateData(anchorset = anchors)
      } else if (input$integration_method == "Harmony") {
        integrated_data <- merge(data_list[[1]], y = data_list[-1])
        integrated_data <- RunPCA(integrated_data)
        integrated_data <- RunHarmony(integrated_data, group.by.vars = "orig.ident")
      }
      
      # Scaling data
      integrated_data <- ScaleData(integrated_data)
      
      # Pass integrated data to the next step
      data_for_dim_reduction(integrated_data)
      
      output$integration_status <- renderText("Integration and scaling completed.")
      shinyjs::enable("run_dim_reduction")
    })
    
    observeEvent(input$run_dim_reduction, {
      data <- data_for_dim_reduction()
      
      if (is.null(data)) return()
      
      output$integration_status <- renderText("Running PCA and UMAP...")
      
      # Run PCA
      data <- RunPCA(data, features = VariableFeatures(object = data))
      
      # Run UMAP or t-SNE for dimensionality reduction
      data <- RunUMAP(data, dims = 1:10)
      
      # Plot UMAP results
      output$dim_reduction_plot <- renderPlot({
        DimPlot(data, reduction = "umap", label = TRUE)
      })
      
      shinyjs::enable("next_step")
    })
  })
}
