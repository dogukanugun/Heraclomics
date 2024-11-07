# Load Necessary Libraries
library(shiny)
library(monocle3)
library(SeuratWrappers)
library(plotly)
library(ggplot2)
library(shinyWidgets)
library(shinyalert)
library(dplyr)
library(shinyjs)
library(DT)

# UI Function for Trajectory Analysis Module
trajectoryAnalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyalert(),
    h2("Trajectory Analysis"),
    fluidRow(
      box(
        title = "Trajectory Analysis Overview",
        width = 12,
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        p("This module performs trajectory analysis to determine cellular progressions through pseudotime. 
           You can select specific clusters or samples to analyze and visualize the transition of cells through different biological states."),
        p("The analysis is computationally intensive and may take some time depending on the dataset size.")
      )
    ),
    fluidRow(
      box(
        title = "Settings",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        pickerInput(
          ns("ta_samples"), 
          "Select Sample(s) for Analysis:",
          choices = NULL, 
          multiple = TRUE, 
          options = list(`actions-box` = TRUE, `live-search` = TRUE)
        ),
        pickerInput(
          ns("ta_clusters"), 
          "Select Cluster(s) for Analysis:",
          choices = NULL, 
          multiple = TRUE, 
          options = list(`actions-box` = TRUE, `live-search` = TRUE)
        ),
        numericInput(
          ns("num_cells_ta"), 
          "Number of Cells for Sketching (default: 1000):",
          value = 1000, min = 100, max = 50000, step = 100
        ),
        actionButton(
          ns("run_trajectory"), 
          "Start Trajectory Analysis", 
          icon = icon("play"), 
          class = "btn-success"
        )
      )
    ),
    fluidRow(
      box(
        title = "Pseudotime Visualization",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotlyOutput(ns("trajectory_plot")),
        downloadButton(ns("download_trajectory_plot"), "Download Plot")
      )
    ),
    fluidRow(
      box(
        title = "Gene Expression Along Pseudotime",
        width = 6,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        selectizeInput(
          ns("gene_of_interest"), 
          "Select Gene of Interest:",
          choices = NULL, 
          options = list(
            placeholder = 'Start typing to search for a gene...',
            maxOptions = 1000  # Adjust based on performance
          )
        ),
        plotOutput(ns("gene_expression_plot"))
      )
    )
  )
}

# Server Function for Trajectory Analysis Module
trajectoryAnalysisServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize Shinyalert
    shinyalert::useShinyalert()
    
    # Update Sample and Cluster Choices on App Initialization
    observe({
      req(rv$corrected_data)
      samples <- unique(rv$corrected_data@meta.data$orig.ident)
      clusters <- unique(rv$corrected_data$seurat_clusters)
      
      updatePickerInput(session, "ta_samples", choices = samples)
      updatePickerInput(session, "ta_clusters", choices = clusters)
    })
    
    # Update Gene Choices with Server-Side Selectize for Performance
    observe({
      req(rv$corrected_data)
      all_genes <- rownames(rv$corrected_data)
      
      updateSelectizeInput(
        session, 
        "gene_of_interest", 
        choices = all_genes, 
        server = TRUE,
        options = list(
          placeholder = 'Start typing to search for a gene...',
          maxOptions = 1000
        )
      )
    })
    
    # Perform Trajectory Analysis on Button Click
    observeEvent(input$run_trajectory, {
      req(rv$corrected_data)
      
      tryCatch({
        showModal(modalDialog(
          title = "Trajectory Analysis in Progress",
          "Please wait while the trajectory analysis is being performed...",
          footer = NULL,
          easyClose = FALSE
        ))
        
        withProgress(message = 'Running Trajectory Analysis...', value = 0, {
          # Step 1: Subset Seurat Object Based on Selected Samples and Clusters
          incProgress(0.1, detail = "Subsetting data based on selections...")
          
          cells_to_use <- rv$corrected_data %>%
            subset(
              orig.ident %in% input$ta_samples & 
                seurat_clusters %in% input$ta_clusters
            )
          
          # Diagnostic Check: Ensure Sufficient Cells and Genes
          if (ncol(cells_to_use) < 50) {
            stop("Not enough cells selected for trajectory analysis. Please select more cells.")
          }
          if (nrow(cells_to_use) < 100) {
            stop("Not enough genes selected for trajectory analysis. Please select more genes.")
          }
          
          # Step 2: Convert Seurat Object to Monocle CDS Object
          incProgress(0.2, detail = "Converting Seurat object to Monocle CDS...")
          cds <- as.cell_data_set(cells_to_use)
          print("CDS Conversion Successful:")
          print(cds)
          
          # Step 3: Preprocess Data
          incProgress(0.3, detail = "Preprocessing data with Monocle...")
          cds <- preprocess_cds(cds, num_dim = 30)
          print("Preprocessing Completed:")
          print(dim(reducedDims(cds)$PCA))
          
          # Step 4: Reduce Dimensions using UMAP
          incProgress(0.4, detail = "Reducing dimensions with UMAP...")
          cds <- reduce_dimension(cds)
          print("Dimensionality Reduction Completed:")
          print(dim(reducedDims(cds)$UMAP))
          
          # Step 5: Cluster Cells
          incProgress(0.5, detail = "Clustering cells with Monocle...")
          cds <- cluster_cells(cds, partition_cells = TRUE)
          print("Clustering Completed:")
          print(head(partitions(cds)))
          
          # Step 6: Learn Graph
          incProgress(0.6, detail = "Learning graph structure...")
          cds <- learn_graph(cds)
          print("Graph Learning Completed:")
          
          # Step 7: Order Cells and Calculate Pseudotime
          incProgress(0.7, detail = "Ordering cells and calculating pseudotime...")
          cds <- order_cells(cds)
          print("Ordering Cells Completed:")
          print(head(colData(cds)[, "pseudotime"]))
          
          # Step 8: Store Monocle Object in Reactive Values
          rv$trajectory_data <- cds
          
          # Finalize Progress and Notify Success
          incProgress(1, detail = "Finalizing Trajectory Analysis...")
          removeModal()
          shinyalert("Success", "Trajectory Analysis completed successfully!", type = "success")
        })
        
        # Render Pseudotime Plot
        output$trajectory_plot <- renderPlotly({
          req(rv$trajectory_data)
          cds <- rv$trajectory_data
          
          if (!"pseudotime" %in% colnames(colData(cds))) {
            shinyalert::shinyalert(
              title = "Pseudotime Error",
              text = "Pseudotime data is missing. Please check the analysis steps.",
              type = "error"
            )
            return(NULL)
          }
          
          plot <- plot_cells(cds, color_cells_by = "pseudotime")
          ggplotly(plot)
        })
        
        # Render Gene Expression Plot Along Pseudotime
        output$gene_expression_plot <- renderPlot({
          req(rv$trajectory_data, input$gene_of_interest)
          
          cds <- rv$trajectory_data
          gene <- input$gene_of_interest
          
          if (!(gene %in% rownames(cds))) {
            shinyalert::shinyalert(
              title = "Gene Error",
              text = paste("Gene", gene, "not found in the dataset."),
              type = "error"
            )
            return(NULL)
          }
          
          gene_expr_plot <- plot_genes_in_pseudotime(cds[gene, ]) + 
            ggtitle(paste("Expression of", gene, "Along Pseudotime"))
          
          gene_expr_plot
        })
        
        # Download Pseudotime Plot as HTML
        output$download_trajectory_plot <- downloadHandler(
          filename = function() {
            paste0("trajectory_analysis_", Sys.Date(), ".html")
          },
          content = function(file) {
            req(rv$trajectory_data)
            cds <- rv$trajectory_data
            plot <- plot_cells(cds, color_cells_by = "pseudotime")
            p <- ggplotly(plot)
            htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
          }
        )
        
      })
    })
  })
}
