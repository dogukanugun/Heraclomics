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
library(igraph)  # For graph-related functions

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
          options = list(actionsBox = TRUE, liveSearch = TRUE)
        ),
        pickerInput(
          ns("ta_clusters"), 
          "Select Cluster(s) for Analysis:",
          choices = NULL, 
          multiple = TRUE, 
          options = list(actionsBox = TRUE, liveSearch = TRUE)
        ),
        numericInput(
          ns("min_cells_ta"), 
          "Minimum Number of Cells (default: 50):",
          value = 50, min = 10, max = 1000, step = 10
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
        title = "Trajectory Graph Visualization",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotlyOutput(ns("trajectory_plot"))
      )
    ),
    fluidRow(
      box(
        title = "Cluster Visualization",
        width = 6,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotlyOutput(ns("cluster_plot"))
      ),
      box(
        title = "Partition Visualization",
        width = 6,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotlyOutput(ns("partition_plot"))
      )
    ),
    fluidRow(
      box(
        title = "Gene Expression Along Pseudotime",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        selectizeInput(
          ns("gene_of_interest"), 
          "Select Gene of Interest:",
          choices = NULL, 
          options = list(
            placeholder = 'Start typing to search for a gene...',
            maxOptions = 1000
          )
        ),
        plotlyOutput(ns("gene_expression_plot"))
      )
    )
  )
}

# Server Function for Trajectory Analysis Module
trajectoryAnalysisServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Update Sample and Cluster Choices from the Seurat object metadata
    observe({
      req(rv$corrected_data)
      samples <- unique(rv$corrected_data@meta.data$orig.ident)
      clusters <- unique(rv$corrected_data$seurat_clusters)
      
      updatePickerInput(session, "ta_samples", choices = samples)
      updatePickerInput(session, "ta_clusters", choices = clusters)
    })
    
    # Update Gene Choices for Selectize Input
    observe({
      req(rv$corrected_data)
      all_genes <- rownames(rv$corrected_data)
      
      updateSelectizeInput(
        session, 
        "gene_of_interest", 
        choices = all_genes, 
        server = TRUE,
        options = list(placeholder = 'Start typing to search for a gene...')
      )
    })
    
    # Perform Trajectory Analysis
    observeEvent(input$run_trajectory, {
      req(rv$corrected_data)
      
      showModal(modalDialog(
        title = "Trajectory Analysis in Progress",
        "Please wait while the trajectory analysis is being performed...",
        footer = NULL,
        easyClose = FALSE
      ))
      
      tryCatch({
        withProgress(message = "Running Trajectory Analysis", value = 0, {
          # Step 1: Subset data
          incProgress(0.1, detail = "Subsetting data based on selections...")
          cells_to_use <- subset(
            rv$corrected_data,
            subset = orig.ident %in% input$ta_samples & seurat_clusters %in% input$ta_clusters
          )
          
          num_cells <- ncol(cells_to_use)
          if (num_cells < input$min_cells_ta || nrow(cells_to_use) < 100) {
            removeModal()
            shinyalert("Error", "Not enough cells or genes for trajectory analysis.", type = "error")
            return()
          }
          
          # Step 2: Convert Seurat object to a Monocle cell_data_set (CDS)
          incProgress(0.2, detail = "Converting Seurat object to Monocle CDS...")
          cds <- as.cell_data_set(cells_to_use)
          
          # Remove Seurat cluster labels to allow Monocle3 to recalc its own clustering
          if ("seurat_clusters" %in% colnames(colData(cds))) {
            colData(cds)$seurat_clusters <- NULL
          }
          
          # Step 3: Preprocess the CDS (normalization, variance stabilization, etc.)
          incProgress(0.3, detail = "Preprocessing data with Monocle...")
          cds <- preprocess_cds(cds, num_dim = 50)
          
          # Step 4: Dimensionality reduction using UMAP (explicitly set method)
          incProgress(0.4, detail = "Reducing dimensions using UMAP...")
          cds <- reduce_dimension(cds, reduction_method = "UMAP")
          
          # Step 5: Cluster cells using Monocle3 so that proper partitions are calculated
          incProgress(0.5, detail = "Clustering cells...")
          cds <- cluster_cells(cds, partition_qval = 0.05)
          
          # Verify clusters (using Monocle???s computed clusters)
          cluster_counts <- table(clusters(cds))
          print("Cluster sizes:")
          print(cluster_counts)
          
          if (length(cluster_counts) < 2) {
            removeModal()
            shinyalert("Error", "Clustering did not produce enough clusters. Try including more cells or adjusting parameters.", type = "error")
            return()
          }
          
          # Step 6: Learn the principal graph (trajectory structure)
          incProgress(0.6, detail = "Learning graph structure...")
          cds <- learn_graph(cds, use_partition = TRUE)
          
          # Step 7: Order cells and calculate pseudotime
          incProgress(0.7, detail = "Ordering cells and calculating pseudotime...")
          
          # Use Monocle3's clusters (not seurat_clusters) to compute the starting nodes
          get_earliest_principal_nodes <- function(cds) {
            partitions_list <- unique(partitions(cds))
            root_pr_nodes <- c()
            
            for (partition_id in partitions_list) {
              # Get indices of cells in the current partition
              cells_in_partition <- which(partitions(cds) == partition_id)
              if (length(cells_in_partition) == 0) next
              
              cds_subset <- cds[, cells_in_partition]
              
              # Use the cluster (computed by Monocle3) with the most cells in the partition as the starting cluster
              cluster_counts <- table(colData(cds_subset)$clusters)
              if(length(cluster_counts) == 0) next
              starting_cluster <- names(cluster_counts)[which.max(cluster_counts)]
              
              # Identify the cell indices in that starting cluster
              cell_ids <- which(colData(cds_subset)$clusters == starting_cluster)
              
              # Extract the projected vertices for these cells from the learned principal graph
              closest_vertex <- cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
              closest_vertex <- closest_vertex[cell_ids, , drop = FALSE]
              colnames(closest_vertex) <- "principal_point"
              closest_vertex_df <- as.data.frame(closest_vertex)
              
              counts <- table(closest_vertex_df$principal_point)
              if (length(counts) == 0) next
              root_pr_node_number <- names(counts)[which.max(counts)]
              # Adjust naming as needed; here we prepend "Y_" as per the default naming in Monocle3
              root_pr_node <- paste0("Y_", root_pr_node_number)
              
              root_pr_nodes <- c(root_pr_nodes, root_pr_node)
            }
            
            print("Identified root principal nodes for all partitions:")
            print(root_pr_nodes)
            
            return(unique(root_pr_nodes))
          }
          
          # Get root principal nodes from the CDS
          root_pr_nodes <- get_earliest_principal_nodes(cds)
          
          # Ensure that the identified root nodes exist in the principal graph
          graph_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
          root_pr_nodes <- root_pr_nodes[root_pr_nodes %in% graph_nodes]
          
          if (length(root_pr_nodes) == 0) {
            removeModal()
            shinyalert("Error", "Failed to identify valid root principal nodes in the principal graph.", type = "error")
            return()
          }
          
          # Order cells based on the computed root principal nodes
          cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
          
          # Check if pseudotime was successfully calculated
          if (!"pseudotime" %in% colnames(colData(cds))) {
            removeModal()
            shinyalert("Error", "Pseudotime calculation failed after ordering cells.", type = "error")
            return()
          }
          
          # Save the CDS object (with pseudotime) into reactiveValues for later use
          rv$trajectory_data <- cds
          incProgress(1, detail = "Trajectory analysis completed.")
        })
        
        removeModal()
        shinyalert("Success", "Trajectory analysis completed successfully!", type = "success")
        
      }, error = function(e) {
        removeModal()
        shinyalert("Error", paste("An error occurred:", e$message), type = "error")
      })
    })
    
    # Render Interactive Trajectory Plot
    output$trajectory_plot <- renderPlotly({
      req(rv$trajectory_data)
      plot <- plot_cells(
        rv$trajectory_data,
        color_cells_by = "pseudotime",
        label_cell_groups = FALSE,
        label_leaves = TRUE,
        label_branch_points = TRUE,
        graph_label_size = 5,
        show_trajectory_graph = TRUE
      )
      ggplotly(plot, tooltip = c("cell_id", "pseudotime"))
    })
    
    # Render Interactive Cluster Plot
    output$cluster_plot <- renderPlotly({
      req(rv$trajectory_data)
      plot <- plot_cells(
        rv$trajectory_data,
        color_cells_by = "cluster",
        label_cell_groups = TRUE,
        label_leaves = TRUE,
        label_branch_points = TRUE,
        graph_label_size = 5
      )
      ggplotly(plot, tooltip = c("cell_id", "cluster"))
    })
    
    # Render Interactive Partition Plot
    output$partition_plot <- renderPlotly({
      req(rv$trajectory_data)
      plot <- plot_cells(
        rv$trajectory_data,
        color_cells_by = "partition",
        label_cell_groups = TRUE,
        label_leaves = TRUE,
        label_branch_points = TRUE,
        graph_label_size = 5
      )
      ggplotly(plot, tooltip = c("cell_id", "partition"))
    })
    
    # Render Interactive Gene Expression Plot
    output$gene_expression_plot <- renderPlotly({
      req(rv$trajectory_data, input$gene_of_interest)
      gene <- input$gene_of_interest
      
      if (!(gene %in% rownames(rv$trajectory_data))) {
        shinyalert("Error", paste("Gene", gene, "not found in the dataset."), type = "error")
        return(NULL)
      }
      
      plot <- plot_genes_in_pseudotime(
        rv$trajectory_data[gene, ], 
        color_cells_by = "pseudotime"
      ) +
        ggtitle(paste("Expression of", gene, "Along Pseudotime"))
      
      ggplotly(plot, tooltip = c("pseudotime", "expression"))
    })
    
    # Download Pseudotime Plot
    output$download_trajectory_plot <- downloadHandler(
      filename = function() {
        paste0("trajectory_analysis_", Sys.Date(), ".html")
      },
      content = function(file) {
        plot <- plot_cells(
          rv$trajectory_data,
          color_cells_by = "pseudotime",
          label_cell_groups = FALSE,
          label_leaves = TRUE,
          label_branch_points = TRUE,
          graph_label_size = 5,
          show_trajectory_graph = TRUE
        )
        htmlwidgets::saveWidget(ggplotly(plot), file, selfcontained = TRUE)
      }
    )
  })
}
