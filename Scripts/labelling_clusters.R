# Load necessary libraries
library(shiny)
library(Seurat)
library(SingleR)
library(celldex)
library(shinyWidgets)
library(dplyr)
library(HGNChelper)
library(openxlsx)
library(ggplot2)  # For plotting
library(DT)
library(SingleCellExperiment)  # Added for SingleCellExperiment functions

# Load scType functions from the source
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# UI Function for Labelling Clusters Module
labellingClustersUI <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Cluster Labelling"),
    fluidRow(
      box(
        title = "Cluster Marker Gene Analysis",
        width = 12,
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        p("Analyze marker genes for each cluster and assign cell types using custom labels or reference datasets."),
        selectInput(ns("labelling_method"), "Select Labelling Method:",
                    choices = list("Custom" = "custom", "scType Reference" = "sctype", "SingleR Reference" = "singleR")),
        uiOutput(ns("reference_dataset_selector")),
        actionButton(ns("label_clusters"), "Label Clusters", class = "btn-success")
      )
    ),
    fluidRow(
      box(
        title = "Cluster Labelling Results",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        dataTableOutput(ns("cluster_labels_table")),
        downloadButton(ns("download_labels"), "Download Cluster Labels")
      )
    ),
    fluidRow(
      box(
        title = "UMAP Plot with Cluster Labels",
        width = 12,
        status = "warning",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotOutput(ns("umap_plot"))
      )
    )
  )
}

# Server Function for Labelling Clusters Module
labellingClustersServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Dynamically generate the reference dataset selector based on the selected labelling method
    output$reference_dataset_selector <- renderUI({
      if (input$labelling_method == "sctype") {
        selectizeInput(ns("sctype_ref"), "Select scType Reference Dataset:",
                       choices = c("Brain", "Immune", "Pancreas"),
                       selected = "Brain",
                       options = list(placeholder = 'Select a scType Reference Dataset'))
      } else if (input$labelling_method == "singleR") {
        selectizeInput(ns("singleR_ref"), "Select SingleR Reference Dataset:",
                       choices = list(
                         "Human Primary Cell Atlas" = "HPCAData",
                         "Blueprint Encode" = "BlueprintEncodeData",
                         "Database Immune Cell Expression" = "DatabaseImmuneCellExpressionData",
                         "ImmGen" = "ImmGenData",
                         "Monaco Immune Data" = "MonacoImmuneData",
                         "Novershtern Hematopoietic Data" = "NovershternHematopoieticData",
                         "Mouse RNAseq" = "MouseRNAseqData",
                         "ENCODE and Blueprint Epigenomics" = "ENCODEBlueprintData"
                       ),
                       selected = "HPCAData",
                       options = list(placeholder = 'Select a SingleR Reference Dataset'),
                       multiple = FALSE)
      }
    })
    
    observeEvent(input$label_clusters, {
      # Ensure that corrected_data and assay_used are available
      req(rv$corrected_data, rv$assay_used)
      print("Corrected data and assay are available.")
      
      assay_used <- rv$assay_used
      
      # Check if 'seurat_clusters' exists in metadata
      if (!"seurat_clusters" %in% colnames(rv$corrected_data@meta.data)) {
        showModal(modalDialog(
          title = "Error",
          "The Seurat object does not contain 'seurat_clusters'. Please ensure clustering has been performed.",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        return(NULL)
      }
      
      print("Seurat clusters are available.")
      
      # Set cluster identities
      Idents(rv$corrected_data) <- "seurat_clusters"
      showModal(modalDialog(
        title = "Labelling Clusters",
        "Please wait while clusters are being labelled...",
        footer = NULL
      ))
      
      withProgress(message = 'Labelling Clusters...', value = 0, {
        tryCatch({
          if (input$labelling_method == "singleR") {
            incProgress(0.2, detail = "Using SingleR reference dataset...")
            print("Using SingleR reference dataset...")
            ref_data <- switch(input$singleR_ref,
                               "HPCAData" = celldex::HumanPrimaryCellAtlasData(),
                               "BlueprintEncodeData" = celldex::BlueprintEncodeData(),
                               "DatabaseImmuneCellExpressionData" = celldex::DatabaseImmuneCellExpressionData(),
                               "ImmGenData" = celldex::ImmGenData(),
                               "MonacoImmuneData" = celldex::MonacoImmuneData(),
                               "NovershternHematopoieticData" = celldex::NovershternHematopoieticData(),
                               "MouseRNAseqData" = celldex::MouseRNAseqData(),
                               "ENCODEBlueprintData" = celldex::ENCODEBlueprintData(),
                               stop("Invalid SingleR reference dataset selected."))
            
            # Convert Seurat object to SingleCellExperiment
            sce <- as.SingleCellExperiment(rv$corrected_data, assay = assay_used)
            
            # Ensure gene overlap
            common_genes <- intersect(rownames(ref_data), rownames(sce))
            print(paste("Number of common genes:", length(common_genes)))
            if (length(common_genes) < 100) {
              stop("Insufficient gene overlap between reference and dataset. Please check the datasets.")
            }
            sce <- sce[common_genes, ]
            ref_data <- ref_data[common_genes, ]
            
            # Run SingleR on clusters
            pred_cluster <- SingleR(test = sce, ref = ref_data, labels = ref_data$label.main, clusters = sce$seurat_clusters)
            incProgress(0.3, detail = "SingleR predictions completed.")
            print("SingleR predictions completed.")
            print(pred_cluster$labels)
            
            # Assign cluster labels
            cluster_labels <- pred_cluster$labels
            names(cluster_labels) <- rownames(pred_cluster)
            rv$cluster_labels <- cluster_labels
            print("SingleR cluster labelling completed.")
            
            # Map cluster labels back to cells
            cell_types <- cluster_labels[as.character(rv$corrected_data$seurat_clusters)]
            names(cell_types) <- Cells(rv$corrected_data)
            
            # Update Seurat object metadata
            rv$corrected_data <- AddMetaData(rv$corrected_data, metadata = cell_types, col.name = "cell_type")
            incProgress(0.1, detail = "Metadata assignment completed.")
            print("Metadata assignment completed.")
            
            # Render the cluster labels table
            output$cluster_labels_table <- renderDataTable({
              req(rv$cluster_labels)
              datatable(data.frame(Cluster = names(rv$cluster_labels), Label = rv$cluster_labels, row.names = NULL))
            })
            
            # Render the UMAP plot
            output$umap_plot <- renderPlot({
              req(rv$corrected_data$cell_type)
              if (!"umap" %in% names(rv$corrected_data@reductions)) {
                showModal(modalDialog(
                  title = "Error",
                  "UMAP reduction not found. Please run UMAP before plotting.",
                  easyClose = TRUE,
                  footer = modalButton("Close")
                ))
                return(NULL)
              }
              DimPlot(rv$corrected_data, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
                ggtitle("UMAP Plot Colored by Cell Type Labels")
            })
            
            incProgress(1, detail = "Completed")
            print("Labelling process completed successfully.")
            
          } else if (input$labelling_method == "custom") {
            # Implement your custom labelling logic here
            cluster_ids <- levels(Idents(rv$corrected_data))
            # For example, you can prompt the user to input labels or assign default labels
            cluster_labels <- setNames(rep("Custom_Label", length(cluster_ids)), cluster_ids)
            rv$cluster_labels <- cluster_labels
            incProgress(0.3, detail = "Custom cluster labels assigned.")
            print("Custom cluster labels assigned.")
            
            # Map cluster labels back to cells
            cell_types <- cluster_labels[as.character(rv$corrected_data$seurat_clusters)]
            names(cell_types) <- Cells(rv$corrected_data)
            
            # Update Seurat object metadata
            rv$corrected_data <- AddMetaData(rv$corrected_data, metadata = cell_types, col.name = "cell_type")
            incProgress(0.1, detail = "Metadata assignment completed.")
            print("Metadata assignment completed.")
            
            # Render the cluster labels table
            output$cluster_labels_table <- renderDataTable({
              req(rv$cluster_labels)
              datatable(data.frame(Cluster = names(rv$cluster_labels), Label = rv$cluster_labels, row.names = NULL))
            })
            
            # Render the UMAP plot
            output$umap_plot <- renderPlot({
              req(rv$corrected_data$cell_type)
              if (!"umap" %in% names(rv$corrected_data@reductions)) {
                showModal(modalDialog(
                  title = "Error",
                  "UMAP reduction not found. Please run UMAP before plotting.",
                  easyClose = TRUE,
                  footer = modalButton("Close")
                ))
                return(NULL)
              }
              DimPlot(rv$corrected_data, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
                ggtitle("UMAP Plot Colored by Cell Type Labels")
            })
            
            incProgress(1, detail = "Completed")
            print("Custom labelling process completed successfully.")
            
          } else if (input$labelling_method == "sctype") {
            # Implement scType labelling logic here
            incProgress(1, detail = "scType labelling method not implemented yet.")
            print("scType labelling method not implemented.")
          }
          
        }, error = function(e) {
          removeModal()
          showModal(modalDialog(
            title = "Error",
            paste("An error occurred during labelling:", e$message),
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
          incProgress(1, detail = "Error occurred.")
          print(paste("Error during labelling:", e$message))
        })
      })
      
      removeModal()
      
      if (!is.null(rv$cluster_labels)) {
        showModal(modalDialog(
          title = "Labelling Complete",
          "The clusters have been successfully labelled.",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        print("Clusters have been successfully labelled.")
      }
    })
    
    # Download Handler for Cluster Labels
    output$download_labels <- downloadHandler(
      filename = function() {
        paste("cluster_labels_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        labels_df <- isolate({
          req(rv$cluster_labels)
          data.frame(Cluster = names(rv$cluster_labels), Label = rv$cluster_labels, row.names = NULL)
        })
        write.csv(labels_df, file, row.names = FALSE)
      }
    )
  })
}
