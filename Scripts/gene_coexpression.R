# Load necessary libraries
library(shiny)
library(shinyWidgets)
library(Seurat)
library(MEGENA)
library(plotly)
library(igraph)
library(tibble)
library(DT)
library(clusterProfiler)
library(org.Hs.eg.db)  # Replace with appropriate organism database if not human
library(ggplot2)
library(enrichplot)
library(dplyr)

# UI Function for Gene Co-Expression Networks Module
geneCoexpressionUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    h2("Gene Co-Expression Networks"),
    fluidRow(
      actionButton(ns("startMEGENA"), "Start MEGENA Analysis", class = "btn btn-primary")
    ),
    hr(),
    tabsetPanel(
      id = ns("tabs"),
      tabPanel(
        "Network Plot",
        fluidRow(
          box(
            title = "Network Plot",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput(ns("network_plot"))
          )
        )
      ),
      tabPanel(
        "Module Summary",
        fluidRow(
          box(
            title = "Module Summary Table",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            DT::dataTableOutput(ns("module_summary_table"))
          )
        )
      ),
      tabPanel(
        "Hub Genes",
        fluidRow(
          box(
            title = "Hub Genes Information",
            width = 12,
            status = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            DT::dataTableOutput(ns("hub_genes_table"))
          )
        )
      ),
      tabPanel(
        "Gene Module Membership",
        fluidRow(
          box(
            title = "Gene-Module Membership",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            DT::dataTableOutput(ns("gene_module_table"))
          )
        )
      ),
      tabPanel(
        "Module Visualization",
        fluidRow(
          box(
            title = "Visualize Individual Modules",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            selectInput(ns("module_to_plot"), "Select Module:", choices = NULL),
            plotlyOutput(ns("module_plot"))
          )
        )
      ),
      tabPanel(
        "Enrichment Analysis",
        fluidRow(
          box(
            title = "Functional Enrichment Analysis",
            width = 12,
            status = "danger",
            solidHeader = TRUE,
            collapsible = TRUE,
            selectInput(ns("selected_module"), "Select Module:", choices = NULL),
            DT::dataTableOutput(ns("enrichment_table")),
            plotlyOutput(ns("enrichment_barplot")),
            plotlyOutput(ns("enrichment_dotplot")),
            plotlyOutput(ns("enrichment_cnetplot"))
          )
        )
      ),
      tabPanel(
        "Download Results",
        fluidRow(
          box(
            title = "Downloadable Results",
            width = 12,
            status = "success",
            solidHeader = TRUE,
            collapsible = TRUE,
            downloadButton(ns("download_module_summary"), "Download Module Summary"),
            downloadButton(ns("download_hub_genes"), "Download Hub Genes"),
            downloadButton(ns("download_gene_module"), "Download Gene-Module Membership")
          )
        )
      )
    )
  )
}

# Server Function for Gene Co-Expression Networks
geneCoexpressionServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive expression to store MEGENA results
    megena_result <- reactiveVal(NULL)
    
    # Define a minimum number of cells required per cluster
    min_cells_required <- 10
    
    # Event to start MEGENA process
    observeEvent(input$startMEGENA, {
      cat("Event: Start MEGENA Analysis Initiated\n")
      
      # Check if rv$corrected_data is available before proceeding
      if (is.null(rv$corrected_data) || length(rv$corrected_data) == 0) {
        showNotification("No corrected data available. Please complete data correction before proceeding.", type = "error")
        cat("Error: No corrected data available.\n")
        return(NULL)  # Stop further processing if no data
      }
      
      # Ensure Idents are set to clusters
      if (!"seurat_clusters" %in% colnames(rv$corrected_data@meta.data)) {
        showNotification("Cluster information not found in the Seurat object.", type = "error")
        cat("Error: 'seurat_clusters' not found in metadata.\n")
        return(NULL)
      }
      
      Idents(rv$corrected_data) <- "seurat_clusters"  # Set to your cluster column
      cat("Set Idents to 'seurat_clusters'.\n")
      
      # Show modal dialog for MEGENA settings
      showModal(modalDialog(
        title = strong("MEGENA Settings"),
        selectInput(ns("MEGENA_selectmethod"), "Select input method:",
                    choices = c("Variable Genes", "Differentially Expressed Genes")),
        numericInput(ns("variablegenesMEGENA"), "Pick number of variable genes:", value = 500, min = 1, step = 50),
        selectInput(ns("MEGENA_selectcluster"), "Select cells by:",
                    choices = c("Seurat Clusters", "Labelled Cells")),
        # Use conditionalPanel to show/hide metadata column selection
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Labelled Cells'", ns("MEGENA_selectcluster")),
          selectInput(ns("MEGENA_cluster_column"), "Select Metadata Column:", choices = NULL)
        ),
        selectizeInput(ns("MEGENA_filterclusters"), "Select Clusters:", choices = NULL, 
                       multiple = TRUE, options = list(placeholder = 'Select Clusters')),
        numericInput(ns("MEGENApval"), "Module p-value threshold:", value = 0.05, min = 0, max = 1),
        numericInput(ns("MEGENAmin_size"), "Minimum module size:", value = 10, min = 2, max = 100),
        actionButton(ns("startbuttonMEGENA"), "Apply", class = "btn btn-lg btn-primary"),
        easyClose = TRUE, footer = NULL
      ))
      
      # Update cluster selection options based on the user's input
      observeEvent(input$MEGENA_selectcluster, {
        cat("Event: MEGENA_selectcluster changed to", input$MEGENA_selectcluster, "\n")
        if (input$MEGENA_selectcluster == "Seurat Clusters") {
          # Use "seurat_clusters" as the cluster column
          Idents(rv$corrected_data) <- "seurat_clusters"
          cluster_table <- table(Idents(rv$corrected_data))
          eligible_clusters <- names(cluster_table[cluster_table >= min_cells_required])
          updateSelectizeInput(session, "MEGENA_filterclusters", choices = eligible_clusters, server = TRUE)
          cat("Updated 'MEGENA_filterclusters' with Seurat clusters.\n")
        } else {
          # For Labelled Cells, allow the user to select the metadata column
          metadata_columns <- colnames(rv$corrected_data@meta.data)
          updateSelectInput(session, "MEGENA_cluster_column", choices = metadata_columns)
          cat("Updated 'MEGENA_cluster_column' with metadata columns.\n")
          
          # Update clusters based on selected metadata column
          observeEvent(input$MEGENA_cluster_column, {
            selected_column <- input$MEGENA_cluster_column
            cat("Event: MEGENA_cluster_column changed to", selected_column, "\n")
            cluster_values <- rv$corrected_data@meta.data[[selected_column]]
            # Filter clusters based on minimum cell count
            cluster_table <- table(cluster_values)
            eligible_clusters <- names(cluster_table[cluster_table >= min_cells_required])
            updateSelectizeInput(session, "MEGENA_filterclusters", choices = eligible_clusters, server = TRUE)
            cat("Updated 'MEGENA_filterclusters' with values from", selected_column, "\n")
          }, ignoreInit = TRUE)
        }
      }, ignoreInit = TRUE)
    })
    
    # Event to apply MEGENA settings and run the analysis
    observeEvent(input$startbuttonMEGENA, {
      cat("Event: Apply MEGENA Settings and Start Analysis\n")
      
      # Check again if rv$corrected_data is available when the user applies settings
      if (is.null(rv$corrected_data) || length(rv$corrected_data) == 0) {
        showNotification("No corrected data available for MEGENA analysis.", type = "error")
        cat("Error: No corrected data available for MEGENA analysis.\n")
        return(NULL)  # Stop if no data
      }
      
      # Proceed with MEGENA analysis
      removeModal()  # Close the modal
      cat("Modal closed. Starting MEGENA analysis.\n")
      
      withProgress(message = "Running MEGENA Analysis", value = 0, {
        incProgress(1/10, detail = "Preparing data")
        cat("Step 1/10: Preparing data for MEGENA.\n")
        
        # Subset based on input method
        gene_set <- switch(input$MEGENA_selectmethod,
                           "Variable Genes" = {
                             incProgress(2/10, detail = "Using variable genes")
                             cat("Step 2/10: Selecting Variable Genes.\n")
                             # Ensure variable genes are ordered by variance
                             variable_features <- VariableFeatures(rv$corrected_data)
                             if (length(variable_features) == 0) {
                               showNotification("No variable features found. Please run FindVariableFeatures() first.", type = "error")
                               cat("Error: No variable features found.\n")
                               return(NULL)
                             }
                             gene_data <- GetAssayData(rv$corrected_data, slot = "data")[variable_features, , drop = FALSE]
                             gene_vars <- apply(gene_data, 1, var, na.rm = TRUE)
                             # Handle cases where gene_vars may have NAs
                             gene_vars[is.na(gene_vars)] <- 0
                             top_variable_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(input$variablegenesMEGENA, length(variable_features))]
                             cat("Selected", length(top_variable_genes), "variable genes.\n")
                             showNotification(paste("Selected", length(top_variable_genes), "variable genes."), type = "message")
                             top_variable_genes
                           },
                           "Differentially Expressed Genes" = {
                             if (length(input$MEGENA_filterclusters) < 2) {
                               showNotification("Please select at least two clusters for differential expression analysis.", type = "error")
                               cat("Error: Less than two clusters selected for DEGs.\n")
                               return(NULL)
                             }
                             # Inform the user that only the first two selected clusters will be used
                             if (length(input$MEGENA_filterclusters) > 2) {
                               showNotification("Only the first two selected clusters will be used for differential expression analysis.", type = "warning")
                               cat("Warning: More than two clusters selected. Only the first two will be used for DEGs.\n")
                             }
                             incProgress(2/10, detail = "Selecting DEGs")
                             cat("Step 2/10: Selecting Differentially Expressed Genes.\n")
                             ident.1 <- input$MEGENA_filterclusters[1]
                             ident.2 <- input$MEGENA_filterclusters[2]
                             
                             # Set Idents based on selected cluster column
                             if (input$MEGENA_selectcluster == "Seurat Clusters") {
                               Idents(rv$corrected_data) <- "seurat_clusters"
                             } else {
                               Idents(rv$corrected_data) <- input$MEGENA_cluster_column
                             }
                             cat("Set Idents based on selection:", input$MEGENA_selectcluster, "\n")
                             
                             markers <- FindMarkers(rv$corrected_data, ident.1 = ident.1, ident.2 = ident.2)
                             if (nrow(markers) == 0) {
                               showNotification("No differentially expressed genes found between the selected clusters.", type = "error")
                               cat("Error: No DEGs found between clusters", ident.1, "and", ident.2, ".\n")
                               return(NULL)
                             }
                             cat("Found", nrow(markers), "differentially expressed genes.\n")
                             showNotification(paste("Found", nrow(markers), "differentially expressed genes."), type = "message")
                             rownames(markers)
                           })
        
        if (is.null(gene_set) || length(gene_set) == 0) {
          showNotification("No genes found for the selected input method.", type = "error")
          cat("Error: No genes found for the selected input method.\n")
          return(NULL)  # Stop if no genes found
        }
        
        # Filter expression matrix based on selected clusters and genes
        if (input$MEGENA_selectcluster == "Seurat Clusters") {
          Idents(rv$corrected_data) <- "seurat_clusters"
          selected_clusters <- input$MEGENA_filterclusters
          cat("Selected clusters (Seurat Clusters):", paste(selected_clusters, collapse = ", "), "\n")
        } else {
          Idents(rv$corrected_data) <- input$MEGENA_cluster_column
          selected_clusters <- input$MEGENA_filterclusters
          cat("Selected clusters (Labelled Cells -", input$MEGENA_cluster_column, "):", paste(selected_clusters, collapse = ", "), "\n")
        }
        
        cells_to_use <- WhichCells(rv$corrected_data, idents = selected_clusters)
        cat("Number of cells selected:", length(cells_to_use), "\n")
        
        # Enforce a minimum number of cells
        if (length(cells_to_use) < min_cells_required) {
          showNotification(paste("Not enough cells selected for MEGENA analysis. Please select at least", min_cells_required, "cells."), type = "error")
          cat("Error: Only", length(cells_to_use), "cells selected. Minimum required is", min_cells_required, "\n")
          return(NULL)
        }
        
        expr_matrix <- GetAssayData(rv$corrected_data, slot = "data")[gene_set, cells_to_use, drop = FALSE]
        cat("Expression matrix dimensions:", dim(expr_matrix), "\n")  # Should be genes x samples
        
        # Remove genes with zero variance
        gene_vars <- apply(expr_matrix, 1, var)
        zero_var_genes <- which(gene_vars == 0)
        if (length(zero_var_genes) > 0) {
          expr_matrix <- expr_matrix[-zero_var_genes, , drop = FALSE]
          cat("Removed", length(zero_var_genes), "genes with zero variance.\n")
        }
        
        cat("Number of genes after variance filtering:", nrow(expr_matrix), "\n")
        
        if (nrow(expr_matrix) < 2) {
          showNotification("Not enough genes after filtering for MEGENA analysis.", type = "error")
          cat("Error: Less than two genes after variance filtering.\n")
          return(NULL)
        }
        
        if (ncol(expr_matrix) < 2) {
          showNotification("Not enough cells in the selected clusters for MEGENA analysis.", type = "error")
          cat("Error: Less than two cells in the selected clusters.\n")
          return(NULL)
        }
        
        incProgress(5/10, detail = "Calculating Correlations")
        cat("Step 5/10: Calculating Correlations.\n")
        
        # Set correlation method to Pearson
        correlation_method <- tolower("pearson")
        cat("Correlation method selected:", correlation_method, "\n")
        
        # Ensure the expression matrix is numeric
        expr_matrix_numeric <- as.matrix(expr_matrix)
        cat("Expression matrix dimensions:", dim(expr_matrix_numeric), "\n")  # Should be genes x samples
        
        if (!is.numeric(expr_matrix_numeric)) {
          showNotification("Expression matrix must be numeric.", type = "error")
          cat("Error: Expression matrix is not numeric.\n")
          return(NULL)
        }
        
        # Check for invalid values
        if (any(is.na(expr_matrix_numeric)) || any(is.nan(expr_matrix_numeric)) || any(is.infinite(expr_matrix_numeric))) {
          showNotification("Expression matrix contains NA, NaN, or infinite values.", type = "error")
          cat("Error: Expression matrix contains invalid values.\n")
          return(NULL)
        }
        
        # Remove genes with zero variance (already done above, but double-check)
        gene_vars <- apply(expr_matrix_numeric, 1, var)
        zero_var_genes <- which(gene_vars == 0)
        if (length(zero_var_genes) > 0) {
          expr_matrix_numeric <- expr_matrix_numeric[-zero_var_genes, , drop = FALSE]
          cat("Removed", length(zero_var_genes), "genes with zero variance after numeric conversion.\n")
        }
        
        # Ensure sufficient number of genes and samples
        if (nrow(expr_matrix_numeric) < 2 || ncol(expr_matrix_numeric) < 2) {
          showNotification("Not enough genes or samples after filtering.", type = "error")
          cat("Error: Insufficient data for correlation calculation.\n")
          return(NULL)
        }
        
        # Calculate correlations using calculate.correlation from MEGENA
        rho <- tryCatch({
          calculate.correlation(
            datExpr = expr_matrix_numeric,
            doPerm = 10,  # Use a small number of permutations
            method = correlation_method,
            FDR.cutoff = 0.05,
            n.increment = 100,
            is.signed = FALSE,
            output.permFDR = TRUE,
            output.corTable = TRUE,
            saveto = NULL,
            doPar = FALSE,
            num.cores = 1
          )
        }, error = function(e) {
          showNotification(paste("Error in correlation calculation:", e$message), type = "error")
          cat("Error: Correlation calculation failed -", e$message, "\n")
          return(NULL)
        })
        
        if (is.null(rho) || nrow(rho) == 0) {
          showNotification("Correlation matrix is empty or invalid.", type = "error")
          cat("Error: Correlation matrix is empty or invalid.\n")
          return(NULL)
        }
        
        # Output column names for debugging
        cat("Column names in rho:", paste(colnames(rho), collapse = ", "), "\n")
        
        # Report the number of significant gene pairs
        num_gene_pairs <- nrow(rho)
        showNotification(paste("Calculated", num_gene_pairs, "gene pairs with significant correlations."), type = "message")
        cat("Calculated", num_gene_pairs, "gene pairs with significant correlations.\n")
        
        incProgress(7/10, detail = "Constructing Planar Filtered Network")
        cat("Step 7/10: Constructing Planar Filtered Network.\n")
        
        # Construct Planar Filtered Network using calculate.PFN
        pfn <- tryCatch({
          # Ensure that rho has columns "row", "col", "rho"
          required_cols <- c("row", "col", "rho")
          if (!all(required_cols %in% colnames(rho))) {
            stop("rho does not contain the required columns: row, col, rho")
          }
          
          calculate.PFN(
            rho[, required_cols],
            doPar = FALSE,
            num.cores = 1,
            keep.track = FALSE
          )
        }, error = function(e) {
          showNotification(paste("Error in PFN construction:", e$message), type = "error")
          cat("Error: PFN construction failed -", e$message, "\n")
          return(NULL)
        })
        
        if (is.null(pfn) || nrow(pfn) == 0) {
          showNotification("Planar Filtered Network construction failed.", type = "error")
          cat("Error: PFN construction returned empty or NULL.\n")
          return(NULL)
        }
        
        # Convert adjacency edgelist to igraph object
        cat("Step 8/10: Converting PFN to igraph object.\n")
        pfn_igraph <- tryCatch({
          graph.data.frame(pfn, directed = FALSE)
        }, error = function(e) {
          showNotification(paste("Conversion to igraph object failed:", e$message), type = "error")
          cat("Error: Conversion to igraph object failed -", e$message, "\n")
          return(NULL)
        })
        
        # Validate the igraph object
        if (is.null(pfn_igraph) || !is_igraph(pfn_igraph)) {
          showNotification("Conversion to igraph object failed.", type = "error")
          cat("Error: Invalid igraph object after conversion.\n")
          return(NULL)
        }
        
        # Check if pfn_igraph has any edges
        if (ecount(pfn_igraph) == 0) {
          showNotification("No edges found in the PFN. Cannot proceed with MEGENA analysis.", type = "error")
          cat("Error: igraph object has zero edges.\n")
          return(NULL)
        }
        
        cat("Number of edges in igraph object:", ecount(pfn_igraph), "\n")
        
        incProgress(9/10, detail = "Performing MEGENA Analysis")
        cat("Step 9/10: Performing MEGENA Analysis.\n")
        
        # Run MEGENA with the PFN igraph object
        MEGENA_output <- tryCatch({
          cat("Running do.MEGENA...\n")
          do.MEGENA(
            g = pfn_igraph,
            mod.pval = input$MEGENApval,
            hub.pval = input$MEGENApval,
            min.size = input$MEGENAmin_size,
            doPar = FALSE,
            num.cores = 1,
            n.perm = 100  # Adjust as needed
          )
        }, error = function(e) {
          showNotification(paste("MEGENA analysis failed:", e$message), type = "error")
          cat("Error: MEGENA analysis failed -", e$message, "\n")
          return(NULL)
        })
        
        if (is.null(MEGENA_output)) {
          cat("Error: MEGENA_output is NULL.\n")
          return(NULL)
        }
        
        incProgress(10/10, detail = "Completed")
        cat("Step 10/10: MEGENA Analysis Completed.\n")
        megena_result(list(MEGENA_output = MEGENA_output, pfn_igraph = pfn_igraph))
        
        showNotification("MEGENA analysis completed!", type = "message")
        cat("Notification: MEGENA analysis completed successfully.\n")
        
        # Update module choices for visualization and enrichment
        modules <- MEGENA_output$module.output$modules
        updateSelectInput(session, "module_to_plot", choices = names(modules))
        updateSelectInput(session, "selected_module", choices = names(modules))
      })
    })
    
    # Render the network plot
    output$network_plot <- renderPlotly({
      req(megena_result())
      cat("Rendering network plot.\n")
      
      # Extract the PFN igraph object
      graph <- megena_result()$pfn_igraph
      
      # Ensure the graph object is valid
      if (is.null(graph) || !is_igraph(graph)) {
        showNotification("Graph object is empty or invalid.", type = "error")
        cat("Error: Invalid graph object for plotting.\n")
        return(NULL)
      }
      
      # Generate layout
      cat("Generating layout for the graph.\n")
      layout_fr <- tryCatch({
        layout_with_fr(graph)
      }, error = function(e) {
        showNotification(paste("Error in layout generation:", e$message), type = "error")
        cat("Error: Layout generation failed -", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(layout_fr)) {
        cat("Error: Layout generation returned NULL.\n")
        return(NULL)
      }
      
      # Convert layout to data frame
      layout_df <- as.data.frame(layout_fr)
      colnames(layout_df) <- c("x", "y")
      layout_df$gene <- V(graph)$name
      
      cat("Number of nodes:", nrow(layout_df), "\n")
      
      # Create edge list with coordinates
      edges <- as.data.frame(as_edgelist(graph, names = TRUE))
      colnames(edges) <- c("from", "to")
      
      # If edge weights exist, add them
      if ("weight" %in% edge_attr_names(graph)) {
        edges$weight <- E(graph)$weight
      } else {
        edges$weight <- 1  # Default weight
      }
      
      cat("Number of edges:", nrow(edges), "\n")
      
      # Add coordinates for edges
      edges$x <- layout_df$x[match(edges$from, layout_df$gene)]
      edges$y <- layout_df$y[match(edges$from, layout_df$gene)]
      edges$xend <- layout_df$x[match(edges$to, layout_df$gene)]
      edges$yend <- layout_df$y[match(edges$to, layout_df$gene)]
      
      # Plot using plotly
      cat("Plotting network using Plotly.\n")
      p <- plot_ly() %>%
        add_segments(
          data = edges,
          x = ~x,
          y = ~y,
          xend = ~xend,
          yend = ~yend,
          line = list(color = 'grey', width = 1),
          showlegend = FALSE,
          hoverinfo = "none"
        ) %>%
        add_markers(
          data = layout_df,
          x = ~x,
          y = ~y,
          type = 'scatter',
          mode = 'markers',
          text = ~gene,
          hoverinfo = 'text',
          marker = list(size = 5, color = 'blue')  # Adjust size for clarity
        ) %>%
        layout(title = "Gene Co-Expression Network", showlegend = FALSE)
      
      cat("Network plot rendered successfully.\n")
      p  # Return the plotly object
    })
    
    # Module Summary Table
    output$module_summary_table <- DT::renderDataTable({
      req(megena_result())
      MEGENA_output <- megena_result()$MEGENA_output
      
      # Generate module summary
      module_summary <- MEGENA.ModuleSummary(
        MEGENA_output,
        mod.pvalue = input$MEGENApval,
        hub.pvalue = input$MEGENApval,
        min.size = input$MEGENAmin_size,
        max.size = 10000,
        output.sig = TRUE
      )
      
      # Extract module table
      module_table <- module_summary$module.table
      
      # Display as a DataTable
      DT::datatable(module_table, options = list(pageLength = 10))
    })
    
    # Hub Genes Information
    output$hub_genes_table <- DT::renderDataTable({
      req(megena_result())
      MEGENA_output <- megena_result()$MEGENA_output
      
      # Extract hub genes
      hub_genes <- get.hub.genes(MEGENA_output)
      
      cat("Number of hub genes extracted:", nrow(hub_genes), "\n")
      
      if (nrow(hub_genes) == 0) {
        showNotification("No hub genes found.", type = "warning")
      }
      
      # Display as a DataTable
      DT::datatable(hub_genes, options = list(pageLength = 10))
    })
    
    # Helper function to extract hub genes
    get.hub.genes <- function(MEGENA_output) {
      if (is.null(MEGENA_output$hub.genes)) {
        return(data.frame(Module = character(), Gene = character(), stringsAsFactors = FALSE))
      }
      
      hub_list <- MEGENA_output$hub.genes
      
      if (length(hub_list) == 0) {
        return(data.frame(Module = character(), Gene = character(), stringsAsFactors = FALSE))
      }
      
      # Convert list to data frame
      hub_df <- do.call(rbind, lapply(names(hub_list), function(module_id) {
        data.frame(
          Module = module_id,
          Gene = hub_list[[module_id]],
          stringsAsFactors = FALSE
        )
      }))
      return(hub_df)
    }
    
    # Gene Module Membership
    output$gene_module_table <- DT::renderDataTable({
      req(megena_result())
      MEGENA_output <- megena_result()$MEGENA_output
      
      # Extract modules
      modules <- MEGENA_output$module.output$modules
      
      # Create a data frame of gene-module assignments
      gene_module_df <- data.frame(
        Gene = unlist(modules),
        Module = rep(names(modules), sapply(modules, length)),
        stringsAsFactors = FALSE
      )
      
      # Display as a DataTable
      DT::datatable(gene_module_df, options = list(pageLength = 10))
    })
    
    # Module Visualization
    observeEvent(megena_result(), {
      modules <- megena_result()$MEGENA_output$module.output$modules
      updateSelectInput(session, "module_to_plot", choices = names(modules))
    })
    
    output$module_plot <- renderPlotly({
      req(input$module_to_plot)
      module_genes <- megena_result()$MEGENA_output$module.output$modules[[input$module_to_plot]]
      
      # Subgraph for the selected module
      subgraph <- induced_subgraph(megena_result()$pfn_igraph, vids = module_genes)
      
      # Generate layout
      layout_sub <- layout_with_fr(subgraph)
      layout_df <- as.data.frame(layout_sub)
      colnames(layout_df) <- c("x", "y")
      layout_df$gene <- V(subgraph)$name
      
      # Create edge list with coordinates
      edges <- as.data.frame(as_edgelist(subgraph, names = TRUE))
      colnames(edges) <- c("from", "to")
      
      # Add coordinates for edges
      edges$x <- layout_df$x[match(edges$from, layout_df$gene)]
      edges$y <- layout_df$y[match(edges$from, layout_df$gene)]
      edges$xend <- layout_df$x[match(edges$to, layout_df$gene)]
      edges$yend <- layout_df$y[match(edges$to, layout_df$gene)]
      
      # Plot using plotly
      plot_ly() %>%
        add_segments(
          data = edges,
          x = ~x,
          y = ~y,
          xend = ~xend,
          yend = ~yend,
          line = list(color = 'grey', width = 1),
          showlegend = FALSE,
          hoverinfo = "none"
        ) %>%
        add_markers(
          data = layout_df,
          x = ~x,
          y = ~y,
          type = 'scatter',
          mode = 'markers',
          text = ~gene,
          hoverinfo = 'text',
          marker = list(size = 10, color = 'red')  # Highlight module genes
        ) %>%
        layout(title = paste("Module", input$module_to_plot), showlegend = FALSE)
    })
    
    # Enrichment Analysis
    observeEvent(megena_result(), {
      modules <- megena_result()$MEGENA_output$module.output$modules
      updateSelectInput(session, "selected_module", choices = names(modules))
    })
    
    output$enrichment_table <- DT::renderDataTable({
      req(input$selected_module)
      modules <- megena_result()$MEGENA_output$module.output$modules
      genes <- modules[[input$selected_module]]
      
      # Convert gene symbols to Entrez IDs
      entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)$ENTREZID
      
      # Perform GO enrichment analysis
      ego <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      # Check if enrichment results are available
      if (is.null(ego) || nrow(ego) == 0) {
        showNotification("No enrichment results found for this module.", type = "warning")
        return(NULL)
      }
      
      # Display results
      DT::datatable(as.data.frame(ego), options = list(pageLength = 10))
    })
    
    # Bar Plot of Top Enriched GO Terms
    output$enrichment_barplot <- renderPlotly({
      req(input$selected_module)
      modules <- megena_result()$MEGENA_output$module.output$modules
      genes <- modules[[input$selected_module]]
      
      # Convert gene symbols to Entrez IDs
      entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)$ENTREZID
      
      # Perform GO enrichment analysis
      ego <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      # Check if enrichment results are available
      if (is.null(ego) || nrow(ego) == 0) {
        showNotification("No enrichment results found for this module.", type = "warning")
        return(NULL)
      }
      
      # Select top 10 enriched GO terms
      top_ego <- ego@result %>% 
        arrange(p.adjust) %>% 
        head(10)
      
      # Create bar plot using ggplot2
      p <- ggplot(top_ego, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        xlab("GO Term") +
        ylab("-log10 Adjusted P-value") +
        ggtitle(paste("Top 10 Enriched GO Terms for Module", input$selected_module)) +
        theme_minimal()
      
      # Convert to plotly
      ggplotly(p)
    })
    
    # Dot Plot of Enriched GO Terms
    output$enrichment_dotplot <- renderPlotly({
      req(input$selected_module)
      modules <- megena_result()$MEGENA_output$module.output$modules
      genes <- modules[[input$selected_module]]
      
      # Convert gene symbols to Entrez IDs
      entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)$ENTREZID
      
      # Perform GO enrichment analysis
      ego <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      # Check if enrichment results are available
      if (is.null(ego) || nrow(ego) == 0) {
        showNotification("No enrichment results found for this module.", type = "warning")
        return(NULL)
      }
      
      # Select top 10 enriched GO terms
      top_ego <- ego@result %>% 
        arrange(p.adjust) %>% 
        head(10)
      
      # Create dot plot using ggplot2
      p <- ggplot(top_ego, aes(x = Count, y = reorder(Description, -p.adjust), size = Count, color = p.adjust)) +
        geom_point() +
        scale_color_continuous(low = "red", high = "blue") +
        xlab("Gene Count") +
        ylab("GO Term") +
        ggtitle(paste("Enriched GO Terms for Module", input$selected_module)) +
        theme_minimal()
      
      # Convert to plotly
      ggplotly(p)
    })
    
    # Cnet Plot (Enrichment Map)
    output$enrichment_cnetplot <- renderPlotly({
      req(input$selected_module)
      modules <- megena_result()$MEGENA_output$module.output$modules
      genes <- modules[[input$selected_module]]
      
      # Convert gene symbols to Entrez IDs
      entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)$ENTREZID
      
      # Perform GO enrichment analysis
      ego <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      # Check if enrichment results are available
      if (is.null(ego) || nrow(ego) == 0) {
        showNotification("No enrichment results found for this module.", type = "warning")
        return(NULL)
      }
      
      # Create Cnet plot using enrichplot
      p <- cnetplot(ego, showCategory = 5, foldChange = NULL, circular = FALSE, colorEdge = TRUE)
      
      # Convert to plotly
      ggplotly(p)
    })
    
    # Downloadable Results
    output$download_module_summary <- downloadHandler(
      filename = function() {
        paste0("module_summary_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(megena_result())
        MEGENA_output <- megena_result()$MEGENA_output
        
        # Generate module summary
        module_summary <- MEGENA.ModuleSummary(
          MEGENA_output,
          mod.pvalue = input$MEGENApval,
          hub.pvalue = input$MEGENApval,
          min.size = input$MEGENAmin_size,
          max.size = 10000,
          output.sig = TRUE
        )
        
        # Extract module table
        module_table <- module_summary$module.table
        
        write.csv(module_table, file, row.names = FALSE)
      }
    )
    
    output$download_hub_genes <- downloadHandler(
      filename = function() {
        paste0("hub_genes_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(megena_result())
        hub_genes <- get.hub.genes(megena_result()$MEGENA_output)
        
        if (nrow(hub_genes) == 0) {
          write.csv(data.frame(Message = "No hub genes found."), file, row.names = FALSE)
        } else {
          write.csv(hub_genes, file, row.names = FALSE)
        }
      }
    )
    
    output$download_gene_module <- downloadHandler(
      filename = function() {
        paste0("gene_module_membership_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(megena_result())
        modules <- megena_result()$MEGENA_output$module.output$modules
        gene_module_df <- data.frame(
          Gene = unlist(modules),
          Module = rep(names(modules), sapply(modules, length)),
          stringsAsFactors = FALSE
        )
        write.csv(gene_module_df, file, row.names = FALSE)
      }
    )
    
    # Helper function to extract hub genes
    get.hub.genes <- function(MEGENA_output) {
      if (is.null(MEGENA_output$hub.genes)) {
        return(data.frame(Module = character(), Gene = character(), stringsAsFactors = FALSE))
      }
      
      hub_list <- MEGENA_output$hub.genes
      
      if (length(hub_list) == 0) {
        return(data.frame(Module = character(), Gene = character(), stringsAsFactors = FALSE))
      }
      
      # Convert list to data frame
      hub_df <- do.call(rbind, lapply(names(hub_list), function(module_id) {
        data.frame(
          Module = module_id,
          Gene = hub_list[[module_id]],
          stringsAsFactors = FALSE
        )
      }))
      return(hub_df)
    }
  })
}
