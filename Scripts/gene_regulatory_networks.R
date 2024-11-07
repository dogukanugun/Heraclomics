# Load required libraries within the module file GNR2
library(shiny)
library(shinydashboard)
library(shiny)
library(DT)
library(GENIE3)
library(RTN)
library(AUCell)
library(Seurat)
library(igraph)
library(visNetwork)
library(plotly)
library(pheatmap)
library(shinyalert)
library(uwot)
library(shinyjs)
library(tidyverse)
library(tidyr)
library(msigdbr)

# UI Function for Gene Regulatory Networks Module
geneRegulatoryNetworkUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyalert(),   # For alert messages
    h2("Gene Regulatory Networks (GRN) Analysis"),
    fluidRow(
      box(
        title = "GRN Analysis Controls",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        actionButton(
          ns("start_grn"),
          "Start GRN Analysis",
          icon = icon("play"),
          class = "btn-success"
        ),
        br(),
        br(),
        verbatimTextOutput(ns("grn_log"))
      )
    ),
    # Outputs for plots and tables
    fluidRow(
      box(
        title = "Regulon Activity UMAP",
        width = 6,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotlyOutput(ns("regulon_umap"))
      ),
      box(
        title = "Regulon Activity Heatmap",
        width = 6,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotOutput(ns("regulon_heatmap"))
      )
    ),
    fluidRow(
      box(
        title = "Interactive Regulatory Network",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        visNetworkOutput(ns("interactive_network"))
      )
    ),
    fluidRow(
      box(
        title = "Regulon Table",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        DTOutput(ns("regulon_table"))
      )
    )
  )
}

# Server Function for Gene Regulatory Networks Module
geneRegulatoryNetworkServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize Shinyalert
    shinyalert::useShinyalert()
    
    # Reactive values to store GRN results
    grn_results <- reactiveValues(
      genie3_links = NULL,
      regulons = NULL,
      auc_scores = NULL,
      umap_coords = NULL
    )
    
    # Pre-defined TF list loaded manually or from alternative sources
    predefined_tf_list <- c("TP53", "NFYA", "SP1", "E2F1", "MYC", "STAT3", "AP1", "CREB1", "SMAD4", "REL")  # Example list
    
    # Observe the Start GRN Analysis button
    observeEvent(input$start_grn, {
      print("Starting GRN analysis...")
      print(is.null(rv$corrected_data))  # Check if rv$corrected_data is available
      print(str(rv$corrected_data))  # Inspect its structure
      
      req(rv$corrected_data)  # Ensure corrected_data is available
      
      # Show a modal dialog to get GRN settings
      showModal(modalDialog(
        title = "GRN Analysis Settings",
        radioButtons(
          ns("grn_tf_selection"),
          "Select Transcription Factor (TF) Selection Method:",
          choices = c("Use Pre-defined TF List from MSigDB", "Specify TFs Manually"),
          selected = "Use Pre-defined TF List from MSigDB"
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Specify TFs Manually'", ns("grn_tf_selection")),
          textAreaInput(
            ns("grn_tfs"),
            "Enter Transcription Factors (one per line):",
            placeholder = "e.g.,\nTF1\nTF2\nTF3"
          )
        ),
        actionButton(
          ns("confirm_grn"),
          "Run GRN Analysis",
          class = "btn-primary"
        ),
        footer = modalButton("Cancel"),
        easyClose = TRUE
      ))
    })
    
    # Confirm GRN Analysis settings and start analysis
    observeEvent(input$confirm_grn, {
      print("Confirm GRN event triggered")  # Debugging
      print(input$grn_tf_selection)
      
      # Initialize tf_input
      tf_input <- NULL
      
      # Validate TF selection
      if (input$grn_tf_selection == "Specify TFs Manually") {
        tf_input <- strsplit(input$grn_tfs, "\n")[[1]]
        tf_input <- trimws(tf_input)
        tf_input <- tf_input[tf_input != ""]
        if (length(tf_input) == 0) {
          shinyalert::shinyalert(
            title = "Input Error",
            text = "Please enter at least one transcription factor.",
            type = "error"
          )
          return(NULL)
        }
      } else if (input$grn_tf_selection == "Use Pre-defined TF List from MSigDB") {
        # Use the manually loaded predefined TF list
        tf_input <- predefined_tf_list
        
        if (length(tf_input) == 0) {
          shinyalert::shinyalert(
            title = "TF List Error",
            text = "No transcription factors available in the predefined TF list.",
            type = "error"
          )
          return(NULL)
        }
      }
      
      # Remove the modal dialog
      removeModal()
      
      # Proceed with GRN analysis
      perform_grn_analysis(tf_input)
    })
    
    # Function to perform GRN analysis
    perform_grn_analysis <- function(tf_input) {
      showModal(modalDialog(
        title = "GRN Analysis in Progress",
        "Please wait while GRN analysis is being performed...",
        footer = NULL,
        easyClose = FALSE
      ))
      
      tryCatch({
        withProgress(message = 'Running GRN Analysis...', value = 0, {
          # Step 1: Prepare expression matrix
          incProgress(0.05, detail = "Preparing expression data...")
          
          # Filter to highly variable genes
          rv$corrected_data <- FindVariableFeatures(rv$corrected_data, selection.method = "vst", nfeatures = 2000)
          variable_genes <- VariableFeatures(rv$corrected_data)
          
          # Extract the expression matrix for variable genes
          exprMat <- as.matrix(GetAssayData(rv$corrected_data, slot = "data")[variable_genes, ])
          
          # Diagnostic Check: Expression Matrix Dimensions
          print("Expression Matrix Dimensions:")
          print(dim(exprMat))
          
          if (is.null(exprMat) || length(dim(exprMat)) < 2 || dim(exprMat)[1] < 2 || dim(exprMat)[2] < 2) {
            stop("Expression matrix must have at least two genes and two cells.")
          }
          
          # Optional: Subset to a reasonable number of cells to save memory
          # Uncomment and adjust as needed
          # set.seed(123)
          # sampled_cells <- sample(colnames(exprMat), size = min(2000, ncol(exprMat)))
          # exprMat <- exprMat[, sampled_cells]
          
          # Diagnostic Check: After optional subsetting
          # print("Expression Matrix Dimensions after subsetting:")
          # print(dim(exprMat))
          
          # Convert to sparse matrix to save memory (if compatible with downstream functions)
          # exprMat <- Matrix(exprMat, sparse = TRUE)
          
          # Step 2: Determine TFs based on user selection
          incProgress(0.1, detail = "Selecting transcription factors...")
          if (!is.null(tf_input)) {
            tfs <- intersect(tf_input, rownames(exprMat))
            if (length(tfs) == 0) {
              stop("None of the specified TFs are present in the dataset.")
            }
            print(paste("Number of TFs after intersection:", length(tfs)))
            print("List of Selected TFs:")
            print(tfs)
            
            # Optional: Select top TFs based on average expression to reduce memory usage
            tf_expression <- rowMeans(exprMat[tfs, ])
            top_n <- min(100, length(tfs))  # Adjust '100' based on your requirements
            top_tfs <- names(sort(tf_expression, decreasing = TRUE))[1:top_n]
            tfs <- top_tfs
            print(paste("Number of TFs after selecting top", top_n, ":", length(tfs)))
            print("Top TFs:")
            print(tfs)
          } else {
            stop("No transcription factors provided.")
          }
          
          # Step 3: Run GENIE3
          incProgress(0.2, detail = "Running GENIE3 to infer regulatory links...")
          
          # Adjust GENIE3 parameters to optimize memory usage
          weightMatrix <- GENIE3(exprMat, regulators = tfs, nTrees = 50, verbose = TRUE, nCores = max(1, parallel::detectCores() - 1))
          
          # Diagnostic Check: Weight Matrix Dimensions
          print("Weight Matrix Dimensions:")
          print(dim(weightMatrix))
          
          if (is.null(weightMatrix) || length(dim(weightMatrix)) < 2) {
            stop("GENIE3 weight matrix is not two-dimensional.")
          }
          
          # Extract regulatory links with a higher threshold to reduce the number of links
          linkList <- getLinkList(weightMatrix, threshold = 0.05)  # Example threshold; adjust as needed
          
          # Diagnostic Check: Regulatory Links
          print("Regulatory Links Sample:")
          print(head(linkList))
          
          if (nrow(linkList) == 0) {
            stop("GENIE3 did not identify any regulatory links. Consider adjusting the threshold or verifying the TF list.")
          }
          
          grn_results$genie3_links <- linkList
          
          # Free up memory
          rm(weightMatrix)
          gc()
          
          # Step 4: Run RTN to construct regulons
          incProgress(0.3, detail = "Running RTN to construct regulons...")
          
          tni <- tni.constructor(expData = exprMat, regulatoryElements = tfs)
          tni <- tni.permutation(tni, steps = 100)    # Increase steps if needed
          tni <- tni.bootstrap(tni, steps = 100)      # Increase steps if needed
          tni <- tni.dpi.filter(tni, dpi_threshold = 0.8)  # Adjust threshold as needed
          regulons <- tni.get(tni, what = "regulons")
          
          # Diagnostic Check: Regulons
          print(paste("Number of Regulons Identified:", length(regulons)))
          if (length(regulons) == 0) {
            stop("RTN did not identify any regulons. Check the input data and TF list.")
          }
          print("First Regulon:")
          print(regulons[[1]])
          
          grn_results$regulons <- regulons
          
          # Free up memory
          rm(tni)
          gc()
          
          # Step 5: Calculate Regulon Activity with AUCell
          incProgress(0.4, detail = "Calculating regulon activity with AUCell...")
          
          geneSets <- regulons
          
          # Create gene expression rankings for AUCell
          cells_rankings <- AUCell_buildRankings(exprMat, nCores = max(1, parallel::detectCores() - 1), plotStats = FALSE)
          
          # Calculate AUCell scores
          aucScores <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(exprMat)))
          
          # Diagnostic Check: AUCell Scores
          print("AUCell Scores Dimensions:")
          print(dim(aucScores))
          
          if (is.null(aucScores) || length(dim(aucScores)) < 2) {
            stop("AUCell failed to calculate regulon activity scores.")
          }
          
          grn_results$auc_scores <- aucScores
          
          # Free up memory
          rm(cells_rankings)
          gc()
          
          # Step 6: Perform UMAP on AUCell AUC scores
          incProgress(0.5, detail = "Performing dimensionality reduction (UMAP)...")
          
          aucMatrix <- getAUC(aucScores)
          
          # Diagnostic Check: AUCell Matrix Dimensions
          print("AUCell Matrix Dimensions:")
          print(dim(aucMatrix))
          
          if (is.null(aucMatrix) || length(dim(aucMatrix)) < 2 || dim(aucMatrix)[1] < 2 || dim(aucMatrix)[2] < 2) {
            stop("AUCell matrix must have at least two regulons and two cells for UMAP.")
          }
          
          set.seed(123)  # For reproducibility
          umap_result <- uwot::umap(t(aucMatrix), n_neighbors = 30, min_dist = 0.3, metric = "cosine")
          rownames(umap_result) <- colnames(aucMatrix)
          grn_results$umap_coords <- umap_result
          
          # Diagnostic Check: UMAP Result
          print("UMAP Result Dimensions:")
          print(dim(umap_result))
          
          # Free up memory
          rm(aucMatrix)
          gc()
          
          incProgress(0.6, detail = "Finalizing GRN Analysis...")
          
          # Update reactive values
          rv$genie3_links <- linkList
          rv$regulons <- regulons
          rv$auc_scores <- aucScores
          
          incProgress(0.7, detail = "Rendering Plots...")
          
          # Additional steps like plotting can be done here
          
          # Finalize progress
          incProgress(1, detail = "GRN Analysis Completed!")
          removeModal()
          shinyalert::shinyalert(
            title = "GRN Analysis Completed",
            text = "Gene Regulatory Network analysis has been completed successfully!",
            type = "success"
          )
        })
      }, error = function(e) {
        removeModal()
        # Capture the stack trace
        error_message <- paste("An error occurred during GRN analysis:", e$message)
        stack_trace <- paste(as.character(sys.calls()), collapse = "\n")
        shinyalert::shinyalert(
          title = "GRN Analysis Error",
          text = paste(error_message, "\n\nStack Trace:\n", stack_trace),
          type = "error"
        )
      })
    }
    
    # Output: GRN Analysis Log
    output$grn_log <- renderPrint({
      if (is.null(grn_results$genie3_links)) {
        cat("GRN analysis not yet performed.")
      } else {
        cat("GRN analysis completed.\n")
        cat("Number of regulatory links:", nrow(grn_results$genie3_links), "\n")
        cat("Number of regulons identified:", length(grn_results$regulons), "\n")
      }
    })
    
    # Output: Regulon Activity UMAP Plot
    output$regulon_umap <- renderPlotly({
      req(grn_results$umap_coords)
      plot_data <- data.frame(
        UMAP1 = grn_results$umap_coords[, 1],
        UMAP2 = grn_results$umap_coords[, 2],
        Cell = rownames(grn_results$umap_coords)
      )
      
      # Add cluster information from Seurat object, if available
      if ("seurat_clusters" %in% colnames(rv$corrected_data@meta.data)) {
        plot_data$Cluster <- rv$corrected_data@meta.data$seurat_clusters[plot_data$Cell]
      } else {
        plot_data$Cluster <- "NA"
      }
      
      # Plot using plotly for interactivity
      plot_ly(
        data = plot_data,
        x = ~UMAP1,
        y = ~UMAP2,
        color = ~Cluster,
        text = ~Cell,
        type = 'scatter',
        mode = 'markers',
        marker = list(size = 5)
      ) %>%
        layout(title = "UMAP of Regulon Activity",
               xaxis = list(title = "UMAP 1"),
               yaxis = list(title = "UMAP 2"))
    })
    
    # Output: Regulon Activity Heatmap
    output$regulon_heatmap <- renderPlot({
      req(grn_results$auc_scores)
      aucMatrix <- getAUC(grn_results$auc_scores)
      
      # Select top regulons based on average activity
      top_regulons <- head(order(rowMeans(aucMatrix), decreasing = TRUE), 50)
      selected_auc <- aucMatrix[top_regulons, ]
      
      # Diagnostic Check: Heatmap Input Dimensions
      print("Heatmap Input Dimensions:")
      print(dim(selected_auc))
      
      # Plot heatmap
      pheatmap::pheatmap(
        selected_auc,
        show_rownames = TRUE,
        show_colnames = FALSE,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = "row",
        main = "Regulon Activity Heatmap (Top 50)"
      )
    })
    
    # Output: Interactive Regulatory Network
    output$interactive_network <- renderVisNetwork({
      req(grn_results$genie3_links)
      nodes <- data.frame(
        id = unique(c(grn_results$genie3_links$regulatoryGene, grn_results$genie3_links$targetGene)),
        label = unique(c(grn_results$genie3_links$regulatoryGene, grn_results$genie3_links$targetGene)),
        stringsAsFactors = FALSE
      )
      
      edges <- grn_results$genie3_links[, c("regulatoryGene", "targetGene", "weight")]
      colnames(edges) <- c("from", "to", "weight")
      
      # Apply a weight threshold to filter edges
      weight_threshold <- quantile(edges$weight, 0.99)
      filtered_edges <- edges[edges$weight >= weight_threshold, ]
      
      if (nrow(filtered_edges) == 0) {
        shinyalert::shinyalert(
          title = "No Significant Links",
          text = "No regulatory links met the threshold criteria.",
          type = "info"
        )
        return(NULL)
      }
      
      filtered_edges$color <- "gray"
      filtered_edges$arrows <- "to"
      
      nodes$group <- ifelse(nodes$id %in% grn_results$genie3_links$regulatoryGene, "TF", "Target")
      
      visNetwork::visNetwork(nodes, filtered_edges) %>%
        visNetwork::visNodes(color = list(border = "black", background = "lightgray")) %>%
        visNetwork::visGroups(groupname = "TF", color = "orange") %>%
        visNetwork::visGroups(groupname = "Target", color = "lightblue") %>%
        visNetwork::visEdges(arrows = "to") %>%
        visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
    })
    
    # Output: Regulon Table
    output$regulon_table <- renderDT({
      req(grn_results$regulons)
      regulons_df <- data.frame(
        Regulon = names(grn_results$regulons),
        Genes = sapply(grn_results$regulons, paste, collapse = ", "),
        stringsAsFactors = FALSE
      )
      
      DT::datatable(regulons_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })
  })
}

