# Load required libraries
library(shiny)
library(shinydashboard)
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
library(msigdbr)
library(Matrix)
library(future)
library(logger)
library(SCENIC)      # For SCENIC analysis
library(RcisTarget)  # For motif enrichment
library(doParallel)  # For parallel processing in SCENIC

# UI Function for Gene Regulatory Networks Module
geneRegulatoryNetworkUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyalert(),
    useShinyjs(),
    
    h2("Gene Regulatory Networks (GRN) Analysis"),
    
    fluidRow(
      box(
        title = "GRN Analysis Controls",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(6,
                 radioButtons(
                   ns("analysis_method"),
                   "Select Analysis Method:",
                   choices = c("GENIE3 + RTN" = "genie3",
                               "SCENIC" = "scenic"),
                   selected = "genie3"
                 ),
                 actionButton(
                   ns("start_grn"),
                   "Start GRN Analysis",
                   icon = icon("play"),
                   class = "btn-success"
                 ),
                 actionButton(
                   ns("save_results"),
                   "Save Results",
                   icon = icon("save"),
                   class = "btn-info"
                 )
          ),
          column(6,
                 conditionalPanel(
                   condition = sprintf("input['%s'] == 'genie3'", ns("analysis_method")),
                   numericInput(
                     ns("n_trees"),
                     "Number of Trees (GENIE3):",
                     value = 200,
                     min = 50,
                     max = 1000
                   )
                 ),
                 numericInput(
                   ns("auc_max_rank"),
                   "AUCell Max Rank (% of genes):",
                   value = 1,
                   min = 0.1,
                   max = 10,
                   step = 0.1
                 )
          )
        ),
        br(),
        verbatimTextOutput(ns("grn_log"))
      )
    ),
    
    # Analysis Settings Box
    fluidRow(
      box(
        title = "Network Settings",
        width = 12,
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        conditionalPanel(
          condition = sprintf("input['%s'] == 'genie3'", ns("analysis_method")),
          sliderInput(
            ns("edge_weight_threshold"),
            "Edge Weight Threshold Percentile:",
            min = 90,
            max = 99.9,
            value = 95,
            step = 0.1
          ),
          numericInput(
            ns("n_permutations"),
            "Number of Permutations (RTN):",
            value = 100,
            min = 50,
            max = 1000
          ),
          numericInput(
            ns("n_bootstraps"),
            "Number of Bootstraps (RTN):",
            value = 100,
            min = 50,
            max = 1000
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'scenic'", ns("analysis_method")),
          numericInput(
            ns("scenic_n_cores"),
            "Number of cores for SCENIC:",
            value = 4,
            min = 1,
            max = parallel::detectCores()
          ),
          selectInput(
            ns("scenic_species"),
            "Species:",
            choices = c("Human" = "hgnc",
                        "Mouse" = "mgi"),
            selected = "hgnc"
          ),
          numericInput(
            ns("scenic_motif_threshold"),
            "Motif similarity threshold:",
            value = 0.75,
            min = 0,
            max = 1,
            step = 0.05
          )
        ),
        numericInput(
          ns("max_regulons"),
          "Maximum Regulons to Display:",
          value = 50,
          min = 10,
          max = 200
        ),
        selectInput(
          ns("umap_metric"),
          "UMAP Distance Metric:",
          choices = c("cosine", "euclidean", "manhattan", "pearson"),
          selected = "cosine"
        ),
        checkboxInput(
          ns("parallel_processing"),
          "Enable Parallel Processing",
          value = TRUE
        )
      )
    ),
    
    # Outputs for plots and tables
    tabsetPanel(
      tabPanel("Visualizations",
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
               )
      ),
      tabPanel("Results",
               fluidRow(
                 box(
                   title = "Regulon Table",
                   width = 12,
                   status = "primary",
                   solidHeader = TRUE,
                   collapsible = TRUE,
                   DTOutput(ns("regulon_table"))
                 )
               ),
               conditionalPanel(
                 condition = sprintf("input['%s'] == 'scenic'", ns("analysis_method")),
                 fluidRow(
                   box(
                     title = "SCENIC Regulon Specificity Scores",
                     width = 12,
                     status = "primary",
                     solidHeader = TRUE,
                     collapsible = TRUE,
                     plotOutput(ns("regulon_specificity"))
                   )
                 )
               ),
               fluidRow(
                 box(
                   title = "Top Regulatory Links",
                   width = 12,
                   status = "primary",
                   solidHeader = TRUE,
                   collapsible = TRUE,
                   DTOutput(ns("links_table"))
                 )
               )
      ),
      tabPanel("Summary",
               fluidRow(
                 box(
                   title = "GRN Summary Statistics",
                   width = 12,
                   status = "primary",
                   solidHeader = TRUE,
                   collapsible = TRUE,
                   verbatimTextOutput(ns("summary_stats"))
                 )
               )
      )
    )
  )
}

# Server Function for Gene Regulatory Networks Module
geneRegulatoryNetworkServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize logger
    log_threshold(INFO)
    log_layout(layout_glue_colors)
    
    # Reactive values to store GRN results
    grn_results <- reactiveValues(
      genie3_links = NULL,
      regulons = NULL,
      auc_scores = NULL,
      umap_coords = NULL,
      tf_list = NULL,
      analysis_complete = FALSE,
      scenic_results = NULL
    )
    
    # Get transcription factors
    get_tf_list <- function() {
      log_info("Fetching transcription factors from MSigDB")
      tfs <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
        filter(gs_subcat %in% c("TFT:GTRD", "TFT:JASPAR")) %>% 
        distinct(gene_symbol) %>% 
        pull(gene_symbol)
      
      log_info("Found {length(tfs)} transcription factors in MSigDB")
      return(tfs)
    }
    
    # Run GENIE3 + RTN analysis
    run_genie3_analysis <- function(exprMat, tfs) {
      withProgress(message = 'Running GENIE3 + RTN Analysis...', value = 0, {
        # Step 1: Run GENIE3 (requires a dense numeric matrix)
        incProgress(0.2, detail = "Running GENIE3 to infer regulatory links...")
        weightMatrix <- GENIE3(
          exprMat,
          regulators = tfs,
          nTrees = input$n_trees,
          verbose = TRUE,
          nCores = ifelse(input$parallel_processing, availableCores() - 1, 1)
        )
        
        linkList <- getLinkList(weightMatrix, threshold = 0.01)
        grn_results$genie3_links <- linkList
        rm(weightMatrix)
        gc()
        
        # Step 2: Run RTN to construct regulons
        incProgress(0.5, detail = "Running RTN to construct regulons...")
        tni <- tni.constructor(
          expData = exprMat,
          regulatoryElements = tfs
        )
        
        tni <- tni.permutation(
          tni,
          nPermutations = input$n_permutations
        )
        
        tni <- tni.bootstrap(
          tni,
          nBootstraps = input$n_bootstraps
        )
        
        tni <- tni.dpi.filter(tni)
        regulons <- tni.get(tni, what = "regulons")
        grn_results$regulons <- regulons
        rm(tni)
        gc()
        
        return(regulons)
      })
    }
    
    # Run SCENIC analysis
    run_scenic_analysis <- function(exprMat) {
      withProgress(message = 'Running SCENIC Analysis...', value = 0, {
        # Set up SCENIC
        incProgress(0.1, detail = "Initializing SCENIC...")
        
        # Prepare expression matrix
        exprMat <- as.matrix(exprMat)
        
        # Initialize SCENIC settings with correct database paths
        org <- switch(input$scenic_species,
                      "hgnc" = "hgnc",
                      "mgi" = "mgi")
        
        # Point to your downloaded databases
        dbDir <- "C:" # CHANGE THIS TO YOUR ACTUAL PATH
        
        # Select appropriate database based on species
        if(org == "hgnc") {
          dbs <- c("10kb" = "Users/USER78/Documents/Heraclomics/Scripts/cisTarget_databases/hgnc/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                   "500bp" = "Users/USER78/Documents/Heraclomics/Scripts/cisTarget_databases/hgnc/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
        } else { # mgi
          dbs <- c("10kb" = "Users/USER78/Documents/Heraclomics/Scripts/cisTarget_databases/mgi/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
                   "500bp" = "Users/USER78/Documents/Heraclomics/Scripts/cisTarget_databases/mgi/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
        }
        
        scenicOptions <- initializeScenic(
          org = org,
          dbDir = dbDir,
          dbs = dbs,
          datasetTitle = "SCENIC Analysis",
          nCores = input$scenic_n_cores
        )
        
        # Co-expression network
        incProgress(0.2, detail = "Inferring co-expression network...")
        genesKept <- geneFiltering(exprMat, scenicOptions)
        exprMat_filtered <- exprMat[genesKept, ]
        runCorrelation(exprMat_filtered, scenicOptions)
        
        # GENIE3 (run within SCENIC)
        incProgress(0.3, detail = "Running GENIE3...")
        runGenie3(exprMat_filtered, scenicOptions)
        
        # Build and score GRN
        incProgress(0.4, detail = "Building and scoring GRN...")
        runSCENIC_1_coexNetwork2modules(scenicOptions)
        runSCENIC_2_createRegulons(scenicOptions, 
                                   minMotifScore = input$scenic_motif_threshold)
        
        # Calculate AUC
        incProgress(0.6, detail = "Calculating regulon activity...")
        aucellApp <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
        
        # Save SCENIC results
        grn_results$scenic_results <- list(
          regulons = loadInt(scenicOptions, "regulons"),
          aucell = aucellApp,
          scenicOptions = scenicOptions
        )
        
        # Format regulons for consistency with GENIE3 output
        scenic_regulons <- lapply(grn_results$scenic_results$regulons, function(x) x$targets)
        names(scenic_regulons) <- names(grn_results$scenic_results$regulons)
        grn_results$regulons <- scenic_regulons
        
        # Clean up temporary SCENIC directory
        unlink(scenic_dir, recursive = TRUE)
        
        return(scenic_regulons)
      })
    }
    
    # Calculate AUC scores (common for both methods)
    calculate_auc_scores <- function(exprMat, regulons) {
      withProgress(message = 'Calculating AUC scores...', value = 0.7, {
        cells_rankings <- AUCell_buildRankings(
          exprMat,
          nCores = ifelse(input$parallel_processing, availableCores() - 1, 1),
          plotStats = FALSE
        )
        
        aucScores <- AUCell_calcAUC(
          regulons,
          cells_rankings,
          aucMaxRank = ceiling(input$auc_max_rank/100 * nrow(exprMat))
        )
        
        grn_results$auc_scores <- aucScores
        rm(cells_rankings)
        gc()
        
        return(aucScores)
      })
    }
    
    # Perform UMAP (common for both methods)
    perform_umap <- function(aucScores) {
      withProgress(message = 'Performing UMAP...', value = 0.9, {
        aucMatrix <- getAUC(aucScores)
        set.seed(123)
        
        umap_result <- uwot::umap(
          t(aucMatrix),
          n_neighbors = 30,
          min_dist = 0.3,
          metric = input$umap_metric
        )
        
        rownames(umap_result) <- colnames(aucMatrix)
        grn_results$umap_coords <- umap_result
        
        return(umap_result)
      })
    }
    
    # Observe the Start GRN Analysis button
    observeEvent(input$start_grn, {
      req(rv$corrected_data)
      
      log_info("Starting GRN analysis")
      
      # Show settings modal
      showModal(
        modalDialog(
          title = "GRN Analysis Settings",
          size = "l",
          fluidRow(
            column(6,
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
                       placeholder = "e.g.,\nTP53\nMYC\nSTAT3",
                       rows = 10
                     )
                   )
            ),
            column(6,
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'genie3'", ns("analysis_method")),
                     helpText("GENIE3 + RTN settings shown")
                   ),
                   conditionalPanel(
                     condition = sprintf("input['%s'] == 'scenic'", ns("analysis_method")),
                     helpText("SCENIC settings shown")
                   )
            )
          ),
          footer = tagList(
            actionButton(ns("confirm_grn"), "Run Analysis", class = "btn-primary"),
            modalButton("Cancel")
          )
        )
      )
    })
    
    # Confirm GRN Analysis settings and start analysis
    observeEvent(input$confirm_grn, {
      removeModal()
      
      # Initialize parallel processing if enabled
      if (input$parallel_processing) {
        plan(multisession)
        log_info("Parallel processing enabled")
      } else {
        plan(sequential)
      }
      
      tryCatch({
        withProgress(message = 'Running GRN Analysis...', value = 0, {
          # Prepare expression matrix
          incProgress(0.05, detail = "Preparing expression data...")
          
          # Filter to highly variable genes
          rv$corrected_data <- FindVariableFeatures(
            rv$corrected_data, 
            selection.method = "vst", 
            nfeatures = 2000
          )
          
          variable_genes <- VariableFeatures(rv$corrected_data)
          exprMat <- as.matrix(GetAssayData(rv$corrected_data, slot = "data")[variable_genes, ])
          
          # For GENIE3, do NOT convert to a sparse matrix because GENIE3 requires a dense matrix.
          if (input$analysis_method == "scenic") {
            exprMat <- Matrix(exprMat, sparse = TRUE)
          }
          
          if (input$analysis_method == "genie3") {
            # Step 1: Prepare TF list
            incProgress(0.1, detail = "Preparing transcription factors...")
            
            if (input$grn_tf_selection == "Specify TFs Manually") {
              tf_input <- strsplit(input$grn_tfs, "\n")[[1]] %>% 
                trimws() %>% 
                .[. != ""]
              
              if (length(tf_input) == 0) {
                stop("No valid transcription factors entered")
              }
            } else {
              tf_input <- get_tf_list()
            }
            
            # Filter TFs to those present in the data
            tfs <- intersect(tf_input, rownames(exprMat))
            if (length(tfs) == 0) {
              stop("None of the specified TFs are present in the dataset")
            }
            
            # Limit to top TFs by expression to manage computation
            tf_expr <- rowMeans(exprMat[tfs, ])
            top_n <- min(200, length(tfs))
            tfs <- names(sort(tf_expr, decreasing = TRUE))[1:top_n]
            grn_results$tf_list <- tfs
            
            # Run GENIE3 + RTN
            regulons <- run_genie3_analysis(exprMat, tfs)
            
          } else if (input$analysis_method == "scenic") {
            # Run SCENIC
            regulons <- run_scenic_analysis(exprMat)
          }
          
          # Common steps for both methods
          aucScores <- calculate_auc_scores(exprMat, regulons)
          umap_result <- perform_umap(aucScores)
          
          # Mark analysis as complete
          grn_results$analysis_complete <- TRUE
          
          incProgress(1, detail = "Analysis complete!")
          shinyalert(
            title = "Success",
            text = "GRN analysis completed successfully!",
            type = "success"
          )
        })
      }, error = function(e) {
        shinyalert(
          title = "Error",
          text = paste("GRN analysis failed:", e$message),
          type = "error"
        )
        log_error("GRN analysis failed: {e$message}")
      }, finally = {
        plan(sequential) # Reset to sequential processing
      })
    })
    
    # SCENIC regulon specificity plot
    output$regulon_specificity <- renderPlot({
      req(grn_results$scenic_results, rv$corrected_data)
      
      # Get cell annotations
      cellInfo <- data.frame(
        Cluster = Idents(rv$corrected_data),
        row.names = colnames(rv$corrected_data)
      )
      
      # Get regulon AUC matrix
      aucellResults <- grn_results$scenic_results$aucell
      regulonAUC <- aucellResults@assays@data@listData$AUC
      
      # Calculate specificity scores
      regulonActivity_byCluster <- sapply(split(rownames(cellInfo), cellInfo$Cluster),
                                          function(cells) rowMeans(regulonAUC[, cells]))
      
      # Normalize by row (regulon)
      regulonActivity_byCluster <- t(scale(t(regulonActivity_byCluster),
                                           center = TRUE,
                                           scale = TRUE))
      
      # Convert to long format for plotting
      plotData <- regulonActivity_byCluster %>%
        as.data.frame() %>%
        rownames_to_column("Regulon") %>%
        pivot_longer(cols = -Regulon, names_to = "Cluster", values_to = "Score")
      
      # Plot
      ggplot(plotData, aes(x = Cluster, y = Regulon, fill = Score)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Regulon Specificity by Cluster",
             x = "Cluster",
             y = "Regulon",
             fill = "Z-score")
    })
    
    # Output: GRN Analysis Log
    output$grn_log <- renderPrint({
      if (!grn_results$analysis_complete) {
        cat("GRN analysis not yet performed.\n")
      } else {
        cat("=== GRN Analysis Summary ===\n")
        cat("Method:", ifelse(input$analysis_method == "genie3", "GENIE3 + RTN", "SCENIC"), "\n")
        cat("Number of TFs used:", length(grn_results$tf_list), "\n")
        cat("Number of regulatory links:", 
            ifelse(!is.null(grn_results$genie3_links), nrow(grn_results$genie3_links), "N/A (SCENIC)"), "\n")
        cat("Number of regulons identified:", length(grn_results$regulons), "\n")
        cat("Analysis completed successfully!\n")
      }
    })
    
    # Output: Summary Statistics
    output$summary_stats <- renderPrint({
      req(grn_results$analysis_complete)
      
      cat("=== GRN Network Statistics ===\n\n")
      
      if (input$analysis_method == "genie3") {
        # Regulatory links summary
        links <- grn_results$genie3_links
        cat("Regulatory Links:\n")
        cat("  Total links:", nrow(links), "\n")
        cat("  Mean weight:", mean(links$weight), "\n")
        cat("  Median weight:", median(links$weight), "\n")
        cat("  Max weight:", max(links$weight), "\n\n")
      }
      
      # Regulon summary
      regulons <- grn_results$regulons
      regulon_sizes <- lengths(regulons)
      cat("Regulons:\n")
      cat("  Total regulons:", length(regulons), "\n")
      cat("  Mean targets per regulon:", mean(regulon_sizes), "\n")
      cat("  Median targets per regulon:", median(regulon_sizes), "\n")
      cat("  Largest regulon:", max(regulon_sizes), "targets\n")
      cat("  Smallest regulon:", min(regulon_sizes), "targets\n\n")
      
      if (input$analysis_method == "genie3") {
        # Top TFs by number of targets
        tf_degrees <- grn_results$genie3_links %>% 
          count(regulatoryGene, name = "degree") %>% 
          arrange(desc(degree))
        
        cat("Top 10 TFs by number of targets:\n")
        print(head(tf_degrees, 10))
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
      
      # Add metadata if available
      if (!is.null(rv$corrected_data)) {
        meta_data <- rv$corrected_data@meta.data
        common_cells <- intersect(plot_data$Cell, rownames(meta_data))
        
        if (length(common_cells) > 0) {
          plot_data <- plot_data %>% 
            left_join(
              meta_data %>% 
                rownames_to_column("Cell") %>% 
                select(Cell, any_of(c("seurat_clusters", "cell_type"))),
              by = "Cell"
            )
        }
      }
      
      # Create plot
      p <- plot_ly(
        data = plot_data,
        x = ~UMAP1,
        y = ~UMAP2,
        
        text = ~paste("Cell:", Cell),
        hoverinfo = "text",
        type = 'scatter',
        mode = 'markers',
        marker = list(size = 5, opacity = 0.7)
      ) %>% 
        layout(
          title = "UMAP of Regulon Activity",
          xaxis = list(title = "UMAP 1"),
          yaxis = list(title = "UMAP 2"),
          showlegend = !is.null(plot_data$seurat_clusters)
        )
      
      return(p)
    })
    
    # Output: Regulon Activity Heatmap
    output$regulon_heatmap <- renderPlot({
      req(grn_results$auc_scores)
      
      aucMatrix <- getAUC(grn_results$auc_scores)
      n_regulons <- min(input$max_regulons, nrow(aucMatrix))
      
      # Select top regulons by variance
      regulon_vars <- apply(aucMatrix, 1, var)
      top_regulons <- names(sort(regulon_vars, decreasing = TRUE))[1:n_regulons]
      selected_auc <- aucMatrix[top_regulons, ]
      
      # Cluster cells if reasonable number
      cluster_cols <- ncol(selected_auc) <= 1000
      
      pheatmap(
        selected_auc,
        show_rownames = TRUE,
        show_colnames = FALSE,
        cluster_rows = TRUE,
        cluster_cols = cluster_cols,
        scale = "row",
        color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
        main = paste("Top", n_regulons, "Variable Regulons")
      )
    })
    
    # Output: Interactive Regulatory Network
    output$interactive_network <- renderVisNetwork({
      if (input$analysis_method == "genie3") {
        req(grn_results$genie3_links)
        
        links <- grn_results$genie3_links
        threshold <- quantile(links$weight, probs = input$edge_weight_threshold/100)
        filtered_links <- links[links$weight >= threshold, ]
        
        if (nrow(filtered_links) == 0) {
          showNotification("No edges meet the current threshold. Try lowering the threshold.", type = "warning")
          return(NULL)
        }
        
        # Prepare nodes for GENIE3 branch
        all_nodes <- unique(c(filtered_links$regulatoryGene, filtered_links$targetGene))
        nodes <- data.frame(
          id = all_nodes,
          label = all_nodes,
          group = ifelse(all_nodes %in% filtered_links$regulatoryGene, "TF", "Target"),
          stringsAsFactors = FALSE
        )
        
        # Prepare edges for GENIE3 branch
        edges <- filtered_links %>% 
          rename(from = regulatoryGene, to = targetGene, value = weight) %>% 
          mutate(
            color = "gray",
            arrows = "to",
            smooth = TRUE
          )
        
        # Create network for GENIE3 branch
        visNetwork(nodes, edges) %>%
          visNodes(
            shape = "dot",
            shadow = TRUE,
            size = 20
          ) %>%
          visGroups(
            groupname = "TF",
            color = list(
              background = "orange",
              border = "darkorange",
              highlight = "red"
            ),
            shape = "diamond"
          ) %>%
          visGroups(
            groupname = "Target",
            color = list(
              background = "lightblue",
              border = "darkblue",
              highlight = "blue"
            )
          ) %>%
          visLegend() %>%
          visOptions(
            highlightNearest = list(enabled = TRUE, degree = 1),
            nodesIdSelection = TRUE,
            selectedBy = "group"
          ) %>%
          visPhysics(
            solver = "forceAtlas2Based",
            forceAtlas2Based = list(gravitationalConstant = -50)
          ) %>%
          visInteraction(
            navigationButtons = TRUE,
            keyboard = TRUE
          )
      } else {
        # For SCENIC, show the top regulons and their targets
        req(grn_results$regulons)
        
        # Get top regulons by size
        top_regulons <- grn_results$regulons %>%
          enframe(name = "TF", value = "targets") %>%
          mutate(size = map_int(targets, length)) %>%
          arrange(desc(size)) %>%
          head(input$max_regulons)
        
        # Prepare nodes and edges
        nodes <- data.frame(
          id = unique(c(top_regulons$TF, unlist(top_regulons$targets))),
          label = unique(c(top_regulons$TF, unlist(top_regulons$targets))),
          group = ifelse(unique(c(top_regulons$TF, unlist(top_regulons$targets))) %in% top_regulons$TF,
                         "TF", "Target"),
          stringsAsFactors = FALSE
        )
        
        edges <- top_regulons %>%
          select(TF, targets) %>%
          unnest(targets) %>%
          rename(from = TF, to = targets) %>%
          mutate(
            color = "gray",
            arrows = "to",
            smooth = TRUE
          )
        
        # Create network for SCENIC branch
        visNetwork(nodes, edges) %>%
          visNodes(
            shape = "dot",
            shadow = TRUE,
            size = 20
          ) %>%
          visGroups(
            groupname = "TF",
            color = list(
              background = "orange",
              border = "darkorange",
              highlight = "red"
            ),
            shape = "diamond"
          ) %>%
          visGroups(
            groupname = "Target",
            color = list(
              background = "lightblue",
              border = "darkblue",
              highlight = "blue"
            )
          ) %>%
          visLegend() %>%
          visOptions(
            highlightNearest = list(enabled = TRUE, degree = 1),
            nodesIdSelection = TRUE,
            selectedBy = "group"
          ) %>%
          visPhysics(
            solver = "forceAtlas2Based",
            forceAtlas2Based = list(gravitationalConstant = -50)
          ) %>%
          visInteraction(
            navigationButtons = TRUE,
            keyboard = TRUE
          )
      }
    })
    
    # Output: Regulon Table
    output$regulon_table <- renderDT({
      req(grn_results$regulons)
      
      regulons_df <- enframe(grn_results$regulons, name = "Regulon", value = "Targets") %>% 
        mutate(
          Size = map_int(Targets, length),
          Targets = map_chr(Targets, ~paste(.x, collapse = ", "))
        ) %>% 
        arrange(desc(Size))
      
      datatable(
        regulons_df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        rownames = FALSE,
        extensions = 'Buttons',
        selection = 'single'
      )
    })
    
    # Output: Top Regulatory Links Table
    output$links_table <- renderDT({
      if (input$analysis_method == "genie3") {
        req(grn_results$genie3_links)
        
        links <- grn_results$genie3_links %>% 
          arrange(desc(weight)) %>% 
          rename(
            "TF" = regulatoryGene,
            "Target" = targetGene,
            "Weight" = weight
          )
        
        datatable(
          links,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel')
          ),
          rownames = FALSE,
          extensions = 'Buttons'
        )
      } else {
        # For SCENIC, show regulons and their sizes
        req(grn_results$regulons)
        
        regulons_df <- enframe(grn_results$regulons, name = "Regulon", value = "Targets") %>% 
          mutate(
            Size = map_int(Targets, length),
            Targets = map_chr(Targets, ~paste(.x, collapse = ", "))
          ) %>% 
          arrange(desc(Size))
        
        datatable(
          regulons_df,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel')
          ),
          rownames = FALSE,
          extensions = 'Buttons'
        )
      }
    })
    
    # Save results handler
    observeEvent(input$save_results, {
      req(grn_results$analysis_complete)
      
      showModal(
        modalDialog(
          title = "Save GRN Results",
          textInput(ns("save_name"), "Name for this analysis:"),
          selectInput(
            ns("save_format"),
            "Format:",
            choices = c("RDS", "CSV (selected components)")
          ),
          footer = tagList(
            modalButton("Cancel"),
            downloadButton(ns("download_results"), "Save")
          )
        )
      )
    })
    
    # Download handler for saving results
    output$download_results <- downloadHandler(
      filename = function() {
        paste0(input$save_name, "_GRN_results.", 
               ifelse(input$save_format == "RDS", "rds", "zip"))
      },
      content = function(file) {
        if (input$save_format == "RDS") {
          saveRDS(reactiveValuesToList(grn_results), file)
        } else {
          # Save selected components as CSV
          tmpdir <- tempdir()
          
          if (input$analysis_method == "genie3") {
            write_csv(grn_results$genie3_links, file.path(tmpdir, "regulatory_links.csv"))
          }
          
          write_csv(
            enframe(grn_results$regulons) %>% 
              mutate(targets = map_chr(value, ~paste(.x, collapse = ","))) %>% 
              select(-value),
            file.path(tmpdir, "regulons.csv")
          )
          
          write_csv(
            as.data.frame(t(getAUC(grn_results$auc_scores))),
            file.path(tmpdir, "auc_scores.csv")
          )
          
          if (input$analysis_method == "scenic") {
            # Save SCENIC-specific results
            scenic_regulon_scores <- grn_results$scenic_results$aucell@assays@data@listData$AUC
            write_csv(
              as.data.frame(scenic_regulon_scores) %>% rownames_to_column("Regulon"),
              file.path(tmpdir, "scenic_regulon_scores.csv")
            )
          }
          
          zip(file, files = list.files(tmpdir, pattern = "\\.csv$"), flags = "-j")
        }
      }
    )
  })
}
