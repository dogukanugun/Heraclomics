# Load Necessary Libraries
library(shiny)
library(CellChat)
library(Seurat)
library(shinyWidgets)
library(shinyalert)
library(ggplot2)
library(dplyr)
library(DT)
library(igraph)
library(NMF)
library(ggalluvial)
library(visNetwork)  # For interactive network plots

# UI Function for Cell-Cell Communication Module
cellCommunicationUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyalert(),
    
    h2("Cell-Cell Communication Analysis"),
    
    fluidRow(
      box(
        title = "Overview",
        width = 12,
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        p("This module employs the CellChat tool to quantitatively infer and analyze intercellular communication networks from single-cell RNA-seq data."),
        p("CellChat uses a comprehensive database of ligand-receptor pairs and applies social network analysis tools to identify significant cell-cell interactions."),
        p("Note: Cluster labels in plots may differ from original labels.")
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
          ns("cc_cell_types"),
          "Select Cell Types/Clusters:",
          choices = NULL,
          multiple = TRUE,
          options = list(actionsBox = TRUE, liveSearch = TRUE)
        ),
        
        actionButton(
          ns("run_cell_comm"),
          "Run Cell-Cell Communication Analysis",
          icon = icon("play"),
          class = "btn-success"
        )
      )
    ),
    
    # New Input: Select Signaling Pathways
    fluidRow(
      box(
        title = "Select Signaling Pathways",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        
        pickerInput(
          ns("select_pathways"),
          "Select Signaling Pathways to Visualize:",
          choices = NULL,
          multiple = TRUE,
          options = list(actionsBox = TRUE, liveSearch = TRUE)
        )
      )
    ),
    
    # Communication Network Visualization
    fluidRow(
      box(
        title = "Communication Network Visualization",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        # Use visNetworkOutput for interactive network plots
        visNetworkOutput(ns("cell_comm_network_plot")),
        downloadButton(ns("download_comm_network_plot"), "Download Plot")
      )
    ),
    
    # Visualization of Specific Signaling Pathway
    fluidRow(
      box(
        title = "Specific Signaling Pathway Visualization",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        # Static plot for pathway visualization
        plotOutput(ns("pathway_network_plot")),
        downloadButton(ns("download_pathway_network_plot"), "Download Plot")
      )
    ),
    
    # Communication Patterns Heatmap
    fluidRow(
      box(
        title = "Communication Patterns Heatmap",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        # Static plot for heatmap
        plotOutput(ns("communication_heatmap")),
        downloadButton(ns("download_communication_heatmap"), "Download Heatmap")
      )
    ),
    
    # Bubble Plot of Interactions
    fluidRow(
      box(
        title = "Bubble Plot of Interactions",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        # Static plot for bubble plot
        plotOutput(ns("bubble_plot")),
        downloadButton(ns("download_bubble_plot"), "Download Plot")
      )
    ),
    
    # Significant Ligand-Receptor Interactions Table
    fluidRow(
      box(
        title = "Significant Ligand-Receptor Interactions",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        DT::dataTableOutput(ns("signif_lig_recep_table")),
        downloadButton(ns("download_signif_lig_recep"), "Download Table")
      )
    )
  )
}

# Server Function for Cell-Cell Communication Module
cellCommunicationServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize Shinyalert
    shinyalert::useShinyalert()
    
    # Reactive Value to Store CellChat Object
    cellchat_reactive <- reactiveVal(NULL)
    
    # Update Cell Types/Clusters Choices After Data Correction Completes
    observeEvent(rv$corrected_data, {
      req(rv$corrected_data)
      local_corrected_data <- rv$corrected_data
      
      original_cluster_labels <- as.character(local_corrected_data$seurat_clusters)
      adjusted_cluster_labels <- as.numeric(as.character(local_corrected_data$seurat_clusters)) + 1
      local_corrected_data$adjusted_clusters <- as.character(adjusted_cluster_labels)
      Idents(local_corrected_data) <- local_corrected_data$adjusted_clusters
      
      cell_types <- unique(original_cluster_labels)
      cell_types <- sort(as.character(cell_types))
      updatePickerInput(session, "cc_cell_types", choices = cell_types, selected = cell_types)
    }, ignoreInit = TRUE)
    
    # Perform Cell-Cell Communication Analysis on Button Click
    observeEvent(input$run_cell_comm, {
      req(rv$corrected_data)
      req(input$cc_cell_types)
      
      if (length(input$cc_cell_types) < 2) {
        shinyalert("Input Error", "Please select at least two cell types/clusters for analysis.", type = "error")
        return(NULL)
      }
      
      tryCatch({
        showModal(modalDialog(
          title = "Cell-Cell Communication Analysis in Progress",
          "Please wait while the analysis is being performed...",
          footer = NULL,
          easyClose = FALSE
        ))
        
        withProgress(message = 'Performing Cell-Cell Communication Analysis...', value = 0, {
          local_corrected_data <- rv$corrected_data
          original_cluster_labels <- as.character(local_corrected_data$seurat_clusters)
          adjusted_cluster_labels <- as.numeric(as.character(local_corrected_data$seurat_clusters)) + 1
          local_corrected_data$adjusted_clusters <- as.character(adjusted_cluster_labels)
          Idents(local_corrected_data) <- local_corrected_data$adjusted_clusters
          
          cluster_mapping <- unique(data.frame(
            original = original_cluster_labels,
            adjusted = local_corrected_data$adjusted_clusters,
            stringsAsFactors = FALSE
          ))
          
          selected_clusters <- input$cc_cell_types
          selected_adjusted_clusters <- cluster_mapping$adjusted[cluster_mapping$original %in% selected_clusters]
          
          seurat_subset <- subset(local_corrected_data, idents = selected_adjusted_clusters)
          Idents(seurat_subset) <- seurat_subset$adjusted_clusters
          rownames(seurat_subset) <- toupper(rownames(seurat_subset))
          
          if (any(duplicated(rownames(seurat_subset)))) {
            removeModal()
            shinyalert("Data Error", "Duplicate gene names found after converting to uppercase. Please resolve duplicates.", type = "error")
            return(NULL)
          }
          
          cellchat <- createCellChat(object = seurat_subset, group.by = "adjusted_clusters")
          CellChatDB <- CellChatDB.human
          cellchat@DB <- CellChatDB
          
          signaling_genes <- unique(c(CellChatDB$interaction$ligand, CellChatDB$interaction$receptor))
          overlapping_genes <- intersect(rownames(cellchat@data), signaling_genes)
          
          if (length(overlapping_genes) == 0) {
            removeModal()
            shinyalert("Data Error", "No overlapping genes between your data and the CellChat database. Please check gene names.", type = "error")
            return(NULL)
          }
          
          cellchat <- subsetData(cellchat)
          
          if (dim(cellchat@data.signaling)[1] == 0) {
            removeModal()
            shinyalert("Data Error", "data.signaling is empty after subsetting. Please ensure that your data contains signaling genes.", type = "error")
            return(NULL)
          }
          
          incProgress(0.4, detail = "Identifying overexpressed ligands and receptors...")
          cellchat <- identifyOverExpressedGenes(cellchat)
          cellchat <- identifyOverExpressedInteractions(cellchat)
          
          incProgress(0.5, detail = "Computing communication probabilities...")
          cellchat <- computeCommunProb(cellchat)
          
          incProgress(0.6, detail = "Filtering communication probabilities...")
          cellchat <- filterCommunication(cellchat, min.cells = 10)
          
          incProgress(0.7, detail = "Computing communication probabilities at signaling pathway level...")
          cellchat <- computeCommunProbPathway(cellchat)
          cellchat <- aggregateNet(cellchat)
          
          incProgress(0.8, detail = "Computing network centrality scores...")
          cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
          
          cellchat_reactive(list(cellchat = cellchat, cluster_mapping = cluster_mapping))
          incProgress(1, detail = "Finalizing analysis...")
          removeModal()
          shinyalert("Success", "Cell-Cell Communication Analysis completed successfully!", type = "success")
        })
      }, error = function(e) {
        removeModal()
        shinyalert("Analysis Error", paste("An error occurred:", e$message), type = "error")
      })
    })
    
    # Provide Choices for Signaling Pathways After Analysis
    observeEvent(cellchat_reactive(), {
      req(cellchat_reactive())
      cellchat_obj <- cellchat_reactive()$cellchat
      pathways <- cellchat_obj@netP$pathways
      updatePickerInput(session, "select_pathways", choices = pathways, selected = pathways[1])
    })
    
    # Helper function to map adjusted cluster labels back to original
    map_clusters <- function(cluster_ids, cluster_mapping) {
      sapply(cluster_ids, function(x) {
        original <- cluster_mapping$original[cluster_mapping$adjusted == x]
        if (length(original) == 0) {
          return(NA)
        } else {
          return(original)
        }
      })
    }
    
    # Render Communication Network Plot (Interactive using visNetwork)
    output$cell_comm_network_plot <- renderVisNetwork({
      req(cellchat_reactive())
      cellchat_obj <- cellchat_reactive()$cellchat
      cluster_mapping <- cellchat_reactive()$cluster_mapping
      
      tryCatch({
        adjacency_matrix <- cellchat_obj@net$count
        
        # Convert the adjacency matrix to an edge list
        edge_list <- as.data.frame(as.table(adjacency_matrix))
        colnames(edge_list) <- c("from", "to", "weight")
        edge_list <- edge_list[edge_list$weight > 0, ]
        
        # Map clusters back to original labels
        edge_list$from <- map_clusters(edge_list$from, cluster_mapping)
        edge_list$to <- map_clusters(edge_list$to, cluster_mapping)
        
        # Create node data frame
        nodes <- data.frame(id = unique(c(edge_list$from, edge_list$to)))
        nodes$label <- nodes$id
        
        # Create visNetwork object
        visNetwork(nodes, edge_list) %>%
          visEdges(arrows = 'to') %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
        
      }, error = function(e) {
        shinyalert("Plot Error", paste("An error occurred while rendering the interactive plot:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # Download Communication Network Plot as PNG
    output$download_comm_network_plot <- downloadHandler(
      filename = function() {
        paste0("cell_comm_network_", Sys.Date(), ".html")
      },
      content = function(file) {
        req(cellchat_reactive())
        cellchat_obj <- cellchat_reactive()$cellchat
        cluster_mapping <- cellchat_reactive()$cluster_mapping
        tryCatch({
          adjacency_matrix <- cellchat_obj@net$count
          
          # Convert the adjacency matrix to an edge list
          edge_list <- as.data.frame(as.table(adjacency_matrix))
          colnames(edge_list) <- c("from", "to", "weight")
          edge_list <- edge_list[edge_list$weight > 0, ]
          
          # Map clusters back to original labels
          edge_list$from <- map_clusters(edge_list$from, cluster_mapping)
          edge_list$to <- map_clusters(edge_list$to, cluster_mapping)
          
          # Create node data frame
          nodes <- data.frame(id = unique(c(edge_list$from, edge_list$to)))
          nodes$label <- nodes$id
          
          # Create visNetwork object
          network_plot <- visNetwork(nodes, edge_list) %>%
            visEdges(arrows = 'to') %>%
            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
          
          # Save the interactive plot as an HTML file
          visSave(network_plot, file = file)
        }, error = function(e) {
          shinyalert("Download Error", paste("An error occurred while downloading the plot:", e$message), type = "error")
          return(NULL)
        })
      }
    )
    
    # Visualize Specific Signaling Pathways (Static)
    output$pathway_network_plot <- renderPlot({
      req(cellchat_reactive())
      req(input$select_pathways)
      cellchat_obj <- cellchat_reactive()$cellchat
      
      tryCatch({
        netVisual_aggregate(cellchat_obj, signaling = input$select_pathways, layout = "circle")
      }, error = function(e) {
        shinyalert("Plot Error", paste("An error occurred while rendering the pathway plot:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # Download Pathway Network Plot as PNG
    output$download_pathway_network_plot <- downloadHandler(
      filename = function() {
        paste0("pathway_network_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(cellchat_reactive())
        req(input$select_pathways)
        cellchat_obj <- cellchat_reactive()$cellchat
        tryCatch({
          png(file, width = 800, height = 800)
          netVisual_aggregate(cellchat_obj, signaling = input$select_pathways, layout = "circle")
          dev.off()
        }, error = function(e) {
          shinyalert("Download Error", paste("An error occurred while downloading the pathway plot:", e$message), type = "error")
          return(NULL)
        })
      }
    )
    
    # Visualize Communication Patterns Heatmap (Static)
    output$communication_heatmap <- renderPlot({
      req(cellchat_reactive())
      cellchat_obj <- cellchat_reactive()$cellchat
      
      tryCatch({
        netAnalysis_signalingRole_heatmap(cellchat_obj, pattern = "outgoing")
      }, error = function(e) {
        shinyalert("Plot Error", paste("An error occurred while rendering the heatmap:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # Download Communication Patterns Heatmap as PNG
    output$download_communication_heatmap <- downloadHandler(
      filename = function() {
        paste0("communication_heatmap_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(cellchat_reactive())
        cellchat_obj <- cellchat_reactive()$cellchat
        tryCatch({
          png(file, width = 800, height = 600)
          netAnalysis_signalingRole_heatmap(cellchat_obj, pattern = "outgoing")
          dev.off()
        }, error = function(e) {
          shinyalert("Download Error", paste("An error occurred while downloading the heatmap:", e$message), type = "error")
          return(NULL)
        })
      }
    )
    
    # Visualize Bubble Plot of Interactions (Static)
    output$bubble_plot <- renderPlot({
      req(cellchat_reactive())
      cellchat_obj <- cellchat_reactive()$cellchat
      
      tryCatch({
        netVisual_bubble(cellchat_obj, remove.isolate = TRUE)
      }, error = function(e) {
        shinyalert("Plot Error", paste("An error occurred while rendering the bubble plot:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # Download Bubble Plot as PNG
    output$download_bubble_plot <- downloadHandler(
      filename = function() {
        paste0("bubble_plot_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(cellchat_reactive())
        cellchat_obj <- cellchat_reactive()$cellchat
        tryCatch({
          png(file, width = 1000, height = 800)
          netVisual_bubble(cellchat_obj, remove.isolate = TRUE)
          dev.off()
        }, error = function(e) {
          shinyalert("Download Error", paste("An error occurred while downloading the bubble plot:", e$message), type = "error")
          return(NULL)
        })
      }
    )
    
    # Render Significant Ligand-Receptor Interactions Table
    output$signif_lig_recep_table <- DT::renderDataTable({
      req(cellchat_reactive())
      cellchat_obj <- cellchat_reactive()$cellchat
      cluster_mapping <- cellchat_reactive()$cluster_mapping
      
      tryCatch({
        significants <- subsetCommunication(cellchat_obj)
        if (nrow(significants) == 0) {
          shinyalert("Table Error", "No significant ligand-receptor interactions found.", type = "error")
          return(NULL)
        }
        df <- as.data.frame(significants)
        df$source <- map_clusters(df$source, cluster_mapping)
        df$target <- map_clusters(df$target, cluster_mapping)
        DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
      }, error = function(e) {
        shinyalert("Table Error", paste("An error occurred while rendering the table:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # Download Significant Ligand-Receptor Interactions Table as CSV
    output$download_signif_lig_recep <- downloadHandler(
      filename = function() {
        paste0("signif_lig_recep_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(cellchat_reactive())
        cellchat_obj <- cellchat_reactive()$cellchat
        cluster_mapping <- cellchat_reactive()$cluster_mapping
        tryCatch({
          significants <- subsetCommunication(cellchat_obj)
          if (nrow(significants) == 0) {
            shinyalert("Download Error", "No significant interactions to download.", type = "error")
            return(NULL)
          }
          df <- as.data.frame(significants)
          df$source <- map_clusters(df$source, cluster_mapping)
          df$target <- map_clusters(df$target, cluster_mapping)
          write.csv(df, file, row.names = FALSE)
        }, error = function(e) {
          shinyalert("Download Error", paste("An error occurred while downloading the table:", e$message), type = "error")
          return(NULL)
        })
      }
    )
  })
}
