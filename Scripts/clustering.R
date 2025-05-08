# clustering.R

# Load necessary libraries
library(shiny)
library(Seurat)
library(ggplot2)
library(plotly)
library(dplyr)
library(htmlwidgets)
library(RColorBrewer)


# Clustering Module UI Function
clusteringUI <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Clustering"),
    fluidRow(
      box(
        title = "Select Parameters",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        tabsetPanel(
          tabPanel("Nearest-neighbour Graph",
                   numericInput(ns("nn_dims"), "Number of Dimensions:", value = 15, min = 1),
                   numericInput(ns("k_param"), "k.param:", value = 20, min = 1),
                   numericInput(ns("n_trees"), "n.trees:", value = 50, min = 1)
          ),
          tabPanel("Clustering Parameters",
                   numericInput(ns("resolution"), "Resolution:", value = 0.3, min = 0.1, step = 0.1),
                   selectInput(ns("clustering_algorithm"), "Clustering Algorithm:",
                               choices = c("Louvain", "SLM", "Leiden"), selected = "Louvain")
          ),
          tabPanel("UMAP Parameters",
                   numericInput(ns("umap_dims"), "Number of Dimensions:", value = 15, min = 1),
                   numericInput(ns("umap_n_neighbors"), "k-nearest-neighbours:", value = 20, min = 1),
                   numericInput(ns("umap_min_dist"), "min.dist:", value = 0.3, min = 0.0, max = 1.0, step = 0.1)
          ),
          tabPanel("t-SNE Parameters",
                   numericInput(ns("tsne_dims"), "Number of Dimensions:", value = 15, min = 1),
                   numericInput(ns("tsne_perplexity"), "Perplexity:", value = 30, min = 1),
                   numericInput(ns("tsne_max_iter"), "Max iterations:", value = 1000, min = 250)
          )
        ),
        br(),
        actionButton(ns("run_clustering"), "Apply Clustering", icon = icon("play"), class = "btn-success")
      )
    ),
    fluidRow(
      box(
        title = "UMAP/t-SNE Plot",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        selectInput(ns("plot_type"), "Select Plot Type:",
                    choices = c("UMAP", "3D UMAP", "t-SNE", "3D t-SNE"), selected = "UMAP"),
        # Removed "sample" from choices
        selectInput(ns("label_by"), "Label Cells By:",
                    choices = c("seurat_clusters"), selected = "seurat_clusters"),
        plotlyOutput(ns("cluster_plot")),
        downloadButton(ns("download_cluster_plot"), "Download Plot")
      )
    ),
    fluidRow(
      box(
        title = "Cell Composition Summary",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotOutput(ns("cell_composition_plot")),
        downloadButton(ns("download_composition_plot"), "Download Composition Plot")
      )
    ),
    fluidRow(
      box(
        width = 12,
        actionButton(ns("next_step"), "Continue to Next Step", icon = icon("arrow-right"), class = "btn-primary")
      )
    )
  )
}

# Clustering Module Server Function
clusteringServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Observe when Apply Clustering button is clicked
    observeEvent(input$run_clustering, {
      req(rv$seurat_integrated)
      
      # Show a modal dialog to confirm clustering parameters
      showModal(
        modalDialog(
          title = "Confirm Clustering Parameters",
          "Are you sure you want to proceed with clustering using the selected parameters?",
          easyClose = FALSE,
          footer = tagList(
            modalButton("Cancel"),
            actionButton(ns("confirm_clustering"), "Yes, Proceed", class = "btn-primary")
          )
        )
      )
    })
    
    # Observe when the user confirms clustering
    observeEvent(input$confirm_clustering, {
      removeModal()  # Remove the confirmation modal
      
      # Show a modal dialog indicating clustering is running
      showModal(
        modalDialog(
          title = "Clustering in Progress",
          "Please wait while clustering is being performed...",
          footer = NULL
        )
      )
      
      output$clustering_status <- renderText({
        "Clustering started..."
      })
      
      tryCatch({
        ### Nearest-neighbour Graph Construction ###
        dims_used <- 1:input$nn_dims
        k_param <- input$k_param
        n_trees <- input$n_trees
        
        rv$seurat_integrated <- FindNeighbors(
          rv$seurat_integrated,
          dims = dims_used,
          k.param = k_param,
          nn.method = "annoy",
          n.trees = n_trees
        )
        
        ### Clustering ###
        clustering_algorithm <- switch(
          input$clustering_algorithm,
          "Louvain" = 1,
          "SLM" = 3,
          "Leiden" = 4
        )
        
        rv$seurat_integrated <- FindClusters(
          rv$seurat_integrated,
          resolution = input$resolution,
          algorithm = clustering_algorithm
        )
        
        ### UMAP ###
        rv$seurat_integrated <- RunUMAP(
          rv$seurat_integrated,
          dims = 1:input$umap_dims,
          n.neighbors = input$umap_n_neighbors,
          min.dist = input$umap_min_dist,
          n.components = 3
        )
        
        ### t-SNE ###
        rv$seurat_integrated <- RunTSNE(
          rv$seurat_integrated,
          dims = 1:input$tsne_dims,
          perplexity = input$tsne_perplexity,
          max_iter = input$tsne_max_iter,
          dim.embed = 3
        )
        
        output$clustering_status <- renderText({
          "Clustering completed."
        })
        
        # Remove the modal dialog when clustering is complete
        removeModal()
        
        ### Generate UMAP/t-SNE Plot ###
        output$cluster_plot <- renderPlotly({
          req(input$plot_type)
          label_by <- input$label_by
          
          if (input$plot_type == "UMAP") {
            p <- DimPlot(
              rv$seurat_integrated,
              reduction = "umap",
              group.by = label_by,
              label = TRUE
            ) +
              ggtitle("UMAP Plot")
            ggplotly(p)
          } else if (input$plot_type == "3D UMAP") {
            umap_coords <- Embeddings(rv$seurat_integrated, "umap")
            plot_data <- data.frame(
              x = umap_coords[, 1],
              y = umap_coords[, 2],
              z = umap_coords[, 3],
              group = rv$seurat_integrated@meta.data[[label_by]]
            )
            plot_ly(
              data = plot_data,
              x = ~x, y = ~y, z = ~z,
              color = ~group,
              colors = RColorBrewer::brewer.pal(n = 11, name = "Set1"),
              type = "scatter3d",
              mode = "markers",
              marker = list(size = 2)
            ) %>%
              layout(title = "3D UMAP Plot")
          } else if (input$plot_type == "t-SNE") {
            p <- DimPlot(
              rv$seurat_integrated,
              reduction = "tsne",
              group.by = label_by,
              label = TRUE
            ) +
              ggtitle("t-SNE Plot")
            ggplotly(p)
          } else if (input$plot_type == "3D t-SNE") {
            tsne_coords <- Embeddings(rv$seurat_integrated, "tsne")
            plot_data <- data.frame(
              x = tsne_coords[, 1],
              y = tsne_coords[, 2],
              z = tsne_coords[, 3],
              group = rv$seurat_integrated@meta.data[[label_by]]
            )
            plot_ly(
              data = plot_data,
              x = ~x, y = ~y, z = ~z,
              color = ~group,
              colors = RColorBrewer::brewer.pal(n = 11, name = "Set1"),
              type = "scatter3d",
              mode = "markers",
              marker = list(size = 2)
            ) %>%
              layout(title = "3D t-SNE Plot")
          }
        })
        
        ### Cell Composition Summary ###
        output$cell_composition_plot <- renderPlot({
          label_by <- input$label_by
          data <- rv$seurat_integrated@meta.data %>%
            group_by(!!sym(label_by)) %>%
            summarise(Count = n()) %>%
            mutate(Percentage = Count / sum(Count) * 100)
          ggplot(data, aes_string(x = label_by, y = "Percentage", fill = label_by)) +
            geom_bar(stat = "identity") +
            theme_minimal() +
            ylab("Percentage (%)") +
            xlab("") +
            ggtitle("Cell Composition Summary")
        })
        
      }, error = function(e) {
        output$clustering_status <- renderText({
          paste("Error:", e$message)
        })
        removeModal()  # Remove the progress modal in case of error
        # Show an error modal dialog
        showModal(
          modalDialog(
            title = "Error",
            paste("An error occurred during clustering:", e$message),
            easyClose = TRUE,
            footer = modalButton("Dismiss")
          )
        )
      })
    })
    
    ### Download Handlers ###
    # Download UMAP/t-SNE Plot
    output$download_cluster_plot <- downloadHandler(
      filename = function() {
        paste0("cluster_plot_", Sys.Date(), ".html")
      },
      content = function(file) {
        label_by <- input$label_by
        
        if (input$plot_type == "UMAP") {
          p <- DimPlot(
            rv$seurat_integrated,
            reduction = "umap",
            group.by = label_by,
            label = TRUE
          ) +
            ggtitle("UMAP Plot")
          p <- ggplotly(p)
          htmlwidgets::saveWidget(p, file)
        } else if (input$plot_type == "3D UMAP") {
          umap_coords <- Embeddings(rv$seurat_integrated, "umap")
          plot_data <- data.frame(
            x = umap_coords[, 1],
            y = umap_coords[, 2],
            z = umap_coords[, 3],
            group = rv$seurat_integrated@meta.data[[label_by]]
          )
          p <- plot_ly(
            data = plot_data,
            x = ~x, y = ~y, z = ~z,
            color = ~group,
            colors = RColorBrewer::brewer.pal(n = 11, name = "Set1"),
            type = "scatter3d",
            mode = "markers",
            marker = list(size = 2)
          ) %>%
            layout(title = "3D UMAP Plot")
          htmlwidgets::saveWidget(p, file)
        } else if (input$plot_type == "t-SNE") {
          p <- DimPlot(
            rv$seurat_integrated,
            reduction = "tsne",
            group.by = label_by,
            label = TRUE
          ) +
            ggtitle("t-SNE Plot")
          p <- ggplotly(p)
          htmlwidgets::saveWidget(p, file)
        } else if (input$plot_type == "3D t-SNE") {
          tsne_coords <- Embeddings(rv$seurat_integrated, "tsne")
          plot_data <- data.frame(
            x = tsne_coords[, 1],
            y = tsne_coords[, 2],
            z = tsne_coords[, 3],
            group = rv$seurat_integrated@meta.data[[label_by]]
          )
          p <- plot_ly(
            data = plot_data,
            x = ~x, y = ~y, z = ~z,
            color = ~group,
            colors = RColorBrewer::brewer.pal(n = 11, name = "Set1"),
            type = "scatter3d",
            mode = "markers",
            marker = list(size = 2)
          ) %>%
            layout(title = "3D t-SNE Plot")
          htmlwidgets::saveWidget(p, file)
        }
      }
    )
    
    # Download Cell Composition Plot
    output$download_composition_plot <- downloadHandler(
      filename = function() {
        paste0("cell_composition_plot_", Sys.Date(), ".png")
      },
      content = function(file) {
        label_by <- input$label_by
        data <- rv$seurat_integrated@meta.data %>%
          group_by(!!sym(label_by)) %>%
          summarise(Count = n()) %>%
          mutate(Percentage = Count / sum(Count) * 100)
        p <- ggplot(data, aes_string(x = label_by, y = "Percentage", fill = label_by)) +
          geom_bar(stat = "identity") +
          theme_minimal() +
          ylab("Percentage (%)") +
          xlab("") +
          ggtitle("Cell Composition Summary")
        ggsave(file, plot = p, width = 8, height = 6)
      }
    )
    
    ### Next Step Button ###
    observeEvent(input$next_step, {
      rv$proceed_to_next_step <- TRUE  # Set this reactive value
    })
    
    # Return the reactive value to signal the next step
    return(list(proceed_to_next_step = reactive(rv$proceed_to_next_step)))
  })
    
  
}
