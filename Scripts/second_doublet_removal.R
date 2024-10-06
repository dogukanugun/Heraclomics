# Scripts/second_doublet_removal.R

library(shiny)
library(shinyjs)
library(shinyalert)
library(Seurat)
library(DoubletFinder)
library(plotly)
library(parallel)
library(htmlwidgets)
library(future)
library(future.apply)


# Second Doublet Removal Module UI Function
secondDoubletRemovalUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),  # Enable shinyjs
    h2("Second Dataset - Doublet Removal"),
    fluidRow(
      box(
        title = "Doublet Removal Controls",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(
            width = 6,
            actionButton(ns("start_doublet_removal"), "Start Doublet Removal", icon = icon("play"), class = "btn-success"),
            actionButton(ns("skip_doublet_removal"), "SKIP Doublet Removal", icon = icon("forward"), class = "btn-warning")
          ),
          column(
            width = 6,
            shinyWidgets::radioGroupButtons(
              ns("plot_type"),
              "Select Plot Type:",
              choices = c("UMAP" = "UMAP", "3D UMAP" = "3D UMAP", "t-SNE" = "t-SNE", "3D t-SNE" = "3D t-SNE"),
              selected = "UMAP",
              justified = TRUE,
              checkIcon = list(
                yes = icon("ok", lib = "glyphicon"),
                no = icon("remove", lib = "glyphicon")
              )
            )
          )
        )
      )
    ),
    fluidRow(
      box(
        title = "Plot Adjustments",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(
            width = 4,
            sliderInput(
              ns("label_size"),
              "Label Size:",
              min = 5,
              max = 20,
              value = 10
            )
          ),
          column(
            width = 4,
            sliderInput(
              ns("point_size"),
              "Point Size:",
              min = 0.1,
              max = 5,
              value = 1,
              step = 0.1
            )
          ),
          column(
            width = 4,
            sliderInput(
              ns("plot_width"),
              "Plot Width (px):",
              min = 300,
              max = 1200,
              value = 800,
              step = 50
            )
          )
        )
      )
    ),
    fluidRow(
      box(
        title = "Doublet Removal Plots",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        tabsetPanel(
          tabPanel(
            "Plot",
            downloadButton(ns("download_plot"), "Download Plot"),
            plotlyOutput(ns("plot_output"), height = "auto")
          )
        )
      )
    ),
    fluidRow(
      box(
        title = "Doublet Operations",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(
            width = 6,
            actionButton(ns("simulate_doublets"), "Simulate Doublets", icon = icon("dice"), class = "btn-info")
          ),
          column(
            width = 6,
            actionButton(ns("remove_doublets"), "Remove Doublets", icon = icon("trash"), class = "btn-danger")
          )
        )
      )
    ),
    fluidRow(
      box(
        title = "Proceed to Next Step",
        width = 12,
        status = "success",
        solidHeader = TRUE,
        collapsible = TRUE,
        actionButton(ns("continue_to_load_integration"), "Continue to Next Step", icon = icon("arrow-right"))
      )
    )
  )
}

# Second Doublet Removal Module Server Function
secondDoubletRemovalServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values to store processed data
    processed_data <- reactiveVal(NULL)
    proceed_to_load_integration <- reactiveVal(FALSE)
    
    # Disable certain buttons initially
    observe({
      shinyjs::disable(ns("simulate_doublets"))
      shinyjs::disable(ns("remove_doublets"))
      shinyjs::disable(ns("continue_to_load_integration"))
    })
    
    # Start Doublet Removal
    observeEvent(input$start_doublet_removal, {
      req(rv$seurat_second)
      
      showModal(modalDialog("Doublet Removal in progress for the second dataset. Please wait...", footer = NULL))
      
      withProgress(message = 'Preprocessing Data and Detecting Doublets...', value = 0, {
        data <- rv$seurat_second
        
        # Preprocessing Steps
        incProgress(0.2, detail = "Normalizing Data")
        data <- NormalizeData(data)
        
        incProgress(0.2, detail = "Finding Variable Features")
        data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
        
        incProgress(0.1, detail = "Scaling Data")
        data <- ScaleData(data, features = VariableFeatures(object = data))
        
        incProgress(0.1, detail = "Running PCA")
        data <- RunPCA(data, features = VariableFeatures(object = data))
        
        # Verify PCA results
        if (!"pca" %in% names(data@reductions)) {
          removeModal()
          showNotification("PCA results not found in the Seurat object. Please check the PCA step.", type = "error")
          return(NULL)
        }
        
        incProgress(0.1, detail = "Running UMAP")
        data <- RunUMAP(data, dims = 1:10)
        
        incProgress(0.1, detail = "Running t-SNE")
        data <- RunTSNE(data, dims = 1:10)
        
        # Update rv with processed data
        rv$seurat_second <- data
        processed_data(data)
      })
      
      removeModal()
      shinyalert("Success", "Starting Doublet Removal for the second dataset completed successfully! Now you can Simulate Doublets.", type = "success")
      
      # Enable Simulate and Remove Doublets buttons
      shinyjs::enable(ns("simulate_doublets"))
      shinyjs::enable(ns("remove_doublets"))
    })
    
    # Skip Doublet Removal
    observeEvent(input$skip_doublet_removal, {
      req(rv$seurat_second)
      processed_data(rv$seurat_second)
      
      shinyalert("Skipped", "Doublet removal for the second dataset has been skipped.", type = "info")
      
      # Enable Next Step button
      shinyjs::enable(ns("continue_to_load_integration"))
    })
    
    # Simulate Doublets
    observeEvent(input$simulate_doublets, {
      req(processed_data())
      
      showModal(modalDialog("Simulating Doublets for the second dataset. Please wait...", footer = NULL))
      
      withProgress(message = 'Simulating and Detecting Doublets...', value = 0, {
        data <- processed_data()
        
        # DoubletFinder Steps
        incProgress(0.2, detail = "Doublet Detection with DoubletFinder")
        
        # Parameter Sweep
        sweep.res.list <- paramSweep(data, PCs = 1:10, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        
        # Choose optimal pK
        pK_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
        if (length(pK_opt) == 0 || is.na(pK_opt)) {
          pK_opt <- 0.1  # Fallback pK
          showNotification("Optimal pK not found. Using default pK = 0.1.", type = "warning")
        }
        
        incProgress(0.2, detail = "Estimating Doublet Proportion")
        
        # Estimate Homotypic Doublet Proportion
        data <- FindNeighbors(data, dims = 1:10)
        data <- FindClusters(data, resolution = 0.5)
        annotations <- data@meta.data$seurat_clusters
        homotypic.prop <- modelHomotypic(annotations)
        
        # Adjust expected doublets
        nExp_poi <- round(0.075 * ncol(data))  # 7.5% expected doublets
        nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
        
        incProgress(0.4, detail = "Running DoubletFinder")
        
        # Run DoubletFinder
        data <- doubletFinder(data, PCs = 1:10, pN = 0.25, pK = pK_opt,
                              nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
        
        # Extract Doublet Classification Column
        df_class_col <- grep("DF.classifications", colnames(data@meta.data), value = TRUE)
        if (length(df_class_col) == 0) {
          removeModal()
          showNotification("Doublet classification column not found.", type = "error")
          return(NULL)
        }
        
        data$Doublet <- data@meta.data[[df_class_col]]
        
        incProgress(0.2, detail = "Updating Data")
        
        # Update rv with doublet classifications
        rv$seurat_second <- data
        processed_data(data)
      })
      
      removeModal()
      shinyalert("Success", "Doublets have been simulated and detected successfully for the second dataset!", type = "success")
      
      # Enable "Remove Doublets" and "Continue to Next Step" buttons
      shinyjs::enable(ns("remove_doublets"))
      shinyjs::enable(ns("continue_to_load_integration"))
      
      # Generate Plots
      output$plot_output <- renderPlotly({
        req(processed_data())  # Ensure data is available
        
        data <- processed_data()
        plot_type <- input$plot_type  # Get selected plot type (UMAP, 3D UMAP, t-SNE, 3D t-SNE)
        
        # Determine reduction method and dimensionality
        if (plot_type %in% c("UMAP", "3D UMAP")) {
          reduction <- "umap"
        } else {
          reduction <- "tsne"
        }
        
        # Check if the selected reduction exists in the Seurat object
        if (!(reduction %in% names(data@reductions))) {
          showNotification(paste0(reduction, " not found in the Seurat object."), type = "error")
          return(NULL)
        }
        
        # Extract embeddings
        embeddings <- Embeddings(data, reduction = reduction)
        
        # Check dimensionality for 3D plots
        if (plot_type %in% c("3D UMAP", "3D t-SNE") && ncol(embeddings) < 3) {
          showNotification("3D embeddings not available. Please ensure that the reduction has at least 3 dimensions.", type = "error")
          return(NULL)
        }
        
        # Prepare plot data
        plot_data <- data.frame(
          x = embeddings[, 1],
          y = embeddings[, 2],
          group = data@meta.data$Doublet
        )
        
        # Include z-coordinate for 3D plots
        if (plot_type %in% c("3D UMAP", "3D t-SNE")) {
          plot_data$z <- embeddings[, 3]
        }
        
        # Dynamic color palette based on the number of groups
        n_groups <- length(unique(plot_data$group))
        colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_groups)
        
        # Create plots based on plot type
        if (plot_type %in% c("UMAP", "t-SNE")) {
          # 2D Plot
          p <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
            geom_point(size = input$point_size) +
            theme_minimal() +
            theme(text = element_text(size = input$label_size)) +
            labs(title = paste(plot_type, "Plot"), color = "Group")
          
          ggplotly(p)
        } else {
          # 3D Plot
          plot_ly(
            data = plot_data,
            x = ~x, y = ~y, z = ~z,
            color = ~group,
            colors = colors,
            type = "scatter3d",
            mode = "markers",
            marker = list(size = input$point_size)
          ) %>%
            layout(title = paste("3D", plot_type, "Plot"))
        }
      })
    })
    
    # Remove Doublets
    observeEvent(input$remove_doublets, {
      req(processed_data())
      data <- processed_data()
      
      # Filter out doublets
      singlet_data <- subset(data, subset = Doublet == "Singlet")
      
      # Update rv with singlet data
      rv$seurat_second <- singlet_data
      processed_data(singlet_data)
      
      shinyalert("Success", "Doublets have been removed successfully for the second dataset!", type = "success")
      
      # Update plots if needed
      # (Implement plot updates if necessary)
    })
    
    # Continue to Next Step
    observeEvent(input$continue_to_load_integration, {
      shinyalert(
        title = "Proceed to Integration?",
        text = "Are you sure you want to proceed to the integration step?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, proceed",
        cancelButtonText = "No, stay here",
        callbackR = function(value) {
          if (isTRUE(value)) {
            # Debugging print
            print("Confirmed: Proceeding to Integration.")
            # Set reactive value to TRUE
            proceed_to_load_integration(TRUE)
          }
        }
      )
    })
    
    # Return the reactive value
    return(list(proceed_to_load_integration = proceed_to_load_integration))
  })
}
