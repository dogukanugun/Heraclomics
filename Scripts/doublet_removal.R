doublet_removal_server <- function(input, output, session, app_state) {
  ns <- session$ns
  
  cleaned_data <- reactiveVal(NULL)
  doublet_classification_column <- reactiveVal(NULL)
  
  # Disable controls initially
  observe({
    shinyjs::toggleState("skip_doublet_removal", !is.null(app_state$data))
    shinyjs::toggleState("simulate_doublets", !is.null(app_state$data))
    shinyjs::toggleState("remove_doublets", !is.null(cleaned_data()))
    shinyjs::toggleState("next_step", !is.null(cleaned_data()))
  })
  
  observeEvent(input$start_doublet_removal, {
    req(app_state$data)
    withProgress(message = 'Processing data...', value = 0, {
      data <- app_state$data
      
      # Preprocess data: normalization, scaling, PCA, UMAP, t-SNE
      incProgress(0.2, detail = "Normalizing data")
      data <- NormalizeData(data)
      data <- FindVariableFeatures(data)
      data <- ScaleData(data)
      
      # Run PCA
      incProgress(0.2, detail = "Running PCA")
      data <- RunPCA(data)
      
      # Run UMAP and t-SNE
      incProgress(0.3, detail = "Running UMAP and t-SNE")
      data <- RunUMAP(data, dims = 1:10)
      data <- RunTSNE(data, dims = 1:10)
      
      app_state$data <- data  # Update the app state with the modified data
      cleaned_data(data)  # Set initial cleaned data
    })
  })
  
  observeEvent(input$skip_doublet_removal, {
    req(app_state$data)
    cleaned_data(app_state$data)
  })
  
  observeEvent(input$simulate_doublets, {
    shinyjs::disable("simulate_doublets")
    req(cleaned_data())
    
    withProgress(message = 'Simulating doublets...', value = 0, {
      data <- cleaned_data()
      
      incProgress(0.3, detail = "Determining optimal parameters")
      tryCatch({
        sweep.res <- paramSweep(data, PCs = 1:10, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        
        if (nrow(bcmvn) == 0) {
          stop("No pK values were computed.")
        }
        
        # Find optimal pK
        optimal_pK <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
        
        # Validate optimal pK
        if (is.na(optimal_pK) || optimal_pK <= 0 || optimal_pK > 1) {
          optimal_pK <- 0.1
          showNotification("Using fallback pK of 0.1 due to invalid optimal pK.", type = "warning")
        }
        
        # Doublet detection
        homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
        nExp_poi <- round(0.075 * ncol(data))
        nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
        
        data <- DoubletFinder::doubletFinder(data, PCs = 1:10, pN = 0.25, pK = optimal_pK, 
                                             nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
        
        doublet_classification_col <- grep("DF.classifications", colnames(data@meta.data), value = TRUE)
        doublet_classification_column(doublet_classification_col)
        cleaned_data(data)
        
        showNotification(paste("Doublet simulation completed successfully. Optimal pK:", optimal_pK), type = "message")
      }, error = function(e) {
        showNotification(paste("Error in doublet simulation:", e$message), type = "error", duration = NULL)
        message("Error in doublet simulation: ", e$message)
      })
    })
  })
  
  observeEvent(input$remove_doublets, {
    req(cleaned_data())
    
    data <- cleaned_data()
    doublet_classification_col <- doublet_classification_column()
    
    # Filter out doublets
    data_filtered <- subset(data, get(doublet_classification_col) == "Singlet")
    cleaned_data(data_filtered)
  })
  
  output$doublet_plot <- renderPlotly({
    req(cleaned_data())
    
    data <- cleaned_data()
    
    plot_type <- input$plot_type
    doublet_classification_col <- doublet_classification_column()
    
    # Check the plot type and create the plot
    if (plot_type %in% c("UMAP-2D", "t-SNE-2D")) {
      reduction <- ifelse(plot_type == "UMAP-2D", "umap", "tsne")
      plot <- DimPlot(data, reduction = reduction, group.by = doublet_classification_col, pt.size = input$point_size) +
        theme(text = element_text(size = input$label_size))
      ggplotly(plot)
    } else if (plot_type %in% c("UMAP-3D", "t-SNE-3D")) {
      reduction <- ifelse(plot_type == "UMAP-3D", "umap", "tsne")
      
      # Extract the 3D embedding
      embeddings <- data@reductions[[reduction]]@cell.embeddings
      
      # Check that we have at least 3 dimensions
      if (ncol(embeddings) >= 3) {
        plot_ly(data = as.data.frame(embeddings), 
                x = ~V1, y = ~V2, z = ~V3, 
                color = ~data@meta.data[[doublet_classification_col]], 
                colors = c("blue", "red")) %>%
          layout(scene = list(xaxis = list(title = "X1"),
                              yaxis = list(title = "X2"),
                              zaxis = list(title = "X3")))
      } else {
        showNotification("Not enough dimensions for 3D plot.", type = "error")
      }
    }
  })
  
  observeEvent(input$next_step, {
    updateTabItems(session, "tabs", "clustering")
  })
  
  return(cleaned_data)
}
