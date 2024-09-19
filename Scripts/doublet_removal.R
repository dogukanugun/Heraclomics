library(shiny)
library(Seurat)
library(DoubletFinder)
library(shinyjs)

doublet_removal_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 3: Doublet Removal"),
    actionButton(ns("remove_doublets"), "Remove Doublets"),
    plotOutput(ns("doublet_plot")),
    actionButton(ns("next_step"), "Continue to Clustering", disabled = TRUE)
  )
}

doublet_removal_server <- function(id, data_for_doublets) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    cleaned_data <- reactiveVal(NULL)
    
    observeEvent(input$remove_doublets, {
      req(data_for_doublets())
      
      tryCatch({
        data <- data_for_doublets()
        
        # Placeholder for determining optimal pK value for DoubletFinder
        sweep.res <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
        
        # Running DoubletFinder with estimated homotypic proportion
        homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
        nExp_poi <- round(0.075 * ncol(data))  # Adjust doublet rate as needed
        nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
        
        # Identifying doublets
        data <- DoubletFinder::doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
        
        # Filter out doublets
        cleaned_data <- subset(data, subset = DF.classifications_0.25_0.75 == "Singlet")  # Use correct naming from DoubletFinder
        
        # Visualize the results (nFeature_RNA vs nCount_RNA after doublet removal)
        output$doublet_plot <- renderPlot({
          FeatureScatter(cleaned_data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", 
                         pt.size = 0.5) + ggtitle("After Doublet Removal")
        })
        
        # Store cleaned data for the next step
        cleaned_data(cleaned_data)
        shinyjs::enable("next_step")
        
      }, error = function(e) {
        showNotification(paste("Error during doublet removal:", e$message), type = "error")
      })
    })
    
    # Move to next step with cleaned data
    observeEvent(input$next_step, {
      req(cleaned_data())
      data_for_clustering(cleaned_data())  # Pass cleaned data to clustering
    })
  })
}
