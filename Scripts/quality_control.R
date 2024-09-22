# Load required libraries
library(shiny)
library(shinyjs)
library(Seurat)

# UI for Quality Control step
quality_control_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Quality Control"),
    
    # Display scatter plots
    plotOutput(ns("scatter_plot1")),
    plotOutput(ns("scatter_plot2")),
    
    # Display QC summary
    verbatimTextOutput(ns("qc_summary")),
    
    # Button to proceed to the next step
    actionButton(ns("next_step"), "Next Step", disabled = TRUE)
  )
}

# Server logic for Quality Control step
quality_control_server <- function(input, output, session, app_state) {
  ns <- session$ns
  
  observe({
    req(app_state$data)
    
    # Perform quality control
    data_for_qc <- app_state$data
    
    # Calculate percent.mt
    data_for_qc[["percent.mt"]] <- PercentageFeatureSet(data_for_qc, pattern = "^MT-")
    
    # Filter cells based on quality control metrics
    filtered <- subset(data_for_qc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    
    # Generate scatter plots for nCount_RNA vs nFeature_RNA and nCount_RNA vs percent.mt
    output$scatter_plot1 <- renderPlot({
      FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    })
    
    output$scatter_plot2 <- renderPlot({
      FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
    })
    
    # Display a summary of the QC process
    output$qc_summary <- renderPrint({
      cat("QC Summary:\n")
      cat("Total cells before filtering:", ncol(data_for_qc), "\n")
      cat("Total cells after filtering:", ncol(filtered), "\n")
      cat("Cells removed:", ncol(data_for_qc) - ncol(filtered), "\n")
      cat("Genes detected:", nrow(filtered), "\n")
    })
    
    shinyjs::enable("next_step")
  })
  
  # Return the filtered data for the next step
  return(reactive({ filtered }))
}
