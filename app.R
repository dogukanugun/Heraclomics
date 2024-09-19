# Increase file upload size to 500 MB
options(shiny.maxRequestSize = 500 * 1024^2)

# Load required libraries
library(shiny)
library(shinyjs)
library(magrittr)
library(Matrix)
library(Seurat)
library(DT)
library(harmony)
library(SingleR)
library(celldex)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)
library(promises)
library(memoise)

# Set up parallel processing
plan(multisession, workers = 4)  # Adjust number of workers as needed

# Load modules
tryCatch({
  source("Scripts/load_data.R")
  source("Scripts/quality_control.R")
  source("Scripts/doublet_removal.R")
  source("Scripts/clustering.R")
  source("Scripts/integration_dimensionality_reduction.R")
  source("Scripts/second_dataset.R")
  source("Scripts/data_correction.R")
  source("Scripts/cluster_labeling.R")
  source("Scripts/gene_expression.R")
  source("Scripts/gene_coexpression.R")
  source("Scripts/gene_regulatory_networks.R")
  source("Scripts/gwas_associations.R")
  source("Scripts/trajectory_analysis.R")
  source("Scripts/multimodal_analysis.R")
  source("Scripts/differential_expression.R")
  source("Scripts/custom_differential_expression.R")
}, error = function(e) {
  message("Error loading modules: ", e$message)
})

# Define UI for application
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Heracleomics: Single-cell RNA-seq Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      navbarPage("Heracleomics",
                 tabPanel("Home", 
                          h3("Welcome to Heracleomics"),
                          p("A comprehensive tool for single-cell RNA-seq analysis."),
                          p("Choose 'New' to start a new analysis or 'Old' to load a previously saved environment."),
                          actionButton("start_new", "New"),
                          actionButton("continue_old", "Continue")
                 ),
                 tabPanel("Load Data", 
                          uiOutput("load_data_ui")
                 ),
                 # Add other tabPanels as needed
      )
    ),
    mainPanel(
      textOutput("current_step"),
      uiOutput("step_ui"),
      uiOutput("progress"),
      verbatimTextOutput("error_output")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive value for storing the loaded data for QC
  data_for_qc <- reactiveVal(NULL)
  
  # Initialize reactive values to store data for different stages
  data <- reactiveValues(
    for_qc = NULL,
    for_integration = NULL,
    for_clustering = NULL,
    for_correction = NULL,
    for_labeling = NULL,
    for_gene_expression = NULL,
    for_gene_coexpression = NULL,
    for_gene_regulatory_networks = NULL,
    for_gwas_associations = NULL,
    for_trajectory_analysis = NULL,
    for_cell_communication = NULL,
    for_multimodal_analysis = NULL,
    for_differential_expression = NULL,
    for_custom_differential_expression = NULL
  )
  
  # Helper function to ensure Seurat object and metadata compatibility
  ensure_seurat_metadata_compatibility <- function(seurat_obj, new_metadata) {
    cells_in_object <- colnames(seurat_obj)
    cells_in_metadata <- rownames(new_metadata)
    
    common_cells <- intersect(cells_in_object, cells_in_metadata)
    if (length(common_cells) == 0) {
      stop("No overlapping cells found between Seurat object and new metadata")
    }
    
    seurat_obj <- subset(seurat_obj, cells = common_cells)
    new_metadata <- new_metadata[common_cells, ]
    
    seurat_obj <- AddMetaData(seurat_obj, new_metadata)
    return(seurat_obj)
  }
  
  # Caching function for expensive operations
  cached_operation <- memoise(function(seurat_obj, operation) {
    operation(seurat_obj)
  })
  
  # Event to handle starting a new analysis
  observeEvent(input$start_new, {
    tryCatch({
      output$load_data_ui <- renderUI({
        load_data_ui("load_data")
      })
      
      # Assign the loaded data to the data_for_qc reactive value
      load_data_server("load_data", data_for_qc)
      
      output$current_step <- renderText("Step 1: Load Your Data")
    }, error = function(e) {
      output$error_output <- renderPrint({
        paste("Error in start_new:", e$message)
      })
    })
  })
  
  # Event to handle continuing an old analysis
  observeEvent(input$continue_old, {
    # Implement logic for loading previously saved environment
    output$current_step <- renderText("Loading Previous Environment...")
  })
  
  # Error handling output
  output$error_output <- renderPrint({
    if (!is.null(attr(session$output, "error"))) {
      paste("An error occurred:", attr(session$output, "error")$message)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
