library(shiny)
library(shinyjs)
library(Seurat)
library(Matrix)
library(tools)
library(future)
library(promises)

# Set up parallel processing
plan(multisession, workers = 4)  # Adjust number of workers as needed

# UI for the data loading step
load_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 1: Load Your Data"),
    radioButtons(ns("dataset_type"), "Dataset Type", 
                 choices = c("10X Data", "Seurat Object", "Count Matrix", "AnnData", "Loom", "SingleCellExperiment")),
    uiOutput(ns("file_input_ui")),
    actionButton(ns("submit_data"), "Submit"),
    br(),  # Add space
    tableOutput(ns("preview")),
    verbatimTextOutput(ns("status")),
    br(),  # Add space
    actionButton(ns("next_step"), "Continue to Quality Control", disabled = TRUE)
  )
}

# Server logic for data loading
load_data_server <- function(id, data_for_qc, go_to_qc) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    data <- reactiveVal(NULL)
    
    # Dynamically update file input based on selected dataset type
    output$file_input_ui <- renderUI({
      switch(input$dataset_type,
             "10X Data" = tagList(
               fileInput(ns("matrix"), "Choose matrix.mtx file", accept = c(".mtx", ".mtx.gz")),
               fileInput(ns("features"), "Choose features.tsv file", accept = c(".tsv", ".tsv.gz")),
               fileInput(ns("barcodes"), "Choose barcodes.tsv file", accept = c(".tsv", ".tsv.gz"))
             ),
             "Seurat Object" = fileInput(ns("datafile"), "Choose Seurat object file (.rds)", accept = ".rds"),
             "Count Matrix" = fileInput(ns("datafile"), "Choose Count Matrix (csv, txt, tsv)", accept = c(".txt", ".csv", ".tsv")),
             "AnnData" = fileInput(ns("datafile"), "Choose AnnData file (.h5ad)", accept = ".h5ad"),
             "Loom" = fileInput(ns("datafile"), "Choose Loom file (.loom)", accept = ".loom"),
             "SingleCellExperiment" = fileInput(ns("datafile"), "Choose SingleCellExperiment object (.rds)", accept = ".rds")
      )
    })
    
    # Handle data submission and processing
    observeEvent(input$submit_data, {
      shinyjs::disable("submit_data")
      output$status <- renderText("Loading data... Please wait.")
      
      # Use reactive values inside isolate to avoid dependency issues
      dataset_type <- isolate(input$dataset_type)
      
      # Prepare files based on dataset type
      files <- isolate({
        if (dataset_type == "10X Data") {
          req(input$matrix, input$features, input$barcodes)
          list(matrix = input$matrix, features = input$features, barcodes = input$barcodes)
        } else {
          req(input$datafile)
          input$datafile
        }
      })
      
      # Asynchronous data loading
      future_promise({
        tryCatch({
          loaded_data <- switch(dataset_type,
                                "10X Data" = read10x_files(files),
                                "Seurat Object" = readRDS(files$datapath),
                                "Count Matrix" = read_count_matrix(files$datapath, file_ext(files$name)),
                                "AnnData" = read_anndata(files$datapath),
                                "Loom" = read_loom(files$datapath),
                                "SingleCellExperiment" = readRDS(files$datapath)
          )
          
          if (!inherits(loaded_data, "Seurat")) {
            loaded_data <- CreateSeuratObject(counts = loaded_data)
          }
          
          list(success = TRUE, data = loaded_data)
        }, error = function(e) {
          list(success = FALSE, error = e$message)
        })
      }) %...>% 
        (function(result) {
          if (result$success) {
            data(result$data)
            data_for_qc(result$data)  # Update the reactive value
            
            output$preview <- renderTable({
              head(GetAssayData(result$data, slot = "counts")[1:10, 1:6])
            })
            
            output$status <- renderText("Data loaded successfully")
            shinyjs::enable("next_step")
          } else {
            output$status <- renderText(paste("Error loading data:", result$error))
          }
          
          shinyjs::enable("submit_data")
        })
    })
    
    # Navigate to the Quality Control step
    observeEvent(input$next_step, {
      go_to_qc()  # Call the function to move to the QC step
    })
    
    return(data)
  })
}
# Helper function for loading 10X data from individual files
# Helper function for loading 10X data from individual files
read10x_files <- function(files) {
  matrix_path <- files$matrix$datapath
  features_path <- files$features$datapath
  barcodes_path <- files$barcodes$datapath
  
  # Read the matrix
  matrix <- readMM(if(endsWith(matrix_path, ".gz")) gzfile(matrix_path) else matrix_path)
  
  # Read the features
  features <- read.table(
    if(endsWith(features_path, ".gz")) gzfile(features_path) else features_path,
    header = FALSE, stringsAsFactors = FALSE
  )
  
  # Clean feature names
  features$V2 <- gsub("_", "-", features$V2)  # Replace underscores with dashes
  
  # Read the barcodes
  barcodes <- read.table(
    if(endsWith(barcodes_path, ".gz")) gzfile(barcodes_path) else barcodes_path,
    header = FALSE, stringsAsFactors = FALSE
  )
  
  # Set row and column names
  rownames(matrix) <- make.unique(features$V2)
  colnames(matrix) <- barcodes$V1
  
  # Create the Seurat object
  CreateSeuratObject(
    counts = matrix,
    project = "scRNA_seq",
    min.cells = 3,
    min.features = 200
  )
}

# Helper function for reading count matrices
read_count_matrix <- function(file_path, ext) {
  data <- switch(tolower(ext),
                 "txt" = , "tsv" = read.delim(file_path, row.names = 1),
                 "csv" = read.csv(file_path, row.names = 1),
                 stop("Unsupported file format")
  )
  
  # Clean row names
  rownames(data) <- gsub("_", "-", rownames(data))  # Replace underscores with dashes
  as(as.matrix(data), "dgCMatrix")
}


# Helper function for reading AnnData files
read_anndata <- function(file_path) {
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("Package 'anndata' is required to read AnnData files. Please install it.")
  }
  adata <- anndata::read_h5ad(file_path)
  as(adata$X, "dgCMatrix")
}

# Helper function for reading Loom files
read_loom <- function(file_path) {
  if (!requireNamespace("loomR", quietly = TRUE)) {
    stop("Package 'loomR' is required to read Loom files. Please install it.")
  }
  loom <- loomR::connect(filename = file_path, mode = "r")
  counts <- t(loom$matrix[,])
  loom$close_all()
  as(counts, "dgCMatrix")
}
