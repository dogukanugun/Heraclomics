# UI for the data loading step
load_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 1: Load Your Data"),
    radioButtons(ns("dataset_type"), "Dataset Type", choices = c("10X Data", "Seurat Object", "Count Matrix")),
    uiOutput(ns("file_input_ui")), # Dynamic UI for file input based on dataset type
    actionButton(ns("submit_data"), "Submit"),
    tableOutput(ns("preview")), # Preview of the first few rows of the data
    verbatimTextOutput(ns("status")), # Status message for data loading
    actionButton(ns("next_step"), "Continue to Quality Control", disabled = TRUE)
  )
}

# Server logic for data loading
load_data_server <- function(id, data_for_qc) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    data <- reactiveVal(NULL)
    
    # Dynamically update file input based on selected dataset type
    observeEvent(input$dataset_type, {
      output$file_input_ui <- renderUI({
        if (input$dataset_type == "10X Data") {
          fileInput(ns("datafile"), "Choose 10X data folder (filtered feature matrices)", 
                    multiple = TRUE, accept = c(".gz", ".tsv", ".mtx"))
        } else if (input$dataset_type == "Seurat Object") {
          fileInput(ns("datafile"), "Choose Seurat object file (.rds)", 
                    multiple = FALSE, accept = ".rds")
        } else {
          fileInput(ns("datafile"), "Choose Count Matrix (csv, txt, tsv)", 
                    multiple = FALSE, accept = c(".txt", ".csv", ".tsv"))
        }
      })
    })
    
    # Handle data submission and processing
    observeEvent(input$submit_data, {
      file <- input$datafile
      if (is.null(file)) {
        output$status <- renderText("No file uploaded")
        return()
      }
      
      tryCatch({
        # Process data based on dataset type
        if (input$dataset_type == "10X Data") {
          # Handle 10X Genomics data folder input
          data_dir <- dirname(file$datapath[1])
          data(read10x(data_dir))
        } else if (input$dataset_type == "Seurat Object") {
          # Handle Seurat object file
          data(readRDS(file$datapath))
        } else {
          # Handle Count Matrix (txt, csv, tsv)
          ext <- tools::file_ext(file$name)
          if (ext == "txt" || ext == "tsv") {
            data(read.delim(file$datapath))
          } else if (ext == "csv") {
            data(read.csv(file$datapath))
          } else {
            stop("Unsupported file format")
          }
        }
        
        # Save data for quality control
        data_for_qc(data())
        
        # Preview the data
        output$preview <- renderTable({
          head(data(), 10)
        })
        
        # Success message
        output$status <- renderText("Data loaded successfully")
        
        # Enable the "next step" button
        shinyjs::enable("next_step")
        
      }, error = function(e) {
        # Handle errors gracefully
        output$status <- renderText(paste("Error loading data:", e$message))
      })
    })
  })
}

# Helper function for loading 10X data using Seurat
read10x <- function(data_dir) {
  Seurat::Read10X(data.dir = data_dir)
}
