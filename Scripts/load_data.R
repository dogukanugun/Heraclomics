# Load required libraries
library(shiny)
library(shinyjs)
library(Seurat)

# UI for Load Data step
load_data_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Load Your Data"),
    
    # Selection for single or multiple samples
    radioButtons(ns("sample_type"), "Select sample type:",
                 choices = c("Single Sample", "Multiple Samples")),
    
    # Data input methods
    radioButtons(ns("input_method"), "Choose your data input method:",
                 choices = c("10X Data Folder", "Seurat Object (RDS)", 
                             "Count Matrix (txt)")),
    
    # Upload files
    conditionalPanel(
      condition = sprintf("input['%s'] == '10X Data Folder'", ns("input_method")),
      fileInput(ns("data_folder"), "Upload 10X Data Folder", multiple = TRUE)
    ),
    
    conditionalPanel(
      condition = sprintf("input['%s'] == 'Seurat Object (RDS)'", ns("input_method")),
      fileInput(ns("seurat_file"), "Upload Seurat RDS file")
    ),
    
    conditionalPanel(
      condition = sprintf("input['%s'] == 'Count Matrix (txt)'", ns("input_method")),
      fileInput(ns("count_matrix"), "Upload Count Matrix (txt) file")
    ),
    
    # Label experiment and filter criteria
    textInput(ns("sample_name"), "Sample Name:"),
    numericInput(ns("min_cells"), 
                 "Include genes detected in at least this many cells:", 3),
    numericInput(ns("min_genes"), 
                 "Include cells where at least this many genes are detected:", 
                 200),
    
    # Preview data
    actionButton(ns("preview"), "Preview Data"),
    tableOutput(ns("data_preview")),
    
    # Submit button
    actionButton(ns("submit_data"), "Submit Data")
  )
}

# Server logic for Load Data step
load_data_server <- function(input, output, session, app_state) {
  ns <- session$ns
  
  # Reactive value to store the loaded data
  loaded_data <- reactiveVal(NULL)
  
  # Preview data
  observeEvent(input$preview, {
    if (input$input_method == "10X Data Folder") {
      showNotification("Previewing 10X Data Folder...", type = "message")
      
      # Ensure the folder is uploaded
      req(input$data_folder)
      
      # Create a temporary directory to store the uploaded files
      temp_dir <- tempdir()
      
      # Copy the uploaded files to the temporary directory
      for (i in 1:nrow(input$data_folder)) {
        file.copy(input$data_folder$datapath[i], file.path(temp_dir, input$data_folder$name[i]))
      }
      
      # Check if the required files are present
      required_files <- c("barcodes.tsv", "genes.tsv", "matrix.mtx")
      uploaded_files <- list.files(temp_dir)
      if (!all(required_files %in% uploaded_files)) {
        showNotification("Missing required 10X files.", type = "error")
        return(NULL)
      }
      
      # Load the 10X data
      seurat_obj <- Read10X(data.dir = temp_dir)
      
      # Replace underscores with dashes in feature names
      rownames(seurat_obj) <- gsub("_", "-", rownames(seurat_obj))
      
      # Create a Seurat object
      seurat_obj <- CreateSeuratObject(counts = seurat_obj, 
                                       project = input$sample_name, 
                                       min.cells = input$min_cells, 
                                       min.features = input$min_genes)
      
      # Store the loaded data
      loaded_data(seurat_obj)
      
      # Render a preview of the data
      output$data_preview <- renderTable({
        head(seurat_obj@meta.data)
      })
    } else if (input$input_method == "Seurat Object (RDS)") {
      showNotification("Previewing Seurat Object...", type = "message")
      req(input$seurat_file)
      seurat_obj <- readRDS(input$seurat_file$datapath)
      loaded_data(seurat_obj)
      output$data_preview <- renderTable({
        head(seurat_obj@meta.data)
      })
    } else if (input$input_method == "Count Matrix (txt)") {
      showNotification("Previewing Count Matrix...", type = "message")
      req(input$count_matrix)
      count_matrix <- read.table(input$count_matrix$datapath, 
                                 header = TRUE, row.names = 1)
      rownames(count_matrix) <- gsub("_", "-", rownames(count_matrix))
      loaded_data(CreateSeuratObject(counts = count_matrix, 
                                     project = input$sample_name, 
                                     min.cells = input$min_cells, 
                                     min.features = input$min_genes))
      output$data_preview <- renderTable({
        head(count_matrix)
      })
    }
  })
  
  # Submit data
  observeEvent(input$submit_data, {
    showNotification("Data submitted successfully!", type = "message")
    app_state$data <- loaded_data()  # Store the loaded data in app_state
    shinyjs::enable("next_step")  # Enable the next step button
  })
}
