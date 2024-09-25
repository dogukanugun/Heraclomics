library(shiny)
library(shinyjs)
library(Seurat)

second_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    h3("Step 4: Load Second Dataset"),
    
    fluidRow(
      column(12, radioButtons(ns("data_input_method"), "Select Data Input Method", choices = c("Tab Delimited Table", "10X Data Folder", "Saved Seurat Object"))),
      column(12, conditionalPanel(
        condition = sprintf("input['%s'] == 'Tab Delimited Table'", ns("data_input_method")),
        fileInput(ns("tab_delimited_file"), "Upload Tab Delimited Table", accept = ".txt")
      )),
      column(12, conditionalPanel(
        condition = sprintf("input['%s'] == '10X Data Folder'", ns("data_input_method")),
        fileInput(ns("data_folder"), "Upload 10X Data Folder", multiple = TRUE)
      )),
      column(12, conditionalPanel(
        condition = sprintf("input['%s'] == 'Saved Seurat Object'", ns("data_input_method")),
        fileInput(ns("seurat_object_file"), "Upload Seurat Object (RDS file)", accept = ".rds")
      )),
      column(6, actionButton(ns("start_upload"), "Start Upload")),
      column(6, actionButton(ns("skip_upload"), "Skip Upload"))
    )
  )
}

second_dataset_server <- function(input, output, session, app_state) {
  ns <- session$ns
  
  observeEvent(input$start_upload, {
    req(input$data_input_method)  # Ensure input is available
    
    withProgress(message = 'Loading second dataset...', value = 0, {
      tryCatch({
        # Initialize second_dataset
        second_dataset <- NULL
        
        if (input$data_input_method == "Tab Delimited Table") {
          req(input$tab_delimited_file)
          incProgress(0.5, detail = "Reading tab-delimited file...")
          data <- read.table(input$tab_delimited_file$datapath, header = TRUE, sep = "\t", row.names = 1)
          second_dataset <- CreateSeuratObject(counts = data)
          
        } else if (input$data_input_method == "10X Data Folder") {
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
          incProgress(0.5, detail = "Loading 10X data...")
          seurat_obj <- Read10X(data.dir = temp_dir)
          
          # Replace underscores with dashes in feature names
          rownames(seurat_obj) <- gsub("_", "-", rownames(seurat_obj))
          
          # Create a Seurat object
          second_dataset <- CreateSeuratObject(counts = seurat_obj, 
                                               project = "Second Dataset",  # You can customize this
                                               min.cells = 3,  # Example threshold; modify as needed
                                               min.features = 200)  # Example threshold; modify as needed
          
        } else if (input$data_input_method == "Saved Seurat Object") {
          req(input$seurat_object_file)
          incProgress(0.5, detail = "Reading saved Seurat object...")
          second_dataset <- readRDS(input$seurat_object_file$datapath)
        }
        
        # Check if second_dataset was created successfully
        req(second_dataset)  # Ensure second_dataset is not NULL
        
        # Integrate the second dataset with the first dataset
        incProgress(0.5, detail = "Integrating datasets...")
        integrated_data <- merge(app_state$data, y = second_dataset)
        app_state$data <- integrated_data
        
        showNotification("Second dataset loaded and integrated successfully.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
  
  observeEvent(input$skip_upload, {
    showNotification("Skipped loading second dataset.", type = "message")
  })
}
