# Load necessary libraries
library(shiny)
library(shinyjs)
library(shinyalert)
library(Seurat)
library(DT)
library(tensorflow)
library(keras)

# Load Data Module UI Function
loadDataUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),     # Enable shinyjs
    # No need for useShinyalert()
    h2("Load Your Data"),
    fluidRow(
      box(
        title = "Data Input Method",
        width = 4,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        radioButtons(
          ns("data_method"),
          "Choose your load method:",
          choices = c("10X Data" = "tenx", "Seurat Object (RDS file)" = "rds", "Count Matrix (txt)" = "txt")
        ),
        conditionalPanel(
          condition = paste0("input['", ns("data_method"), "'] == 'tenx'"),
          fileInput(
            ns("tenx_files"),
            "Upload 10X Data Files (barcodes.tsv, features.tsv/genes.tsv, matrix.mtx):",
            multiple = TRUE,
            accept = c(".gz", ".tsv", ".mtx")
          )
        ),
        conditionalPanel(
          condition = paste0("input['", ns("data_method"), "'] == 'rds'"),
          fileInput(ns("rds_file"), "Upload Seurat Object (RDS file):", accept = ".rds")
        ),
        conditionalPanel(
          condition = paste0("input['", ns("data_method"), "'] == 'txt'"),
          fileInput(ns("count_matrix"), "Upload Count Matrix (txt):", accept = c(".txt", ".tsv"))
        ),
        actionButton(ns("submit_data"), "Submit", icon = icon("upload"))
      ),
      box(
        title = "Sample Information",
        width = 4,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        textInput(ns("sample_name"), "Sample Name:", value = "Sample1"),
        selectInput(
          ns("sample_type"),
          "Select single or multiple samples:",
          choices = c("Single Sample", "Multiple Samples")
        ),
        numericInput(
          ns("min_cells"),
          "Include genes detected in at least this many cells:",
          value = 3,
          min = 1
        ),
        numericInput(
          ns("min_features"),
          "Include cells where at least this many genes are detected:",
          value = 200,
          min = 1
        )
      ),
      box(
        title = "Data Preview",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        DT::dataTableOutput(ns("data_preview"))
      ),
      box(
        title = "Continue to Next Step",
        width = 12,
        status = "success",
        solidHeader = TRUE,
        collapsible = FALSE,
        collapsed = FALSE,
        
        actionButton(ns("continue_to_qc"), "Continue to Next Step", icon = icon("arrow-right")
        )
      )
    )
  )
}

# Load Data Module Server Function
loadDataServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    proceed_to_qc <- reactiveVal(FALSE)
    
    # Reactive value to store the Seurat object
    seurat_obj <- reactiveVal(NULL)
    
    # Observe when the submit button is clicked
    observeEvent(input$submit_data, {
      print("Submit button clicked")
      tryCatch({
        # Check which data method is selected
        if (input$data_method == "tenx") {
          print("Data method: 10X")
          # Ensure all files are uploaded
          if (!is.null(input$tenx_files)) {
            files <- input$tenx_files$datapath
            filenames <- input$tenx_files$name
            print(paste("Uploaded files:", filenames))
            
            # Identify files by name
            barcodes_idx <- grep("barcodes", filenames, ignore.case = TRUE)
            features_idx <- grep("features|genes", filenames, ignore.case = TRUE)
            matrix_idx <- grep("matrix", filenames, ignore.case = TRUE)
            
            if (length(barcodes_idx) == 1 && length(features_idx) == 1 && length(matrix_idx) == 1) {
              barcodes_file <- files[barcodes_idx]
              features_file <- files[features_idx]
              matrix_file <- files[matrix_idx]
              
              # Read the data using ReadMtx
              counts <- ReadMtx(
                mtx = matrix_file,
                features = features_file,
                cells = barcodes_file,
                feature.column = 2,  # Adjust based on your features file
                cell.column = 1,
                mtx.transpose = FALSE
              )
              print("Data read using ReadMtx")
              
              # Create Seurat object
              seurat_obj(
                CreateSeuratObject(
                  counts = counts,
                  project = input$sample_name,
                  min.cells = input$min_cells,
                  min.features = input$min_features
                )
              )
              print("Seurat object created")
              rv$data_loaded <- TRUE
              print("rv$data_loaded set to TRUE")
            } else {
              shinyalert(
                title = "Error",
                text = "Please upload barcodes.tsv, features.tsv or genes.tsv, and matrix.mtx files.",
                type = "error"
              )
            }
          } else {
            shinyalert(
              title = "Error",
              text = "Please upload the required files.",
              type = "error"
            )
          }
        } else if (input$data_method == "rds") {
          print("Data method: RDS")
          if (!is.null(input$rds_file)) {
            seurat_obj(readRDS(input$rds_file$datapath))
            rv$data_loaded <- TRUE
            print("Seurat object loaded from RDS")
          } else {
            shinyalert(
              title = "Error",
              text = "Please upload an RDS file.",
              type = "error"
            )
          }
        } else if (input$data_method == "txt") {
          print("Data method: TXT")
          if (!is.null(input$count_matrix)) {
            counts <- read.table(input$count_matrix$datapath, header = TRUE, row.names = 1)
            seurat_obj(
              CreateSeuratObject(
                counts = counts,
                project = input$sample_name,
                min.cells = input$min_cells,
                min.features = input$min_features
              )
            )
            rv$data_loaded <- TRUE
            print("Seurat object created from count matrix")
          } else {
            shinyalert(
              title = "Error",
              text = "Please upload a count matrix file.",
              type = "error"
            )
          }
        }
        
        # If data is loaded successfully, update the data preview
        if (isTRUE(rv$data_loaded)) {
          # Update the reactive value in rv
          rv$seurat_object <- seurat_obj()
          print("rv$seurat_object updated")
          
          # Calculate percent.mt and percent.ribo
          rv$seurat_object[["percent.mt"]] <- PercentageFeatureSet(rv$seurat_object, pattern = "^MT-")
          rv$seurat_object[["percent.ribo"]] <- PercentageFeatureSet(rv$seurat_object, pattern = "^RPS|^RPL")
          print("QC metrics calculated")
          
          # Display a preview of the data
          output$data_preview <- DT::renderDataTable({
            # Use GetAssayData to access the counts matrix
            counts_data <- GetAssayData(rv$seurat_object, assay = "RNA", slot = "counts")
            
            # Subset the counts data
            counts_data_subset <- counts_data[1:min(500, nrow(counts_data)), 1:min(10, ncol(counts_data))]
            
            # Convert to a dense matrix for display purposes
            counts_matrix <- as.matrix(counts_data_subset)
            
            # Convert to a data frame and add gene names as a column
            counts_df <- data.frame(
              Gene = rownames(counts_matrix),
              counts_matrix,
              check.names = FALSE
            )
            
            # Render the data table
            DT::datatable(
              counts_df,
              options = list(scrollX = TRUE)
            )
          })
          print("Data preview rendered")
          
          # Enable the "Continue to Next Step" button
          shinyjs::enable(ns("continue_to_qc"))
          print("Continue button enabled")
          
          # Show a success alert using shinyalert
          shinyalert(
            title = "Success",
            text = "Data loaded successfully! You can now proceed to the next step.",
            type = "success",
            timer = 3000,
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            showConfirmButton = FALSE
          )
        }
      }, error = function(e) {
        print(paste("Error in observeEvent:", e$message))
        # Show an error alert using shinyalert
        shinyalert(
          title = "Error",
          text = paste("Error:", e$message),
          type = "error"
        )
      })
    })
    
    # Observe the "Continue to Next Step" button
    observeEvent(input$continue_to_qc, {
      shinyalert(
        title = "Proceed to Quality Control?",
        text = "Are you sure you want to proceed to the next step?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, proceed",
        cancelButtonText = "No, stay here",
        callbackR = function(value) {
          if (value) {
            # User confirmed, set proceed_to_qc to TRUE
            proceed_to_qc(TRUE)
          }
        }
      )
    })
    
    # Return the reactive value
    return(list(proceed_to_qc = proceed_to_qc))
  })
}