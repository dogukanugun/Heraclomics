# Scripts/load_second_data.R

# Load necessary libraries
library(shiny)
library(shinyjs)
library(shinyalert)
library(Seurat)
library(DT)
library(shinyWidgets)

# Load Second Dataset Module UI Function
loadSecondDataUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),  # Enable shinyjs
    h2("Load Second Dataset"),
    fluidRow(
      box(
        title = "Data Input",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        # Input controls for data upload
        selectInput(
          ns("data_method2"),
          "Select Data Input Method:",
          choices = c("10X Genomics" = "tenx", "Seurat Object (RDS file)" = "rds", "Count Matrix (TXT)" = "txt")
        ),
        conditionalPanel(
          condition = "input.data_method2 == 'tenx'",
          ns = ns,
          fileInput(
            ns("tenx_files"),
            "Upload 10X Files (barcodes.tsv, features.tsv/genes.tsv, matrix.mtx):",
            multiple = TRUE,
            accept = c(".mtx", ".tsv", ".gz")
          )
        ),
        conditionalPanel(
          condition = "input.data_method2 == 'rds'",
          ns = ns,
          fileInput(
            ns("rds_file"),
            "Upload Seurat Object (RDS file):",
            accept = c(".rds")
          )
        ),
        conditionalPanel(
          condition = "input.data_method2 == 'txt'",
          ns = ns,
          fileInput(
            ns("count_matrix"),
            "Upload Count Matrix (TXT):",
            accept = c(".txt", ".csv")
          )
        ),
        textInput(ns("sample_name"), "Sample Name:", value = "Sample2"),
        numericInput(ns("min_cells"), "Filter genes detected in at least this many cells:", value = 3, min = 0),
        numericInput(ns("min_features"), "Filter cells with at least this many genes detected:", value = 200, min = 0),
        actionButton(ns("submit_data2"), "Submit", icon = icon("upload"))
      )
    ),
    fluidRow(
      box(
        title = "Data Preview",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        DT::dataTableOutput(ns("data_preview"))
      )
    ),
    fluidRow(
      box(
        title = "Proceed to Next Step",
        width = 12,
        status = "success",
        solidHeader = TRUE,
        collapsible = FALSE,
        actionButton(ns("skip_step"), "Skip this Step", icon = icon("forward")),
        actionButton(ns("continue_to_integration"), "Continue to Next Step", icon = icon("arrow-right"))
      )
    )
  )
}

# Load Second Dataset Module Server Function
loadSecondDataServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    proceed_to_integration <- reactiveVal(FALSE)
    
    # Reactive value to store the second Seurat object
    seurat_second <- reactiveVal(NULL)
    
    # Observe when the submit button is clicked
    observeEvent(input$submit_data2, {
      print("Submit button clicked for second dataset")
      
      tryCatch({
        # Ensure that the selected data method has valid input
        print(paste("Selected data method:", input$data_method2))
        
        # Handling 10X data
        if (input$data_method2 == "tenx") {
          print("Data method: 10X for second dataset")
          req(input$tenx_files)
          
          # Check if all required files are uploaded
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
              feature.column = 2,
              cell.column = 1,
              mtx.transpose = FALSE
            )
            print("Data read using ReadMtx for second dataset")
            
            # Create Seurat object
            seurat_second(
              CreateSeuratObject(
                counts = counts,
                project = input$sample_name,
                min.cells = input$min_cells,
                min.features = input$min_features
              )
            )
            print("Seurat object created for second dataset")
            rv$seurat_second <- seurat_second()
            rv$second_data_loaded <- TRUE
          } else {
            shinyalert(
              title = "Error",
              text = "Please upload barcodes.tsv, features.tsv or genes.tsv, and matrix.mtx files.",
              type = "error"
            )
          }
          
          # Handling RDS files
        } else if (input$data_method2 == "rds") {
          print("Data method: RDS for second dataset")
          req(input$rds_file)
          seurat_second(readRDS(input$rds_file$datapath))
          print("Seurat object loaded from RDS for second dataset")
          rv$seurat_second <- seurat_second()
          rv$second_data_loaded <- TRUE
          
          # Handling plain count matrix
        } else if (input$data_method2 == "txt") {
          print("Data method: TXT for second dataset")
          req(input$count_matrix)
          counts <- read.table(input$count_matrix$datapath, header = TRUE, row.names = 1)
          seurat_second(
            CreateSeuratObject(
              counts = counts,
              project = input$sample_name,
              min.cells = input$min_cells,
              min.features = input$min_features
            )
          )
          print("Seurat object created from count matrix for second dataset")
          rv$seurat_second <- seurat_second()
          rv$second_data_loaded <- TRUE
        } else {
          shinyalert(
            title = "Error",
            text = "Please upload a valid file.",
            type = "error"
          )
        }
        
        # If the second dataset is loaded successfully, update the data preview
        if (isTRUE(rv$second_data_loaded)) {
          print("rv$seurat_second updated")
          
          # Calculate percent.mt and percent.ribo
          rv$seurat_second[["percent.mt"]] <- PercentageFeatureSet(rv$seurat_second, pattern = "^MT-")
          rv$seurat_second[["percent.ribo"]] <- PercentageFeatureSet(rv$seurat_second, pattern = "^RPS|^RPL")
          print("QC metrics calculated")
          
          # Display a preview of the data
          output$data_preview <- DT::renderDataTable({
            # Use GetAssayData to access the counts matrix
            counts_data <- GetAssayData(rv$seurat_second, assay = "RNA", slot = "counts")
            
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
          print("Data preview rendered for second dataset")
          
          # Enable the "Continue to Integration" button
          shinyjs::enable(ns("continue_to_integration"))
          print("Continue button enabled")
          
          # Show a success alert using shinyalert
          shinyalert(
            title = "Success",
            text = "Second dataset loaded successfully! You can now proceed to the next step.",
            type = "success",
            timer = 3000,
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            showConfirmButton = FALSE
          )
        }
        
      }, error = function(e) {
        print(paste("Error in observeEvent for second dataset:", e$message))
        shinyalert(
          title = "Error",
          text = paste("Error:", e$message),
          type = "error"
        )
      })
    })
    
    # Observe the "Skip this Step" button
    observeEvent(input$skip_step, {
      print("Skip button clicked for second dataset")
      shinyalert(
        title = "Skipping Step",
        text = "You have chosen to skip the second dataset. Proceeding to the next step.",
        type = "info",
        timer = 2000,
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        showConfirmButton = FALSE
      )
      # Ensure rv$seurat_second remains NULL
      rv$seurat_second <- NULL
      proceed_to_integration(TRUE)
    })
    
    # Observe the "Continue to Next Step" button for the second dataset
    observeEvent(input$continue_to_integration, {
      shinyalert(
        title = "Proceed to Next Step?",
        text = "Are you sure you want to proceed to the next step?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, proceed",
        cancelButtonText = "No, stay here",
        callbackR = function(value) {
          if (isTRUE(value)) {
            # Proceed to the next step
            proceed_to_integration(TRUE)
            print("Proceeding to next step after second dataset.")
          }
        }
      )
    })
    
    return(list(proceed_to_integration = proceed_to_integration))
  })
}
