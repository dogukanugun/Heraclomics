# second_dataset.R

library(shiny)
library(Seurat)

# UI for loading a second dataset
load_second_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 4: Load Second Dataset (Optional)"),
    p("Select a method to load your second dataset, or skip this step."),
    radioButtons(ns("second_dataset_type"), "Dataset Type", 
                 choices = c("Tab-delimited Table", "10X Data", "Saved Seurat Object")),
    uiOutput(ns("second_file_input_ui")), # Dynamic file input based on dataset type
    actionButton(ns("submit_second_data"), "Submit Second Dataset"),
    actionButton(ns("skip_step"), "Skip Step"),
    tableOutput(ns("second_preview")),
    verbatimTextOutput(ns("second_status")),
    actionButton(ns("next_step"), "Continue to Integration", disabled = TRUE)
  )
}

# Server logic for loading a second dataset
load_second_dataset_server <- function(id, data_for_integration) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    second_data <- reactiveVal(NULL)
    
    # Dynamically change file input based on selected dataset type
    observeEvent(input$second_dataset_type, {
      output$second_file_input_ui <- renderUI({
        if (input$second_dataset_type == "Tab-delimited Table") {
          fileInput(ns("second_datafile"), "Choose tab-delimited table", multiple = FALSE, accept = c(".txt", ".tsv"))
        } else if (input$second_dataset_type == "10X Data") {
          fileInput(ns("second_datafile"), "Choose 10X data folder", multiple = FALSE, accept = c(".gz", ".tsv", ".mtx"))
        } else if (input$second_dataset_type == "Saved Seurat Object") {
          fileInput(ns("second_datafile"), "Choose Seurat object file", multiple = FALSE, accept = c(".rds"))
        }
      })
    })
    
    # Handle file upload based on selected dataset type
    observeEvent(input$submit_second_data, {
      file <- input$second_datafile
      if (is.null(file)) {
        output$second_status <- renderText("No file uploaded")
        return()
      }
      
      # Load data based on dataset type
      if (input$second_dataset_type == "Tab-delimited Table") {
        second_data(read.delim(file$datapath))
      } else if (input$second_dataset_type == "10X Data") {
        data_dir <- dirname(file$datapath)
        second_data(Seurat::Read10X(data.dir = data_dir))
      } else {
        second_data(readRDS(file$datapath))
      }
      
      data_for_integration(second_data()) # Update reactive value
      
      # Preview the uploaded second dataset
      output$second_preview <- renderTable({
        head(second_data(), 10)
      })
      
      output$second_status <- renderText("Second dataset loaded successfully")
      shinyjs::enable("next_step")
    })
    
    # Skip step functionality
    observeEvent(input$skip_step, {
      output$second_status <- renderText("Skipping second dataset loading")
      shinyjs::enable("next_step")
    })
  })
}
