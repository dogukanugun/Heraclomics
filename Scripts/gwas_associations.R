# gwas_associations.R

# Load necessary libraries
library(shiny)
library(MAGMA.Celltyping)  # Example library for GWAS and MAGMA
library(DT)
library(dplyr)
library(shinyalert)

# UI Function for GWAS Associations Module
gwasAssociationsUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyalert(),
    h2("GWAS Associations"),
    fluidRow(
      box(
        title = "GWAS Associations Overview",
        width = 12,
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        p("Genome-wide association studies have identified several loci in various genes associated with traits of interest."),
        p("The expression profiles of cell clusters can be compared against these GWAS loci to identify potential causal cell types underlying complex traits."),
        p("This module employs the MAGMA methodology to identify increased linear associations between cluster-derived gene sets and GWAS traits."),
        p("You may upload your own GWAS summary statistics or select from curated datasets.")
      )
    ),
    fluidRow(
      box(
        title = "GWAS Input",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fileInput(ns("gwas_file"), "Upload GWAS Summary Statistics", accept = ".txt"),
        selectInput(ns("gwas_dataset"), "Select GWAS Trait", choices = NULL, multiple = FALSE),
        helpText("You can either upload a custom GWAS file or choose from available curated datasets."),
        actionButton(ns("run_analysis"), "Run GWAS Analysis", icon = icon("play"), class = "btn-success")
      )
    ),
    fluidRow(
      box(
        title = "GWAS Results",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        DTOutput(ns("gwas_results")),
        downloadButton(ns("download_gwas_results"), "Download Results")
      )
    )
  )
}

# Server Function for GWAS Associations Module
gwasAssociationsServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Placeholder function for loading GWAS datasets
    load_gwas_datasets <- function() {
      # Replace with real datasets
      return(c("Alzheimer's Disease", "Parkinson's Disease", "Schizophrenia"))
    }
    
    # Update available GWAS datasets on startup
    observe({
      available_datasets <- load_gwas_datasets()
      updateSelectInput(session, "gwas_dataset", choices = available_datasets)
    })
    
    # Placeholder for running the GWAS analysis with MAGMA
    run_gwas_analysis <- function(gwas_data, seurat_obj) {
      # Replace with MAGMA methodology, include necessary data manipulations
      return(data.frame(
        Gene = c("GeneA", "GeneB", "GeneC"),
        P_Value = c(0.001, 0.05, 0.01),
        Association = c("Positive", "Negative", "Positive")
      ))
    }
    
    # Action to run GWAS analysis
    observeEvent(input$run_analysis, {
      req(rv$corrected_data)
      shinyalert::shinyalert(
        title = "Confirm GWAS Analysis",
        text = "Proceed with the GWAS analysis?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, Proceed",
        cancelButtonText = "Cancel",
        callbackR = function(x) {
          if (x) {
            # Example processing, replace with actual file read and MAGMA analysis
            gwas_file_data <- NULL
            if (!is.null(input$gwas_file)) {
              gwas_file_data <- read.table(input$gwas_file$datapath, header = TRUE)
            }
            
            results <- run_gwas_analysis(gwas_file_data, rv$corrected_data)
            output$gwas_results <- renderDT({ results })
          }
        }
      )
    })
    
    # Allow users to download GWAS results
    output$download_gwas_results <- downloadHandler(
      filename = function() { paste0("GWAS_Results_", Sys.Date(), ".csv") },
      content = function(file) {
        results <- run_gwas_analysis(NULL, rv$corrected_data)
        write.csv(results, file, row.names = FALSE)
      }
    )
  })
}
