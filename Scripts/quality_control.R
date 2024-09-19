library(shiny)
library(Seurat)
library(ggplot2)
library(shinyjs)

quality_control_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 2: Quality Control"),
    fluidRow(
      column(3, numericInput(ns("min_nFeature_RNA"), "Min nFeature_RNA:", min = 0, value = 200)),
      column(3, numericInput(ns("max_nFeature_RNA"), "Max nFeature_RNA:", min = 0, value = 6000)),
      column(3, numericInput(ns("min_nCount_RNA"), "Min nCount_RNA:", min = 0, value = 500)),
      column(3, numericInput(ns("max_nCount_RNA"), "Max nCount_RNA:", min = 0, value = 50000))
    ),
    fluidRow(
      column(3, numericInput(ns("min_percent_mt"), "Min percent.mt:", min = 0, max = 100, value = 0, step = 0.1)),
      column(3, numericInput(ns("max_percent_mt"), "Max percent.mt:", min = 0, max = 100, value = 10, step = 0.1)),
      column(3, numericInput(ns("min_percent_ribo"), "Min percent.ribo:", min = 0, max = 100, value = 0, step = 0.1)),
      column(3, numericInput(ns("max_percent_ribo"), "Max percent.ribo:", min = 0, max = 100, value = 50, step = 0.1))
    ),
    actionButton(ns("apply_filters"), "Apply Filters"),
    plotOutput(ns("violin_plot")),
    plotOutput(ns("scatter_plot1")),
    plotOutput(ns("scatter_plot2")),
    verbatimTextOutput(ns("qc_summary")),
    actionButton(ns("next_step"), "Continue to Doublet Removal", disabled = TRUE)
  )
}

quality_control_server <- function(id, data_for_qc) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    filtered_data <- reactiveVal(NULL)
    
    # Apply filters based on user input and update the Seurat object
    observeEvent(input$apply_filters, {
      req(data_for_qc())
      
      tryCatch({
        data <- data_for_qc()
        
        # Calculate percent.mt and percent.ribo if not present
        if (!"percent.mt" %in% colnames(data@meta.data)) {
          data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
        }
        if (!"percent.ribo" %in% colnames(data@meta.data)) {
          data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
        }
        
        # Filter cells based on QC thresholds
        filtered_data(subset(data, 
                             nFeature_RNA > input$min_nFeature_RNA & nFeature_RNA < input$max_nFeature_RNA &
                               nCount_RNA > input$min_nCount_RNA & nCount_RNA < input$max_nCount_RNA &
                               percent.mt > input$min_percent_mt & percent.mt < input$max_percent_mt &
                               percent.ribo > input$min_percent_ribo & percent.ribo < input$max_percent_ribo
        ))
        
        filtered <- filtered_data()
        
        # Generate violin plot for quality control metrics
        output$violin_plot <- renderPlot({
          VlnPlot(filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
                  ncol = 4, pt.size = 0.1) + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        })
        
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
          cat("Total cells before filtering:", ncol(data_for_qc()), "\n")
          cat("Total cells after filtering:", ncol(filtered), "\n")
          cat("Cells removed:", ncol(data_for_qc()) - ncol(filtered), "\n")
          cat("Genes detected:", nrow(filtered), "\n")
        })
        
        shinyjs::enable("next_step")
      }, error = function(e) {
        showNotification(paste("Error in quality control:", e$message), type = "error")
      })
    })
    
    # Return the filtered data for the next step
    return(reactive({ filtered_data() }))
  })
}
