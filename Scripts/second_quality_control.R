# Scripts/second_quality_control.R

# Load necessary libraries
library(shiny)
library(shinyjs)
library(shinyalert)
library(Seurat)
library(ggplot2)
library(plotly)
library(shinyWidgets)
library(cowplot)
library(htmlwidgets)  # For saving interactive plots as HTML

# Second Quality Control Module UI Function
secondQualityControlUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyjs(),  # Enable shinyjs
    h2("Second Dataset - Quality Control"),
    fluidRow(
      box(
        title = "QC Parameters",
        width = 4,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        # Input controls for QC thresholds
        numericInput(
          ns("min_features"),
          "Minimum Features (Genes) per Cell:",
          value = 200,
          min = 0
        ),
        numericInput(
          ns("max_features"),
          "Maximum Features (Genes) per Cell:",
          value = 3500,
          min = 0
        ),
        numericInput(
          ns("min_counts"),
          "Minimum Counts per Cell:",
          value = 500,
          min = 0
        ),
        numericInput(
          ns("max_counts"),
          "Maximum Counts per Cell:",
          value = 10000,
          min = 0,
          max = 100000
        ),
        numericInput(
          ns("max_mt"),
          "Maximum Percentage of Mitochondrial Genes:",
          value = 10,
          min = 0,
          max = 100
        ),
        numericInput(
          ns("max_ribo"),
          "Maximum Percentage of Ribosomal Genes:",
          value = 50,
          min = 0,
          max = 100
        ),
        actionButton(ns("apply_filters"), "Apply", icon = icon("filter"))
      ),
      box(
        title = "Plot Adjustments",
        width = 4,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        sliderInput(
          ns("label_size"),
          "Label Size:",
          min = 6,
          max = 20,
          value = 10
        ),
        sliderInput(
          ns("point_size"),
          "Point Size:",
          min = 0.1,
          max = 5,
          value = 1,
          step = 0.1
        ),
        sliderInput(
          ns("plot_height"),
          "Plot Height (px):",
          min = 100,
          max = 1000,
          value = 300,
          step = 50
        )
      ),
      box(
        title = "Select Violin Plots to Display",
        width = 4,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        checkboxGroupInput(
          ns("selected_violin"),
          "Choose Violin Plots:",
          choices = list(
            "nFeature_RNA" = "nFeature_RNA",
            "nCount_RNA" = "nCount_RNA",
            "percent.mt" = "percent.mt",
            "percent.ribo" = "percent.ribo"
          ),
          selected = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
        )
      )
    ),
    fluidRow(
      box(
        title = "QC Plots",
        width = 12,
        height ="auto",
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        tabsetPanel(
          tabPanel(
            "Violin Plots",
            downloadButton(ns("download_violin"), "Download Selected Violin Plots"),
            uiOutput(ns("violin_plots_ui"))
          ),
          tabPanel(
            "Scatter Plots",
            downloadButton(ns("download_scatter"), "Download Scatter Plots"),
            plotlyOutput(ns("scatter_plots"), height = "auto")
          ),
          tabPanel(
            "Histograms",
            downloadButton(ns("download_histogram"), "Download Histograms"),
            plotlyOutput(ns("histograms"), height = "auto")
          )
        )
      )
    ),
    fluidRow(
      box(
        title = "Proceed to Doublet Removal",
        width = 12,
        status = "success",
        solidHeader = TRUE,
        collapsible = TRUE,
        actionButton(ns("continue_to_doublet"), "Continue to Next Step", icon = icon("arrow-right"))
      )
    )
  )
}

# Second Quality Control Module Server Function
secondQualityControlServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    proceed_to_second_doublet <- reactiveVal(FALSE)
    
    # Reactive values to store plots
    violin_plots <- reactiveValues()
    scatter_plots <- reactiveVal(NULL)
    histogram_plot <- reactiveVal(NULL)
    
    # Generate Violin Plots UI based on selection
    output$violin_plots_ui <- renderUI({
      req(rv$seurat_second)
      selected <- input$selected_violin
      plot_output_list <- lapply(selected, function(feature) {
        plotlyOutput(ns(paste0("violin_plot_", feature)))
      })
      do.call(tagList, plot_output_list)
    })
    
    # Observe when the "Apply" button is clicked
    observeEvent(input$apply_filters, {
      req(rv$seurat_second)
      seurat_obj <- rv$seurat_second
      
      # Ensure 'percent.ribo' is calculated
      if (!"percent.ribo" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
      }
      
      # Define filters based on user input
      filters <- WhichCells(seurat_obj, expression = 
                              nFeature_RNA > input$min_features & 
                              nFeature_RNA < input$max_features & 
                              nCount_RNA > input$min_counts & 
                              nCount_RNA < input$max_counts & 
                              percent.mt < input$max_mt &
                              percent.ribo < input$max_ribo)
      
      # Subset the Seurat object
      seurat_obj_filtered <- subset(seurat_obj, cells = filters)
      
      # Update the reactive value
      rv$seurat_second <- seurat_obj_filtered
      
      # Generate Individual Violin Plots based on selection
      observe({
        req(input$selected_violin)
        lapply(input$selected_violin, function(feature) {
          output[[paste0("violin_plot_", feature)]] <- renderPlotly({
            p <- VlnPlot(
              seurat_obj_filtered,
              features = feature,
              pt.size = input$point_size
            ) + theme(text = element_text(size = input$label_size))
            ggplotly(p, height = input$plot_height)
          })
        })
      })
      
      # Generate Scatter Plots
      plot1 <- FeatureScatter(
        seurat_obj_filtered,
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA",
        pt.size = input$point_size
      ) + theme(text = element_text(size = input$label_size))
      
      plot2 <- FeatureScatter(
        seurat_obj_filtered,
        feature1 = "nCount_RNA",
        feature2 = "percent.mt",
        pt.size = input$point_size
      ) + theme(text = element_text(size = input$label_size))
      
      plot1_plotly <- ggplotly(plot1)
      plot2_plotly <- ggplotly(plot2)
      
      scatter_plots(plot1_plotly)
      scatter_plots(plot2_plotly)
      
      output$scatter_plots <- renderPlotly({
        subplot(plot1_plotly, plot2_plotly, nrows = 1, margin = 0.05) %>%
          layout(height = input$plot_height)
      })
      
      # Generate Histograms
      hist_features <- c("nFeature_RNA", "nCount_RNA")
      histogram_list <- lapply(hist_features, function(feature) {
        p <- ggplot(seurat_obj_filtered@meta.data, aes(x = .data[[feature]])) +
          geom_histogram(binwidth = 30, fill = "steelblue", color = "black") +
          theme_minimal() +
          theme(text = element_text(size = input$label_size)) +
          labs(x = feature, y = "Count")
        ggplotly(p)
      })
      
      histogram_plot(subplot(histogram_list, nrows = 2, margin = 0.05) %>%
                       layout(height = input$plot_height))
      
      output$histograms <- renderPlotly({
        histogram_plot()
      })
      
      # Show success message
      shinyalert(
        title = "Success",
        text = "Filters applied and plots generated successfully for the second dataset!",
        type = "success",
        timer = 2000,
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        showConfirmButton = FALSE
      )
      
      # Enable the "Continue to Next Step" button
      shinyjs::enable(ns("continue_to_doublet"))
    })
    
    # Download handlers for Violin Plots
    # (Implement as needed)
    
    # Observe the "Continue to Next Step" button
    observeEvent(input$continue_to_doublet, {
      shinyalert(
        title = "Proceed to Doublet Removal?",
        text = "Are you sure you want to proceed to the next step?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, proceed",
        cancelButtonText = "No, stay here",
        callbackR = function(value) {
          if (value) {
            # User confirmed, set proceed_to_second_doublet to TRUE
            proceed_to_second_doublet(TRUE)
          }
        }
      )
    })
    
    # Return the reactive value
    return(list(proceed_to_second_doublet = proceed_to_second_doublet))
  })
}
