# Load required libraries
library(shiny)
library(shinyjs)
library(Seurat)
library(colourpicker)  # Load the colourpicker package

# UI for Quality Control step
quality_control_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Quality Control"),
    
    # Adjusters for label size, point size, plot height, and point color
    fluidRow(
      column(3, sliderInput(ns("label_size"), "Label Size", min = 8, max = 20, value = 12)),
      column(3, sliderInput(ns("point_size"), "Point Size", min = 1, max = 10, value = 3)),
      column(3, sliderInput(ns("plot_height"), "Plot Height", min = 300, max = 1000, value = 500)),
      column(3, colourInput(ns("point_color"), "Point Color", value = "black"))
    ),
    
    # Filtering thresholds
    fluidRow(
      column(4, numericInput(ns("min_features"), "Min Features", value = 200, min = 0)),
      column(4, numericInput(ns("max_features"), "Max Features", value = 2500, min = 0)),
      column(4, numericInput(ns("max_percent_mt"), "Max Percent MT", value = 5, min = 0, max = 100))
    ),
    
    # Display scatter plots
    fluidRow(
      column(6, plotOutput(ns("scatter_plot1")), downloadButton(ns("download_scatter_plot1"), "Download Scatter Plot 1")),
      column(6, plotOutput(ns("scatter_plot2")), downloadButton(ns("download_scatter_plot2"), "Download Scatter Plot 2"))
    ),
    
    # Display violin plots
    fluidRow(
      column(4, plotOutput(ns("violin_plot1")), downloadButton(ns("download_violin_plot1"), "Download Violin Plot 1")),
      column(4, plotOutput(ns("violin_plot2")), downloadButton(ns("download_violin_plot2"), "Download Violin Plot 2")),
      column(4, plotOutput(ns("violin_plot3")), downloadButton(ns("download_violin_plot3"), "Download Violin Plot 3"))
    ),
    
    # Display histogram
    fluidRow(
      column(12, plotOutput(ns("histogram")), downloadButton(ns("download_histogram"), "Download Histogram"))
    ),
    
    # Display QC summary
    fluidRow(
      column(12, verbatimTextOutput(ns("qc_summary")))
    ),
    
    # Button to proceed to the next step
    fluidRow(
      column(12, actionButton(ns("next_step"), "Next Step", disabled = TRUE))
    )
  )
}

# Server logic for Quality Control step
quality_control_server <- function(input, output, session, app_state) {
  ns <- session$ns
  
  observe({
    req(app_state$data)
    
    # Perform quality control
    data_for_qc <- app_state$data
    
    # Calculate percent.mt
    data_for_qc[["percent.mt"]] <- PercentageFeatureSet(data_for_qc, pattern = "^MT-")
    
    # Filter cells based on user-defined quality control metrics
    filtered <- subset(data_for_qc, subset = nFeature_RNA > input$min_features & nFeature_RNA < input$max_features & percent.mt < input$max_percent_mt)
    
    # Generate scatter plots for nCount_RNA vs nFeature_RNA and nCount_RNA vs percent.mt
    output$scatter_plot1 <- renderPlot({
      FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        theme(text = element_text(size = input$label_size)) +
        geom_point(size = input$point_size, color = input$point_color)
    }, height = function() { input$plot_height })
    
    output$scatter_plot2 <- renderPlot({
      FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") +
        theme(text = element_text(size = input$label_size)) +
        geom_point(size = input$point_size, color = input$point_color)
    }, height = function() { input$plot_height })
    
    # Generate violin plots for nFeature_RNA, nCount_RNA, and percent.mt
    output$violin_plot1 <- renderPlot({
      VlnPlot(filtered, features = "nFeature_RNA") +
        theme(text = element_text(size = input$label_size))
    }, height = function() { input$plot_height })
    
    output$violin_plot2 <- renderPlot({
      VlnPlot(filtered, features = "nCount_RNA") +
        theme(text = element_text(size = input$label_size))
    }, height = function() { input$plot_height })
    
    output$violin_plot3 <- renderPlot({
      VlnPlot(filtered, features = "percent.mt") +
        theme(text = element_text(size = input$label_size))
    }, height = function() { input$plot_height })
    
    # Generate histogram for percent.mt
    output$histogram <- renderPlot({
      hist(filtered[["percent.mt"]], main = "Histogram of percent.mt", xlab = "percent.mt", breaks = 50)
    }, height = function() { input$plot_height })
    
    # Display a summary of the QC process
    output$qc_summary <- renderPrint({
      cat("QC Summary:\n")
      cat("Total cells before filtering:", ncol(data_for_qc), "\n")
      cat("Total cells after filtering:", ncol(filtered), "\n")
      cat("Cells removed:", ncol(data_for_qc) - ncol(filtered), "\n")
      cat("Genes detected:", nrow(filtered), "\n")
    })
    
    shinyjs::enable("next_step")
  })
  
  # Download handlers for plots
  output$download_scatter_plot1 <- downloadHandler(
    filename = function() { "scatter_plot1.png" },
    content = function(file) {
      png(file)
      print(FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
              theme(text = element_text(size = input$label_size)) +
              geom_point(size = input$point_size, color = input$point_color))
      dev.off()
    }
  )
  
  output$download_scatter_plot2 <- downloadHandler(
    filename = function() { "scatter_plot2.png" },
    content = function(file) {
      png(file)
      print(FeatureScatter(filtered, feature1 = "nCount_RNA", feature2 = "percent.mt") +
              theme(text = element_text(size = input$label_size)) +
              geom_point(size = input$point_size, color = input$point_color))
      dev.off()
    }
  )
  
  output$download_violin_plot1 <- downloadHandler(
    filename = function() { "violin_plot1.png" },
    content = function(file) {
      png(file)
      print(VlnPlot(filtered, features = "nFeature_RNA") +
              theme(text = element_text(size = input$label_size)))
      dev.off()
    }
  )
  
  output$download_violin_plot2 <- downloadHandler(
    filename = function() { "violin_plot2.png" },
    content = function(file) {
      png(file)
      print(VlnPlot(filtered, features = "nCount_RNA") +
              theme(text = element_text(size = input$label_size)))
      dev.off()
    }
  )
  
  output$download_violin_plot3 <- downloadHandler(
    filename = function() { "violin_plot3.png" },
    content = function(file) {
      png(file)
      print(VlnPlot(filtered, features = "percent.mt") +
              theme(text = element_text(size = input$label_size)))
      dev.off()
    }
  )
  
  output$download_histogram <- downloadHandler(
    filename = function() { "histogram.png" },
    content = function(file) {
      png(file)
      hist(filtered[["percent.mt"]], main = "Histogram of percent.mt", xlab = "percent.mt", breaks = 50)
      dev.off()
    }
  )
  
  # Return the filtered data for the next step
  return(reactive({ filtered }))
}
