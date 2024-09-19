# gene_expression.R

library(shiny)
library(Seurat)

gene_expression_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 9: Gene Expression"),
    textInput(ns("gene_name"), "Enter Gene Name:"),
    actionButton(ns("plot_gene_expression"), "Plot Gene Expression"),
    plotOutput(ns("gene_expression_plot")),
    actionButton(ns("next_step"), "Continue to Gene Co-expression", disabled = TRUE)
  )
}

gene_expression_server <- function(id, data_for_gene_expression) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$plot_gene_expression, {
      data <- data_for_gene_expression()
      
      if (is.null(data) || input$gene_name == "") return()
      
      # Plot gene expression
      output$gene_expression_plot <- renderPlot({
        VlnPlot(data, features = input$gene_name, pt.size = 0.1) +
          ggtitle(paste("Expression of", input$gene_name))
      })
      
      shinyjs::enable("next_step")
    })
  })
}
