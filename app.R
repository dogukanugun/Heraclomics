# app.R

# Increase file upload size to 500 MB
options(shiny.maxRequestSize = 500 * 1024^2)

# Load required libraries
library(shiny)
library(shinyjs)
library(magrittr)  # Load magrittr for pipe operator
library(Matrix)
library(Seurat)
library(DT)
library(harmony)
library(SingleR)
library(celldex)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(ggplot2)

# Load modules
source("Scripts/load_data.R")
source("Scripts/quality_control.R")
source("Scripts/doublet_removal.R")
source("Scripts/clustering.R")
source("Scripts/integration_dimensionality_reduction.R")
source("Scripts/second_dataset.R")
source("Scripts/data_correction.R")
source("Scripts/cluster_labeling.R")
source("Scripts/gene_expression.R")
source("Scripts/gene_coexpression.R")
source("Scripts/gene_regulatory_networks.R")
source("Scripts/gwas_associations.R")
source("Scripts/trajectory_analysis.R")
#source("Scripts/cell_communication.R")
source("Scripts/multimodal_analysis.R")
source("Scripts/differential_expression.R")
source("Scripts/custom_differential_expression.R")

# Define UI for application
ui <- fluidPage(
  useShinyjs(), # Include shinyjs
  titlePanel("Heracleomics: Single-cell RNA-seq Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      navbarPage("Heracleomics",
                 tabPanel("Home", 
                          h3("Welcome to Heracleomics"),
                          p("A comprehensive tool for single-cell RNA-seq analysis."),
                          p("Choose 'New' to start a new analysis or 'Old' to load a previously saved environment."),
                          actionButton("start_new", "New"),
                          actionButton("continue_old", "Continue")
                 ),
                 tabPanel("Load Data", 
                          uiOutput("load_data_ui")
                 ),
                 tabPanel("Quality Control", 
                          uiOutput("quality_control_ui")
                 ),
                 tabPanel("Integration and Dimensionality Reduction", 
                          uiOutput("integration_dimensionality_reduction_ui")
                 ),
                 tabPanel("Clustering", 
                          uiOutput("clustering_ui")
                 ),
                 tabPanel("Data Correction", 
                          uiOutput("data_correction_ui")
                 ),
                 tabPanel("Cluster Labeling", 
                          uiOutput("cluster_labeling_ui")
                 ),
                 tabPanel("Gene Expression", 
                          uiOutput("gene_expression_ui")
                 ),
                 tabPanel("Gene Co-expression", 
                          uiOutput("gene_coexpression_ui")
                 ),
                 tabPanel("Gene Regulatory Networks", 
                          uiOutput("gene_regulatory_networks_ui")
                 ),
                 tabPanel("GWAS Associations", 
                          uiOutput("gwas_associations_ui")
                 ),
                 tabPanel("Trajectory Analysis", 
                          uiOutput("trajectory_analysis_ui")
                 ),
                 tabPanel("Cell Communication", 
                          uiOutput("cell_communication_ui")
                 ),
                 tabPanel("Multimodal Analysis", 
                          uiOutput("multimodal_analysis_ui")
                 ),
                 tabPanel("Differential Expression", 
                          uiOutput("differential_expression_ui")
                 ),
                 tabPanel("Custom Differential Expression", 
                          uiOutput("custom_differential_expression_ui")
                 )
      )
    ),
    mainPanel(
      textOutput("current_step"),
      uiOutput("step_ui")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Initialize data reactive values
  data <- reactiveValues(
    for_qc = NULL,
    for_integration = NULL,
    for_clustering = NULL,
    for_correction = NULL,
    for_labeling = NULL,
    for_gene_expression = NULL,
    for_gene_coexpression = NULL,
    for_gene_regulatory_networks = NULL,
    for_gwas_associations = NULL,
    for_trajectory_analysis = NULL,
    for_cell_communication = NULL,
    for_multimodal_analysis = NULL,
    for_differential_expression = NULL,
    for_custom_differential_expression = NULL
  )
  
  # Helper function to ensure Seurat object and metadata compatibility
  ensure_seurat_metadata_compatibility <- function(seurat_obj, new_metadata) {
    cells_in_object <- colnames(seurat_obj)
    cells_in_metadata <- rownames(new_metadata)
    
    common_cells <- intersect(cells_in_object, cells_in_metadata)
    if (length(common_cells) == 0) {
      stop("No overlapping cells found between Seurat object and new metadata")
    }
    
    seurat_obj <- subset(seurat_obj, cells = common_cells)
    new_metadata <- new_metadata[common_cells, ]
    
    seurat_obj <- AddMetaData(seurat_obj, new_metadata)
    return(seurat_obj)
  }
  
  observeEvent(input$start_new, {
    output$load_data_ui <- renderUI({
      load_data_ui("load_data")
    })
    data$for_qc <- load_data_server("load_data", reactive(data$for_qc))
    output$current_step <- renderText("Step 1: Load Your Data")
  })
  
  observeEvent(input$continue_old, {
    # Implement loading of previously saved environment
    output$current_step <- renderText("Loading Previous Environment...")
  })
  
  observeEvent(data$for_qc(), {
    req(data$for_qc())
    output$quality_control_ui <- renderUI({
      quality_control_ui("quality_control", data_for_qc = reactive(data$for_qc()))
    })
    output$step_ui <- renderUI({
      quality_control_server("quality_control", data_for_qc = reactive(data$for_qc()))
    })
    output$current_step <- renderText("Step 2: Quality Control")
  })
  
  observeEvent(input$next_step_quality_control, {
    req(data$for_qc())
    data$for_integration <- data$for_qc()
    output$integration_dimensionality_reduction_ui <- renderUI({
      integration_dimensionality_reduction_ui("integration_dimensionality_reduction", data_for_integration = reactive(data$for_integration))
    })
    output$step_ui <- renderUI({
      integration_dimensionality_reduction_server("integration_dimensionality_reduction", data_for_integration = reactive(data$for_integration))
    })
    output$current_step <- renderText("Step 3: Integration and Dimensionality Reduction")
  })
  
  observeEvent(input$next_step_integration, {
    req(data$for_integration)
    data$for_clustering <- data$for_integration
    output$clustering_ui <- renderUI({
      clustering_ui("clustering", data_for_clustering = reactive(data$for_clustering))
    })
    output$step_ui <- renderUI({
      clustering_server("clustering", data_for_clustering = reactive(data$for_clustering))
    })
    output$current_step <- renderText("Step 4: Clustering")
  })
  
  observeEvent(input$next_step_clustering, {
    req(data$for_clustering)
    data$for_correction <- data$for_clustering
    output$data_correction_ui <- renderUI({
      data_correction_ui("data_correction", data_for_correction = reactive(data$for_correction))
    })
    output$step_ui <- renderUI({
      data_correction_server("data_correction", data_for_correction = reactive(data$for_correction))
    })
    output$current_step <- renderText("Step 5: Data Correction")
  })
  
  observeEvent(input$next_step_data_correction, {
    req(data$for_correction)
    data$for_labeling <- data$for_correction
    output$cluster_labeling_ui <- renderUI({
      cluster_labeling_ui("cluster_labeling", data_for_labeling = reactive(data$for_labeling))
    })
    output$step_ui <- renderUI({
      cluster_labeling_server("cluster_labeling", data_for_labeling = reactive(data$for_labeling))
    })
    output$current_step <- renderText("Step 6: Cluster Labeling")
  })
  
  observeEvent(input$next_step_labeling, {
    req(data$for_labeling)
    data$for_gene_expression <- data$for_labeling
    output$gene_expression_ui <- renderUI({
      gene_expression_ui("gene_expression", data_for_gene_expression = reactive(data$for_gene_expression))
    })
    output$step_ui <- renderUI({
      gene_expression_server("gene_expression", data_for_gene_expression = reactive(data$for_gene_expression))
    })
    output$current_step <- renderText("Step 7: Gene Expression")
  })
  
  observeEvent(input$next_step_gene_expression, {
    req(data$for_gene_expression)
    data$for_gene_coexpression <- data$for_gene_expression
    output$gene_coexpression_ui <- renderUI({
      gene_coexpression_ui("gene_coexpression", data_for_gene_coexpression = reactive(data$for_gene_coexpression))
    })
    output$step_ui <- renderUI({
      gene_coexpression_server("gene_coexpression", data_for_gene_coexpression = reactive(data$for_gene_coexpression))
    })
    output$current_step <- renderText("Step 8: Gene Co-expression")
  })
  
  observeEvent(input$next_step_gene_coexpression, {
    req(data$for_gene_coexpression)
    data$for_gene_regulatory_networks <- data$for_gene_coexpression
    output$gene_regulatory_networks_ui <- renderUI({
      gene_regulatory_networks_ui("gene_regulatory_networks", data_for_gene_regulatory_networks = reactive(data$for_gene_regulatory_networks))
    })
    output$step_ui <- renderUI({
      gene_regulatory_networks_server("gene_regulatory_networks", data_for_gene_regulatory_networks = reactive(data$for_gene_regulatory_networks))
    })
    output$current_step <- renderText("Step 9: Gene Regulatory Networks")
  })
  
  observeEvent(input$next_step_gene_regulatory_networks, {
    req(data$for_gene_regulatory_networks)
    data$for_gwas_associations <- data$for_gene_regulatory_networks
    output$gwas_associations_ui <- renderUI({
      gwas_associations_ui("gwas_associations", data_for_gwas_associations = reactive(data$for_gwas_associations))
    })
    output$step_ui <- renderUI({
      gwas_associations_server("gwas_associations", data_for_gwas_associations = reactive(data$for_gwas_associations))
    })
    output$current_step <- renderText("Step 10: GWAS Associations")
  })
  
  observeEvent(input$next_step_gwas_associations, {
    req(data$for_gwas_associations)
    data$for_trajectory_analysis <- data$for_gwas_associations
    output$trajectory_analysis_ui <- renderUI({
      trajectory_analysis_ui("trajectory_analysis", data_for_trajectory_analysis = reactive(data$for_trajectory_analysis))
    })
    output$step_ui <- renderUI({
      trajectory_analysis_server("trajectory_analysis", data_for_trajectory_analysis = reactive(data$for_trajectory_analysis))
    })
    output$current_step <- renderText("Step 11: Trajectory Analysis")
  })
  
  #observeEvent(input$next_step_trajectory_analysis, {
  #  req(data$for_trajectory_analysis)
 #   data$for_cell_communication <- data$for_trajectory_analysis
  #  output$cell_communication_ui <- renderUI({
   #   cell_communication_ui("cell_communication", data_for_cell_communication = reactive(data$for_cell_communication))
   # })
  #  output$step_ui <- renderUI({
  #   cell_communication_server("cell_communication", data_for_cell_communication = reactive(data$for_cell_communication))
    #})
    #output$current_step <- renderText("Step 12: Cell Communication")
  #})
  
  observeEvent(input$next_step_trajectory_analysis, {
    req(data$for_trajectory_analysis)
    data$for_multimodal_analysis <- data$for_trajectory_analysis
    output$multimodal_analysis_ui <- renderUI({
      multimodal_analysis_ui("multimodal_analysis", data_for_multimodal_analysis = reactive(data$for_multimodal_analysis))
    })
    output$step_ui <- renderUI({
      multimodal_analysis_server("multimodal_analysis", data_for_multimodal_analysis = reactive(data$for_multimodal_analysis))
    })
    output$current_step <- renderText("Step 13: Multimodal Analysis")
  })
  
  observeEvent(input$next_step_multimodal_analysis, {
    req(data$for_multimodal_analysis)
    data$for_differential_expression <- data$for_multimodal_analysis
    output$differential_expression_ui <- renderUI({
      differential_expression_ui("differential_expression", data_for_differential_expression = reactive(data$for_differential_expression))
    })
    output$step_ui <- renderUI({
      differential_expression_server("differential_expression", data_for_differential_expression = reactive(data$for_differential_expression))
    })
    output$current_step <- renderText("Step 14: Differential Expression")
  })
  
  observeEvent(input$next_step_differential_expression, {
    req(data$for_differential_expression)
    data$for_custom_differential_expression <- data$for_differential_expression
    output$custom_differential_expression_ui <- renderUI({
      custom_differential_expression_ui("custom_differential_expression", data_for_custom_differential_expression = reactive(data$for_custom_differential_expression))
    })
    output$step_ui <- renderUI({
      custom_differential_expression_server("custom_differential_expression", data_for_custom_differential_expression = reactive(data$for_custom_differential_expression))
    })
    output$current_step <- renderText("Step 15: Custom Differential Expression")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
