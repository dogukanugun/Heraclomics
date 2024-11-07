options(future.globals.maxSize = 64 * 1024^3)  # Increase to 64 GiB
options(shiny.maxRequestSize = 5000 * 1024^2)

# Load necessary libraries
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyalert)
library(Seurat)
library(DT)
library(ggplot2)
library(plotly)
library(cowplot)
library(DoubletFinder)
library(shinyWidgets)
library(tensorflow)
library(keras)
library(harmony)
library(ggrepel)
library(MEGENA)

# Source all the module scripts
source("Scripts/introduction.R")
source("Scripts/load_data.R")
source("Scripts/quality_control.R")
source("Scripts/doublet_removal.R")
source("Scripts/load_second_data.R")
source("Scripts/second_quality_control.R")
source("Scripts/second_doublet_removal.R")
source("Scripts/integration_dimensionality_reduction.R")
source("Scripts/clustering.R")
source("Scripts/data_correction.R")
source("Scripts/labelling_clusters.R")
source("Scripts/gene_expression.R")
source("Scripts/gene_coexpression.R")
source("Scripts/gene_regulatory_networks.R")
source("Scripts/trajectory_analysis.R")
source("Scripts/cell_communication.R")

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Heraclomics"),
  dashboardSidebar(
    sidebarMenu(id = "sidebarMenu",
                menuItem("Introduction", tabName = "introduction", icon = icon("home")),
                menuItem("Load Data", tabName = "load_data", icon = icon("upload")),
                menuItem("Quality Control", tabName = "qc", icon = icon("filter")),
                menuItem("Doublet Removal", tabName = "doublet_removal", icon = icon("trash")),
                menuItem("Load Second Dataset", tabName = "load_second_data", icon = icon("upload")),
                menuItem("Second Dataset QC", tabName = "second_qc", icon = icon("filter")),
                menuItem("Second Dataset Doublet Removal", tabName = "doublet_removal_second", icon = icon("trash")),
                menuItem("Integration and Dimensionality Reduction", tabName = "dim_reduction", icon = icon("project-diagram")),
                menuItem("Clustering", tabName = "clustering", icon = icon("sitemap")),
                menuItem("Data Correction", tabName = "data_correction", icon = icon("magic")),
                menuItem("Further Analysis", tabName = "further_analysis", icon = icon("tasks"), startExpanded = TRUE,
                         menuSubItem("Labelling Clusters", tabName = "labelling_clusters", icon = icon("tags")),
                         menuSubItem("Gene Expression", tabName = "gene_expression", icon = icon("chart-bar")),
                         menuSubItem("Gene Co-Expression Networks", tabName = "gene_coexpression"),
                         menuSubItem("Gene Regulatory Networks", tabName = "gene_regulatory_network"),
                         menuSubItem("Trajectory Analysis", tabName = "trajectory_analysis"),
                         menuSubItem("Cell Communication", tabName = "cell_comm")
                         
                         
                )
    )
  ),
  dashboardBody(
    useShinyjs(),
    useShinyalert(),
    tabItems(
      tabItem(tabName = "introduction", introductionUI("introduction_ui")),
      tabItem(tabName = "load_data", loadDataUI("load_data")),
      tabItem(tabName = "qc", qualityControlUI("quality_control")),
      tabItem(tabName = "doublet_removal", doubletRemovalUI("doublet_removal")),
      tabItem(tabName = "load_second_data", loadSecondDataUI("load_second_data")),
      tabItem(tabName = "second_qc", secondQualityControlUI("second_quality_control")),
      tabItem(tabName = "doublet_removal_second", secondDoubletRemovalUI("second_doublet_removal")),
      tabItem(tabName = "dim_reduction", integrationDimReductionUI("dim_reduction_ui")),
      tabItem(tabName = "clustering", clusteringUI("clustering_ui")),
      tabItem(tabName = "data_correction", dataCorrectionUI("data_correction_ui")),
      # Further Analysis Tab
      tabItem(tabName = "labelling_clusters", labellingClustersUI("labelling_clusters_ui")),
      tabItem(tabName = "gene_expression", geneExpressionUI("gene_expression")), 
      tabItem(tabName = "gene_coexpression", geneCoexpressionUI("gene_coexpression")),
      tabItem(tabName = "gene_regulatory_network", geneRegulatoryNetworkUI("gene_regulatory_network_ui")),
      tabItem(tabName = "trajectory_analysis", trajectoryAnalysisUI("trajectory_analysis_ui")),
      tabItem(tabName = "cell_comm", cellCommunicationUI("cell_comm_ui"))
     
      
    )
  )
)

# Define Server
server <- function(input, output, session) {
  # Initialize reactive values
  rv <- reactiveValues(
    data_loaded = FALSE,
    qc_done = FALSE,
    doublet_removed = FALSE,
    second_data_loaded = FALSE,
    second_qc_done = FALSE,
    second_doublet_removed = FALSE,
    integration_done = FALSE,
    seurat_object = NULL,
    seurat_second = NULL,
    seurat_integrated = NULL,
    proceed_to_clustering = FALSE,
    proceed_to_further_analysis = FALSE,  # New reactive value for further analysis,
    cluster_labels = NULL,
    trajectory_data = NULL
  )
  
  # Hide all tabs except Introduction at startup
  observe({
    hideTab(inputId = "sidebarMenu", target = "load_data")
    hideTab(inputId = "sidebarMenu", target = "qc")
    hideTab(inputId = "sidebarMenu", target = "doublet_removal")
    hideTab(inputId = "sidebarMenu", target = "load_second_data")
    hideTab(inputId = "sidebarMenu", target = "second_qc")
    hideTab(inputId = "sidebarMenu", target = "doublet_removal_second")
    hideTab(inputId = "sidebarMenu", target = "dim_reduction")
    hideTab(inputId = "sidebarMenu", target = "clustering")
    hideTab(inputId = "sidebarMenu", target = "data_correction")
    hideTab(inputId = "sidebarMenu", target = "further_analysis")
    updateTabItems(session, "sidebarMenu", "introduction")
  })
  
  # Call the Introduction module and capture its return value
  intro_return <- introductionServer("introduction_ui", rv)
  
  # Observe when the user chooses to start a new analysis
  observeEvent(intro_return$start_new_analysis(), {
    if (isTRUE(intro_return$start_new_analysis())) {
      # Show Load Data tab and switch to it
      showTab(inputId = "sidebarMenu", target = "load_data")
      updateTabItems(session, "sidebarMenu", "load_data")
      # Reset the reactive value
      intro_return$start_new_analysis(FALSE)
    }
  })
  
  # Call the Load Data module and capture its return value
  load_data_return <- loadDataServer("load_data", rv)
  
  # Observe when to proceed to QC
  observeEvent(load_data_return$proceed_to_qc(), {
    if (isTRUE(load_data_return$proceed_to_qc())) {
      # Enable QC tab and switch to it
      showTab(inputId = "sidebarMenu", target = "qc")
      updateTabItems(session, "sidebarMenu", "qc")
      # Reset the reactive value
      load_data_return$proceed_to_qc(FALSE)
    }
  })
  
  # Call the Quality Control module and capture its return value
  qc_return <- qualityControlServer("quality_control", rv)
  
  # Observe when to proceed to Doublet Removal
  observeEvent(qc_return$proceed_to_doublet(), {
    if (isTRUE(qc_return$proceed_to_doublet())) {
      # Enable Doublet Removal tab and switch to it
      showTab(inputId = "sidebarMenu", target = "doublet_removal")
      updateTabItems(session, "sidebarMenu", "doublet_removal")
      # Reset the reactive value
      qc_return$proceed_to_doublet(FALSE)
    }
  })
  
  # Call the Doublet Removal module and capture its return value
  doublet_return <- doubletRemovalServer("doublet_removal", rv)
  
  # Observe when to proceed to Load Second Dataset
  observeEvent(doublet_return$proceed_to_load_second_data(), {
    if (isTRUE(doublet_return$proceed_to_load_second_data())) {
      # Enable Load Second Dataset tab and switch to it
      showTab(inputId = "sidebarMenu", target = "load_second_data")
      updateTabItems(session, "sidebarMenu", "load_second_data")
      # Reset the reactive value
      doublet_return$proceed_to_load_second_data(FALSE)
    }
  })
  
  # Call the Load Second Data module and capture its return value
  load_second_return <- loadSecondDataServer("load_second_data", rv)
  
  # Observe when to proceed to Second QC or Integration
  observeEvent(load_second_return$proceed_to_integration(), {
    if (isTRUE(load_second_return$proceed_to_integration())) {
      if (!is.null(rv$seurat_second)) {
        # Second dataset is loaded, proceed to Second QC
        showTab(inputId = "sidebarMenu", target = "second_qc")
        updateTabItems(session, "sidebarMenu", "second_qc")
        
        # Call the Second Quality Control module and capture its return value
        second_qc_return <- secondQualityControlServer("second_quality_control", rv)
        
        # Observe when to proceed to Second Doublet Removal
        observeEvent(second_qc_return$proceed_to_second_doublet(), {
          if (isTRUE(second_qc_return$proceed_to_second_doublet())) {
            # Enable Second Doublet Removal tab and switch to it
            showTab(inputId = "sidebarMenu", target = "doublet_removal_second")
            updateTabItems(session, "sidebarMenu", "doublet_removal_second")
            # Reset the reactive value
            second_qc_return$proceed_to_second_doublet(FALSE)
            
            # Call the Second Doublet Removal module
            second_doublet_return <- secondDoubletRemovalServer("second_doublet_removal", rv)
            
            # Observe when to proceed to Integration
            observeEvent(second_doublet_return$proceed_to_load_integration(), {
              if (isTRUE(second_doublet_return$proceed_to_load_integration())) {
                # Enable Integration tab and switch to it
                showTab(inputId = "sidebarMenu", target = "dim_reduction")
                updateTabItems(session, "sidebarMenu", "dim_reduction")
                # Reset the reactive value
                second_doublet_return$proceed_to_load_integration(FALSE)
              }
            })
          }
        })
      } else {
        # Second dataset is not loaded, proceed directly to Integration
        showTab(inputId = "sidebarMenu", target = "dim_reduction")
        updateTabItems(session, "sidebarMenu", "dim_reduction")
      }
      # Reset the reactive value
      load_second_return$proceed_to_integration(FALSE)
    }
  })
  
  # Call the Integration and Dimensionality Reduction module
  integrationDimReductionServer("dim_reduction_ui", rv)
  
  # Observe when to proceed to Clustering
  observeEvent(rv$proceed_to_clustering, {
    if (isTRUE(rv$proceed_to_clustering)) {
      # Enable Clustering tab and switch to it
      showTab(inputId = "sidebarMenu", target = "clustering")
      updateTabItems(session, "sidebarMenu", "clustering")
      # Reset the reactive value
      rv$proceed_to_clustering <- FALSE
    }
  })
  
  # Call the Clustering module
  clustering_return <- clusteringServer("clustering_ui", rv)
  
  # Observe when to proceed to Data Correction
  observeEvent(clustering_return$proceed_to_next_step(), {
    if (isTRUE(clustering_return$proceed_to_next_step())) {
      # Enable Data Correction tab and switch to it
      showTab(inputId = "sidebarMenu", target = "data_correction")
      updateTabItems(session, "sidebarMenu", "data_correction")
      # Reset the reactive value
      rv$proceed_to_next_step <- FALSE  # Ensure this is reset in rv
    }
  })
  
  # Call the Data Correction module and capture its return value
  
  data_correction_return <- dataCorrectionServer("data_correction_ui", rv)
  
  # Observe when to proceed to Further Analysis after Data Correction
  observeEvent(data_correction_return$proceed_to_next_step, {
    if (isTRUE(data_correction_return$proceed_to_next_step)) {
      # Enable Further Analysis tab and switch to it
      showTab(inputId = "sidebarMenu", target = "further_analysis")
      updateTabItems(session, "sidebarMenu", "labelling_clusters")  # Show Labelling Clusters first
      # Reset the reactive value
      data_correction_return$proceed_to_next_step <- FALSE
    }
  })
  
  # Call the Labelling Clusters module (currently only this module is enabled under Further Analysis)
  labellingClustersServer("labelling_clusters_ui", rv)
  # Call the Gene Expression module
  geneExpressionServer("gene_expression", rv)
  
  geneCoexpressionServer("gene_coexpression", rv)
  
  
  gene_regulatory_network_return <- geneRegulatoryNetworkServer("gene_regulatory_network_ui", rv)
  
  trajectoryAnalysisServer("trajectory_analysis_ui", rv)
  
  cellCommunicationServer("cell_comm_ui", rv)
  
  
}

# Run the application
shinyApp(ui = ui, server = server)   
