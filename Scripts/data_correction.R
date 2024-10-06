# data_correction.R

# Load necessary libraries
library(shiny)
library(Seurat)
library(ggplot2)
library(plotly)
library(dplyr)
library(shinyWidgets)
library(msigdbr)  # For gene pathways
library(htmlwidgets)
library(shinyalert)  # Added
library(Matrix)       # Added
library(viridis)      # For scalable color palettes

# Data Correction Module UI Function
dataCorrectionUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyalert(),  # Initialize shinyalert
    h2("Data Correction"),
    fluidRow(
      box(
        title = "Correction Parameters",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        checkboxGroupInput(
          ns("correction_options"), 
          "Select Effects to Regress:",
          choices = list(
            "Confounding Effects" = "confounders",
            "Gene(s) of Interest" = "genes_of_interest",
            "Gene Pathway" = "gene_pathway"
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] && input['%s'].includes('confounders')", 
                              ns("correction_options"), ns("correction_options")),
          checkboxGroupInput(
            ns("confounders"), 
            "Select Confounding Effects to Regress:",
            choices = list(
              "Cell Cycle Phase (S.Score)" = "S.Score",
              "Cell Cycle Phase (G2M.Score)" = "G2M.Score",  # Added
              "nFeature_RNA" = "nFeature_RNA",
              "nCount_RNA" = "nCount_RNA",
              "Percent Ribo" = "percent.ribo",
              "Percent MT" = "percent.mt"
            ),
            selected = c("S.Score", "G2M.Score", "nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt")
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] && input['%s'].includes('genes_of_interest')", 
                              ns("correction_options"), ns("correction_options")),
          pickerInput(
            ns("genes"), 
            "Select Gene(s) of Interest:",
            choices = NULL,  # To be populated dynamically
            options = list(
              `actions-box` = TRUE, 
              `live-search` = TRUE  # Enables search functionality
            ),
            multiple = TRUE
          ),
          helpText("You can select up to 100 genes.")  # Informative text
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] && input['%s'].includes('gene_pathway')", 
                              ns("correction_options"), ns("correction_options")),
          selectInput(
            ns("pathway"), 
            "Select Gene Pathway:",
            choices = NULL,  # To be populated dynamically
            selectize = TRUE  # Enables search functionality
          )
        ),
        br(),
        actionButton(
          ns("run_correction"), 
          "Apply Data Correction", 
          icon = icon("play"), 
          class = "btn-success"
        )
      )
    ),
    fluidRow(
      box(
        title = "UMAP/t-SNE Plot After Correction",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        selectInput(
          ns("post_correction_plot_type"), 
          "Select Plot Type:",
          choices = c("UMAP", "3D UMAP", "t-SNE", "3D t-SNE"), 
          selected = "UMAP"
        ),
        selectInput(
          ns("post_label_by"), 
          "Label Cells By:",
          choices = list(
            "Seurat Clusters" = "seurat_clusters",
            "Sample" = "Sample",
            "Cell Cycle Phase (S.Score)" = "S.Score",
            "Cell Cycle Phase (G2M.Score)" = "G2M.Score",
            "Number of Features (nFeature_RNA)" = "nFeature_RNA",
            "Number of Counts (nCount_RNA)" = "nCount_RNA",
            "Percent Ribo" = "percent.ribo",
            "Percent MT" = "percent.mt"
          ),
          selected = "seurat_clusters"
        ),
        plotlyOutput(ns("post_correction_plot")),
        downloadButton(ns("download_post_plot"), "Download Plot")
      )
    ),
    fluidRow(
      box(
        width = 12,
        actionButton(
          ns("next_step"), 
          "Continue to Next Step", 
          icon = icon("arrow-right"), 
          class = "btn-primary"
        )
      )
    )
  )
}
# Data Correction Module Server Function
dataCorrectionServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Populate gene and pathway selections
    observe({
      req(rv$seurat_integrated)
      all_genes <- rownames(rv$seurat_integrated)
      updatePickerInput(session, "genes", choices = all_genes)
      
      # Populate pathways from MSigDB
      pathways <- msigdbr(species = "Homo sapiens", category = "C2") %>%
        dplyr::select(gs_name) %>%
        distinct() %>%
        arrange(gs_name) %>%
        pull(gs_name)
      updateSelectInput(session, "pathway", choices = pathways, selected = pathways[1])  # Default selection
    })
    
    # Observe when Apply Data Correction button is clicked
    observeEvent(input$run_correction, {
      req(rv$seurat_integrated)
      
      # Show a modal dialog to confirm correction parameters
      shinyalert::shinyalert(
        title = "Confirm Data Correction",
        text = "Are you sure you want to proceed with data correction using the selected parameters?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, Proceed",
        cancelButtonText = "Cancel",
        callback = function(x) {
          if (x) {
            # Proceed with correction
            perform_data_correction()
          }
        }
      )
    })
    
    # Function to perform data correction
    perform_data_correction <- function() {
      showModal(
        modalDialog(
          title = "Data Correction in Progress",
          "Please wait while data correction is being performed...",
          footer = NULL
        )
      )
      
      tryCatch({
        # Initialize variables to regress
        vars_to_regress <- c()
        
        # Confounding Effects
        if ("confounders" %in% input$correction_options) {
          selected_confounders <- input$confounders
          vars_to_regress <- c(vars_to_regress, selected_confounders)
        }
        
        # Gene(s) of Interest
        if ("genes_of_interest" %in% input$correction_options) {
          selected_genes <- input$genes
          if (length(selected_genes) > 0) {
            avg_gene_expression <- Matrix::colMeans(
              GetAssayData(rv$seurat_integrated, slot = "data")[selected_genes, , drop = FALSE]
            )
            rv$seurat_integrated$avg_gene_expr <- avg_gene_expression
            vars_to_regress <- c(vars_to_regress, "avg_gene_expr")
          }
        }
        
        # Gene Pathway
        if ("gene_pathway" %in% input$correction_options) {
          selected_pathway <- input$pathway
          pathway_genes <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
            filter(gs_name == selected_pathway) %>%
            pull(gene_symbol)
          pathway_genes <- intersect(pathway_genes, rownames(rv$seurat_integrated))
          if (length(pathway_genes) > 0) {
            avg_pathway_expression <- Matrix::colMeans(
              GetAssayData(rv$seurat_integrated, slot = "data")[pathway_genes, , drop = FALSE]
            )
            rv$seurat_integrated$avg_pathway_expr <- avg_pathway_expression
            vars_to_regress <- c(vars_to_regress, "avg_pathway_expr")
          }
        }
        
        # Perform Scaling with Regression
        withProgress(message = 'Scaling Data...', value = 0, {
          rv$seurat_integrated <- ScaleData(
            rv$seurat_integrated,
            vars.to.regress = vars_to_regress,
            verbose = FALSE
          )
          incProgress(1)
        })
        
        # Remove the modal dialog when correction is complete
        removeModal()
        
        # Notify the user
        shinyalert::shinyalert(
          title = "Success", 
          text = "Data correction has been applied successfully!", 
          type = "success"
        )
        
        # Generate Post-Correction Plot
        output$post_correction_plot <- renderPlotly({
          post_correction_plot_reactive()
        })
        
      }, error = function(e) {
        removeModal()  # Remove the progress modal in case of error
        shinyalert::shinyalert(
          title = "Error", 
          text = paste("An error occurred during data correction:", e$message, "Please check your inputs and try again."), 
          type = "error"
        )
      })
    }
    
    ### Reactive Expression for Post-Correction Plot ###
    post_correction_plot_reactive <- reactive({
      req(rv$seurat_integrated)
      req(input$post_correction_plot_type)
      label_by <- input$post_label_by
      
      if (input$post_correction_plot_type == "UMAP") {
        p <- DimPlot(
          rv$seurat_integrated,
          reduction = "umap",
          group.by = label_by,
          label = TRUE
        ) +
          ggtitle("UMAP Plot After Correction")
        ggplotly(p)
      } else if (input$post_correction_plot_type == "3D UMAP") {
        umap_coords <- Embeddings(rv$seurat_integrated, "umap")
        if (ncol(umap_coords) < 3) {
          shinyalert::shinyalert(
            title = "Error",
            text = "3D UMAP embeddings not available. Please compute 3D UMAP first.",
            type = "error"
          )
          return(NULL)
        }
        plot_data <- data.frame(
          x = umap_coords[, 1],
          y = umap_coords[, 2],
          z = umap_coords[, 3],
          group = rv$seurat_integrated@meta.data[[label_by]]
        )
        plot_ly(
          data = plot_data,
          x = ~x, y = ~y, z = ~z,
          color = ~group,
          colors = viridis::viridis(length(unique(plot_data$group))),
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 2)
        ) %>%
          layout(title = "3D UMAP Plot After Correction")
      } else if (input$post_correction_plot_type == "t-SNE") {
        p <- DimPlot(
          rv$seurat_integrated,
          reduction = "tsne",
          group.by = label_by,
          label = TRUE
        ) +
          ggtitle("t-SNE Plot After Correction")
        ggplotly(p)
      } else if (input$post_correction_plot_type == "3D t-SNE") {
        tsne_coords <- Embeddings(rv$seurat_integrated, "tsne")
        if (ncol(tsne_coords) < 3) {
          shinyalert::shinyalert(
            title = "Error",
            text = "3D t-SNE embeddings not available. Please compute 3D t-SNE first.",
            type = "error"
          )
          return(NULL)
        }
        plot_data <- data.frame(
          x = tsne_coords[, 1],
          y = tsne_coords[, 2],
          z = tsne_coords[, 3],
          group = rv$seurat_integrated@meta.data[[label_by]]
        )
        plot_ly(
          data = plot_data,
          x = ~x, y = ~y, z = ~z,
          color = ~group,
          colors = viridis::viridis(length(unique(plot_data$group))),
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 2)
        ) %>%
          layout(title = "3D t-SNE Plot After Correction")
      }
    })
    
    ### Download Handlers ###
    # Download Post-Correction Plot
    output$download_post_plot <- downloadHandler(
      filename = function() {
        paste0("post_correction_plot_", Sys.Date(), ".html")
      },
      content = function(file) {
        plot <- post_correction_plot_reactive()
        req(plot)
        htmlwidgets::saveWidget(plot, file)
      }
    )
    
    ### Next Step Button ###
    observeEvent(input$next_step, {
      # Signal to proceed to the next step
      rv$proceed_to_next_step <- TRUE
    })
    
    ### Selection Limit Observer ###
    observeEvent(input$genes, {
      max_genes <- 100
      selected_genes <- input$genes
      
      if (length(selected_genes) > max_genes) {
        # Remove the extra selections
        updatePickerInput(
          session, 
          "genes", 
          selected = selected_genes[1:max_genes]
        )
        
        # Notify the user
        shinyalert::shinyalert(
          title = "Selection Limit Exceeded",
          text = paste("You can select up to", max_genes, "genes. Extra selections have been removed."),
          type = "warning"
        )
      }
    })
    
  })
}
