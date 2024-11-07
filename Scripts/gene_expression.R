library(shiny)
library(Seurat)
library(GSVA)
library(msigdbr)
library(ggplot2)
library(plotly)
library(dplyr)
library(RColorBrewer)
library(tidyr)

# UI Function for Gene Expression Module
geneExpressionUI <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Gene Expression and Pathway Analysis"),
    
    # Gene and Pathway Selection Section
    fluidRow(
      column(
        width = 12,
        box(
          title = "Gene and Pathway Selection",
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          
          # Select Input for Plot Type (2D or 3D UMAP/t-SNE)
          selectInput(
            ns("plot_type"),
            label = "Select Plot Type:",
            choices = c("UMAP", "3D UMAP", "t-SNE", "3D t-SNE"),
            selected = "UMAP"
          ),
          
          # Gene Selection Input with selectizeInput for placeholder support
          selectizeInput(
            ns("selected_gene"),
            label = "Select Gene(s):",
            choices = NULL,  # Will be populated by server
            multiple = TRUE,
            options = list(placeholder = "Choose genes to visualize")
          ),
          
          # Pathway Selection Input with Correct Categories
          selectInput(
            ns("selected_pathway"),
            label = "Select Gene Pathway (Optional):",
            choices = c(
              "None",
              "KEGG" = "CP:KEGG",
              "Reactome" = "CP:REACTOME",
              "WikiPathways" = "CP:WIKIPATHWAYS",
              "GO Biological Process" = "GO:BP",
              "GO Molecular Function" = "GO:MF"
            ),
            selected = "None"
          ),
          
          # Action Button to Run Expression/Pathway Analysis
          actionButton(
            ns("run_gene_expression"),
            "Run Analysis",
            class = "btn-success",
            style = "margin-top: 10px"
          )
        )
      )
    ),
    
    # Plot Results Section
    fluidRow(
      column(
        width = 12,
        box(
          title = "Visualization of Gene Expression or Pathway Scores",
          width = 12,
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          
          # Filter next to the plots for Seurat Clusters or Cell Cycle Phase
          radioButtons(
            ns("group_by"),
            label = "Group by:",
            choices = c("Seurat Clusters" = "seurat_clusters", "Cell Cycle Phases" = "cell_cycle_phase"),
            selected = "seurat_clusters",
            inline = TRUE
          ),
          
          # Add the plot output for UMAP/t-SNE (display only the selected one)
          uiOutput(ns("plot_ui"))
        )
      )
    ),
    
    # Violin Plot and Dot Plot Tabs
    tabsetPanel(
      tabPanel("Violin Plot", plotlyOutput(ns("violin_plot"))),
      tabPanel("Dot Plot", plotlyOutput(ns("dot_plot")))
    ),
    
    # Download Data Section
    fluidRow(
      column(
        width = 12,
        box(
          title = "Download Gene Expression or Pathway Data",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          downloadButton(ns("download_gene_expression_data"), "Download Data")
        )
      )
    )
  )
}

# Server Function for Gene Expression Module
geneExpressionServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Populate gene choices based on available features in corrected data
    observe({
      req(rv$corrected_data)
      updateSelectizeInput(
        session,
        "selected_gene",
        choices = rownames(rv$corrected_data),
        server = TRUE
      )
    })
    
    # Observe the button to run the gene expression analysis
    observeEvent(input$run_gene_expression, {
      req(rv$corrected_data)
      
      # Confirm that the user wants to run the analysis
      showModal(modalDialog(
        title = "Confirm Gene Expression or Pathway Analysis",
        "Do you want to proceed with the selected gene expression or pathway analysis?",
        easyClose = FALSE,
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("confirm_gene_expression"), "Yes, Proceed", class = "btn-primary")
        )
      ))
    })
    
    # Handle the confirmation of gene expression or pathway analysis
    observeEvent(input$confirm_gene_expression, {
      removeModal()
      
      # Show progress modal
      showModal(modalDialog(
        title = "Analysis in Progress",
        "Please wait while the analysis is running...",
        footer = NULL
      ))
      
      # Try-Catch block to handle errors
      tryCatch({
        selected_genes <- input$selected_gene
        selected_pathway <- input$selected_pathway
        
        # Determine if the user is plotting genes or pathway scores
        is_pathway <- selected_pathway != "None"
        plot_data <- NULL
        
        if (is_pathway) {
          # Split the selected_pathway into category and subcategory
          pathway_split <- strsplit(selected_pathway, ":")[[1]]
          pathway_category <- pathway_split[1]
          pathway_subcategory <- paste(pathway_split[1], pathway_split[2], sep = ":")
          
          # Get the pathway genes using msigdbr
          msigdb_data <- msigdbr(species = "Homo sapiens", category = pathway_category, subcategory = pathway_subcategory)
          pathway_genes <- msigdb_data$gene_symbol %>% unique()
          
          # Filter to only include pathway genes that are in the dataset
          pathway_genes_in_data <- pathway_genes[pathway_genes %in% rownames(rv$corrected_data)]
          
          # Check if there are enough genes in the pathway
          if (length(pathway_genes_in_data) < 5) {  # Minimum of 5 genes for a reliable GSVA score
            showModal(modalDialog(
              title = "Insufficient Pathway Genes",
              paste("Not enough genes from the selected pathway are present in the dataset."),
              easyClose = TRUE
            ))
            return()
          }
          
          # Filter pathway genes in the dataset and run GSVA
          expr_matrix <- as.matrix(rv$corrected_data@assays$RNA@data[pathway_genes_in_data, ])
          gsva_results <- gsva(expr_matrix, list(pathway = pathway_genes_in_data), method = "gsva", verbose = FALSE)
          
          # Combine GSVA results with meta data
          plot_data <- rv$corrected_data@meta.data
          plot_data$pathway_score <- gsva_results[1, ]
          rv$gene_expr_results <- plot_data
        } else {
          # Ensure at least one gene is selected
          if (length(selected_genes) == 0) {
            showModal(modalDialog(
              title = "Error",
              "Please select at least one gene or a pathway.",
              easyClose = TRUE
            ))
            return()
          }
          
          # Fetch data for selected genes
          valid_genes <- selected_genes[selected_genes %in% rownames(rv$corrected_data)]
          if (length(valid_genes) == 0) {
            showModal(modalDialog(
              title = "Error",
              "None of the selected genes are available in the dataset.",
              easyClose = TRUE
            ))
            return()
          }
          plot_data <- FetchData(rv$corrected_data, vars = valid_genes)
          rv$gene_expr_results <- plot_data
        }
        
        # Determine grouping column
        label_by <- input$group_by  # This should match the actual meta data column name
        
        # UMAP and t-SNE Plotting Logic
        if (input$plot_type == "UMAP") {
          # 2D UMAP Plot
          p <- DimPlot(
            rv$corrected_data,
            reduction = "umap",
            group.by = label_by,
            label = TRUE
          ) +
            ggtitle("UMAP Plot")
          output$umap_plot <- renderPlotly(ggplotly(p))
          
        } else if (input$plot_type == "3D UMAP") {
          # 3D UMAP Plot
          if (!"umap" %in% Reductions(rv$corrected_data)) {
            # Compute 3D UMAP if not present
            rv$corrected_data <- RunUMAP(rv$corrected_data, reduction = "pca", dims = 1:30, n.components = 3)
          }
          umap_coords <- Embeddings(rv$corrected_data, "umap")
          plot_data_plotly <- data.frame(
            x = umap_coords[, 1],
            y = umap_coords[, 2],
            z = umap_coords[, 3],
            group = rv$corrected_data@meta.data[[label_by]]
          )
          output$umap_plot <- renderPlotly({
            plot_ly(
              data = plot_data_plotly,
              x = ~x, y = ~y, z = ~z,
              color = ~group,
              colors = RColorBrewer::brewer.pal(n = 11, name = "Set1"),
              type = "scatter3d",
              mode = "markers",
              marker = list(size = 2)
            ) %>%
              layout(title = "3D UMAP Plot")
          })
          
        } else if (input$plot_type == "t-SNE") {
          # 2D t-SNE Plot
          p <- DimPlot(
            rv$corrected_data,
            reduction = "tsne",
            group.by = label_by,
            label = TRUE
          ) +
            ggtitle("t-SNE Plot")
          output$tsne_plot <- renderPlotly(ggplotly(p))
          
        } else if (input$plot_type == "3D t-SNE") {
          # 3D t-SNE Plot
          if (!"tsne" %in% Reductions(rv$corrected_data)) {
            # Compute 3D t-SNE if not present
            rv$corrected_data <- RunTSNE(rv$corrected_data, dims = 1:30, n.components = 3)
          }
          tsne_coords <- Embeddings(rv$corrected_data, "tsne")
          plot_data_plotly <- data.frame(
            x = tsne_coords[, 1],
            y = tsne_coords[, 2],
            z = tsne_coords[, 3],
            group = rv$corrected_data@meta.data[[label_by]]
          )
          output$tsne_plot <- renderPlotly({
            plot_ly(
              data = plot_data_plotly,
              x = ~x, y = ~y, z = ~z,
              color = ~group,
              colors = RColorBrewer::brewer.pal(n = 11, name = "Set1"),
              type = "scatter3d",
              mode = "markers",
              marker = list(size = 2)
            ) %>%
              layout(title = "3D t-SNE Plot")
          })
        }
        
        # Render the appropriate plot UI
        output$plot_ui <- renderUI({
          if (input$plot_type %in% c("UMAP", "3D UMAP")) {
            plotlyOutput(ns("umap_plot"))
          } else if (input$plot_type %in% c("t-SNE", "3D t-SNE")) {
            plotlyOutput(ns("tsne_plot"))
          }
        })
        
        # Violin Plot
        output$violin_plot <- renderPlotly({
          req(rv$gene_expr_results)
          
          if (is_pathway) {
            p <- ggplot(rv$gene_expr_results, aes(x = .data[[input$group_by]], y = pathway_score, fill = .data[[input$group_by]])) +
              geom_violin() +
              geom_jitter(width = 0.2, size = 1.5) +
              theme_minimal() +
              ggtitle("Violin Plot for Pathway Scores by Group")
          } else {
            plot_long <- rv$gene_expr_results %>%
              pivot_longer(cols = everything(), names_to = "gene", values_to = "expression")
            
            p <- ggplot(plot_long, aes(x = .data[[input$group_by]], y = expression, fill = .data[[input$group_by]])) +
              geom_violin() +
              geom_jitter(width = 0.2, size = 1.5) +
              theme_minimal() +
              ggtitle("Violin Plot for Selected Genes by Group")
          }
          
          ggplotly(p)
        })
        
        # Dot Plot
        output$dot_plot <- renderPlotly({
          req(rv$gene_expr_results)
          
          if (is_pathway) {
            # Dot plot i??in yolak skorlar?? uygun olmayabilir. Alternatif bir g??rselle??tirme d??????nebilirsiniz.
            showNotification("Dot Plot is not applicable for pathway scores.", type = "warning")
            return(NULL)
          } else {
            p <- DotPlot(rv$corrected_data, features = input$selected_gene, group.by = input$group_by) +
              ggtitle("Dot Plot for Selected Genes")
            ggplotly(p)
          }
        })
        
        # Hide progress modal
        removeModal()
        
      }, error = function(e) {
        # Handle errors by showing an error modal
        removeModal()
        showModal(modalDialog(
          title = "Error",
          paste("An error occurred during the analysis:", e$message),
          easyClose = TRUE,
          footer = modalButton("Dismiss")
        ))
      })
    })
    
    # Download Gene Expression Data
    output$download_gene_expression_data <- downloadHandler(
      filename = function() {
        paste("gene_expression_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(rv$gene_expr_results)
        write.csv(rv$gene_expr_results, file)
      }
    )
  })
}
