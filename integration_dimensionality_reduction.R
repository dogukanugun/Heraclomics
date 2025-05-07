library(shiny)
library(shinyjs)
library(shinyalert)
library(Seurat)
library(ggplot2)
library(plotly)
library(cowplot)
library(harmony)
library(ggrepel)
library(dplyr)
# Harmony Integration and Dimensionality Reduction Module UI Function
integrationDimReductionUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    useShinyalert(),
    h2("Harmony Integration and Dimensionality Reduction"),
    
    # Section: Normalization & Scaling
    fluidRow(
      box(
        title = "Normalization & Scaling",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        
        # Sample Names Input (for batch information)
        textInput(
          ns("sample_name_1"),
          "Sample Name 1:",
          value = "Sample1"
        ),
        textInput(
          ns("sample_name_2"),
          "Sample Name 2:",
          value = "Sample2"
        ),
        
        # Normalization Method Selection
        selectInput(
          ns("normalization_method"),
          "Select Normalization Method",
          choices = c("SCTransform" = "SCT", "LogNormalize" = "LogNormalize"),
          selected = "SCT"
        ),
        
        # Conditional Scale Factor Input
        conditionalPanel(
          condition = sprintf("input['%s'] == 'LogNormalize'", ns("normalization_method")),
          ns = ns,  
          numericInput(
            ns("scale_factor"),
            "Scale Factor:",
            value = 10000,
            min = 1,
            step = 1000
          )
        )
      )
    ),
    
    # Section: Start Integration and Dimensionality Reduction
    fluidRow(
      box(
        title = "Integration and Dimensionality Reduction",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        
        # Single Action Button
        actionButton(
          ns("run_analysis"),
          "Run Harmony Integration and Dimensionality Reduction",
          icon = icon("play"),
          class = "btn-success"
        ),
        
        br(), br(),
        
        # Progress Bar
        uiOutput(ns("progress_ui")),
        
        # Status Output
        verbatimTextOutput(ns("analysis_status"))
      )
    ),
    
    # Section: Plot Adjustments
    fluidRow(
      box(
        title = "Plot Adjustments",
        width = 12,
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
        numericInput(
          ns("plot_height"),
          "Plot Height (px):",
          min = 300,
          max = 1000,
          value = 600,
          step = 50
        )
      )
    ),
    
    # Section: Results
    h2("Results"),
    
    # General UMAP Plot Box
    fluidRow(
      box(
        title = "General UMAP Plot",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(6,
                 plotOutput(ns("general_umap_plot"), height = "auto")
          ),
          column(6,
                 # Only "batch" is available
                 selectInput(
                   ns("umap_color_by"),
                   label = "Color UMAP By:",
                   choices = "batch",
                   selected = "batch"
                 ),
                 downloadButton(ns("download_general_umap_plot"), "Download General UMAP Plot")
          )
        )
      )
    ),
    
    # General t-SNE Plot Box
    fluidRow(
      box(
        title = "General t-SNE Plot",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(6,
                 plotOutput(ns("general_tsne_plot"), height = "auto")
          ),
          column(6,
                 # Only "batch" is available
                 selectInput(
                   ns("tsne_color_by"),
                   label = "Color t-SNE By:",
                   choices = "batch",
                   selected = "batch"
                 ),
                 downloadButton(ns("download_general_tsne_plot"), "Download General t-SNE Plot")
          )
        )
      )
    ),
    
    # UMAP Gene Expression Plot Box
    fluidRow(
      box(
        title = "UMAP Gene Expression",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(6,
                 plotOutput(ns("umap_gene_plot"), height = "auto")
          ),
          column(6,
                 selectInput(
                   ns("selected_gene"),
                   label = "Select Gene for UMAP Plot",
                   choices = NULL  # Options will be updated in server
                 ),
                 downloadButton(ns("download_umap_gene_plot"), "Download UMAP Gene Plot")
          )
        )
      )
    ),
    
    # t-SNE Gene Expression Plot Box
    fluidRow(
      box(
        title = "t-SNE Gene Expression",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(6,
                 plotOutput(ns("tsne_gene_plot"), height = "auto")
          ),
          column(6,
                 selectInput(
                   ns("selected_gene_tsne"),
                   label = "Select Gene for t-SNE Plot",
                   choices = NULL  # Options will be updated in server
                 ),
                 downloadButton(ns("download_tsne_gene_plot"), "Download t-SNE Gene Plot")
          )
        )
      )
    ),
    
    # Dimensionality Reduction Heatmap Box
    fluidRow(
      box(
        title = "Dimensionality Reduction Heatmap",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(6,
                 plotOutput(ns("dimension_reduction_heatmap"), height = "auto")
          ),
          column(6,
                 selectInput(ns("dimension_select"),
                             label = "Select Dimensions for Heatmap",
                             choices = 1:30,
                             multiple = TRUE,
                             selected = 1:10),
                 numericInput(ns("heatmap_height"),
                              label = "Heatmap Height (px):",
                              value = 600,
                              min = 400,
                              step = 100),
                 numericInput(ns("heatmap_width"),
                              label = "Heatmap Width (px):",
                              value = 800,
                              min = 400,
                              step = 100),
                 downloadButton(ns("download_dimension_reduction_heatmap"), "Download Heatmap")
          )
        )
      )
    ),
    
    # Elbow Plot Box
    fluidRow(
      box(
        title = "Elbow Plot",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(6,
                 plotOutput(ns("elbow_plot"), height = "auto")
          ),
          column(6,
                 downloadButton(ns("download_elbow_plot"), "Download Elbow Plot")
          )
        )
      )
    ),
    
    # PCA Loadings Plot Box
    fluidRow(
      box(
        title = "PCA Loadings Plot",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        fluidRow(
          column(6,
                 plotOutput(ns("loadings_plot"), height = "auto")
          ),
          column(6,
                 downloadButton(ns("download_loadings_plot"), "Download PCA Loadings Plot")
          )
        )
      )
    ),
    
    # Next Step Button
    fluidRow(
      box(
        title = "Proceed to Clustering",
        width = 12,
        status = "success",
        solidHeader = TRUE,
        actionButton(
          ns("next_step"),
          "Proceed to Clustering",
          icon = icon("arrow-right"),
          class = "btn-primary"
        )
      )
    )
  )
}
# Harmony Integration and Dimensionality Reduction Module Server Function
integrationDimReductionServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    
    ### Integration and Dimensionality Reduction Process ###
    observeEvent(input$run_analysis, {
      req(rv$seurat_object)
      
      # Disable the Next Step button until analysis is complete
      shinyjs::disable("next_step")
      
      # Show a modal dialog indicating that the analysis has started
      showModal(modalDialog(
        title = "Analysis Started",
        "Harmony Integration and Dimensionality Reduction is in progress. Please wait...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Check if a second dataset is provided
      has_second <- !is.null(rv$seurat_second) && input$sample_name_2 != ""
      print("Second dataset check complete")
      
      # Assign batch information
      rv$seurat_object$batch <- input$sample_name_1
      if (has_second) {
        rv$seurat_second$batch <- input$sample_name_2
      }
      print("Batch assignment complete")
      
      output$analysis_status <- renderText({
        "Analysis started..."
      })
      
      tryCatch({
        withProgress(message = 'Running Harmony Integration and Dimensionality Reduction', value = 0, {
          
          ### Step 1: Normalize Data ###
          incProgress(0.1, detail = "Normalizing data...")
          print("Normalization process started")
          
          if (input$normalization_method == "SCT") {
            if (has_second) {
              incProgress(0.05, detail = "Merging datasets...")
              print("Merging datasets...")
              merged_data <- merge(x = rv$seurat_object, y = rv$seurat_second, add.cell.ids = c(input$sample_name_1, input$sample_name_2))
              print("Datasets merged successfully")
              incProgress(0.1, detail = "Applying SCTransform...")
              merged_data <- SCTransform(merged_data, verbose = FALSE)
              print("SCTransform applied")
            } else {
              merged_data <- SCTransform(rv$seurat_object, verbose = FALSE)
              print("SCTransform applied to single dataset")
            }
            rv$assay_used <- "SCT"  # Store the assay used
          } else {
            if (has_second) {
              incProgress(0.05, detail = "Normalizing Sample 1...")
              print("Normalizing Sample 1...")
              rv$seurat_object <- NormalizeData(rv$seurat_object, normalization.method = "LogNormalize", scale.factor = input$scale_factor, verbose = FALSE)
              rv$seurat_object <- FindVariableFeatures(rv$seurat_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
              
              incProgress(0.1, detail = "Normalizing Sample 2...")
              print("Normalizing Sample 2...")
              rv$seurat_second <- NormalizeData(rv$seurat_second, normalization.method = "LogNormalize", scale.factor = input$scale_factor, verbose = FALSE)
              rv$seurat_second <- FindVariableFeatures(rv$seurat_second, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
              
              incProgress(0.15, detail = "Merging datasets...")
              print("Merging datasets after LogNormalize...")
              merged_data <- merge(x = rv$seurat_object, y = rv$seurat_second, add.cell.ids = c(input$sample_name_1, input$sample_name_2))
              print("Datasets merged after LogNormalize")
            } else {
              rv$seurat_integrated <- NormalizeData(rv$seurat_object, normalization.method = "LogNormalize", scale.factor = input$scale_factor, verbose = FALSE)
              rv$seurat_integrated <- FindVariableFeatures(rv$seurat_integrated, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
              merged_data <- rv$seurat_integrated
              print("Single dataset normalized and variable features found")
            }
            rv$assay_used <- "RNA"  # Store the assay used
          }
          
          print(paste("Assay used:", rv$assay_used))
          
          ### Step 2: Set Default Assay ###
          incProgress(0.2, detail = "Setting default assay...")
          print("Setting default assay...")
          DefaultAssay(merged_data) <- rv$assay_used
          print(paste("Default assay set to", rv$assay_used))
          
          ### Step 3: Scale Data ###
          incProgress(0.25, detail = "Scaling data...")
          print("Scaling data...")
          merged_data <- ScaleData(object = merged_data, verbose = FALSE)
          print("Data scaling complete")
          
          ### Step 4: PCA before Harmony ###
          incProgress(0.3, detail = "Running PCA before Harmony...")
          print("Running PCA before Harmony...")
          merged_data <- RunPCA(object = merged_data, npcs = 30, verbose = FALSE)
          print("PCA complete before Harmony")
          
          ### Step 5: Harmony Integration ###
          incProgress(0.4, detail = "Running Harmony integration...")
          if (has_second) {
            print("Running Harmony...")
            if (!"batch" %in% colnames(merged_data@meta.data)) {
              stop("Error: 'batch' variable not found in metadata.")
            }
            
            # Check if the batch variable has at least two levels
            if (length(unique(merged_data$batch)) > 1) {
              merged_data <- RunHarmony(object = merged_data, group.by.vars = "batch", assay.use = rv$assay_used, plot_convergence = FALSE, verbose = FALSE)
              print("Harmony integration complete")
            } else {
              stop("Error: 'batch' must have at least two levels for Harmony integration.")
            }
          } else {
            print("Skipping Harmony integration due to single dataset")
          }
          
          ### Step 6: PCA after Harmony ###
          incProgress(0.55, detail = "Running PCA after Harmony...")
          print("Running PCA...")
          merged_data <- RunPCA(object = merged_data, npcs = 30, verbose = FALSE)
          print("PCA complete")
          
          ### Step 7: Run UMAP and t-SNE ###
          incProgress(0.65, detail = "Running UMAP and t-SNE...")
          print("Running UMAP and t-SNE...")
          merged_data <- RunUMAP(object = merged_data, reduction = "pca", dims = 1:30, verbose = FALSE)
          merged_data <- RunTSNE(object = merged_data, reduction = "pca", dims = 1:30, verbose = FALSE)
          print("UMAP and t-SNE complete")
          
          ### Step 8: Assign Integrated Data ###
          incProgress(0.75, detail = "Assigning integrated data...")
          rv$seurat_integrated <- merged_data
          print("Integrated data assigned")
          
          ### Step 9: Update Input Options ###
          updateSelectInput(
            session = session,
            inputId = "selected_gene",
            choices = rownames(rv$seurat_integrated),
            selected = NULL
          )
          
          updateSelectInput(
            session = session,
            inputId = "selected_gene_tsne",
            choices = rownames(rv$seurat_integrated),
            selected = NULL
          )
          
          # Only "batch" option is available
          updateSelectInput(
            session = session,
            inputId = "umap_color_by",
            choices = "batch",
            selected = "batch"
          )
          updateSelectInput(
            session = session,
            inputId = "tsne_color_by",
            choices = "batch",
            selected = "batch"
          )
          print("Input options updated with 'batch'")
          
          ### Step 10: Generate Plots ###
          incProgress(0.8, detail = "Generating plots...")
          print("Generating plots...")
          
          # General UMAP Plot
          output$general_umap_plot <- renderPlot({
            req(rv$seurat_integrated)
            DimPlot(
              rv$seurat_integrated,
              reduction = "umap",
              group.by = "batch",
              pt.size = input$point_size,
              label = TRUE,
              label.size = input$label_size / 3
            ) +
              ggtitle("UMAP Plot Colored by Batch") +
              theme(
                plot.title = element_text(size = input$label_size)
              )
          }, height = input$plot_height)
          print("General UMAP plot generated")
          
          # General t-SNE Plot
          output$general_tsne_plot <- renderPlot({
            req(rv$seurat_integrated)
            DimPlot(
              rv$seurat_integrated,
              reduction = "tsne",
              group.by = "batch",
              pt.size = input$point_size,
              label = TRUE,
              label.size = input$label_size / 3
            ) +
              ggtitle("t-SNE Plot Colored by Batch") +
              theme(
                plot.title = element_text(size = input$label_size)
              )
          }, height = input$plot_height)
          print("General t-SNE plot generated")
          
          # UMAP Plot for Selected Gene
          output$umap_gene_plot <- renderPlot({
            req(rv$seurat_integrated)
            req(input$selected_gene)
            
            # Check if the selected gene exists in the dataset
            if (!(input$selected_gene %in% rownames(rv$seurat_integrated))) {
              showNotification("Selected gene not found in the dataset.", type = "error")
              return(NULL)
            }
            
            FeaturePlot(
              object = rv$seurat_integrated,
              features = input$selected_gene,
              reduction = "umap",
              cols = c("lightgrey", "red"),
              label = FALSE,
              pt.size = input$point_size
            ) + 
              ggtitle(paste("UMAP Plot -", input$selected_gene, "Expression")) +
              theme(
                plot.title = element_text(size = input$label_size),
                axis.title = element_text(size = input$label_size),
                axis.text = element_text(size = input$label_size)
              )
          }, height = input$plot_height)
          print("UMAP Gene plot generated")
          
          # t-SNE Plot for Selected Gene
          output$tsne_gene_plot <- renderPlot({
            req(rv$seurat_integrated)
            req(input$selected_gene_tsne)
            
            # Check if the selected gene exists in the dataset
            if (!(input$selected_gene_tsne %in% rownames(rv$seurat_integrated))) {
              showNotification("Selected gene not found in the dataset.", type = "error")
              return(NULL)
            }
            
            FeaturePlot(
              object = rv$seurat_integrated,
              features = input$selected_gene_tsne,
              reduction = "tsne",
              cols = c("lightgrey", "red"),
              label = FALSE,
              pt.size = input$point_size
            ) +
              ggtitle(paste("t-SNE Plot -", input$selected_gene_tsne, "Expression")) +
              theme(
                plot.title = element_text(size = input$label_size),
                axis.title = element_text(size = input$label_size),
                axis.text = element_text(size = input$label_size)
              )
          }, height = input$plot_height)
          print("t-SNE Gene plot generated")
          
          # Dimensionality Reduction Heatmap
          output$dimension_reduction_heatmap <- renderPlot({
            req(rv$seurat_integrated)
            selected_dims <- as.numeric(input$dimension_select)
            
            DimHeatmap(
              object = rv$seurat_integrated,
              dims = selected_dims,
              reduction = "pca",
              cells = 500,
              balanced = TRUE
            )
          }, height = input$heatmap_height, width = input$heatmap_width)
          print("Dimensionality Reduction Heatmap generated")
          
          # Elbow Plot
          output$elbow_plot <- renderPlot({
            req(rv$seurat_integrated)
            ElbowPlot(
              object = rv$seurat_integrated,
              ndims = 30
            ) + 
              ggtitle("Elbow Plot") +
              theme(
                plot.title = element_text(size = input$label_size),
                axis.title = element_text(size = input$label_size),
                axis.text = element_text(size = input$label_size)
              )
          }, height = input$plot_height)
          print("Elbow plot generated")
          
          # PCA Loadings Plot
          output$loadings_plot <- renderPlot({
            req(rv$seurat_integrated)
            loadings <- Loadings(rv$seurat_integrated, reduction = "pca")
            loadings_df <- as.data.frame(loadings[, 1:2])
            colnames(loadings_df) <- c("PC_1", "PC_2")
            loadings_df$gene <- rownames(loadings_df)
            top_loadings <- loadings_df %>%
              arrange(desc(abs(PC_1))) %>%
              slice(1:10)
            
            ggplot(loadings_df, aes(x = PC_1, y = PC_2)) +
              geom_point(color = "steelblue", size = input$point_size) +
              geom_text_repel(data = top_loadings, aes(label = gene), size = input$label_size / 3, max.overlaps = Inf) + 
              theme_minimal() +
              labs(title = "Top PCA Loadings", x = "PC 1", y = "PC 2") +
              theme(
                plot.title = element_text(size = input$label_size),
                axis.title = element_text(size = input$label_size),
                axis.text = element_text(size = input$label_size)
              )
          }, height = input$plot_height)
          print("PCA Loadings plot generated")
          
          incProgress(0.9, detail = "Plots generated.")
          
          ### Step 11: Finalize Analysis ###
          incProgress(1, detail = "Analysis complete.")
          # UMAP Plot Download Handler
          output$download_general_umap_plot <- downloadHandler(
            filename = function() {
              paste("general_umap_plot_", Sys.Date(), ".png", sep = "")
            },
            content = function(file) {
              png(file)
              plot(DimPlot(
                rv$seurat_integrated,
                reduction = "umap",
                group.by = "batch",
                pt.size = input$point_size,
                label = TRUE,
                label.size = input$label_size / 3
              ))
              dev.off()
            }
          )
          
          # t-SNE Plot Download Handler
          output$download_general_tsne_plot <- downloadHandler(
            filename = function() {
              paste("general_tsne_plot_", Sys.Date(), ".png", sep = "")
            },
            content = function(file) {
              png(file)
              plot(DimPlot(
                rv$seurat_integrated,
                reduction = "tsne",
                group.by = "batch",
                pt.size = input$point_size,
                label = TRUE,
                label.size = input$label_size / 3
              ))
              dev.off()
            }
          )
          
          # UMAP Gene Expression Plot Download Handler
          output$download_umap_gene_plot <- downloadHandler(
            filename = function() {
              paste("umap_gene_plot_", input$selected_gene, "_", Sys.Date(), ".png", sep = "")
            },
            content = function(file) {
              png(file)
              plot(FeaturePlot(
                object = rv$seurat_integrated,
                features = input$selected_gene,
                reduction = "umap",
                cols = c("lightgrey", "red"),
                label = FALSE,
                pt.size = input$point_size
              ) + ggtitle(paste("UMAP Plot -", input$selected_gene, "Expression")))
              dev.off()
            }
          )
          
          # t-SNE Gene Expression Plot Download Handler
          output$download_tsne_gene_plot <- downloadHandler(
            filename = function() {
              paste("tsne_gene_plot_", input$selected_gene_tsne, "_", Sys.Date(), ".png", sep = "")
            },
            content = function(file) {
              png(file)
              plot(FeaturePlot(
                object = rv$seurat_integrated,
                features = input$selected_gene_tsne,
                reduction = "tsne",
                cols = c("lightgrey", "red"),
                label = FALSE,
                pt.size = input$point_size
              ) + ggtitle(paste("t-SNE Plot -", input$selected_gene_tsne, "Expression")))
              dev.off()
            }
          )
          
          # Dimensionality Reduction Heatmap Download Handler
          output$download_dimension_reduction_heatmap <- downloadHandler(
            filename = function() {
              paste("dimension_reduction_heatmap_", Sys.Date(), ".png", sep = "")
            },
            content = function(file) {
              png(file, width = input$heatmap_width, height = input$heatmap_height)
              selected_dims <- as.numeric(input$dimension_select)
              plot(DimHeatmap(
                object = rv$seurat_integrated,
                dims = selected_dims,
                reduction = "pca",
                cells = 500,
                balanced = TRUE
              ))
              dev.off()
            }
          )
          
          # Elbow Plot Download Handler
          output$download_elbow_plot <- downloadHandler(
            filename = function() {
              paste("elbow_plot_", Sys.Date(), ".png", sep = "")
            },
            content = function(file) {
              png(file)
              plot(ElbowPlot(
                object = rv$seurat_integrated,
                ndims = 30
              ) + ggtitle("Elbow Plot"))
              dev.off()
            }
          )
          
          # PCA Loadings Plot Download Handler
          output$download_loadings_plot <- downloadHandler(
            filename = function() {
              paste("pca_loadings_plot_", Sys.Date(), ".png", sep = "")
            },
            content = function(file) {
              png(file)
              loadings <- Loadings(rv$seurat_integrated, reduction = "pca")
              loadings_df <- as.data.frame(loadings[, 1:2])
              colnames(loadings_df) <- c("PC_1", "PC_2")
              loadings_df$gene <- rownames(loadings_df)
              top_loadings <- loadings_df %>%
                arrange(desc(abs(PC_1))) %>%
                slice(1:10)
              
              plot(ggplot(loadings_df, aes(x = PC_1, y = PC_2)) +
                     geom_point(color = "steelblue", size = input$point_size) +
                     geom_text_repel(data = top_loadings, aes(label = gene), size = input$label_size / 3, max.overlaps = Inf) + 
                     theme_minimal() +
                     labs(title = "Top PCA Loadings", x = "PC 1", y = "PC 2"))
              dev.off()
            }
          )
          
          # Close the modal dialog
          removeModal()
          
          # Show a new modal dialog indicating that the analysis is complete
          showModal(modalDialog(
            title = "Analysis Complete",
            "Harmony Integration and Dimensionality Reduction is complete. You can now proceed to clustering.",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
          
          output$analysis_status <- renderText({
            "Harmony Integration and Dimensionality Reduction complete. You can now proceed to clustering."
          })
          
          # Enable the Next Step button
          shinyjs::enable("next_step")
        })
        
      }, error = function(e) {
        # Close the modal dialog in case of error
        removeModal()
        
        # Show an error message in a modal dialog
        showModal(modalDialog(
          title = "Error",
          paste("An error occurred:", e$message),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        
        output$analysis_status <- renderText({
          paste("Error:", e$message)
        })
        print(paste("Error: ", e$message))
        # Ensure the Next Step button remains disabled in case of error
        shinyjs::disable("next_step")
      })
    })
    
    ### Next Step Button - Proceed to Clustering ###
    observeEvent(input$next_step, {
      # Check if the integrated data is available
      req(rv$seurat_integrated)
      
      shinyalert(
        title = "Proceed to Clustering?",
        text = "Are you sure you want to proceed to clustering?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, proceed",
        cancelButtonText = "No, stay here",
        callbackR = function(value) {
          if (isTRUE(value)) {
            # Debugging print
            print("Confirmed: Proceeding to Clustering.")
            # Send a message to proceed to clustering
            session$sendCustomMessage(type = "proceed_to_clustering", message = list())
            
            # Set reactive value to TRUE
            rv$proceed_to_clustering <- TRUE
            
            # Here you can add the logic to perform clustering
            # For example:
            # perform_clustering(rv$seurat_integrated)
          } else {
            print("Cancelled: Staying on the clustering step.")
          }
        }
      )
    })
    
    
  })
}
            