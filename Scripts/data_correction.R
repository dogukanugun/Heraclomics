# Loading Necessary Libraries
library(shiny)
library(Seurat)
library(ggplot2)
library(plotly)
library(dplyr)
library(shinyWidgets)
library(msigdbr)
library(htmlwidgets)
library(shinyalert)
library(Matrix)
library(viridis)
library(shinyBS)

# UI Function for Data Correction Module
dataCorrectionUI <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyalert(),
    h2("Data Correction"),
    fluidRow(
      box(
        title = "Data Correction Overview",
        width = 12,
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        p("Data correction allows for the removal of potential technical and biological confounders that could mask biological signals, such as cell cycle heterogeneity."),
        p("For the removal of cell cycle effects, a cell cycle score is first assigned to each cell based on its expression of G2/M and S phase markers."),
        p("If the cell shows low or no expression of these markers, then it is likely in the G1 phase. Each cell cycle score is modeled against highly variable genes, and a corrected expression matrix is generated for dimensionality reduction and clustering."),
        p("Alternatively, you can choose to remove the effects of user-selected genes or an entire gene set (MsigDB gene sets).")
      )
    ),
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
        bsTooltip(
          id = ns("correction_options"),
          title = "Select the types of effects you want to regress out from the data. Confounding effects may include cell cycle scores and other technical variables. You can also choose to remove the influence of specific genes or entire gene sets.",
          placement = "right",
          trigger = "hover"
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] && input['%s'].indexOf('confounders') !== -1", ns("correction_options"), ns("correction_options")),
          checkboxGroupInput(
            ns("confounders"),
            "Select Confounding Effects to Regress:",
            choices = list(
              "Cell Cycle Phase (S.Score)" = "S.Score",
              "Cell Cycle Phase (G2M.Score)" = "G2M.Score",
              "nFeature_RNA" = "nFeature_RNA",
              "nCount_RNA" = "nCount_RNA",
              "Percent Ribo" = "percent.ribo",
              "Percent MT" = "percent.mt"
            ),
            selected = c("S.Score", "G2M.Score", "nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt")
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] && input['%s'].indexOf('genes_of_interest') !== -1", ns("correction_options"), ns("correction_options")),
          pickerInput(
            ns("genes"), 
            "Select Gene(s) of Interest:",
            choices = NULL,
            options = list(`actions-box` = TRUE, `live-search` = TRUE),
            multiple = TRUE
          ),
          helpText("You can select up to 100 genes.")
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] && input['%s'].indexOf('gene_pathway') !== -1", ns("correction_options"), ns("correction_options")),
          selectInput(
            ns("pathway"), 
            "Select Gene Pathway:",
            choices = NULL,
            selectize = TRUE
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
        title = "Violin Plot",
        width = 6,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotOutput(ns("violin_plot"))
      ),
      box(
        title = "Ridge Plot",
        width = 6,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        plotOutput(ns("ridge_plot"))
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

## Server Function for Data Correction Module
dataCorrectionServer <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # H??cre d??ng??s?? skorlar??n?? hesaplama fonksiyonu
    calculate_cell_cycle_scores <- function(seurat_integrated_obj) {
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
      
      seurat_integrated_obj <- CellCycleScoring(
        seurat_integrated_obj,
        s.features = s.genes,
        g2m.features = g2m.genes,
        set.ident = TRUE
      )
      return(seurat_integrated_obj)
    }
    
    
    # Ortama ekspresyon hesaplama fonksiyonu
    calculate_avg_expression <- function(seurat_integrated_obj, genes, new_col_name) {
      avg_expr <- Matrix::colMeans(GetAssayData(seurat_integrated_obj, slot = "data")[genes, , drop = FALSE])
      seurat_integrated_obj[[new_col_name]] <- avg_expr
      return(seurat_integrated_obj)
    }
    
    # Veri d??zeltme i??lemi
    perform_data_correction <- function() {
      showModal(modalDialog(
        title = "Data Correction in Progress",
        "Please wait while data correction is being performed...",
        footer = NULL
      ))
      
      tryCatch({
        vars_to_regress <- c()
        
        # H??cre d??ng??s?? skorlar?? hesaplanmad??ysa hesapla
        if (!"S.Score" %in% colnames(rv$seurat_integrated@meta.data) || 
            !"G2M.Score" %in% colnames(rv$seurat_integrated@meta.data)) {
          rv$seurat_integrated <- calculate_cell_cycle_scores(rv$seurat_integrated)
        }
        
        if ("confounders" %in% input$correction_options) {
          selected_confounders <- input$confounders
          missing_confounders <- setdiff(selected_confounders, colnames(rv$seurat_integrated@meta.data))
          if (length(missing_confounders) > 0) {
            removeModal()
            shinyalert::shinyalert(
              title = "Error",
              text = paste("The following confounders are not found in the dataset:", paste(missing_confounders, collapse = ", ")),
              type = "error"
            )
            return(NULL)
          }
          vars_to_regress <- c(vars_to_regress, selected_confounders)
        }
        
        if ("genes_of_interest" %in% input$correction_options) {
          selected_genes <- input$genes
          if (length(selected_genes) > 0) {
            missing_genes <- setdiff(selected_genes, rownames(rv$seurat_integrated))
            if (length(missing_genes) > 0) {
              removeModal()
              shinyalert::shinyalert(
                title = "Error",
                text = paste("The following selected genes are not found in the dataset:", paste(missing_genes, collapse = ", ")),
                type = "error"
              )
              return(NULL)
            }
            rv$seurat_integrated <- calculate_avg_expression(rv$seurat_integrated, selected_genes, "avg_gene_expr")
            vars_to_regress <- c(vars_to_regress, "avg_gene_expr")
          }
        }
        
        if ("gene_pathway" %in% input$correction_options) {
          selected_pathway <- input$pathway
          pathway_genes <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
            filter(gs_name == selected_pathway) %>%
            pull(gene_symbol)
          pathway_genes <- intersect(pathway_genes, rownames(rv$seurat_integrated))
          
          if (length(pathway_genes) > 0) {
            rv$seurat_integrated <- calculate_avg_expression(rv$seurat_integrated, pathway_genes, "avg_pathway_expr")
            vars_to_regress <- c(vars_to_regress, "avg_pathway_expr")
          } else {
            removeModal()
            shinyalert::shinyalert(
              title = "Error",
              text = "No genes found for the selected pathway in the dataset.",
              type = "error"
            )
            return(NULL)
          }
        }
        
        withProgress(message = 'Applying Data Correction...', value = 0, {
          incProgress(0.3, detail = "Scaling data...")
          rv$seurat_integrated <- ScaleData(rv$seurat_integrated, vars.to.regress = vars_to_regress, verbose = FALSE)
          
          incProgress(0.6, detail = "Running PCA...")
          rv$seurat_integrated <- RunPCA(rv$seurat_integrated, verbose = FALSE)
          
          incProgress(1, detail = "Running UMAP...")
          rv$seurat_integrated <- RunUMAP(rv$seurat_integrated, dims = 1:30, n.components = 3)
        })
        
        # Data Correction i??lemi sonras?? d??zeltilmi?? veriyi saklama
        rv$corrected_data <- rv$seurat_integrated
        
        removeModal()
        shinyalert::shinyalert(
          title = "Success",
          text = "Data correction has been applied successfully! You can now proceed to visualize the corrected data or move to the next step.",
          type = "success"
        )
        
        output$post_correction_plot <- renderPlotly({
          post_correction_plot_reactive()
        })
        
      }, error = function(e) {
        log_error <- function(error_message) {
          write(paste(Sys.time(), "-", error_message), file = "error_log.txt", append = TRUE)
        }
        log_error(paste("Data Correction Error: ", e$message))
        log_error(paste("Traceback: ", paste(capture.output(traceback()), collapse = "\n")))
        
        removeModal()
        showModal(modalDialog(
          title = "Data Correction Error",
          paste("An error occurred during data correction:", e$message),
          "Traceback:", paste(capture.output(traceback()), collapse = "\n"),
          easyClose = TRUE
        ))
      })
    }
    
    observe({
      req(rv$seurat_integrated)
      all_genes <- rownames(rv$seurat_integrated)
      updatePickerInput(session, "genes", choices = all_genes)
      
      pathways <- msigdbr(species = "Homo sapiens", category = "C2") %>%
        dplyr::select(gs_name) %>%
        distinct() %>%
        arrange(gs_name) %>%
        pull(gs_name)
      updateSelectInput(session, "pathway", choices = pathways, selected = pathways[1])
    })
    
    observeEvent(input$genes, {
      max_genes <- 100
      selected_genes <- input$genes
      
      if (length(selected_genes) > max_genes) {
        updatePickerInput(session, "genes", selected = selected_genes[1:max_genes])
        shinyalert::shinyalert(
          title = "Selection Limit Exceeded",
          text = paste("You can select up to", max_genes, "genes. Extra selections have been removed."),
          type = "warning"
        )
      }
    })
    
    observeEvent(input$run_correction, {
      req(rv$seurat_integrated)
      shinyalert::shinyalert(
        title = "Confirm Data Correction",
        text = "Are you sure you want to proceed with data correction using the selected parameters?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, Proceed",
        cancelButtonText = "Cancel",
        callbackR = function(x) {
          if (x) {
            perform_data_correction()
          }
        }
      )
    })
    
    post_correction_plot_reactive <- reactive({
      req(rv$seurat_integrated)
      req(input$post_correction_plot_type)
      req(input$post_label_by)
      
      if (!input$post_label_by %in% colnames(rv$seurat_integrated@meta.data)) {
        shinyalert::shinyalert(
          title = "Error",
          text = paste("The selected label", input$post_label_by, "is not found in the dataset."),
          type = "error"
        )
        return(NULL)
      }
      
      if (input$post_correction_plot_type == "UMAP") {
        p <- DimPlot(rv$seurat_integrated, reduction = "umap", group.by = input$post_label_by, label = TRUE) +
          ggtitle("UMAP Plot After Correction")
        ggplotly(p)
      } else if (input$post_correction_plot_type == "3D UMAP") {
        umap_coords <- Embeddings(rv$seurat_integrated, "umap")
        if (ncol(umap_coords) < 3) {
          shinyalert::shinyalert(
            title = "Error",
            text = "3D UMAP embeddings are not available. Please compute 3D UMAP first.",
            type = "error"
          )
          return(NULL)
        }
        plot_data <- data.frame(
          x = umap_coords[, 1],
          y = umap_coords[, 2],
          z = umap_coords[, 3],
          group = rv$seurat_integrated@meta.data[[input$post_label_by]]
        )
        plot_ly(
          data = plot_data,
          x = ~x, y = ~y, z = ~z,
          color = ~group,
          colors = viridis::viridis(length(unique(plot_data$group))),
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 2)
        ) %>% layout(title = "3D UMAP Plot After Correction")
      } else if (input$post_correction_plot_type == "t-SNE") {
        p <- DimPlot(rv$seurat_integrated, reduction = "tsne", group.by = input$post_label_by, label = TRUE) +
          ggtitle("t-SNE Plot After Correction")
        ggplotly(p)
      } else if (input$post_correction_plot_type == "3D t-SNE") {
        tsne_coords <- Embeddings(rv$seurat_integrated, "tsne")
        if (ncol(tsne_coords) < 3) {
          shinyalert::shinyalert(
            title = "Error",
            text = "3D t-SNE embeddings are not available. Please compute 3D t-SNE first.",
            type = "error"
          )
          return(NULL)
        }
        plot_data <- data.frame(
          x = tsne_coords[, 1],
          y = tsne_coords[, 2],
          z = tsne_coords[, 3],
          group = rv$seurat_integrated@meta.data[[input$post_label_by]]
        )
        plot_ly(
          data = plot_data,
          x = ~x, y = ~y, z = ~z,
          color = ~group,
          colors = viridis::viridis(length(unique(plot_data$group))),
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 2)
        ) %>% layout(title = "3D t-SNE Plot After Correction")
      }
    })
    
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
    
    # Violin Plot: Gen ekspresyonunu g??rselle??tirmek i??in
    output$violin_plot <- renderPlot({
      req(rv$seurat_integrated)
      req(input$genes)
      
      missing_genes <- setdiff(input$genes, rownames(rv$seurat_integrated))
      
      if (length(missing_genes) > 0) {
        shinyalert::shinyalert(
          title = "Missing Genes",
          text = paste("The following genes are not found in the dataset:", paste(missing_genes, collapse = ", ")),
          type = "error"
        )
        return(NULL) # E??er genler yoksa grafik ??izilmesin
      }
      
      # Violin plot'u cluster baz??nda ??izme (seurat_clusters kullan??l??yor)
      VlnPlot(rv$seurat_integrated, features = input$genes, group.by = "seurat_clusters", pt.size = 0) +
        ggtitle("Gene Expression Violin Plot by Cluster")
    })
    
    # Ridge Plot: H??cre d??ng??s?? faz??na g??re da????l??m?? g??rmek i??in
    output$ridge_plot <- renderPlot({
      req(rv$seurat_integrated)
      
      RidgePlot(rv$seurat_integrated, features = c("S.Score", "G2M.Score"), group.by = "Phase") +
        ggtitle("Cell Cycle Phase Distribution Ridge Plot")
    })
    
    observeEvent(input$next_step, {
      rv$proceed_to_next_step <- TRUE
      rv$corrected_data <- rv$corrected_data
    })
  })
}
