library(shiny)
library(Seurat)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(DT)
library(plotly)
library(shinyalert)
library(pheatmap)
library(reshape2)
library(EnhancedVolcano)

# UI for Differential Expression Analysis
differential_expression_ui <- function(id) {
  ns <- NS(id)
  tagList(
    useShinyalert(),
    h2("Differential Expression Analysis"),
    fluidRow(
      box(
        title = "Analysis Parameters",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        selectInput(ns("method"), "Select Method:",
                    choices = c("Seurat" = "seurat",
                                "DESeq2" = "deseq2", 
                                "edgeR" = "edger"),
                    selected = "seurat"),
        selectInput(ns("design_factor"), "Select Design Factor:",
                    choices = NULL),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'deseq2'", ns("method")),
          checkboxInput(ns("use_poscounts"),
                        "Use alternative size factor estimation (handles zero counts)",
                        value = FALSE),
          numericInput(ns("min_counts"), "Minimum count threshold:",
                       value = 10, min = 0, step = 1),
          helpText("Genes with counts below this in all samples will be filtered.")
        ),
        fluidRow(
          column(6,
                 numericInput(ns("pval_cutoff"), "Adjusted p-value cutoff:",
                              value = 0.05, min = 0, max = 1, step = 0.01)
          ),
          column(6,
                 numericInput(ns("logfc_cutoff"), "Log2 FC cutoff:",
                              value = 0.5, min = 0, step = 0.1)
          )
        ),
        numericInput(ns("top_n"), "Number of top genes to display:",
                     value = 20, min = 1, step = 1),
        actionButton(ns("run_analysis"), "Run Analysis",
                     icon = icon("play"), class = "btn-success"),
        downloadButton(ns("download_results"), "Download Full Results")
      ),
      box(
        title = "Summary Statistics",
        width = 12,
        status = "info",
        solidHeader = TRUE,
        verbatimTextOutput(ns("summary_stats"))
      )
    ),
    tabBox(
      width = 12,
      tabPanel(
        "Results Table",
        DT::dataTableOutput(ns("de_results_table"))
      ),
      tabPanel(
        "Volcano Plot",
        fluidRow(
          column(3,
                 selectInput(ns("volcano_style"), "Plot Style:",
                             choices = c("EnhancedVolcano", "ggplot2"),
                             selected = "EnhancedVolcano"),
                 numericInput(ns("volcano_label_n"), "Number of genes to label:",
                              value = 10, min = 0)
          ),
          plotOutput(ns("static_volcano")),
          plotlyOutput(ns("interactive_volcano")) %>% 
            shinycssloaders::withSpinner()
        )
      ),
      tabPanel(
        "Heatmap",
        fluidRow(
          column(3,
                 selectInput(ns("heatmap_scale"), "Scale:",
                             choices = c("row", "column", "none"),
                             selected = "row"),
                 checkboxInput(ns("show_rownames"), "Show gene names", value = TRUE)
          )
        ),
        plotOutput(ns("heatmap_plot")) %>% 
          shinycssloaders::withSpinner()
      ),
      tabPanel(
        "MA Plot",
        plotOutput(ns("ma_plot"))
      )
    )
  )
}


# Server for Differential Expression Analysis
differential_expression_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values to store results
    de_results <- reactiveValues(
      raw = NULL,
      filtered = NULL,
      stats = NULL,
      metadata = NULL
    )
    
    # Update design factor choices when data changes
    observe({
      req(rv$corrected_data)
      metadata <- rv$corrected_data@meta.data
      factors <- colnames(metadata)[sapply(metadata, function(x) 
        is.factor(x) | is.character(x))]
      updateSelectInput(session, "design_factor", choices = factors)
    })
    
    # Main analysis function
    observeEvent(input$run_analysis, {
      req(rv$corrected_data, input$design_factor)
      
      tryCatch({
        showModal(modalDialog(
          title = "Running Differential Expression Analysis",
          "This may take several minutes...",
          footer = NULL,
          easyClose = FALSE
        ))
        
        # Get group information
        metadata <- rv$corrected_data@meta.data
        groups <- unique(metadata[[input$design_factor]])
        
        if (length(groups) < 2) {
          stop("Selected design factor must have at least 2 groups")
        }
        
        # Run appropriate DE method
        if (input$method == "seurat") {
          Idents(rv$corrected_data) <- input$design_factor
          
          if (length(groups) > 2) {
            # Handle multi-group comparisons
            de_results$raw <- FindAllMarkers(
              rv$corrected_data,
              logfc.threshold = input$logfc_cutoff,
              only.pos = FALSE,
              min.pct = 0.1
            ) %>%
              rename(
                gene = "gene",
                log2FoldChange = "avg_log2FC",
                pvalue = "p_val",
                padj = "p_val_adj"
              )
          } else {
            # Standard two-group comparison
            de_results$raw <- FindMarkers(
              rv$corrected_data,
              ident.1 = groups[1],
              ident.2 = groups[2],
              logfc.threshold = input$logfc_cutoff
            ) %>%
              tibble::rownames_to_column("gene") %>%
              rename(
                log2FoldChange = "avg_log2FC",
                pvalue = "p_val",
                padj = "p_val_adj"
              )
          }
          
          de_results$stats <- list(
            method = "Seurat",
            groups = paste(groups, collapse = " vs "),
            genes_tested = nrow(de_results$raw),
            genes_removed = 0
          )
          
        } else if (input$method == "deseq2") {
          # DESeq2 analysis
          counts <- as.matrix(GetAssayData(rv$corrected_data, slot = "counts"))
          metadata <- rv$corrected_data@meta.data
          
          # Filter low count genes
          keep <- rowSums(counts >= input$min_counts) >= ceiling(ncol(counts)*0.1)
          counts_filtered <- counts[keep, ]
          genes_removed <- sum(!keep)
          
          if (nrow(counts_filtered) == 0) {
            stop("No genes passed count filtering. Adjust minimum count threshold.")
          }
          
          # Create DESeqDataSet
          dds <- DESeqDataSetFromMatrix(
            countData = counts_filtered,
            colData = metadata,
            design = as.formula(paste("~", input$design_factor))
          )
          
          # Size factor estimation
          if (input$use_poscounts) {
            dds <- estimateSizeFactors(dds, type = "poscounts")
          } else {
            dds <- estimateSizeFactors(dds)
          }
          
          # Differential expression
          dds <- DESeq(dds)
          res <- results(dds, alpha = input$pval_cutoff)
          
          de_results$raw <- as.data.frame(res) %>%
            tibble::rownames_to_column("gene") %>%
            mutate(padj = ifelse(is.na(padj), 1, padj))
          
          de_results$stats <- list(
            method = "DESeq2",
            groups = paste(groups, collapse = " vs "),
            genes_tested = nrow(de_results$raw),
            genes_removed = genes_removed
          )
        } else if (input$method == "edger") {
          # edgeR analysis
          counts <- as.matrix(GetAssayData(rv$corrected_data, slot = "counts"))
          metadata <- rv$corrected_data@meta.data
          
          # Filter low count genes
          keep <- rowSums(counts > 0) >= ceiling(ncol(counts)*0.1)
          counts_filtered <- counts[keep, ]
          genes_removed <- sum(!keep)
          
          # Create DGEList
          dge <- DGEList(counts = counts_filtered,
                         group = metadata[[input$design_factor]])
          dge <- calcNormFactors(dge)
          
          # Design matrix
          design <- model.matrix(~metadata[[input$design_factor]])
          
          # Estimate dispersion
          dge <- estimateDisp(dge, design)
          
          # GLM fit and test
          fit <- glmQLFit(dge, design)
          qlf <- glmQLFTest(fit)
          
          # Get results
          de_results$raw <- topTags(qlf, n = Inf)$table %>%
            tibble::rownames_to_column("gene") %>%
            rename(
              log2FoldChange = "logFC",
              pvalue = "PValue",
              padj = "FDR"
            )
          
          de_results$stats <- list(
            method = "edgeR",
            groups = paste(groups, collapse = " vs "),
            genes_tested = nrow(de_results$raw),
            genes_removed = genes_removed
          )
        }
        
        # Filter results based on thresholds
        de_results$filtered <- de_results$raw %>%
          filter(padj < input$pval_cutoff,
                 abs(log2FoldChange) > input$logfc_cutoff) %>%
          arrange(padj)
        
        removeModal()
        shinyalert("Success", "Analysis completed successfully!", type = "success")
        
      }, error = function(e) {
        removeModal()
        shinyalert("Error", paste("Analysis failed:", e$message), type = "error")
      })
    })
    
    # Summary statistics output
    output$summary_stats <- renderPrint({
      req(de_results$stats)
      
      cat("Differential Expression Analysis Summary\n")
      cat("======================================\n")
      cat("Method:", de_results$stats$method, "\n")
      cat("Comparison:", de_results$stats$groups, "\n")
      cat("Genes tested:", de_results$stats$genes_tested, "\n")
      cat("Genes filtered:", de_results$stats$genes_removed, "\n")
      if (!is.null(de_results$filtered)) {
        cat("Significant genes:", nrow(de_results$filtered), "\n")
        cat("Upregulated:", sum(de_results$filtered$log2FoldChange > 0), "\n")
        cat("Downregulated:", sum(de_results$filtered$log2FoldChange < 0), "\n")
      }
    })
    
    # Interactive results table
    output$de_results_table <- DT::renderDataTable({
      req(de_results$filtered)
      
      DT::datatable(
        de_results$filtered,
        extensions = c('Buttons', 'Scroller'),
        options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'colvis'),
          scrollX = TRUE,
          scrollY = "500px",
          scroller = TRUE,
          pageLength = 10
        ),
        rownames = FALSE,
        selection = 'single'
      ) %>%
        formatRound(columns = c("log2FoldChange", "pvalue", "padj"), digits = 4)
    })
    
    # Static volcano plot
    output$static_volcano <- renderPlot({
      req(de_results$raw, input$volcano_style)
      
      if (input$volcano_style == "EnhancedVolcano") {
        EnhancedVolcano(
          de_results$raw,
          lab = de_results$raw$gene,
          x = 'log2FoldChange',
          y = 'padj',
          pCutoff = input$pval_cutoff,
          FCcutoff = input$logfc_cutoff,
          pointSize = 2.0,
          labSize = 4.0,
          title = "Volcano Plot",
          subtitle = paste("Comparison:", de_results$stats$groups)
        )
      } else {
        ggplot(de_results$raw, aes(x = log2FoldChange, y = -log10(padj))) +
          geom_point(aes(color = ifelse(padj < input$pval_cutoff & 
                                          abs(log2FoldChange) > input$logfc_cutoff,
                                        "Significant", "Not significant")), 
                     alpha = 0.6) +
          scale_color_manual(values = c("Significant" = "red", 
                                        "Not significant" = "gray")) +
          geom_hline(yintercept = -log10(input$pval_cutoff), 
                     linetype = "dashed", color = "blue") +
          geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff), 
                     linetype = "dashed", color = "blue") +
          labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
               color = "Significance") +
          theme_minimal() +
          theme(legend.position = "bottom")
      }
    })
    
    # Interactive volcano plot
    output$interactive_volcano <- renderPlotly({
      req(de_results$raw)
      
      plot_data <- de_results$raw %>%
        mutate(
          significance = case_when(
            padj < input$pval_cutoff & log2FoldChange > input$logfc_cutoff ~ "Up",
            padj < input$pval_cutoff & log2FoldChange < -input$logfc_cutoff ~ "Down",
            TRUE ~ "NS"
          ),
          label = ifelse(padj < input$pval_cutoff & 
                           abs(log2FoldChange) > input$logfc_cutoff,
                         gene, "")
        )
      
      p <- ggplot(plot_data, 
                  aes(x = log2FoldChange, y = -log10(padj),
                      text = paste0(
                        "Gene: ", gene, "\n",
                        "Log2FC: ", round(log2FoldChange, 2), "\n",
                        "P.adj: ", format.pval(padj, digits = 2), "\n",
                        "Significance: ", significance)
                  )) +
        geom_point(aes(color = significance, size = significance), alpha = 0.7) +
        scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
        scale_size_manual(values = c("Up" = 2, "Down" = 2, "NS" = 1)) +
        geom_text(aes(label = label), vjust = 1, hjust = 1, size = 3) +
        theme_minimal()
      
      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "h", y = -0.2))
    })
    
    # Heatmap plot
    output$heatmap_plot <- renderPlot({
      req(de_results$filtered, rv$corrected_data)
      
      top_genes <- de_results$filtered %>%
        arrange(padj) %>%
        head(input$top_n) %>%
        pull(gene)
      
      if (length(top_genes) == 0) {
        return(ggplot() + 
                 annotate("text", x = 1, y = 1, 
                          label = "No significant genes at current thresholds") +
                 theme_void())
      }
      
      # Get normalized expression data
      if (input$method == "seurat") {
        mat <- GetAssayData(rv$corrected_data, slot = "data")[top_genes, ]
      } else {
        mat <- GetAssayData(rv$corrected_data, slot = "counts")[top_genes, ]
        mat <- log1p(mat) # Log transform count data
      }
      
      # Prepare annotation
      annotation_col <- data.frame(
        Group = rv$corrected_data@meta.data[[input$design_factor]]
      )
      rownames(annotation_col) <- colnames(mat)
      
      # Plot heatmap
      pheatmap(
        mat,
        scale = input$heatmap_scale,
        clustering_method = "complete",
        show_colnames = FALSE,
        show_rownames = input$show_rownames,
        annotation_col = annotation_col,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        main = paste("Top", length(top_genes), "Differentially Expressed Genes"),
        fontsize_row = ifelse(input$show_rownames, 8, 0),
        silent = TRUE
      )
    })
    
    # MA plot
    output$ma_plot <- renderPlot({
      req(de_results$raw)
      
      plot_data <- de_results$raw %>%
        mutate(mean_expr = (log2FoldChange + mean(log2FoldChange))/2,
               significance = ifelse(padj < input$pval_cutoff, "Significant", "NS"))
      
      ggplot(plot_data, aes(x = mean_expr, y = log2FoldChange)) +
        geom_point(aes(color = significance), alpha = 0.6) +
        scale_color_manual(values = c("Significant" = "red", "NS" = "gray")) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = c(-input$logfc_cutoff, input$logfc_cutoff), 
                   linetype = "dotted", color = "blue") +
        labs(x = "Mean Expression", y = "Log2 Fold Change",
             title = "MA Plot") +
        theme_minimal()
    })
    
    # Download handler
    output$download_results <- downloadHandler(
      filename = function() {
        paste("DE_results_", input$method, "_", 
              gsub(" ", "_", de_results$stats$groups), "_", 
              Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(de_results$raw, file, row.names = FALSE)
      }
    )
  })
}