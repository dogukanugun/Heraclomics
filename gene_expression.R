# gene_expression.R
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
library(shinycssloaders)  # for withSpinner()

# UI for Differential Expression Analysis 
differentialExpressionUI <- function(id) {
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
                       value = 10, min = 0, step = 1)
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
        actionButton(ns("run_analysis"), "Run Analysis",
                     icon = icon("play"), class = "btn-success"),
        downloadButton(ns("download_results"), "Download Full Results")
      )
    ),
    fluidRow(
      box(
        title = "Summary Statistics",
        width = 12,
        status = "info",
        solidHeader = TRUE,
        verbatimTextOutput(ns("summary_stats"))
      )
    ),
    fluidRow(
      box(
        title = "Results",
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        tabsetPanel(
          tabPanel(
            "Table",
            DT::dataTableOutput(ns("de_results_table")) %>% 
              shinycssloaders::withSpinner()
          ),
          tabPanel(
            "Volcano Plot",
            fluidRow(
              column(4,
                     selectInput(ns("volcano_style"), "Plot Style:",
                                 choices = c("EnhancedVolcano", "ggplot2"),
                                 selected = "EnhancedVolcano")
              ),
              column(4,
                     numericInput(ns("volcano_label_n"), "Top genes to label:",
                                  value = 10, min = 1, max = 50)
              ),
              column(4,
                     checkboxInput(ns("show_all_labels"), "Show all sig. labels", 
                                   value = FALSE)
              )
            ),
            plotOutput(ns("static_volcano"), height = "600px"),
            plotlyOutput(ns("interactive_volcano"), height = "600px")
          ),
          tabPanel(
            "Heatmap",
            fluidRow(
              column(4,
                     numericInput(ns("heatmap_genes"), "Number of genes:",
                                  value = 50, min = 1, max = 200)
              ),
              column(4,
                     selectInput(ns("heatmap_scale"), "Scale:",
                                 choices = c("row", "column", "none"),
                                 selected = "row")
              ),
              column(4,
                     checkboxInput(ns("show_rownames"), "Show gene names", 
                                   value = TRUE)
              )
            ),
            plotOutput(ns("heatmap_plot"), height = "800px") %>% 
              shinycssloaders::withSpinner()
          ),
          tabPanel(
            "MA Plot",
            plotOutput(ns("ma_plot"), height = "600px")
          )
        )
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
      stats = NULL
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
                gene = gene,
                log2FoldChange = avg_log2FC,
                pvalue = p_val,
                padj = p_val_adj
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
                log2FoldChange = avg_log2FC,
                pvalue = p_val,
                padj = p_val_adj
              )
          }
          
          de_results$stats <- list(
            method = "Seurat",
            groups = paste(groups, collapse = " vs "),
            genes_tested = nrow(de_results$raw)
          )
          
        } else if (input$method == "deseq2") {
          # DESeq2 analysis
          counts <- as.matrix(GetAssayData(rv$corrected_data, slot = "counts"))
          metadata <- rv$corrected_data@meta.data
          
          # Filter low count genes
          keep <- rowSums(counts >= input$min_counts) >= ceiling(ncol(counts)*0.1)
          counts_filtered <- counts[keep, ]
          genes_removed <- sum(!keep)
          
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
          keep <- rowSums(counts > 0) >= ceiling(ncol(counts) * 0.1)
          counts_filtered <- counts[keep, ]
          genes_removed <- sum(!keep)
          
          # Create DGEList
          dge <- DGEList(counts = counts_filtered,
                         group = metadata[[input$design_factor]])
          dge <- calcNormFactors(dge)
          
          # Design matrix
          design <- model.matrix(~ metadata[[input$design_factor]])
          
          # Estimate dispersion
          dge <- estimateDisp(dge, design)
          
          # GLM fit and test
          fit <- glmQLFit(dge, design)
          qlf <- glmQLFTest(fit)
          
          # Get results
          de_results$raw <- topTags(qlf, n = Inf)$table %>%
            tibble::rownames_to_column("gene") %>%
            rename(
              log2FoldChange = logFC,
              pvalue = PValue,
              padj = FDR
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
      if (!is.null(de_results$stats$genes_removed)) {
        cat("Genes filtered:", de_results$stats$genes_removed, "\n")
      }
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
        formatRound(columns = c("log2FoldChange", "pvalue", "padj"), digits = 15)
    })
    
    # Static volcano plot with controlled labeling
    output$static_volcano <- renderPlot({
      req(de_results$raw, input$volcano_style)
      
      plot_data <- de_results$raw %>%
        mutate(
          significance = case_when(
            padj < input$pval_cutoff & log2FoldChange > input$logfc_cutoff ~ "Up",
            padj < input$pval_cutoff & log2FoldChange < -input$logfc_cutoff ~ "Down",
            TRUE ~ "NS"
          )
        )
      
      # Get top genes to label
      top_up <- plot_data %>%
        filter(significance == "Up") %>%
        arrange(padj) %>%
        head(input$volcano_label_n)
      
      top_down <- plot_data %>%
        filter(significance == "Down") %>%
        arrange(padj) %>%
        head(input$volcano_label_n)
      
      if (input$volcano_style == "EnhancedVolcano") {
        EnhancedVolcano(
          plot_data,
          lab = ifelse(plot_data$gene %in% c(top_up$gene, top_down$gene) | input$show_all_labels, 
                       plot_data$gene, ""),
          x = 'log2FoldChange',
          y = 'padj',
          pCutoff = input$pval_cutoff,
          FCcutoff = input$logfc_cutoff,
          pointSize = 2.0,
          labSize = 4.0,
          title = "Volcano Plot",
          subtitle = paste("Comparison:", de_results$stats$groups),
          selectLab = c(top_up$gene, top_down$gene),
          drawConnectors = TRUE,
          max.overlaps = 20
        )
      } else {
        ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
          geom_point(aes(color = significance), alpha = 0.6) +
          scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
          geom_hline(yintercept = -log10(input$pval_cutoff), linetype = "dashed", color = "black") +
          geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff), linetype = "dashed", color = "black") +
          geom_text(data = top_up, aes(label = gene), vjust = -0.5, hjust = 0.5, size = 3) +
          geom_text(data = top_down, aes(label = gene), vjust = -0.5, hjust = 0.5, size = 3) +
          labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", color = "Significance") +
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
          label = ifelse((padj < input$pval_cutoff & abs(log2FoldChange) > input$logfc_cutoff) &
                           (gene %in% c(
                             head(arrange(filter(., significance == "Up"), padj)$gene, input$volcano_label_n),
                             head(arrange(filter(., significance == "Down"), padj)$gene, input$volcano_label_n)
                           ) | input$show_all_labels),
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
        geom_text(aes(label = label), vjust = -0.5, hjust = 0.5, size = 3) +
        theme_minimal()
      
      ggplotly(p, tooltip = "text") %>%
        layout(legend = list(orientation = "h", y = -0.2),
               height = 600)
    })
    
    # Heatmap plot
    output$heatmap_plot <- renderPlot({
      req(de_results$filtered, rv$corrected_data, input$design_factor)
      
      tryCatch({
        # Get top significant genes
        top_genes <- de_results$filtered %>%
          arrange(padj) %>%
          head(input$heatmap_genes) %>%
          pull(gene)
        
        if (length(top_genes) == 0) {
          return(ggplot() + 
                   annotate("text", x = 1, y = 1, 
                            label = "No significant genes at current thresholds") +
                   theme_void())
        }
        
        # Get normalized expression data - ensure genes exist in data
        available_genes <- top_genes[top_genes %in% rownames(rv$corrected_data)]
        if (length(available_genes) == 0) {
          return(ggplot() + 
                   annotate("text", x = 1, y = 1, 
                            label = "No significant genes available in the dataset") +
                   theme_void())
        }
        
        mat <- GetAssayData(rv$corrected_data, slot = "data")[available_genes, , drop = FALSE]
        
        # Remove any NA/Inf values
        mat[is.na(mat) | is.infinite(mat)] <- 0
        
        # Prepare annotation - ensure it matches column names
        annotation_col <- data.frame(
          Group = factor(rv$corrected_data@meta.data[[input$design_factor]])
        )
        rownames(annotation_col) <- colnames(mat)
        
        # Handle scaling properly
        if (input$heatmap_scale == "row") {
          mat <- t(scale(t(mat)))
          mat[is.na(mat)] <- 0  # Handle any NA from scaling
        } else if (input$heatmap_scale == "column") {
          mat <- scale(mat)
          mat[is.na(mat)] <- 0
        }
        
        # Set a reasonable font size based on number of genes
        fontsize_row <- ifelse(input$show_rownames, 
                               min(8, 300/length(available_genes)), 
                               0)
        
        # Plot heatmap
        pheatmap(
          mat,
          scale = "none",  # We handle scaling manually above
          clustering_method = "complete",
          show_colnames = FALSE,
          show_rownames = input$show_rownames,
          annotation_col = annotation_col,
          color = colorRampPalette(c("blue", "white", "red"))(100),
          main = paste("Top", length(available_genes), "Differentially Expressed Genes"),
          fontsize_row = fontsize_row,
          fontsize_col = 8,
          silent = TRUE
        )
        
      }, error = function(e) {
        ggplot() + 
          annotate("text", x = 1, y = 1, 
                   label = paste("Error generating heatmap:", e$message)) +
          theme_void()
      })
    })
    
    # MA plot
    output$ma_plot <- renderPlot({
      req(de_results$raw)
      
      tryCatch({
        # Calculate proper mean expression (average across all samples)
        counts <- GetAssayData(rv$corrected_data, slot = "counts")
        mean_expr <- rowMeans(counts)
        
        plot_data <- de_results$raw %>%
          mutate(
            mean_expr = mean_expr[gene],
            significance = case_when(
              padj < input$pval_cutoff & log2FoldChange > input$logfc_cutoff ~ "Up",
              padj < input$pval_cutoff & log2FoldChange < -input$logfc_cutoff ~ "Down",
              TRUE ~ "NS"
            )
          ) %>%
          filter(!is.na(mean_expr), !is.infinite(log2FoldChange))
        
        # Cap extreme values for better visualization
        plot_data <- plot_data %>%
          mutate(
            log2FoldChange = pmin(pmax(log2FoldChange, -5), 5),
            mean_expr = log10(mean_expr + 1)  # Log-transform mean expression
          )
        
        ggplot(plot_data, aes(x = mean_expr, y = log2FoldChange)) +
          geom_point(aes(color = significance), alpha = 0.6, size = 1) +
          scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
          geom_hline(yintercept = c(-input$logfc_cutoff, input$logfc_cutoff), 
                     linetype = "dotted", color = "black") +
          labs(x = "Log10(Mean Expression + 1)", y = "Log2 Fold Change",
               title = "MA Plot",
               subtitle = paste("Comparison:", de_results$stats$groups)) +
          theme_minimal() +
          theme(legend.position = "bottom")
        
      }, error = function(e) {
        ggplot() + 
          annotate("text", x = 1, y = 1, 
                   label = paste("Error generating MA plot:", e$message)) +
          theme_void()
      })
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
