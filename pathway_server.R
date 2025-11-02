# ============================================================
# 🧬 EVOscope Pathway Diversity Server (Final Refined Version)
# ============================================================

pathway_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ─────────────────────────────
    # Reactive storage
    # ─────────────────────────────
    rv <- reactiveValues(
      expr = NULL,
      scores = NULL,
      cluster_means = NULL,
      top30_overall = NULL,
      top10_top30 = NULL,
      heatmap_plot = NULL,
      violin_plot = NULL,
      bar_plot = NULL
    )
    
    # ─────────────────────────────
    # Run Pathway Analysis
    # ─────────────────────────────
    observeEvent(input$run_pathway, {
      output$status_pathway <- renderText("⏳ Initializing pathway analysis…")
      
      tryCatch({
        obj <- get("global_dataset", envir = .GlobalEnv)
        req(obj)
        
        withProgress(message = "Running pathway analysis…", value = 0, {
          incProgress(0.1, detail = "Extracting expression matrix…")
          expr <- get_expr_matrix(obj)
          rv$expr <- expr
          
          if (is.null(gene_sets) || length(gene_sets) < 5)
            stop("No valid GMT gene sets found — check global.R path.")
          
          incProgress(0.4, detail = "Building AUCell rankings…")
          rankings <- AUCell_buildRankings(expr, plotStats = FALSE)
          
          incProgress(0.7, detail = "Calculating AUC per pathway…")
          auc <- AUCell_calcAUC(gene_sets, rankings)
          scores <- as.data.frame(t(getAUC(auc)))
          
          # Add cluster info
          if ("seurat_clusters" %in% colnames(obj@meta.data)) {
            scores$Cluster <- obj$seurat_clusters
          } else {
            scores$Cluster <- as.character(Idents(obj))
          }
          rownames(scores) <- colnames(expr)
          rv$scores <- scores
          
          incProgress(0.85, detail = "Aggregating cluster means…")
          cluster_means <- scores %>%
            dplyr::group_by(Cluster) %>%
            dplyr::summarise(dplyr::across(dplyr::starts_with("HALLMARK"),
                                           \(x) mean(x, na.rm = TRUE))) %>%
            as.data.frame()
          rv$cluster_means <- cluster_means
          
          incProgress(0.9, detail = "Extracting top genes…")
          
          # Helper: high-AUC cell IDs
          get_high_auc_cells <- function(p, df, top_frac = 0.1) {
            if (!p %in% colnames(df)) return(character(0))
            vals <- df[[p]]; names(vals) <- rownames(df)
            vals <- vals[!is.na(vals)]
            if (!length(vals)) return(character(0))
            cutoff <- quantile(vals, 1 - top_frac, na.rm = TRUE)
            names(vals[vals >= cutoff])
          }
          
          # Helper: top genes per pathway (with filtering)
          get_pathway_genes <- function(p, df, expr, top_frac = 0.1, top_n = 30,
                                        drop_genes = NULL) {
            # Define default dropped genes if not supplied
            if (is.null(drop_genes)) {
              drop_genes <- unique(c("MALAT1", "NEAT1",
                                     grep("^MT-", rownames(expr), value = TRUE),
                                     grep("^RPL", rownames(expr), value = TRUE),
                                     grep("^RPS", rownames(expr), value = TRUE)))
            }
            cells <- get_high_auc_cells(p, df, top_frac)
            cells <- intersect(cells, colnames(expr))
            if (!length(cells)) return(NA_character_)
            sub <- expr[, cells, drop = FALSE]
            gmeans <- Matrix::rowMeans(sub)
            genes <- names(sort(gmeans, decreasing = TRUE))
            genes <- setdiff(genes, drop_genes)
            paste(head(genes, top_n), collapse = ", ")
          }
          
          # Numeric columns only
          score_num <- scores[, sapply(scores, is.numeric), drop = FALSE]
          
          all_tbl <- data.frame(
            Pathway = colnames(score_num),
            TopGenes = sapply(colnames(score_num), function(p)
              get_pathway_genes(p, score_num, expr, top_frac = 0.1, top_n = 30))
          )
          rv$top30_overall <- all_tbl
          
          mean_auc <- colMeans(score_num, na.rm = TRUE)
          top10 <- names(sort(mean_auc, decreasing = TRUE))[1:10]
          top10_tbl <- data.frame(
            Pathway = top10,
            TopGenes = sapply(top10, function(p)
              get_pathway_genes(p, score_num, expr, top_frac = 0.1, top_n = 30))
          )
          rv$top10_top30 <- top10_tbl
          
          incProgress(1, detail = "Done.")
        })
        
        output$status_pathway <- renderText("✅ Pathway analysis completed successfully.")
        
        output$summary_pathway <- renderUI({
          nP <- ncol(rv$scores) - 1
          nC <- nrow(rv$scores)
          nCl <- length(unique(rv$scores$Cluster))
          tags$p(strong("Summary: "),
                 sprintf("%d pathways × %d cells across %d clusters.", nP, nC, nCl))
        })
        
      }, error = function(e) {
        output$status_pathway <- renderText(paste("❌ Error in pathway analysis:", e$message))
      })
    })
    
    # ─────────────────────────────
    # Download Handlers
    # ─────────────────────────────
    output$dl_full_scores <- downloadHandler(
      filename = function() "pathway_activity_scores.csv",
      content = function(file) write.csv(rv$scores, file, row.names = TRUE)
    )
    output$dl_cluster_means <- downloadHandler(
      filename = function() "pathway_activity_cluster_means.csv",
      content = function(file) write.csv(rv$cluster_means, file, row.names = FALSE)
    )
    output$dl_top30_overall <- downloadHandler(
      filename = function() "pathway_top30_genes_overall.csv",
      content = function(file) write.csv(rv$top30_overall, file, row.names = FALSE)
    )
    output$dl_top10_top30 <- downloadHandler(
      filename = function() "top10_pathways_top30_genes.csv",
      content = function(file) write.csv(rv$top10_top30, file, row.names = FALSE)
    )
    
    # ─────────────────────────────
    # Heatmap Plot
    # ─────────────────────────────
    observeEvent(input$plot_heatmap, {
      req(rv$scores)
      df <- rv$scores
      
      if (!"Cluster" %in% colnames(df)) {
        showNotification("❌ No 'Cluster' column found.", type = "error")
        return(NULL)
      }
      
      df$Cluster <- as.factor(as.character(df$Cluster))
      cluster_means <- df %>%
        dplyr::group_by(Cluster) %>%
        dplyr::summarise(dplyr::across(dplyr::starts_with("HALLMARK"),
                                       \(x) mean(x, na.rm = TRUE))) %>%
        as.data.frame()
      
      rownames(cluster_means) <- cluster_means$Cluster
      cluster_means$Cluster <- NULL
      if (nrow(cluster_means) == 0 || ncol(cluster_means) == 0) {
        showNotification("⚠️ No valid data to plot heatmap.", type = "warning")
        return(NULL)
      }
      
      mat <- t(as.matrix(cluster_means))
      cluster_names <- colnames(mat)
      if (all(grepl("^[0-9]+$", cluster_names))) {
        mat <- mat[, order(as.numeric(cluster_names)), drop = FALSE]
      } else {
        mat <- mat[, order(cluster_names), drop = FALSE]
      }
      
      vars <- apply(mat, 1, var, na.rm = TRUE)
      vars <- vars[!is.na(vars)]
      if (length(vars) == 0) {
        showNotification("⚠️ All pathway variances are NA.", type = "warning")
        return(NULL)
      }
      
      top_pathways <- names(sort(vars, decreasing = TRUE))[1:min(20, length(vars))]
      mat_top <- mat[top_pathways, , drop = FALSE]
      
      rv$heatmap_plot <- pheatmap::pheatmap(
        mat_top,
        color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        border_color = NA,
        fontsize_row = 9,
        fontsize_col = 10,
        main = "Top 20 Variable Pathways Across Clusters (Ordered)"
      )
      output$heatmap <- renderPlot(rv$heatmap_plot)
    })
    
    # ─────────────────────────────
    # Violin Plot
    # ─────────────────────────────
    observeEvent(input$plot_violin, {
      req(rv$scores)
      df <- rv$scores
      num_cols <- sapply(df, is.numeric)
      pvars <- apply(df[, num_cols, drop = FALSE], 2, var, na.rm = TRUE)
      top5 <- names(sort(pvars, decreasing = TRUE))[1:5]
      m <- reshape2::melt(df[, c(top5, "Cluster")], id.vars = "Cluster")
      rv$violin_plot <- ggplot(m, aes(x = Cluster, y = value, fill = Cluster)) +
        geom_violin(scale = "width") +
        facet_wrap(~ variable, ncol = 2, scales = "free_y") +
        theme_minimal(base_size = 14) +
        theme(panel.grid = element_blank(),
              strip.text = element_text(face = "bold")) +
        labs(title = "Top 5 Diverse Pathways", y = "AUC", x = "Cluster")
      output$violin <- renderPlot(rv$violin_plot)
    })
    
    # ─────────────────────────────
    # Bar Plot (Top 10 pathways × 5 genes, clean)
    # ─────────────────────────────
    observeEvent(input$plot_bar, {
      req(rv$top10_top30)
      df <- rv$top10_top30
      df$PathwayShort <- gsub("^HALLMARK_", "", df$Pathway)
      df$MeanAUC <- colMeans(rv$scores[, df$Pathway, drop = FALSE], na.rm = TRUE)
      
      # Keep top 5 genes only
      df$TopGenesShort <- sapply(df$TopGenes, function(x) {
        genes <- unlist(strsplit(as.character(x), ",\\s*"))
        genes <- genes[!is.na(genes) & genes != ""]
        label <- paste(head(genes, 5), collapse = ", ")
        stringr::str_wrap(label, width = 50)
      })
      
      rv$bar_plot <- ggplot(df, aes(x = reorder(PathwayShort, MeanAUC), y = MeanAUC, fill = MeanAUC)) +
        geom_col(width = 0.6) +
        coord_flip(clip = "off", expand = TRUE) +
        geom_text(aes(label = TopGenesShort),
                  hjust = -0.05, size = 3.2, color = "black") +
        scale_fill_gradient(low = "#a8caff", high = "#08306b") +
        labs(title = "Top 10 Pathways with Top 5 Genes",
             x = "Pathway", y = "Mean AUC") +
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          legend.position = "none",
          plot.margin = margin(10, 120, 10, 10)
        )
      output$barplot <- renderPlot(rv$bar_plot)
    })
    
    # ─────────────────────────────
    # PNG Download Handlers
    # ─────────────────────────────
    output$dl_heatmap_png <- downloadHandler(
      filename = function() "pathway_heatmap.png",
      content = function(file) {
        png(file, width = 2400, height = 2000, res = 300)
        grid::grid.newpage()
        grid::grid.draw(rv$heatmap_plot$gtable)
        dev.off()
      }
    )
    output$dl_violin_png <- downloadHandler(
      filename = function() "pathway_violin.png",
      content = function(file) ggsave(file, plot = rv$violin_plot, dpi = 1000, width = 10, height = 10)
    )
    output$dl_bar_png <- downloadHandler(
      filename = function() "pathway_bar.png",
      content = function(file) ggsave(file, plot = rv$bar_plot, dpi = 1000, width = 14, height = 10)
    )
  })
}
