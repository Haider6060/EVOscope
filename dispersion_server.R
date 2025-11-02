dispersion_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    disp_df_rv    <- reactiveVal(NULL)
    expr_rv       <- reactiveVal(NULL)
    bar_plot_obj  <- reactiveVal(NULL)
    umap_plot_obj <- reactiveVal(NULL)
    
    # ─────────── Universal Expression Extractor ───────────
    get_expr_matrix <- function(obj) {
      assay <- DefaultAssay(obj)
      aobj  <- obj[[assay]]
      
      if ("layers" %in% slotNames(aobj)) {
        data_layers <- grep("^data", names(aobj@layers), value = TRUE)
        if (length(data_layers) == 0) {
          obj <- NormalizeData(obj, verbose = FALSE)
          aobj <- obj[[assay]]
          data_layers <- grep("^data", names(aobj@layers), value = TRUE)
        }
        expr_list <- lapply(data_layers, function(x) as.matrix(aobj@layers[[x]]))
        expr <- do.call(cbind, expr_list)
        colnames(expr) <- Cells(obj)
      } else {
        expr <- GetAssayData(obj, slot = "data")
      }
      expr
    }
    
    # ─────────── Dispersion Computation ───────────
    compute_dispersion <- function(obj, k = 30, ndims = 30) {
      if (!"pca" %in% names(obj@reductions)) {
        if (!"scale.data" %in% names(obj[[DefaultAssay(obj)]]@layers)) {
          obj <- ScaleData(obj, verbose = FALSE)
        }
        obj <- RunPCA(obj, npcs = ndims, verbose = FALSE)
      }
      pcs <- Embeddings(obj, "pca")[, 1:ndims, drop = FALSE]
      knn <- FNN::get.knn(pcs, k = k)
      disp <- sapply(seq_len(nrow(pcs)), function(i) {
        nbr <- knn$nn.index[i, ]
        centroid <- colMeans(pcs[nbr, , drop = FALSE])
        mean(rowSums((t(pcs[i, , drop = FALSE]) - centroid)^2))
      })
      list(values = disp, cells = rownames(pcs))
    }
    
    # ─────────── Run Dispersion ───────────
    observeEvent(input$run_dispersion, {
      output$dispersion_status <- renderText("⏳ Running dispersion analysis ...")
      
      tryCatch({
        obj <- get("global_dataset", envir = .GlobalEnv)
        req(obj)
        
        withProgress(message = "Dispersion Analysis", value = 0, {
          incProgress(0.2, detail = "Extracting expression matrix ...")
          expr <- get_expr_matrix(obj)
          expr_rv(expr)
          
          incProgress(0.6, detail = "Computing dispersion values ...")
          res <- compute_dispersion(obj)
          df <- data.frame(Cell = res$cells, Dispersion = res$values)
          df$Cluster <- Idents(obj)[df$Cell]
          
          high_cells <- df %>% dplyr::filter(Dispersion >= 2) %>% dplyr::pull(Cell)
          high_obj <- if (length(high_cells)) subset(obj, cells = high_cells) else obj
          high_obj <- FindVariableFeatures(high_obj, selection.method = "vst", nfeatures = 50)
          high_genes <- VariableFeatures(high_obj)
          if (length(high_genes) == 0) high_genes <- NA_character_
          
          df$TopGenes <- ifelse(df$Cell %in% high_cells,
                                paste(high_genes, collapse = ", "),
                                NA_character_)
          
          disp_df_rv(df)
          
          output$dispersion_status <- renderText(
            paste0("✅ Dispersion analysis completed successfully. Found ",
                   length(high_cells), " high-dispersion cells (Dᵢ ≥ 2).")
          )
          output$high_dispersion_summary <- renderUI({
            tags$p(style="font-weight:600; color:#2575fc;",
                   paste0("📊 High-dispersion cells: ", length(high_cells),
                          " | Top variable genes identified: ",
                          ifelse(all(is.na(high_genes)), 0, length(high_genes))))
          })
        })
      }, error = function(e) {
        output$dispersion_status <- renderText(paste("❌ Error:", e$message))
      })
    })
    
    # ─────────── BAR PLOT (Top 20 Cells × Top 5 Genes) ───────────
    observeEvent(input$plot_bar_dispersion, {
      df <- disp_df_rv()
      expr <- expr_rv()
      req(df, expr)
      
      obj <- get("global_dataset", envir = .GlobalEnv)
      pcs <- Embeddings(obj, "pca")[, 1:30, drop = FALSE]
      knn <- FNN::get.knn(pcs, k = 30)
      
      # Store knn in a temporary reactive value
      knn_rv <- reactiveVal(knn)
      
      # Select top 20 cells
      top20 <- df %>%
        arrange(desc(Dispersion)) %>%
        slice_head(n = 20)
      top20_idx <- match(top20$Cell, colnames(expr))
      
      # ---- Compute top 5 genes per cell (same as your working code) ----
      get_top_genes <- function(cell_idx, knn_obj, expr_matrix, top_n = 5) {
        if (is.na(cell_idx)) return(NA_character_)
        neighbors <- knn_obj$nn.index[cell_idx, ]
        local_expr <- expr_matrix[, neighbors, drop = FALSE]
        vars <- apply(local_expr, 1, var)
        top_genes <- names(sort(vars, decreasing = TRUE))[1:top_n]
        paste(top_genes, collapse = ", ")
      }
      
      # Make sure knn and expr are visible
      local_knn <- knn_rv()
      local_expr <- expr
      
      top20$Top5Genes <- sapply(top20_idx, function(i)
        get_top_genes(i, local_knn, local_expr, top_n = 5))
      
      # Clean and wrap labels
      top20$Top5Genes <- ifelse(is.na(top20$Top5Genes), "N/A",
                                stringr::str_wrap(top20$Top5Genes, width = 35))
      
      top20$CellShort <- gsub("^GSM[0-9]+_[^_]+_", "", top20$Cell)
      top20$CellShort <- gsub(".*_filtered_feature_bc_matrix_", "", top20$Cell)
      max_d <- max(top20$Dispersion, na.rm = TRUE)
      top20$label_y <- top20$Dispersion + 0.05 * max_d
      
      # ---- Plot ----
      p <- ggplot(top20, aes(x = reorder(CellShort, Dispersion),
                             y = Dispersion, fill = Dispersion)) +
        geom_col(width = 0.7) +
        geom_text(aes(y = label_y, label = Top5Genes),
                  hjust = 0, size = 3.2, color = "black", lineheight = 0.9) +
        scale_fill_gradient(low = "#a8caff", high = "#08306b") +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.4))) +
        coord_flip(clip = "off") +
        labs(title = "Top 20 Cells with Highest Dispersion",
             x = "Cell ID", y = "Dispersion Value", fill = "Dispersion") +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 9),
          plot.margin = margin(10, 20, 10, 10)
        )
      
      bar_plot_obj(p)
      output$bar_dispersion_plot <- renderPlot({ p }, height = 750)
    })
    
    # ─────────── UMAP Plot ───────────
    observeEvent(input$plot_umap_dispersion, {
      df <- disp_df_rv(); req(df)
      obj <- get("global_dataset", envir = .GlobalEnv)
      validate(need("umap" %in% names(obj@reductions),
                    "UMAP not found in dataset."))
      
      umap_df <- as.data.frame(Embeddings(obj, "umap"))
      colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
      obj$Dispersion <- df$Dispersion[match(colnames(obj), df$Cell)]
      umap_df$Dispersion <- obj$Dispersion
      
      p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Dispersion)) +
        geom_point(size = 0.5) +
        scale_color_gradient(low = "#a8caff", high = "#08306b") +
        labs(title = "UMAP Colored by Cell Dispersion (Dᵢ)", color = "Dispersion") +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
        )
      
      umap_plot_obj(p)
      output$umap_dispersion_plot <- renderPlot({ p }, height = 650)
    })
    
    # ─────────── Download Handlers ───────────
    output$download_all_dispersion <- downloadHandler(
      filename = function() "dispersion_all.csv",
      content = function(file) write.csv(disp_df_rv(), file, row.names = FALSE)
    )
    
    output$download_high_dispersion <- downloadHandler(
      filename = function() "dispersion_high.csv",
      content = function(file) {
        df <- disp_df_rv() %>% dplyr::filter(Dispersion >= 2)
        write.csv(df, file, row.names = FALSE)
      }
    )
    
    output$download_bar_plot <- downloadHandler(
      filename = function() "dispersion_barplot.png",
      content = function(file)
        ggsave(file, plot = bar_plot_obj(), width = 12, height = 7, dpi = 1000, bg = "white")
    )
    
    output$download_umap_plot <- downloadHandler(
      filename = function() "dispersion_umap.png",
      content = function(file)
        ggsave(file, plot = umap_plot_obj(), width = 10, height = 6, dpi = 1000, bg = "white")
    )
  })
}
