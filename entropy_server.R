# ─────────────────────────────────────────────
# 📊 Local Entropy Calculation - FINAL STABLE UNIVERSAL VERSION
# Works for Seurat v3–v5, merged objects, and demo datasets
# ─────────────────────────────────────────────

entropy_server <- function(input, output, session) {
  
  entropy_vals  <- reactiveVal(NULL)
  entropy_table <- reactiveVal(NULL)
  
  # ───────────── Helper: universal expression extractor ─────────────
  get_expr_matrix <- function(obj) {
    assay_name <- DefaultAssay(obj)
    assay_obj  <- obj[[assay_name]]
    
    if ("layers" %in% slotNames(assay_obj)) {
      layer_names <- names(assay_obj@layers)
      norm_layers <- grep("^data", layer_names, value = TRUE)
      
      if (length(norm_layers) > 0) {
        message("Detected Seurat v5 layers: ", paste(layer_names, collapse = ", "))
        expr_list <- lapply(norm_layers, function(x) assay_obj@layers[[x]])
        expr <- do.call(cbind, expr_list)
        message("✅ Using normalized layer: ", norm_layers[1])
      } else if ("data" %in% layer_names) {
        expr <- assay_obj@layers[["data"]]
      } else if ("counts" %in% layer_names) {
        expr <- assay_obj@layers[["counts"]]
        message("⚠️ Using raw counts (no normalized layer detected).")
      } else stop("No usable layer found in assay.")
      
    } else {
      expr <- GetAssayData(obj, slot = "data")
    }
    
    # Fix missing gene names if necessary
    if (is.null(rownames(expr)) || length(rownames(expr)) == 0) {
      rownames(expr) <- rownames(obj[[DefaultAssay(obj)]])
    }
    return(expr)
  }
  
  # ───────────── Helper: Shannon entropy for one cell ─────────────
  shannon_entropy <- function(cell_idx, knn_obj, expr_matrix) {
    neighbors <- knn_obj$nn.index[cell_idx, ]
    local_expr <- expr_matrix[, neighbors, drop = FALSE]
    vars <- apply(local_expr, 1, var)
    p <- vars / sum(vars)
    p <- p[p > 0]
    H <- -sum(p * log2(p))
    return(H)
  }
  
  # ───────────── Main Entropy Calculation ─────────────
  observeEvent(input$run_entropy, {
    req(input$datafile)
    
    withProgress(message = "Running Entropy Calculation...", value = 0, {
      tryCatch({
        start_time <- Sys.time()
        
        # Step 1: Load Seurat object
        obj <- readRDS(input$datafile$datapath)
        incProgress(0.05, detail = "Seurat object loaded")
        
        # Step 2: Ensure normalization
        assay_obj <- obj[[DefaultAssay(obj)]]
        if ("layers" %in% slotNames(assay_obj)) {
          layer_names <- names(assay_obj@layers)
          has_data <- any(grepl("^data", layer_names))
        } else {
          has_data <- "data" %in% slotNames(assay_obj)
        }
        if (!has_data) {
          message("⚙️ Running NormalizeData() ...")
          obj <- NormalizeData(obj, verbose = FALSE)
        }
        incProgress(0.10, detail = "Normalization checked")
        
        # Step 3: Ensure variable features exist
        if (length(VariableFeatures(obj)) == 0) {
          message("⚙️ No variable features found — running FindVariableFeatures() ...")
          obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
          message("✅ Variable features computed.")
        }
        incProgress(0.20, detail = "Variable features verified")
        
        # Step 4: Extract expression matrix
        expr <- get_expr_matrix(obj)
        incProgress(0.30, detail = "Extracted normalized expression matrix")
        
        # Step 5: Ensure PCA exists
        if (!"pca" %in% names(obj@reductions)) {
          message("⚙️ PCA not found — running ScaleData() + RunPCA() ...")
          obj <- ScaleData(obj, verbose = FALSE)
          obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
          message("✅ PCA computed successfully.")
        }
        incProgress(0.40, detail = "PCA verified or computed")
        
        # Step 6: Variable genes for entropy
        gene_vars <- apply(expr, 1, var)
        top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:2000]
        expr_var <- as.matrix(expr[top_genes, ])
        expr_var <- log1p(expr_var)
        expr_var <- t(scale(t(expr_var), center = TRUE, scale = FALSE))
        incProgress(0.50, detail = "Selected top 2000 variable genes")
        
        # Step 7: Build KNN graph
        library(FNN)
        pca_coords <- Embeddings(obj, reduction = "pca")[, 1:30]
        knn <- FNN::get.knn(pca_coords, k = 30)
        incProgress(0.65, detail = "Constructed KNN graph (k = 30)")
        
        # Step 8: Compute entropy per cell
        n_cells <- ncol(expr_var)
        entropy_values <- numeric(n_cells)
        batch_size <- 500
        batches <- split(1:n_cells, ceiling(seq_along(1:n_cells) / batch_size))
        
        for (b in seq_along(batches)) {
          idx <- batches[[b]]
          entropy_values[idx] <- sapply(idx, shannon_entropy, 
                                        knn_obj = knn, expr_matrix = expr_var)
          incProgress(0.65 + (0.30 * b / length(batches)),
                      detail = paste("Computing batch", b, "of", length(batches)))
        }
        
        # Step 9: Save & expose results
        result_df <- data.frame(Cell = colnames(obj), Entropy = entropy_values)
        entropy_vals(entropy_values)
        entropy_table(result_df)
        
        elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2)
        incProgress(1, detail = paste("✅ Completed in", elapsed, "min"))
        output$entropy_status <- renderText(
          paste0("✅ Entropy computed successfully in ", elapsed, " minutes!")
        )
        
      }, error = function(e) {
        output$entropy_status <- renderText(paste0("❌ Error: ", e$message))
        message("Entropy error: ", e$message)
      })
    })
  })
  
  # ───────────── CSV Downloads ─────────────
  output$download_entropy_csv <- downloadHandler(
    filename = function() "entropy_scores.csv",
    content = function(file) {
      req(entropy_table())
      write.csv(entropy_table(), file, row.names = FALSE)
    }
  )
  
  output$download_top50_csv <- downloadHandler(
    filename = function() "top50_high_entropy_cells.csv",
    content = function(file) {
      req(entropy_table())
      top50 <- entropy_table() %>% 
        arrange(desc(Entropy)) %>% 
        head(50)
      write.csv(top50, file, row.names = FALSE)
    }
  )
}
