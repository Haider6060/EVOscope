# ─────────────────────────────────────────────
# 🧩 Clustering Module - SERVER LOGIC (universal + timer)
# ─────────────────────────────────────────────

clustering_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ───────────── Reactive values ─────────────
    clustering_in_progress <- reactiveVal(FALSE)
    start_time <- reactiveVal(NULL)
    
    # ───────────── Timer text ─────────────
    output$clustering_status <- renderText("")
    
    # ───────────── Reactive timer ─────────────
    autoInvalidate <- reactiveTimer(1000)  # refresh every second
    
    observe({
      autoInvalidate()
      if (clustering_in_progress()) {
        elapsed <- as.numeric(Sys.time() - start_time(), units = "secs")
        output$clustering_status <- renderText(
          paste0("⏳ Running clustering... (", round(elapsed), " sec elapsed)")
        )
      }
    })
    
    # ───────────── Run Clustering Button ─────────────
    observeEvent(input$run_clustering, {
      start_time(Sys.time())
      clustering_in_progress(TRUE)
      output$clustering_status <- renderText("⏳ Initializing clustering... Please wait.")
      
      tryCatch({
        # --- Step 1: Load Seurat object ---
        obj <- get("global_dataset", envir = .GlobalEnv)
        req(obj)
        
        # --- Step 2: Ensure preprocessing done ---
        message("🔧 Checking Seurat object structure...")
        obj <- prepare_seurat_object(obj)  # from global.R
        
        # --- Step 3: Detect assay type ---
        assay_type <- DefaultAssay(obj)
        message("✅ Using assay type: ", assay_type)
        
        # --- Step 4: Clustering process ---
        resolution <- input$resolution
        message("🚀 Running clustering at resolution = ", resolution)
        
        # Run neighborhood graph & clustering
        obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
        obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
        
        # --- Step 5: Run UMAP if missing ---
        if (!"umap" %in% names(obj@reductions)) {
          message("🧭 Running UMAP...")
          obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
        }
        
        # --- Step 6: Save globally for use in other modules ---
        assign("global_dataset", obj, envir = .GlobalEnv)
        
        # --- Step 7: Render UMAP plot ---
        output$umap_plot <- renderPlot({
          DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
            ggtitle(paste0("UMAP of Clusters (Assay: ", assay_type,
                           ", Resolution: ", resolution, ")")) +
            theme_minimal(base_size = 14) +
            theme(
              panel.grid = element_blank(),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA)
            )
        })
        
        # --- Step 8: Success message ---
        elapsed <- round(as.numeric(Sys.time() - start_time(), units = "secs"))
        output$clustering_status <- renderText(
          paste0("✅ Clustering completed successfully in ", elapsed, " seconds!")
        )
        clustering_in_progress(FALSE)
        
      }, error = function(e) {
        clustering_in_progress(FALSE)
        output$clustering_status <- renderText(paste0("❌ Error: ", e$message))
        message("❌ Error: ", e$message)
      })
    })
    
    # ───────────── Download UMAP PNG ─────────────
    output$download_umap <- downloadHandler(
      filename = function() { paste0("UMAP_Clusters_", Sys.Date(), ".png") },
      content = function(file) {
        obj <- get("global_dataset", envir = .GlobalEnv)
        req(obj)
        assay_type <- DefaultAssay(obj)
        res <- input$resolution
        
        p <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) +
          ggtitle(paste0("UMAP of Clusters (Assay: ", assay_type,
                         ", Resolution: ", res, ")")) +
          theme_minimal(base_size = 14) +
          theme(
            panel.grid = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
          )
        
        ggsave(file, plot = p, width = 8, height = 6, dpi = 1000, bg = "white")
      }
    )
  })
}
