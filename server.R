# ============================================================
# 🧠 EVOscope Main Server (with EPS Visualization + Benchmark)
# ============================================================

server <- function(input, output, session) {
  
  dataset <- reactiveVal(NULL)
  
  # ============================================================
  # 📂 Load and Prepare Uploaded Seurat Object
  # ============================================================
  observeEvent(input$datafile, {
    req(input$datafile)
    file_path <- input$datafile$datapath
    
    tryCatch({
      obj <- readRDS(file_path)
      
      if (!inherits(obj, "Seurat")) {
        showNotification("❌ Uploaded file is not a Seurat object (.rds).", type="error", duration=6)
        return(NULL)
      }
      
      obj <- prepare_seurat_object(obj)
      dataset(obj)
      
      # make globally available for downstream modules
      assign("global_dataset", obj, envir = .GlobalEnv)
      
      showNotification("✅ Seurat object loaded successfully!", type = "message", duration = 5)
      
    }, error = function(e) {
      showNotification(paste("❌ Error loading file:", e$message), type = "error", duration = 6)
    })
  })
  
  # ============================================================
  # 📊 Dataset Summary Panel
  # ============================================================
  output$data_summary <- renderPrint({
    obj <- dataset()
    req(obj)
    cat("📁 File Name:", input$datafile$name, "\n")
    cat("🧬 Object Class:", class(obj), "\n")
    cat("🔹 Number of Cells:", ncol(obj), "\n")
    cat("🔹 Number of Genes:", nrow(obj), "\n")
    cat("⭐ Default Assay:", DefaultAssay(obj), "\n")
  })
  
  # ============================================================
  # ⚙️ Load All Module Servers
  # ============================================================
  clustering_server("cluster_module")
  entropy_server(input, output, session)
  dispersion_server("dispersion_module")
  pathway_server("pathway_module")
  eps_server("eps_module")
  eps_visualization_server("visual_module")
  benchmark_server("benchmark_module")
  
  # ============================================================
  # ✅ Status Messages
  # ============================================================
  cat("✅ EVOscope server initialized.\n")
  cat("📦 Modules loaded: clustering, entropy, dispersion, pathway, EPS, visualization, benchmarking.\n")
}
