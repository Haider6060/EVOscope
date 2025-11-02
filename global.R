# ============================================================
# 🌐 EVOscope Global Configuration (Updated with Benchmark)
# ============================================================

options(shiny.maxRequestSize = 2000 * 1024^2)   # 2 GB upload limit

# ============================================================
# 📦 Load Core Libraries
# ============================================================
library(shiny)
library(shinydashboard)
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(FNN)
library(scales)
library(gridExtra)
library(stringr)
library(GSEABase)
library(AUCell)
library(pheatmap)
library(reshape2)
library(DT)
library(viridis)
library(shinyWidgets)
library(shinyjs)
# ============================================================
# 💾 Increase R memory limit (Windows only)
# ============================================================
if (.Platform$OS.type == "windows") {
  memory.limit(size = 32000)
}

# ============================================================
# 🔧 Helper: Safe Seurat Preprocessing
# ============================================================
prepare_seurat_object <- function(obj) {
  assay_type <- DefaultAssay(obj)
  
  if (!(assay_type %in% c("RNA", "SCT"))) {
    message("⚠️ Unknown assay type detected. Switching to RNA.")
    DefaultAssay(obj) <- "RNA"
    assay_type <- "RNA"
  }
  
  if (!"data" %in% names(obj[[assay_type]]@layers)) {
    message("🧪 Running NormalizeData() ...")
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  
  if (length(VariableFeatures(obj)) == 0) {
    message("🔎 Running FindVariableFeatures() ...")
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  }
  
  if (!"scale.data" %in% names(obj[[assay_type]]@layers)) {
    message("📏 Running ScaleData() ...")
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  }
  
  if (!"pca" %in% names(obj@reductions)) {
    message("📉 Running PCA() ...")
    obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
  }
  
  message("✅ Seurat object prepared successfully.")
  return(obj)
}

# ============================================================
# 🧩 Extract Expression Matrix (Universal v3–v5)
# ============================================================
get_expr_matrix <- function(obj) {
  assay_name <- DefaultAssay(obj)
  assay_obj  <- obj[[assay_name]]
  
  if ("layers" %in% slotNames(assay_obj)) {
    data_layers <- grep("^data", names(assay_obj@layers), value = TRUE)
    if (length(data_layers) == 0) {
      message("⚠️ No 'data' layer found → running NormalizeData()...")
      obj <- NormalizeData(obj, verbose = FALSE)
      expr <- GetAssayData(obj, slot = "data")
    } else {
      message("✅ Using Seurat v5 normalized layers.")
      expr_list <- lapply(data_layers, function(l) as.matrix(assay_obj@layers[[l]]))
      expr <- do.call(cbind, expr_list)
      colnames(expr) <- Cells(obj)
    }
  } else {
    expr <- GetAssayData(obj, slot = "data")
  }
  
  if (is.null(rownames(expr))) rownames(expr) <- rownames(obj[[assay_name]])
  if (is.null(colnames(expr))) colnames(expr) <- Cells(obj)
  
  return(expr)
}

# ============================================================
# 🧬 Load Hallmark Pathways
# ============================================================
gmt_path <- "D:/NEW Tool for 4th Project/Tool/h.all.v2023.1.Hs.symbols.gmt"

if (file.exists(gmt_path)) {
  gene_sets <- getGmt(gmt_path)
  message("✅ Loaded ", length(gene_sets), " hallmark gene sets.")
} else {
  warning("⚠️ GMT file not found: ", gmt_path)
  gene_sets <- NULL
}

# ============================================================
# 🔗 Load All Modules
# ============================================================
source("entropy_ui.R");      source("entropy_server.R")
source("clustering_ui.R");   source("clustering_server.R")
source("dispersion_ui.R");   source("dispersion_server.R")
source("pathway_ui.R");      source("pathway_server.R")
source("eps_ui.R");          source("eps_server.R")
source("eps_visualization_ui.R"); source("eps_visualization_server.R")
source("benchmark_ui.R");    source("benchmark_server.R")

# ============================================================
# ✅ Startup Message
# ============================================================
cat("✅ EVOscope environment initialized.\n")
cat("📦 Loaded modules: entropy, clustering, dispersion, pathway, EPS, visualization, benchmarking\n")
cat("⚙️ File upload limit: 2 GB\n")
