# ─────────────────────────────────────────────
# 🧩 Clustering Module - UI
# ─────────────────────────────────────────────

clustering_ui <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    box(
      width = 12, solidHeader = TRUE, status = "primary",
      title = "🔍 Universal Clustering Module",
      
      p("This module performs clustering on your single-cell dataset using either 
        RNA or SCT assay automatically. You can adjust the clustering resolution 
        and visualize the resulting UMAP. The process includes PCA (if missing), 
        KNN graph construction, and Louvain clustering."),
      
      br(),
      # Resolution control
      sliderInput(
        ns("resolution"),
        label = strong("Select Clustering Resolution:"),
        min = 0.1, max = 1.5, value = 0.5, step = 0.1
      ),
      
      # Run button
      actionButton(ns("run_clustering"), "🚀 Run Clustering", class = "btn-primary"),
      br(), br(),
      
      # Status message
      textOutput(ns("clustering_status")),
      br(),
      
      # UMAP plot display
      plotOutput(ns("umap_plot"), height = "600px"),
      br(),
      
      # Download button
      downloadButton(ns("download_umap"), "📥 Download UMAP")
    )
  )
}
