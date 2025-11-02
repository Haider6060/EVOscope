# ─────────────────────────────────────────────
# 📈 Dispersion Module - UI (Clean Minimal Style)
# ─────────────────────────────────────────────
dispersion_ui <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    box(
      width = 12, solidHeader = TRUE, status = "primary",
      title = "📈 Dispersion Analysis Module",
      
      p("This module quantifies **cell-level dispersion (Dᵢ)** — the degree of deviation 
         of each cell from its local neighborhood centroid in PCA space. 
         High dispersion values indicate locally unstable or heterogeneous cells. 
         A cutoff of **Dᵢ > 2** is used to identify highly dispersed cells."),
      br(),
      
      # Run button
      actionButton(ns("run_dispersion"), "🚀 Run Dispersion Analysis"),
      br(), br(),
      
      # Status text (shows progress and final summary)
      textOutput(ns("dispersion_status")),
      br(),
      
      # Summary of high-dispersion cells (persistent)
      uiOutput(ns("high_dispersion_summary")),
      br(),
      
      # Download and plot buttons
      fluidRow(
        column(3, downloadButton(ns("download_all_dispersion"), "📥 All Cell Values (CSV)")),
        column(3, downloadButton(ns("download_high_dispersion"), "📥 High-Dispersion Cells (CSV)")),
        column(3, actionButton(ns("plot_bar_dispersion"), "📊 Generate Bar Plot")),
        column(3, actionButton(ns("plot_umap_dispersion"), "🧬 Generate UMAP Plot"))
      ),
      br(),
      
      # Plots
      tabsetPanel(
        tabPanel("📊 Bar Plot (Top 20 Cells × Top 5 Genes)",
                 plotOutput(ns("bar_dispersion_plot"), height = "600px"),
                 downloadButton(ns("download_bar_plot"), "📥 Download Bar Plot")
        ),
        tabPanel("🧬 UMAP Colored by Dispersion",
                 plotOutput(ns("umap_dispersion_plot"), height = "600px"),
                 downloadButton(ns("download_umap_plot"), "📥 Download UMAP")
        )
      )
    )
  )
}
