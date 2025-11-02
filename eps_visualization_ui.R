# ============================================================
# 🎨 EPS Visualization UI (final – white buttons, all plots visible)
# ============================================================

eps_visualization_ui <- function(id) {
  ns <- NS(id)
  
  btn_white <- 'background:#fff; color:#333; border:1px solid #ccc;'
  
  tagList(
    h2("🎨 EPS Visualization Dashboard"),
    p("Upload EPS results to visualize UMAP projections and EPS relationships interactively. All buttons use a clean white style for consistency."),
    
    # -------------------------------
    # File Uploads & Actions
    # -------------------------------
    fluidRow(
      box(width = 6, solidHeader = TRUE,
          title = "📄 Upload Full EPS Results (EPS_scores.csv)",
          fileInput(ns("eps_file"), "Select EPS_scores.csv", accept = ".csv"),
          div(
            actionButton(ns("run_umap"), "Generate UMAP", class = "btn", style = btn_white),
            downloadButton(ns("download_umap"), "Download UMAP", class = "btn", style = btn_white),
            style = "display:flex; gap:10px; flex-wrap:wrap; margin-top:10px;"
          )
      ),
      box(width = 6, solidHeader = TRUE,
          title = "⭐ Upload Top 30 EPS Cells (Top30_EPS_cells.csv)",
          fileInput(ns("top_eps_file"), "Select Top30_EPS_cells.csv", accept = ".csv"),
          div(
            actionButton(ns("plot_top_eps"), "Generate Top EPS Combined Plot", class = "btn", style = btn_white),
            downloadButton(ns("download_combined"), "Download Combined Plot", class = "btn", style = btn_white),
            style = "display:flex; gap:10px; flex-wrap:wrap; margin-top:10px;"
          )
      )
    ),
    
    # -------------------------------
    # UMAP
    # -------------------------------
    fluidRow(
      box(width = 12, title = "🗺️ UMAP of EPS Metrics", solidHeader = TRUE,
          plotOutput(ns("umap_plot"), height = "500px")
      )
    ),
    
    # -------------------------------
    # EPS Relationships
    # -------------------------------
    fluidRow(
      box(width = 4, title = "Entropy vs EPS", solidHeader = TRUE,
          plotOutput(ns("entropy_vs_eps"), height = "360px"),
          div(downloadButton(ns("download_entropy"), "Download", class = "btn", style = btn_white),
              style = "margin-top:8px;")
      ),
      box(width = 4, title = "Dispersion vs EPS", solidHeader = TRUE,
          plotOutput(ns("dispersion_vs_eps"), height = "360px"),
          div(downloadButton(ns("download_dispersion"), "Download", class = "btn", style = btn_white),
              style = "margin-top:8px;")
      ),
      box(width = 4, title = "Pathway vs EPS", solidHeader = TRUE,
          plotOutput(ns("pathway_vs_eps"), height = "360px"),
          div(downloadButton(ns("download_pathway"), "Download", class = "btn", style = btn_white),
              style = "margin-top:8px;")
      )
    ),
    
    # -------------------------------
    # Combined EPS Plot (Top 30)
    # -------------------------------
    fluidRow(
      box(width = 12, title = "📊 Combined EPS Metrics (Top 30 Cells)",
          solidHeader = TRUE,
          plotOutput(ns("combined_plot"), height = "500px")
      )
    )
  )
}
