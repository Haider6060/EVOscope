# ============================================================
# 📊 EVOscope EPS Integration — User Interface (Clean White Theme)
# ============================================================

eps_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        6,
        box(
          title = "⚖️ Upload Required Files",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          
          # ---- Upload Panels ----
          fileInput(ns("upload_entropy"), "Upload Entropy Scores CSV", accept = ".csv"),
          tableOutput(ns("preview_entropy")),
          tags$hr(),
          
          fileInput(ns("upload_dispersion"), "Upload Dispersion Scores CSV", accept = ".csv"),
          tableOutput(ns("preview_dispersion")),
          tags$hr(),
          
          fileInput(ns("upload_pathway"), "Upload Pathway Activity CSV", accept = ".csv"),
          tableOutput(ns("preview_pathway")),
          tags$hr(),
          
          # ---- Run EPS ----
          actionButton(ns("run_eps"), "Run EPS Integration", style = "background-color:white; color:black; border:1px solid #ccc;")
        )
      ),
      
      column(
        6,
        box(
          title = "📈 EPS Results Preview",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          
          # ---- Full EPS Results ----
          h4("Full EPS Results"),
          DT::dataTableOutput(ns("eps_table")),
          tags$br(),
          downloadButton(ns("download_eps"), "Download Full EPS Results",
                         style = "background-color:white; color:black; border:1px solid #ccc;"),
          tags$hr(),
          
          # ---- Top 30 ----
          h4("🔥 Top 30 Highly Active Cells"),
          DT::dataTableOutput(ns("top_eps_table")),
          tags$br(),
          downloadButton(ns("download_top_eps"), "Download Top 30 EPS Cells",
                         style = "background-color:white; color:black; border:1px solid #ccc;")
        )
      )
    )
  )
}
