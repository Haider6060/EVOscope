# ============================================================
# 🧬 EVOscope Pathway Diversity - UI (Refined Final)
# ============================================================

pathway_ui <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    box(
      width = 12, solidHeader = TRUE, status = "primary",
      title = "🧬 Pathway Diversity (Hallmark AUCell Scoring)",
      
      # ---- Description ----
      p("This module computes per-cell and per-cluster pathway activities using ",
        strong("AUCell scoring"),
        " against the MSigDB Hallmark gene sets (loaded from GMT).",
        " It helps identify the most transcriptionally diverse and functionally active clusters."
      ),
      tags$hr(),
      
      # ---- Run Button ----
      h4("⚙️ Run Analysis"),
      p("Click below to start the pathway activity scoring. This step may take a few minutes depending on dataset size."),
      actionButton(ns("run_pathway"), "▶ Run Pathway Analysis", class = "btn btn-default"),
      br(), br(),
      textOutput(ns("status_pathway")),
      uiOutput(ns("summary_pathway")),
      tags$hr(),
      
      # ---- Download Section ----
      h4("💾 Download Results"),
      fluidRow(
        column(3, downloadButton(ns("dl_full_scores"), "All Pathway Activity (All Cells)", class = "btn btn-default")),
        column(3, downloadButton(ns("dl_cluster_means"), "Pathway Activity by Cluster", class = "btn btn-default")),
        column(3, downloadButton(ns("dl_top30_overall"), "Top 30 Genes (Overall)", class = "btn btn-default")),
        column(3, downloadButton(ns("dl_top10_top30"), "Top 10 Pathways × Top 30 Genes", class = "btn btn-default"))
      ),
      tags$hr(),
      
      # ---- Visualization Section ----
      h4("📊 Visualization"),
      p("Visualize the most variable pathways and top genes across clusters."),
      fluidRow(
        column(4, actionButton(ns("plot_heatmap"), "🔥 Generate Heatmap (Top Pathways × Clusters)", class = "btn btn-default")),
        column(4, actionButton(ns("plot_violin"), "🎻 Generate Violin Plot (Pathway Variability)", class = "btn btn-default")),
        column(4, actionButton(ns("plot_bar"), "📊 Generate Bar Plot (Top 10 Pathways × Top 5 Genes)", class = "btn btn-default"))
      ),
      
      tags$br(),
      tabsetPanel(
        tabPanel("🔥 Heatmap", 
                 plotOutput(ns("heatmap"), height = "600px"), 
                 downloadButton(ns("dl_heatmap_png"), "⬇️ Download Heatmap", class = "btn btn-default")),
        
        tabPanel("🎻 Violin Plot", 
                 plotOutput(ns("violin"), height = "600px"), 
                 downloadButton(ns("dl_violin_png"), "⬇️ Download Violin", class = "btn btn-default")),
        
        tabPanel("📊 Bar Plot", 
                 plotOutput(ns("barplot"), height = "650px"), 
                 downloadButton(ns("dl_bar_png"), "⬇️ Download Bar Plot", class = "btn btn-default"))
      )
    )
  )
}
