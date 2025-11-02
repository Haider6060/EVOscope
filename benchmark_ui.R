# ============================================================
# 📈 EVOscope Benchmarking UI (Plain White Buttons + Task 3 Added)
# ============================================================

benchmark_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "main-card",
        h2("📈 Benchmarking and Validation"),
        p("This module compares EVOscope’s Evolutionary Potential Score (EPS) 
           with existing single-cell plasticity metrics such as CytoTRACE and scEntropy. 
           Upload your EPS output CSV and run the analysis to evaluate performance, 
           visualize correlations, and download publication-quality figures."),
        hr(),
        
        # --- Upload EPS CSV ---
        fileInput(ns("eps_csv"), "Upload EPS Scores CSV File",
                  accept = c(".csv"), buttonLabel = "Browse..."),
        br(),
        
        # --- Button Style: Plain White ---
        tags$style(HTML("
          .btn-plain {
            background-color: white !important;
            color: black !important;
            border: 1px solid #ccc !important;
            border-radius: 6px !important;
            font-weight: 500 !important;
            box-shadow: none !important;
          }
          .btn-plain:hover {
            background-color: #f9f9f9 !important;
            color: black !important;
          }
        ")),
        
        # --- Run Benchmark ---
        actionButton(ns("run_benchmark"), "🚀 Run Benchmark Analysis", class = "btn btn-plain"),
        br(), br(),
        
        # --- Progress Bar & Status ---
        uiOutput(ns("progress_ui")),
        verbatimTextOutput(ns("correlation_text")),
        hr(),
        
        # --- Plot 1: EPS vs scEntropy ---
        h3("🧠 EPS vs scEntropy Plot"),
        actionButton(ns("plot_entropy"), "Generate Plot", class = "btn btn-plain"),
        downloadButton(ns("download_entropy"), "Download Plot", class = "btn btn-plain"),
        br(), br(),
        plotOutput(ns("entropy_plot"), height = "450px"),
        hr(),
        
        # --- Plot 2: EPS vs CytoTRACE ---
        h3("🔬 EPS vs CytoTRACE Plot"),
        actionButton(ns("plot_cyto"), "Generate Plot", class = "btn btn-plain"),
        downloadButton(ns("download_cyto"), "Download Plot", class = "btn btn-plain"),
        br(), br(),
        plotOutput(ns("cyto_plot"), height = "450px"),
        hr(),
        
        # ============================================================
        # 🧬 Biological Meaning Check (Task 2)
        # ============================================================
        h2("🧬 Biological Meaning Check"),
        p("This step validates the biological significance of the Evolutionary Potential Score (EPS). 
           It identifies marker genes in high-EPS cells, performs pathway enrichment using Enrichr, 
           and visualizes enriched pathways related to proliferation, differentiation, and stress response."),
        hr(),
        
        # --- Run Pathway Analysis ---
        actionButton(ns("run_pathway"), "🧬 Run Pathway Analysis", class = "btn btn-plain"),
        br(), br(),
        uiOutput(ns("pathway_progress")),
        verbatimTextOutput(ns("pathway_message")),
        hr(),
        
        # --- Generate & Download Pathway Plot ---
        h3("📊 EPS Marker Genes and Enriched Pathways"),
        actionButton(ns("generate_pathway_plot"), "Generate Pathway Plot", class = "btn btn-plain"),
        downloadButton(ns("download_pathway_plot"), "Download Plot", class = "btn btn-plain"),
        br(), br(),
        plotOutput(ns("pathway_plot"), height = "550px"),
        hr(),
        
        # ============================================================
        # ♻️ Task 3: Reproducibility Check
        # ============================================================
        h2("♻️ EPS Reproducibility Check"),
        p("This analysis tests the reproducibility and robustness of EVOscope’s EPS values. 
           It randomly samples 70% of the cells, recalculates EPS, and compares the new values 
           to the original EPS. A high correlation (R > 0.9) confirms that EVOscope produces 
           stable and reproducible results across subsets."),
        hr(),
        
        # --- Run Reproducibility Analysis ---
        actionButton(ns("run_reproducibility"), "♻️ Run Reproducibility Analysis", class = "btn btn-plain"),
        br(), br(),
        uiOutput(ns("repro_progress")),
        verbatimTextOutput(ns("repro_message")),
        hr(),
        
        # --- Generate & Download Reproducibility Plot ---
        h3("📈 EPS Reproducibility Scatter Plot"),
        actionButton(ns("generate_repro_plot"), "Generate Plot", class = "btn btn-plain"),
        downloadButton(ns("download_repro_plot"), "Download Plot", class = "btn btn-plain"),
        br(), br(),
        plotOutput(ns("repro_plot"), height = "500px")
    )
  )
}
