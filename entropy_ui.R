# ────────────────────────────────
# 📊 Local Entropy Calculation - UI (Clean Minimal Style)
# ────────────────────────────────

entropy_ui <- function() {
  fluidPage(
    h2("📊 Local Entropy Calculation"),
    p("This module computes the Shannon entropy for each cell using the uploaded Seurat object 
       from the Data Upload step. It measures transcriptional variability and cellular plasticity 
       within the local KNN neighborhood."),
    tags$hr(),
    
    h4("⚙️ Run Analysis"),
    p("Click the button below to start the entropy computation. 
       This may take several minutes depending on dataset size."),
    
    # Run button
    actionButton("run_entropy", "▶ Run Entropy Calculation",
                 style = "margin-bottom:10px;"),
    
    verbatimTextOutput("entropy_status"),
    tags$hr(),
    
    h4("📥 Download Results"),
    p("Once computation is complete, you can download the full entropy results 
       or just the top 50 cells with highest entropy scores."),
    
    # Download buttons (default neutral style)
    downloadButton("download_entropy_csv", "⬇️ Download Full Entropy CSV",
                   style = "margin-right:10px;"),
    downloadButton("download_top50_csv", "⬇️ Download Top 50 High-Entropy Cells"),
    
    tags$br(), tags$br(),
    
    tags$div(style = "color:#666; font-size:13px;",
             "Note: This step uses the dataset uploaded in the Data Upload section.
              Please ensure it is preprocessed and contains PCA embeddings.")
  )
}
