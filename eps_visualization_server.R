# ============================================================
# 🎨 EPS Visualization Server (colored vs EPS plots, final)
# ============================================================

eps_visualization_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      eps = NULL,
      top = NULL,
      umap_df = NULL,
      plot_umap = NULL,
      plot_entropy = NULL,
      plot_dispersion = NULL,
      plot_pathway = NULL,
      plot_combined = NULL
    )
    
    # -------------------------------
    # Uploads
    # -------------------------------
    observeEvent(input$eps_file, {
      req(input$eps_file)
      df <- read.csv(input$eps_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      needed <- c("Cell","Entropy","Dispersion","Pathway","EPS","Cluster")
      keep <- intersect(needed, colnames(df))
      if (length(keep) > 0) df <- df[, keep, drop = FALSE]
      rv$eps <- df
      showNotification("✅ EPS_scores.csv loaded.", type = "message")
      
      # --- Colored scatter plots ---
      if (all(c("Entropy","EPS") %in% colnames(df))) {
        rv$plot_entropy <- ggplot(df, aes(Entropy, EPS)) +
          geom_point(color = "#0072B2", alpha = 0.7, size = 1.8) +  # Blue
          geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
          theme_minimal(base_size = 13) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid = element_blank()) +
          labs(title = "Entropy vs EPS", x = "Entropy", y = "EPS")
      }
      
      if (all(c("Dispersion","EPS") %in% colnames(df))) {
        rv$plot_dispersion <- ggplot(df, aes(Dispersion, EPS)) +
          geom_point(color = "#009E73", alpha = 0.7, size = 1.8) +  # Green
          geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
          theme_minimal(base_size = 13) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid = element_blank()) +
          labs(title = "Dispersion vs EPS", x = "Dispersion", y = "EPS")
      }
      
      if (all(c("Pathway","EPS") %in% colnames(df))) {
        rv$plot_pathway <- ggplot(df, aes(Pathway, EPS)) +
          geom_point(color = "#D55E00", alpha = 0.7, size = 1.8) +  # Red
          geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
          theme_minimal(base_size = 13) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid = element_blank()) +
          labs(title = "Pathway vs EPS", x = "Pathway", y = "EPS")
      }
      
      output$entropy_vs_eps   <- renderPlot({ req(rv$plot_entropy);   rv$plot_entropy })
      output$dispersion_vs_eps<- renderPlot({ req(rv$plot_dispersion); rv$plot_dispersion })
      output$pathway_vs_eps   <- renderPlot({ req(rv$plot_pathway);   rv$plot_pathway })
    })
    
    observeEvent(input$top_eps_file, {
      req(input$top_eps_file)
      rv$top <- read.csv(input$top_eps_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      showNotification("✅ Top30_EPS_cells.csv loaded.", type = "message")
    })
    
    # -------------------------------
    # UMAP
    # -------------------------------
    observeEvent(input$run_umap, {
      req(rv$eps)
      df <- rv$eps
      if (!all(c("Entropy","Dispersion","Pathway","EPS") %in% colnames(df))) {
        showNotification("❌ EPS file must contain Entropy, Dispersion, Pathway, EPS columns.", type = "error")
        return(NULL)
      }
      X <- df[, c("Entropy","Dispersion","Pathway","EPS")]
      X <- as.data.frame(lapply(X, as.numeric))
      X <- X[complete.cases(X), ]
      if (nrow(X) < 10) {
        showNotification("❌ Not enough rows for UMAP.", type = "error")
        return(NULL)
      }
      
      set.seed(42)
      um <- uwot::umap(X, n_neighbors = 15, min_dist = 0.3)
      umap_df <- cbind(X, UMAP1 = um[,1], UMAP2 = um[,2])
      rv$umap_df <- umap_df
      
      rv$plot_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = EPS)) +
        geom_point(size = 1.6, alpha = 0.85) +
        viridis::scale_color_viridis(option = "plasma") +
        theme_minimal(base_size = 14) +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid = element_blank(),
              legend.position = "right") +
        labs(title = "UMAP of EPS Metrics", color = "EPS")
      
      output$umap_plot <- renderPlot({ rv$plot_umap })
    })
    
    # -------------------------------
    # Combined Top-30 EPS Plot
    # -------------------------------
    observeEvent(input$plot_top_eps, {
      req(rv$top)
      top <- rv$top
      req(all(c("Cell","Entropy","Dispersion","Pathway","EPS") %in% colnames(top)))
      
      shortify <- function(x) sub("^.*_", "", x)
      top$ShortCell <- shortify(top$Cell)
      long <- reshape2::melt(top[, c("ShortCell","Entropy","Dispersion","Pathway","EPS")],
                             id.vars = "ShortCell",
                             variable.name = "Metric", value.name = "Value")
      long$Metric <- factor(long$Metric, levels = c("Entropy","Dispersion","Pathway","EPS"))
      
      rv$plot_combined <- ggplot(long, aes(x = ShortCell, y = Value, fill = Metric)) +
        geom_col(position = "dodge", alpha = 0.95) +
        stat_summary(data = subset(long, Metric == "EPS"),
                     aes(group = 1, y = Value),
                     fun = mean, geom = "line", color = "black", linewidth = 0.7) +
        scale_fill_manual(values = c("Entropy"="#0072B2","Dispersion"="#009E73","Pathway"="#D55E00","EPS"="#000000")) +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              panel.background = element_rect(fill = "white"),
              panel.grid = element_blank(),
              legend.position = "right") +
        labs(x = "Cell", y = "Normalized Value", title = "Combined EPS Metrics (Top 30 Cells)")
      
      output$combined_plot <- renderPlot({ rv$plot_combined })
    })
    
    # -------------------------------
    # Downloads (no changes)
    # -------------------------------
    output$download_umap <- downloadHandler(
      filename = function() "UMAP_EPS_1000dpi.png",
      content = function(file) {
        req(rv$plot_umap)
        ggsave(file, plot = rv$plot_umap, dpi = 1000, width = 7, height = 6)
      }
    )
    
    output$download_entropy <- downloadHandler(
      filename = function() "Entropy_vs_EPS_1000dpi.png",
      content = function(file) {
        req(rv$plot_entropy)
        ggsave(file, plot = rv$plot_entropy, dpi = 1000, width = 7, height = 6)
      }
    )
    
    output$download_dispersion <- downloadHandler(
      filename = function() "Dispersion_vs_EPS_1000dpi.png",
      content = function(file) {
        req(rv$plot_dispersion)
        ggsave(file, plot = rv$plot_dispersion, dpi = 1000, width = 7, height = 6)
      }
    )
    
    output$download_pathway <- downloadHandler(
      filename = function() "Pathway_vs_EPS_1000dpi.png",
      content = function(file) {
        req(rv$plot_pathway)
        ggsave(file, plot = rv$plot_pathway, dpi = 1000, width = 7, height = 6)
      }
    )
    
    output$download_combined <- downloadHandler(
      filename = function() "Combined_EPS_1000dpi.png",
      content = function(file) {
        req(rv$plot_combined)
        ggsave(file, plot = rv$plot_combined, dpi = 1000, width = 10, height = 6)
      }
    )
  })
}
