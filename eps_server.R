# ============================================================
# ⚖️ EVOscope EPS Integration Server — Final (Publication-Ready)
# Continuous EPS, Top 30 Export, Histogram Visualization
# ============================================================

eps_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(entropy = NULL, dispersion = NULL, pathway = NULL, eps_data = NULL)
    
    # --------------------------------------------------------
    # 🧩 Safe File Reader
    # --------------------------------------------------------
    safe_read <- function(file) {
      df <- suppressWarnings(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE))
      df[] <- lapply(df, function(x) if (is.character(x)) trimws(x) else x)
      df <- df[, colSums(!is.na(df) & df != "") > 0, drop = FALSE]
      
      # Identify ID column
      if (is.na(colnames(df)[1]) || colnames(df)[1] == "") colnames(df)[1] <- "Cell"
      
      # Assign rownames properly
      if (grepl("GSM|filtered_feature_bc_matrix", df[1, 1])) {
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
      } else if ("Cell" %in% colnames(df)) {
        rownames(df) <- df$Cell
        df$Cell <- NULL
      } else {
        rownames(df) <- make.unique(as.character(seq_len(nrow(df))))
      }
      
      df <- df[rowSums(is.na(df) | df == "") < ncol(df), , drop = FALSE]
      return(df)
    }
    
    # --------------------------------------------------------
    # 📁 Upload Handlers + Previews
    # --------------------------------------------------------
    observeEvent(input$upload_entropy, {
      rv$entropy <- safe_read(input$upload_entropy$datapath)
      showNotification("✅ Entropy file loaded.", type = "message")
      output$preview_entropy <- renderTable(head(data.frame(Cell = rownames(rv$entropy), rv$entropy), 5))
    })
    
    observeEvent(input$upload_dispersion, {
      rv$dispersion <- safe_read(input$upload_dispersion$datapath)
      showNotification("✅ Dispersion file loaded.", type = "message")
      output$preview_dispersion <- renderTable(head(data.frame(Cell = rownames(rv$dispersion), rv$dispersion), 5))
    })
    
    observeEvent(input$upload_pathway, {
      df <- safe_read(input$upload_pathway$datapath)
      first_id <- suppressWarnings(as.numeric(rownames(df)[1]))
      if (!is.na(first_id)) {
        if (!is.null(rv$entropy)) {
          n <- min(nrow(df), nrow(rv$entropy))
          rownames(df)[1:n] <- rownames(rv$entropy)[1:n]
        } else if (!is.null(rv$dispersion)) {
          n <- min(nrow(df), nrow(rv$dispersion))
          rownames(df)[1:n] <- rownames(rv$dispersion)[1:n]
        }
      }
      rv$pathway <- df
      showNotification("✅ Pathway file loaded successfully.", type = "message")
      output$preview_pathway <- renderTable(
        head(data.frame(Cell = rownames(rv$pathway),
                        rv$pathway[, 1:min(6, ncol(rv$pathway)), drop = FALSE]), 5)
      )
    })
    
    # --------------------------------------------------------
    # ⚙️ EPS Integration Logic
    # --------------------------------------------------------
    observeEvent(input$run_eps, {
      tryCatch({
        req(rv$entropy, rv$dispersion, rv$pathway)
        
        # Find overlap
        common <- Reduce(intersect, list(
          rownames(rv$entropy),
          rownames(rv$dispersion),
          rownames(rv$pathway)
        ))
        
        if (length(common) == 0)
          stop("❌ No overlapping Cell IDs found — ensure files share same Cell IDs.")
        
        # Extract metrics
        H <- as.numeric(rv$entropy[common, grep("Entropy", colnames(rv$entropy))[1]])
        D <- as.numeric(rv$dispersion[common, grep("Dispersion", colnames(rv$dispersion))[1]])
        P <- rowMeans(rv$pathway[common, , drop = FALSE], na.rm = TRUE)
        
        # Normalize 0–1
        normalize <- function(x) {
          r <- range(x, na.rm = TRUE)
          if (diff(r) == 0) return(rep(0.5, length(x)))
          (x - r[1]) / diff(r)
        }
        
        Hn <- normalize(H)
        Dn <- normalize(D)
        Pn <- normalize(P)
        
        EPS <- 0.33 * Hn + 0.33 * Dn + 0.34 * Pn
        cluster <- if ("Cluster" %in% colnames(rv$dispersion)) rv$dispersion[common, "Cluster"] else NA
        
        out <- data.frame(
          Cell = common,
          Cluster = cluster,
          Entropy = round(Hn, 4),
          Dispersion = round(Dn, 4),
          Pathway = round(Pn, 4),
          EPS = round(EPS, 4)
        )
        
        dir.create("outputs", showWarnings = FALSE)
        write.csv(out, "outputs/EPS_scores.csv", row.names = FALSE)
        rv$eps_data <- out
        
        # Histogram Visualization
        output$eps_plot <- renderPlot({
          hist(out$EPS, breaks = 50, col = "gray", border = "white",
               main = "EPS Distribution Across Cells", xlab = "EPS (Evolutionary Potential Score)")
          abline(v = quantile(out$EPS, 0.9), col = "red", lwd = 2, lty = 2)
          legend("topright", legend = "Top 10% threshold", col = "red", lwd = 2, lty = 2)
        })
        
        # Table Outputs
        output$eps_table <- DT::renderDataTable(out, options = list(pageLength = 10))
        top30 <- head(out[order(-out$EPS), ], 30)
        output$top_eps_table <- DT::renderDataTable(top30, options = list(pageLength = 10))
        
        showNotification("✅ EPS integration completed successfully.", type = "message")
      },
      error = function(e) {
        showNotification(paste("❌ EPS Integration Error:", e$message), type = "error")
        message("EPS Integration Error:", e$message)
      })
    })
    
    # --------------------------------------------------------
    # 💾 Download Handlers
    # --------------------------------------------------------
    output$download_eps <- downloadHandler(
      filename = function() "EPS_scores.csv",
      content = function(file) {
        req(rv$eps_data)
        write.csv(rv$eps_data, file, row.names = FALSE)
      }
    )
    
    output$download_top_eps <- downloadHandler(
      filename = function() "Top30_EPS_cells.csv",
      content = function(file) {
        req(rv$eps_data)
        top30 <- head(rv$eps_data[order(-rv$eps_data$EPS), ], 30)
        write.csv(top30, file, row.names = FALSE)
      }
    )
  })
}
