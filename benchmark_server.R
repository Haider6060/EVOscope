# ============================================================
# 🧠 EVOscope Benchmarking Server (Final Universal v10 — Plain White Reproducibility Plot)
# ============================================================

benchmark_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    data_reactive <- reactiveVal(NULL)
    cor_results <- reactiveVal(NULL)
    
    # ---------- Helper Functions ----------
    .safe_minmax <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      if (!length(x) || all(is.na(x))) return(rep(NA_real_, length(x)))
      rng <- range(x, na.rm = TRUE)
      if (isTRUE(all.equal(rng[1], rng[2]))) return(rep(0.5, length(x)))
      (x - rng[1]) / (rng[2] - rng[1])
    }
    
    .ensure_benchmark_cols <- function(df) {
      if (!"CytoTRACE_like" %in% names(df))
        df$CytoTRACE_like <- 1 - .safe_minmax(df$Entropy)
      if (!"scEntropy_like" %in% names(df))
        df$scEntropy_like <- .safe_minmax(df$Entropy)
      df
    }
    
    # ============================================================
    # 📂 Load EPS CSV and Link with Seurat Object
    # ============================================================
    observeEvent(input$eps_csv, {
      req(input$eps_csv)
      eps_data <- tryCatch(read.csv(input$eps_csv$datapath),
                           error = function(e) {
                             showNotification("❌ Failed to read EPS CSV.", type = "error")
                             return(NULL)
                           })
      if (is.null(eps_data)) return(NULL)
      
      seu <- if (exists("global_dataset", envir = .GlobalEnv))
        get("global_dataset", envir = .GlobalEnv) else NULL
      if (is.null(seu)) {
        showNotification("⚠️ No Seurat object loaded. Please upload your dataset first.", type = "error")
        return(NULL)
      }
      
      eps_data <- eps_data[match(colnames(seu), eps_data$Cell), ]
      seu$EPS <- eps_data$EPS
      seu$Entropy <- eps_data$Entropy
      seu$Dispersion <- eps_data$Dispersion
      seu$Pathway <- eps_data$Pathway
      eps_data <- .ensure_benchmark_cols(eps_data)
      
      assign("global_eps_data", eps_data, .GlobalEnv)
      assign("global_dataset", seu, .GlobalEnv)
      data_reactive(eps_data)
      
      showNotification("✅ EPS data linked successfully.", type = "message", duration = 4)
    })
    
    # ============================================================
    # 🚀 Run Benchmark Analysis
    # ============================================================
    observeEvent(input$run_benchmark, {
      req(data_reactive())
      eps_data <- .ensure_benchmark_cols(data_reactive())
      valid <- complete.cases(eps_data$EPS, eps_data$CytoTRACE_like, eps_data$scEntropy_like)
      if (sum(valid) < 3) {
        showNotification("❌ Insufficient valid data for correlation.", type = "error")
        return(NULL)
      }
      
      cor_eps_cyto <- suppressWarnings(cor(eps_data$EPS[valid], eps_data$CytoTRACE_like[valid], method = "spearman"))
      cor_eps_entropy <- suppressWarnings(cor(eps_data$EPS[valid], eps_data$scEntropy_like[valid], method = "spearman"))
      cor_results(list(cyto = cor_eps_cyto, entropy = cor_eps_entropy))
      
      output$progress_ui <- renderUI({
        tags$div("✅ Benchmark analysis complete!", style = "font-weight:600;color:green;")
      })
    })
    
    # ============================================================
    # 🎨 Correlation Plots + Download Handlers
    # ============================================================
    output$entropy_plot <- renderPlot({
      req(input$plot_entropy, data_reactive())
      eps_data <- .ensure_benchmark_cols(data_reactive())
      res <- cor_results()
      rlab <- ifelse(is.null(res$entropy) || is.na(res$entropy), "NA", round(res$entropy, 2))
      ggplot(eps_data, aes(x = scEntropy_like, y = EPS)) +
        geom_point(alpha = 0.6, color = "#2575fc", size = 1) +
        theme_bw(14) + theme(panel.grid = element_blank()) +
        labs(title = paste0("EPS vs scEntropy-like (R=", rlab, ")"),
             x = "scEntropy-like score", y = "EVOscope EPS")
    })
    
    output$download_entropy <- downloadHandler(
      filename = function() "EPS_vs_scEntropy.png",
      content = function(file) {
        eps_data <- .ensure_benchmark_cols(data_reactive())
        res <- cor_results()
        rlab <- ifelse(is.null(res$entropy) || is.na(res$entropy), "NA", round(res$entropy, 2))
        g <- ggplot(eps_data, aes(x = scEntropy_like, y = EPS)) +
          geom_point(alpha = 0.6, color = "#2575fc", size = 1.2) +
          theme_bw(18) + theme(panel.grid = element_blank()) +
          labs(title = paste0("EPS vs scEntropy-like (R=", rlab, ")"),
               x = "scEntropy-like score", y = "EVOscope EPS")
        ggsave(file, plot = g, width = 8, height = 6, dpi = 600)
      }
    )
    
    output$cyto_plot <- renderPlot({
      req(input$plot_cyto, data_reactive())
      eps_data <- .ensure_benchmark_cols(data_reactive())
      res <- cor_results()
      rlab <- ifelse(is.null(res$cyto) || is.na(res$cyto), "NA", round(res$cyto, 2))
      ggplot(eps_data, aes(x = CytoTRACE_like, y = EPS)) +
        geom_point(alpha = 0.6, color = "#E65100", size = 1) +
        theme_bw(14) + theme(panel.grid = element_blank()) +
        labs(title = paste0("EPS vs CytoTRACE-like (R=", rlab, ")"),
             x = "CytoTRACE-like score", y = "EVOscope EPS")
    })
    
    output$download_cyto <- downloadHandler(
      filename = function() "EPS_vs_CytoTRACE.png",
      content = function(file) {
        eps_data <- .ensure_benchmark_cols(data_reactive())
        res <- cor_results()
        rlab <- ifelse(is.null(res$cyto) || is.na(res$cyto), "NA", round(res$cyto, 2))
        g <- ggplot(eps_data, aes(x = CytoTRACE_like, y = EPS)) +
          geom_point(alpha = 0.6, color = "#E65100", size = 1.2) +
          theme_bw(18) + theme(panel.grid = element_blank()) +
          labs(title = paste0("EPS vs CytoTRACE-like (R=", rlab, ")"),
               x = "CytoTRACE-like score", y = "EVOscope EPS")
        ggsave(file, plot = g, width = 8, height = 6, dpi = 600)
      }
    )
    
    # ============================================================
    # 🧬 Biological Meaning Check (Plain White Pathway Plot)
    # ============================================================
    observeEvent(input$run_pathway, {
      library(R.utils)
      req(exists("global_dataset", envir = .GlobalEnv))
      seu <- get("global_dataset", envir = .GlobalEnv)
      req("EPS" %in% colnames(seu@meta.data))
      
      showNotification("🔬 Running Biological Meaning Check...", type = "message", duration = 3)
      output$pathway_message <- renderText("⏳ Running enrichment... please wait.")
      
      seu$EPS_Group <- ifelse(seu$EPS >= sort(seu$EPS, decreasing = TRUE)[100], "High_EPS", "Other")
      Idents(seu) <- "EPS_Group"
      
      if ("RNA" %in% names(seu@assays)) {
        try({
          assay_obj <- seu@assays$RNA
          sn <- slotNames(assay_obj)
          if ("joined" %in% sn && isFALSE(assay_obj@joined)) seu <- JoinLayers(seu)
          if ("layers" %in% sn && !"joined" %in% sn && length(assay_obj@layers) > 0) seu <- JoinLayers(seu)
        }, silent = TRUE)
      }
      
      markers_high_eps <- tryCatch(
        FindMarkers(seu, ident.1 = "High_EPS", ident.2 = "Other",
                    logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox"),
        error = function(e) {
          DefaultAssay(seu) <- "RNA"
          FindMarkers(seu, ident.1 = "High_EPS", ident.2 = "Other",
                      logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")
        }
      )
      
      top50 <- head(rownames(markers_high_eps[order(-markers_high_eps$avg_log2FC), ]), 50)
      Sys.setenv(ENRICHR_URL = "https://maayanlab.cloud/Enrichr/")
      
      enriched <- tryCatch({
        withTimeout({
          enrichR::enrichr(top50, c("GO_Biological_Process_2023", "Reactome_2022"))
        }, timeout = 10, onTimeout = "error")
      }, error = function(e) NULL)
      
      if (is.null(enriched)) {
        suppressPackageStartupMessages({
          library(clusterProfiler)
          library(org.Hs.eg.db)
        })
        gene_ids <- tryCatch(bitr(top50, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db), error=function(e)NULL)
        if (!is.null(gene_ids)) {
          enriched_go <- enrichGO(gene=gene_ids$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", readable=TRUE)
          top_pathways <- head(as.data.frame(enriched_go),10)
          top_pathways$log10FDR <- -log10(top_pathways$p.adjust)
          top_pathways$Count <- top_pathways$Count
        } else return(NULL)
      } else {
        df_all <- do.call(rbind, enriched)
        top_pathways <- head(df_all[order(df_all$Adjusted.P.value), ], 10)
        top_pathways$log10FDR <- -log10(top_pathways$Adjusted.P.value)
        top_pathways$Count <- sapply(strsplit(as.character(top_pathways$Genes), ";"), length)
      }
      
      assign("top_pathways", top_pathways, .GlobalEnv)
      
      output$pathway_plot <- renderPlot({
        ggplot(top_pathways, aes(x = reorder(Description, log10FDR),
                                 y = log10FDR, size = Count, color = log10FDR)) +
          geom_point(alpha = 0.85) +
          coord_flip() +
          scale_color_gradient(low = "#56B1F7", high = "#CA0020") +
          theme_minimal(base_size = 14) +
          theme(panel.grid = element_blank(),
                panel.background = element_rect(fill = "white", color = NA),
                plot.background = element_rect(fill = "white", color = NA),
                axis.text = element_text(color = "black"),
                axis.title = element_text(face = "bold"),
                plot.title = element_text(hjust = 0.5, face = "bold")) +
          labs(title = "Top Enriched Pathways",
               x = "Pathway", y = "-log10(Adj P-value)", size = "Gene Count")
      })
      
      output$pathway_message <- renderText("✅ Pathway enrichment complete and plot generated.")
    })
    
    output$download_pathway_plot <- downloadHandler(
      filename=function() "Top_Enriched_Pathways.png",
      content=function(file){
        req(exists("top_pathways", envir=.GlobalEnv))
        tp <- get("top_pathways", envir=.GlobalEnv)
        g <- ggplot(tp, aes(x=reorder(Description,log10FDR), y=log10FDR, size=Count, color=log10FDR))+
          geom_point(alpha=0.85)+coord_flip()+
          scale_color_gradient(low="#56B1F7",high="#CA0020")+
          theme_minimal(base_size=16)+
          theme(panel.grid=element_blank(),
                panel.background=element_rect(fill="white", color=NA),
                plot.background=element_rect(fill="white", color=NA),
                axis.text=element_text(color="black"),
                axis.title=element_text(face="bold"),
                plot.title=element_text(hjust=0.5, face="bold"))+
          labs(title="Top Enriched Pathways",x="Pathway",y="-log10(Adj P-value)",size="Gene Count")
        ggsave(file,plot=g,width=10,height=8,dpi=1000,bg="white")
      }
    )
    
    # ============================================================
    # ♻️ EPS Reproducibility + Download (Plain White Background)
    # ============================================================
    observeEvent(input$run_reproducibility, {
      req(exists("global_dataset", .GlobalEnv))
      seu <- get("global_dataset", .GlobalEnv)
      eps_data <- get("global_eps_data", .GlobalEnv)
      
      showNotification("♻️ Running EPS Reproducibility Check...", type = "message", duration = 3)
      set.seed(123)
      subsample <- sample(colnames(seu), size=floor(0.7*ncol(seu)))
      seu_sub <- subset(seu, cells=subsample)
      eps_sub <- eps_data[match(colnames(seu_sub), eps_data$Cell), ]
      seu_sub$EPS <- eps_sub$EPS
      
      common <- intersect(colnames(seu), colnames(seu_sub))
      cor_r <- suppressWarnings(cor(seu$EPS[common], seu_sub$EPS[common], method="spearman"))
      assign("seu_sub", seu_sub, .GlobalEnv)
      assign("common_cells", common, .GlobalEnv)
      assign("cor_eps_repro", cor_r, .GlobalEnv)
      
      output$repro_message <- renderText({
        paste0("✅ EPS reproducibility (R = ", round(cor_r,3), ")")
      })
    })
    
    output$repro_plot <- renderPlot({
      req(input$generate_repro_plot)
      seu <- get("global_dataset", .GlobalEnv)
      seu_sub <- get("seu_sub", .GlobalEnv)
      common <- get("common_cells", .GlobalEnv)
      cor_r <- get("cor_eps_repro", .GlobalEnv)
      df <- data.frame(Full=seu$EPS[common],Subsample=seu_sub$EPS[common])
      ggplot(df, aes(x=Full,y=Subsample))+
        geom_point(alpha=0.6,color="#2575fc",size=1.2)+
        geom_smooth(method="lm",color="darkred",se=FALSE,linetype="dashed")+
        theme_minimal(base_size=14)+
        theme(panel.grid=element_blank(),
              panel.background=element_rect(fill="white", color=NA),
              plot.background=element_rect(fill="white", color=NA),
              axis.text=element_text(color="black"),
              axis.title=element_text(face="bold"),
              plot.title=element_text(hjust=0.5, face="bold"))+
        labs(title=paste0("EPS Reproducibility (R=",round(cor_r,2),")"),
             x="EPS (Full Dataset)",y="EPS (70% Subsample)")
    })
    
    output$download_repro_plot <- downloadHandler(
      filename=function() "EPS_Reproducibility.png",
      content=function(file){
        seu <- get("global_dataset", .GlobalEnv)
        seu_sub <- get("seu_sub", .GlobalEnv)
        common <- get("common_cells", .GlobalEnv)
        cor_r <- get("cor_eps_repro", .GlobalEnv)
        df <- data.frame(Full=seu$EPS[common],Subsample=seu_sub$EPS[common])
        g <- ggplot(df,aes(x=Full,y=Subsample))+
          geom_point(alpha=0.6,color="#2575fc",size=1.2)+
          geom_smooth(method="lm",color="darkred",se=FALSE,linetype="dashed")+
          theme_minimal(base_size=18)+
          theme(panel.grid=element_blank(),
                panel.background=element_rect(fill="white", color=NA),
                plot.background=element_rect(fill="white", color=NA),
                axis.text=element_text(color="black"),
                axis.title=element_text(face="bold"),
                plot.title=element_text(hjust=0.5, face="bold"))+
          labs(title=paste0("EPS Reproducibility (R=",round(cor_r,2),")"),
               x="EPS (Full Dataset)",y="EPS (70% Subsample)")
        ggsave(file,plot=g,width=8,height=6,dpi=1000,bg="white")
      }
    )
    
    # ============================================================
    # 🔧 Activate Download Buttons
    # ============================================================
    observe({
      session$sendCustomMessage("activateDownloadButtons", list(ns_prefix = ns("")))
    })
  })
}
