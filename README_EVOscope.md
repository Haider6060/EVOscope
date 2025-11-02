# EVOscope: Single-Cell Evolutionary Potential Analysis

## Overview

**EVOscope** is a user‑friendly **Shiny application** and programmatic toolkit for **single‑cell RNA‑seq** (and optional multiome RNA+ATAC) analysis.  
The tool introduces the **Evolutionary Potential Score (EPS)** — a scalar that integrates three orthogonal components to quantify **cellular plasticity** and **adaptive potential**:

- **CellEntropy (H):** local transcriptional disorder (Shannon entropy on neighborhood variance)  
- **Dispersion (D):** positional deviation from the local manifold centroid (PCA/UMAP KNN)  
- **Pathway Diversity (P):** breadth of functional activation (AUCell on MSigDB Hallmark sets)

EVOscope is **dataset‑agnostic** and runs on **any GEO or custom single‑cell dataset**, whether **RNA‑only** or **multi‑omic**, and whether your Seurat object has **layers or not** (the app auto‑detects normalization, variable features, PCA, etc.).  fileciteturn3file0

---

## ✅ Key Features

- **Automated data integrity checks** (Seurat v3–v5): normalize, find variable genes, PCA/UMAP, KNN, clustering  
- **EPS computation (H, D, P → EPS)** with per‑cell and per‑cluster summaries  
- **Interactive visualizations:** UMAP by EPS, entropy/dispersion maps, pathway heatmaps, top‑cells tables  
- **Benchmarking module:** compare EPS with CytoTRACE‑like and scEntropy‑like scores; subsampling reproducibility  
- **Prognostic module (optional):** derive EPS‑associated signatures and run **Cox**, **Kaplan–Meier**, and **time‑dependent ROC** on external cohorts  
- **High‑resolution exports** (PNG/PDF) and **CSV outputs** for downstream use

---

## 📂 Input Requirements

- `.rds` file (Seurat object, v3–v5). RNA‑only or RNA+ATAC are supported.  
- If layers/slots are missing, EVOscope **fills them automatically** (NormalizeData → FindVariableFeatures → ScaleData → RunPCA/UMAP).

### ➤ Preparing Your Data

```r
# Example (from a gene-by-cell matrix):
library(Seurat)
seurat_obj <- CreateSeuratObject(counts = your_matrix)
saveRDS(seurat_obj, file = "your_dataset.rds")
```

---

## 🧪 Datasets Tested (examples)

EVOscope has been successfully tested on multiple real‑world datasets (no accession dependency):

| Dataset Type | Source | Description |
|--------------|--------|-------------|
| Lung cancer | GEO / in‑house | Primary development datasets (tumor + control) |
| Breast cancer | GEO | Validation dataset |
| Pancreatic cancer | GEO | Tumor microenvironment |
| Brain tumor | GEO | Glioblastoma / glioma |
| Non‑cancer lung | GEO | Healthy and disease (non‑malignant) |
| PBMC | Public | Immune reference for cross‑tissue benchmarking |

All processed successfully, confirming **robustness and generalizability** across tissues and conditions.  fileciteturn3file0

---

## 🎯 Quick Start

### Option 1 — Run in RStudio (Shiny)

1. Open `app.R`  
2. Click **Run App**

### Option 2 — R Console

```r
shiny::runApp("path_to_EVOscope_app_folder")
```

### Option 3 — Programmatic EPS (R)

```r
# assuming functions are sourced/packaged
eps <- EVOscope::compute_EPS(seurat_obj, k = 30, alpha = 1, beta = 1, gamma = 1)
EVOscope::plot_eps_umap(seurat_obj, eps)
write.csv(eps, "EPS_scores.csv", row.names = TRUE)
```

---

## 📦 Demo Files Included

| File Name | Description |
|-----------|-------------|
| `demo_lung_seurat.rds` | Lung cancer subsample (≈200 cells) |
| `demo_breast_seurat.rds` | Breast cancer subsample |
| `demo_pancreas_seurat.rds` | Pancreatic cancer subsample |
| `demo_brain_seurat.rds` | Brain tumor subsample |
| `demo_pbmc_seurat.rds` | PBMC subsample (reference) |

> ⚠️ These are small subsamples to meet GitHub size limits. Full datasets used for internal testing are available on request.  fileciteturn3file0

---

## 📊 Outputs

- `EPS_scores.csv` — per‑cell H, D, P, EPS and cluster labels  
- `EPS_cluster_summary.csv` — mean/variance per cluster  
- `EPS_pathway_enrichment.csv` — enriched pathways for high‑EPS genes  
- `plots/` — UMAPs, heatmaps, benchmarking and survival figures (PNG/PDF)

---

## ⚕️ Survival & Prognostic Analysis (Optional)

- Build an EPS‑associated gene signature (top EPS‑correlated genes)  
- **Univariate/Multivariate Cox** models with clinical covariates  
- **Kaplan–Meier** curves for high vs. low risk groups  
- **Time‑dependent ROC** at 1/3/5/8 years

> The prognostic module operates on **external bulk or microarray cohorts** (expression + survival metadata). You can provide your own cohort files; EVOscope includes helpers for model fitting and figure export.

---

## 📄 License

**MIT License** — free academic and commercial use with attribution.

---

## 📬 Contact

- **Lead developers:** Muhammad Zubair, Ali Haider  
- For questions or feature requests, please open a **GitHub Issue** in this repository.
