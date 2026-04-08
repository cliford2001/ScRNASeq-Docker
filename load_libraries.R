# =============================================================================
# load_libraries.R — Required packages for the scRNA-seq pipeline
# =============================================================================
# Usage: source("load_libraries.R")
# =============================================================================

# ── Core Seurat ecosystem ─────────────────────────────────────────────────────
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

# ── Single-cell utilities ─────────────────────────────────────────────────────
library(harmony)
library(DoubletFinder)
library(clustree)
library(monocle3)
library(scater)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(zellkonverter)

# ── Differential expression and GO enrichment ─────────────────────────────────
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)   # Arabidopsis — swap for org.Hs.eg.db / org.Mm.eg.db as needed

# ── File I/O ──────────────────────────────────────────────────────────────────
library(hdf5r)
library(Matrix)

# ── Data wrangling ────────────────────────────────────────────────────────────
library(tidyverse)
library(dplyr)
library(tibble)
library(reshape2)

# ── Visualisation ─────────────────────────────────────────────────────────────
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(cowplot)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(ggvenn)
library(eulerr)
library(UpSetR)

# ── Clustering helpers ────────────────────────────────────────────────────────
library(dynamicTreeCut)
library(WGCNA)

# ── Reporting ─────────────────────────────────────────────────────────────────
library(knitr)
library(kableExtra)

message("All libraries loaded successfully.")
