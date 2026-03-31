# =============================================================================
# Single-Cell RNA-seq Analysis Pipeline
# =============================================================================
# Protocol: Methods in Molecular Biology
#
# Description:
#   End-to-end scRNA-seq pipeline covering:
#     1.  Quality control and ambient-RNA-corrected data loading (CellBender)
#     2.  Doublet detection (DoubletFinder)
#     3.  Data integration and batch correction (Harmony)
#     4.  Clustering and resolution optimization (elbow + clustree)
#     5.  Cell-type annotation (bibliography markers + reference transfer)
#     6.  Differential expression analysis (DESeq2 pseudobulk)
#     7.  GO enrichment analysis (clusterProfiler)
#
# Dependencies:
#   load_libraries.R | ScRNA_Analysis_Functions.R | custom_seurat.R
#
# Usage:
#   Sections 1–13 can be sourced sequentially or run as a whole.
#   Section 11 (Cell-Type Curation) MUST be run interactively — do NOT
#   source the entire file with that section uncommented.
#
# Organism:
#   Default configuration targets Arabidopsis thaliana. To adapt for other
#   organisms, update MT_PATTERN, CP_PATTERN, ORGDB, and KEYTYPE in the
#   CONFIGURATION block below.
#
# Reproducibility:
#   A fixed random seed (1807) is applied at startup. All intermediate objects
#   are written to disk at key checkpoints.
#
# Docker environment:
#   The accompanying Dockerfile (rocker/r-ver:4.5) and docker-compose.yml
#   provision all required R and Python dependencies. The container working
#   directory is /workspace; data should be placed under DATA_DIR before
#   launching the container (see docker-compose.yml volumes section).
# =============================================================================

# =============================================================================
# CONFIGURATION — Edit this block to match your experiment
# =============================================================================

# Path to the directory containing the pipeline helper scripts
# (load_libraries.R, ScRNA_Analysis_Functions.R, custom_seurat.R).
# When running inside the Docker container this is typically /workspace/ScRNASeq-Docker.
PIPELINE_DIR <- "/workspace/ScRNASeq-Docker"

# Root directory for data and results. CellBender output is expected under
# DATA_DIR/cellbender/. All results are written to DATA_DIR/results/.
DATA_DIR     <- "/workspace/ScRNA"
output_dir   <- file.path(DATA_DIR, "results")   # used by save_pdf / save_qc helpers

# ── Sample manifest ───────────────────────────────────────────────────────────
# Each entry requires:
#   file      : path to CellBender-filtered HDF5, relative to DATA_DIR
#   label     : unique sample identifier (used in plots and column names)
#   condition : experimental condition (used for batch correction and DE)
samples <- list(
  list(file = "cellbender/Sample_0N_cellbender_filtered.h5",      label = "0N",      condition = "0N"),
  list(file = "cellbender/Sample_05N_R1_cellbender_filtered.h5",  label = "0.5N_R1", condition = "0.5N"),
  list(file = "cellbender/Sample_05N_2_cellbender_filtered.h5",   label = "0.5N_R2", condition = "0.5N"),
  list(file = "cellbender/Sample_5N_R1_cellbender_filtered.h5",   label = "5N_R1",   condition = "5N"),
  list(file = "cellbender/Sample_5N_2_cellbender_filtered.h5",    label = "5N_R2",   condition = "5N")
)

# ── Color palette (per condition) ─────────────────────────────────────────────
colors <- c(
  "0N"     = "#66c2a5",
  "0.5N_R1" = "#fc8d62", "0.5N_R2" = "#fc8d62",
  "5N_R1"  = "#8da0cb",  "5N_R2"   = "#8da0cb"
)

# ── Organelle gene patterns ────────────────────────────────────────────────────
# Arabidopsis : MT_PATTERN = "^ATMG"  |  CP_PATTERN = "^ATCG"
# Human       : MT_PATTERN = "^MT-"   |  CP_PATTERN = NULL
# Mouse       : MT_PATTERN = "^mt-"   |  CP_PATTERN = NULL
MT_PATTERN <- "^ATMG"
CP_PATTERN <- "^ATCG"

# ── Dimensionality reduction and neighbor graph ────────────────────────────────
DIMS_USE <- 1:30
K_PARAM  <- 30

# ── Clustering resolution sweep ───────────────────────────────────────────────
# Inspect elbow_plot.pdf and clustree.pdf before setting CLUSTER_RESOLUTION.
RESOLUTIONS_TEST <- c(0.15, 0.30, 0.50, 0.8, 1.0)
CLUSTER_RESOLUTION <- 0.3

# ── GO enrichment settings ────────────────────────────────────────────────────
# Organism databases:
#   Arabidopsis thaliana : ORGDB = org.At.tair.db, KEYTYPE = "TAIR"
#   Homo sapiens         : ORGDB = org.Hs.eg.db,   KEYTYPE = "ENSEMBL"
#   Mus musculus         : ORGDB = org.Mm.eg.db,   KEYTYPE = "ENSEMBL"
ORGDB      <- org.At.tair.db
KEYTYPE    <- "TAIR"
GO_ONT     <- "BP"    # Ontology namespace: "BP" | "MF" | "CC"
GO_QVAL    <- 0.05
GO_LEVEL   <- 6       # Maximum GO hierarchy depth for term pruning


# =============================================================================
# INITIALIZATION
# =============================================================================

source(file.path(PIPELINE_DIR, "load_libraries.R"))
source(file.path(PIPELINE_DIR, "custom_seurat.R"))
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))

set.seed(1807)
options(Seurat.allow.s4 = FALSE)

setwd(DATA_DIR)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECTION 1 — DATA LOADING AND PRE-FILTER QC
# =============================================================================
# Each sample is loaded from its CellBender-filtered HDF5 file. CellBender
# removes ambient RNA contamination prior to cell filtering. Mitochondrial
# (percent.mt) and chloroplast (percent.cp) read fractions are computed per
# cell to inform the filtering thresholds in Section 2.

seurat_list_raw <- lapply(samples, load_sample,
                          mt_pattern = MT_PATTERN,
                          cp_pattern = CP_PATTERN)
names(seurat_list_raw) <- sapply(samples, `[[`, "label")

plots_pre <- imap(seurat_list_raw, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))
save_qc(plots_pre, "qc_prefilter.pdf")


# =============================================================================
# SECTION 2 — CELL FILTERING AND DOUBLET DETECTION
# =============================================================================
# Cells failing quality thresholds (< 200 detected genes or > 5%
# mitochondrial reads) are removed. DoubletFinder is subsequently applied
# per sample; putative doublets are discarded before merging.
# Adjust min_features and max_mt based on the distributions in qc_prefilter.pdf.

seurat_list <- lapply(seurat_list_raw, filter_sample, min_features = 200, max_mt = 5)
names(seurat_list) <- sapply(samples, `[[`, "label")

plots_post <- imap(seurat_list, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))
save_qc(plots_post, "qc_postfilter.pdf")


# =============================================================================
# SECTION 3 — MERGE AND INITIAL PREPROCESSING
# =============================================================================
# Filtered samples are merged into a single Seurat object and preprocessed
# through log-normalization, variable feature selection (VST, 2,000 features),
# scaling, PCA (30 PCs), and UMAP. The resulting UMAP visualizes sample-level
# batch effects before integration.

pbmc_harmony <- reduce(seurat_list, merge) %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = DIMS_USE, verbose = FALSE)

pbmc_harmony$orig.ident_uni <- pbmc_harmony$condition

message("Cell counts per condition (pre-integration):")
print(table(pbmc_harmony$condition))
print(table(pbmc_harmony$orig.ident))

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_preharmony.pdf")

# Checkpoint — restore with: pbmc_harmony <- pbmc_harmony.bkp
pbmc_harmony.bkp <- pbmc_harmony


# =============================================================================
# SECTION 4 — HARMONY BATCH CORRECTION
# =============================================================================
# Harmony iteratively adjusts PCA cell embeddings to remove sample-level
# batch effects while preserving biological variation. The corrected
# "harmony" reduction replaces "pca" for all downstream steps.

pbmc_harmony <- pbmc_harmony %>%
  RunHarmony("orig.ident", plot_convergence = FALSE)


# =============================================================================
# SECTION 5 — RESOLUTION OPTIMIZATION
# =============================================================================
# Two complementary diagnostics guide the choice of clustering resolution:
#   (a) K-means elbow plot: inflection in within-cluster sum of squares (WSS)
#       across k = 2–40 suggests the number of biologically distinct clusters.
#   (b) Clustree: tracks cluster membership stability across Leiden resolutions,
#       helping to identify the lowest resolution that separates all major types.

# ── 5a. Elbow plot ────────────────────────────────────────────────────────────
k_range  <- 2:40
pca_data <- Embeddings(pbmc_harmony, "pca")[, DIMS_USE]
wss      <- sapply(k_range, function(k) {
  kmeans(pca_data, centers = k, nstart = 10)$tot.withinss
})

elbow_plot <- ggplot(data.frame(k = k_range, wss = wss), aes(k, wss)) +
  geom_line() + geom_point() +
  labs(x = "Number of clusters (k)", y = "Within-cluster sum of squares") +
  theme_minimal()

save_pdf(elbow_plot, "elbow_plot.pdf", w = 8, h = 6)

# ── 5b. Clustree — resolution sweep ───────────────────────────────────────────
clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = DIMS_USE, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = DIMS_USE,
                k.param = K_PARAM, verbose = FALSE)

for (res in RESOLUTIONS_TEST)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(clustree(clu, prefix = "RNA_snn_res."), "clustree.pdf", w = 14, h = 14)


# =============================================================================
# SECTION 6 — FINAL CLUSTERING
# =============================================================================
# Leiden clustering (algorithm = 4) with the resolution selected after
# inspecting the elbow and clustree diagnostics. All downstream analyses
# use the cluster identities produced here.

pbmc_harmony <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = DIMS_USE, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = DIMS_USE,
                k.param = K_PARAM, verbose = FALSE) %>%
  FindClusters(resolution = CLUSTER_RESOLUTION, algorithm = 4, verbose = FALSE)

message("Cells per cluster:")
print(table(Idents(pbmc_harmony)))

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_postharmony.pdf")

colors_clusters <- sample(colors(distinct = TRUE),
                          length(unique(pbmc_harmony$seurat_clusters)))
Idents(pbmc_harmony) <- "seurat_clusters"
save_pdf(DimPlot(pbmc_harmony, group.by = "seurat_clusters", cols = colors_clusters),
         "umap_seuratclusters.pdf")


# =============================================================================
# SECTION 7 — CELL-TYPE ANNOTATION
# =============================================================================
# Two complementary strategies assign cell-type labels to clusters:
#
#   (a) Bibliography-based markers (biblio_marks.txt): a curated, tab-separated
#       table of known marker genes per cell type (columns: "cell types", "gene").
#       Cluster-level differential genes (FindAllMarkers) are crossed with this
#       table to assign initial labels. R reads "cell types" as "cell.types".
#
#   (b) Reference transfer: a published Seurat object (Arabidopsis leaf atlas,
#       GSE273033; Shahan et al.) serves as reference for label transfer via
#       FindTransferAnchors / TransferData. The resulting prediction is stored
#       in pbmc_harmony$celltype_reference.
#
# biblio_marks.txt format example:
#   cell types         gene
#   guard cell         AT5G26000
#   mesophyll          AT5G38420

marcadores <- read.table(file.path(output_dir, "biblio_marks.txt"),
                         header = TRUE, sep = "\t", quote = "")

# ── 7a. Bibliography-based annotation ─────────────────────────────────────────
markers <- find_markers(pbmc_harmony,
                        output_file = file.path(output_dir, "FindAllMarkers.tsv"))

pbmc_harmony <- annotate_by_markers(pbmc_harmony, markers,
                                    reference_file = file.path(output_dir, "biblio_marks.txt"))

# ── 7b. Reference-based annotation ────────────────────────────────────────────
# The reference RDS file contains a published Seurat object with cell-type
# annotations. It must be present in output_dir prior to running this step.
esp          <- readRDS(file.path(output_dir, "GSE273033_seuratObj_for_publication.rds"))
pbmc_harmony <- annotate_by_reference(pbmc_harmony,
                                      reference_obj = esp,
                                      reference_col = "annotation")
# Annotation stored in: pbmc_harmony$celltype_reference


# =============================================================================
# SECTION 8 — ANNOTATED CLUSTREE
# =============================================================================
# Re-runs the resolution sweep on the annotated object. Each node is labelled
# with the modal cell-type annotation at that resolution, allowing visual
# verification that the chosen resolution cleanly separates known cell types.

Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = DIMS_USE, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = DIMS_USE,
                k.param = K_PARAM, verbose = FALSE)

for (res in RESOLUTIONS_TEST)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(
  clustree(clu, prefix = "RNA_snn_res.",
           node_label = "celltype_reference", node_label_aggr = "Mode"),
  "clustree_annotated.pdf", w = 14, h = 14
)


# =============================================================================
# SECTION 9 — GENE EXPRESSION VISUALIZATION
# =============================================================================
# Violin and feature plots are generated for individual genes and gene sets,
# both across all cell types and within a specific cell type of interest.
# JoinLayers() is required in Seurat 5 before subsetting after merge.

pbmc_harmony     <- JoinLayers(pbmc_harmony)
Idents(pbmc_harmony) <- "celltype_reference"

gene              <- "AT5G26000"
genes_of_interest <- c("AT5G26000", "AT5G54250")
celltype          <- "Guard Cell"   # Must match exact label in celltype_reference

sub_obj <- subset(pbmc_harmony, idents = celltype)

save_vln(VlnPlot(pbmc_harmony, features = gene),                    "vln_gene_all.pdf")
save_pdf(FeaturePlot(pbmc_harmony, features = gene),                "feature_gene_all.pdf")
save_vln(VlnPlot(pbmc_harmony, features = genes_of_interest),       "vln_geneset_all.pdf",
         n = length(genes_of_interest))
save_pdf(FeaturePlot(pbmc_harmony, features = genes_of_interest),   "feature_geneset_all.pdf",
         h = 8 * length(genes_of_interest))

save_vln(VlnPlot(sub_obj, features = gene),                         "vln_gene_celltype.pdf")
save_pdf(FeaturePlot(sub_obj, features = gene),                     "feature_gene_celltype.pdf")
save_vln(VlnPlot(sub_obj, features = genes_of_interest),            "vln_geneset_celltype.pdf",
         n = length(genes_of_interest))
save_pdf(FeaturePlot(sub_obj, features = genes_of_interest),        "feature_geneset_celltype.pdf",
         h = 8 * length(genes_of_interest))


# =============================================================================
# SECTION 10 — CELL-TYPE GROUPING
# =============================================================================
# Fine-grained reference annotations are optionally collapsed into broader
# functional categories for downstream pseudobulk and GO analyses.
# Cell types not listed in 'grouping' retain their original label.

grouping <- c(
  "Companion Cell"    = "Vascular Cell",
  "Cambium"           = "Vascular Cell",
  "Phloem Parenchyma" = "Vascular Cell",
  "Xylem"             = "Vascular Cell",
  "Sieve Element"     = "Vascular Cell",
  "Meristemoid"       = "Stomatal Line"
)

pbmc_harmony$annotation_agrupada <- recode(pbmc_harmony$celltype_reference, !!!grouping)

save_pdf(
  DimPlot(pbmc_harmony, group.by = "annotation_agrupada",
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_annotated.pdf"
)


# =============================================================================
# SECTION 11 — INTERACTIVE CELL-TYPE CURATION
# =============================================================================
# !!! WARNING: Run this section interactively, step by step.
# !!! Do NOT source the entire script while this section is active.
#
# Purpose: manually inspect and correct annotations for populations that
# appear heterogeneous in the UMAP. Subclustering resolves mixed clusters
# that cannot be separated at the global resolution.
#
# Workflow:
#   Step 1. Subcluster cell types that appear heterogeneous.
#   Step 2. Inspect subclusters with DimPlot and FeaturePlot.
#   Step 3. Fill in the reassignment table (reassign) below.
#   Step 4. Apply corrections to the global object.

Idents(pbmc_harmony) <- "annotation_agrupada"

# ── Step 1. Subcluster ────────────────────────────────────────────────────────
meristemoid_umap   <- subclustar_tipo(pbmc_harmony, "Stomatal Line")
pavement_cell_umap <- subclustar_tipo(pbmc_harmony, "Pavement Cell")

# ── Step 2. Inspect ───────────────────────────────────────────────────────────
DimPlot(meristemoid_umap,   group.by = "cluster_subtipo", label = TRUE, raster = FALSE)
DimPlot(pavement_cell_umap, group.by = "cluster_subtipo", label = TRUE, raster = FALSE)

for (i in seq_len(nrow(marcadores)))
  print(FeaturePlot(pavement_cell_umap, features = marcadores$gene[i]) +
        ggtitle(paste(marcadores$cell.types[i], "-", marcadores$gene[i])))

FindAllMarkers(meristemoid_umap,   only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)

FindAllMarkers(pavement_cell_umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)

# ── Step 3. Reassignment table ────────────────────────────────────────────────
# Map subcluster IDs (from $cluster_subtipo, values "0", "1", ...) to final
# cell-type labels. Add one list entry per subclustered cell type.

reassign <- list(
  meristemoid_umap = c(
    "0" = "Stomatal Line",
    "1" = "Stomatal Line",
    "2" = "Pavement Cell",
    "3" = "Stomatal Line",
    "4" = "Stomatal Line"
  ),
  pavement_cell_umap = c(
    "0"      = "Pavement Cell",
    "1"      = "Pavement Cell",
    "2"      = "Pavement Cell",
    "3"      = "Mesophyll",
    "4"      = "Pavement Cell",
    "others" = "Pavement Cell"
  )
)

# ── Step 4. Apply corrections ─────────────────────────────────────────────────
pbmc_harmony$celltype_reference_curated <- pbmc_harmony$annotation_agrupada

for (obj_name in names(reassign)) {
  obj <- get(obj_name)
  pbmc_harmony$celltype_reference_curated[colnames(obj)] <-
    reassign[[obj_name]][obj$cluster_subtipo]
}

save_pdf(
  DimPlot(pbmc_harmony, group.by = "celltype_reference_curated",
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_curada.pdf"
)


# =============================================================================
# SECTION 12 — DOTPLOT: MARKER GENES BY CURATED CELL TYPE
# =============================================================================
# DotPlot of bibliography-derived marker genes across curated cell types.
# Dot size encodes the fraction of expressing cells; color encodes mean
# scaled expression. A near-diagonal expression pattern validates the
# annotation.

cell_order_dotplot <- c(
  "Pavement Cell", "Stomatal lineage", "Guard Cell", "Mesophyll", "Bundle Sheath",
  "Phloem Parenchyma", "Cambium", "Xylem", "Companion Cell",
  "Hydathode", "Cell Cycle: G1-S", "Cell Cycle: G2-M"
)

hacer_dotplot_marcadores(
  pbmc_harmony,
  marcadores,
  annot_col       = "celltype_reference_curated",
  cell_order      = cell_order_dotplot,
  clusters_remove = c("Sieve Element", "Myrosin Idioblast"),
  outfile         = file.path(output_dir, "dotplot_marcadores.pdf"),
  width = 20, height = 10
)


# =============================================================================
# SECTION 13 — EXPORT TO H5AD (Scanpy / Python)
# =============================================================================
# The curated Seurat object is exported to AnnData h5ad format for downstream
# trajectory and RNA velocity analyses in Python. Scanpy, scFates, and
# Palantir are pre-installed in the Docker image.

exportar_para_scanpy(pbmc_harmony,
                     file.path(output_dir, "objs/pbmc_harmony_curated.h5ad"))

# To export a specific cell type subset:
# exportar_para_scanpy(
#   subset(pbmc_harmony, subset = celltype_reference_curated == "Guard Cell"),
#   file.path(output_dir, "objs/GuardCell.h5ad")
# )


# =============================================================================
# PART 2 — PSEUDOBULK DIFFERENTIAL EXPRESSION AND GO ENRICHMENT
# =============================================================================

# =============================================================================
# SECTION 14 — PER-CELL-TYPE SUBSETS
# =============================================================================
# The integrated object is split into one Seurat object per curated cell type
# for independent pseudobulk differential expression analyses.

celular_subsets <- setNames(
  lapply(unique(pbmc_harmony$celltype_reference_curated),
         function(t) subset(pbmc_harmony, subset = celltype_reference_curated == t)),
  gsub("[^[:alnum:]_]", "_", unique(pbmc_harmony$celltype_reference_curated))
)


# =============================================================================
# SECTION 15 — PSEUDO-REPLICATES AND PSEUDOBULK COUNT MATRICES
# =============================================================================
# Because the experiment lacks within-condition biological replicates, cells
# are randomly assigned to n = 3 pseudo-replicates per condition. Counts are
# then aggregated per pseudo-replicate to build the count matrix required by
# DESeq2. Cell types with fewer than two conditions are excluded automatically.

celular_subsets_replicados <- Filter(Negate(is.null),
                                     lapply(celular_subsets, asignar_pseudoreplicados))
pseudobulk_list            <- lapply(celular_subsets_replicados, hacer_pseudobulk)

dir.create(file.path(output_dir, "pseudobulk_replicas"), recursive = TRUE, showWarnings = FALSE)
for (tipo in names(pseudobulk_list))
  write.csv(pseudobulk_list[[tipo]],
            file.path(output_dir, "pseudobulk_replicas", paste0("Pseudobulk_", tipo, ".csv")),
            row.names = TRUE)


# =============================================================================
# SECTION 16 — DESEQ2 DIFFERENTIAL EXPRESSION
# =============================================================================
# Pairwise contrasts are run for all three nitrogen treatment comparisons.
# Results are written as CSV files, one per cell type per comparison, under
# output_dir/deseq2/<tag>/.

comparaciones <- list(
  list(conds = c("0.5N", "5N"),  tag = "05_5"),
  list(conds = c("0N",   "5N"),  tag = "0_5"),
  list(conds = c("0N",   "0.5N"), tag = "0_05")
)

for (tag in sapply(comparaciones, `[[`, "tag"))
  dir.create(file.path(output_dir, "deseq2", tag), recursive = TRUE, showWarnings = FALSE)

for (tipo in names(pseudobulk_list))
  correr_deseq2(as.matrix(pseudobulk_list[[tipo]]), comparaciones,
                output_dir = file.path(output_dir, "deseq2"), tipo = tipo)


# =============================================================================
# SECTION 17 — VOLCANO PLOTS
# =============================================================================
# One volcano plot per cell type per comparison; two plots per PDF page.
# Points are colored by significance category (padj < 0.05, |log2FC| > 1).

dir.create(file.path(output_dir, "volcano_plots"), showWarnings = FALSE)

for (comp in comparaciones) {
  tag       <- comp$tag
  csv_files <- list.files(file.path(output_dir, "deseq2", tag),
                          pattern = "\\.csv$", full.names = TRUE)
  if (!length(csv_files)) { message("No CSV files found for comparison: ", tag); next }

  pdf(file.path(output_dir, "volcano_plots", paste0("VolcanoPlots_", tag, ".pdf")),
      width = 12, height = 6)
  plots <- list()
  for (f in csv_files) {
    plots <- c(plots, list(hacer_volcano(f)))
    if (length(plots) == 2) { grid.arrange(grobs = plots, ncol = 2); plots <- list() }
  }
  if (length(plots)) grid.arrange(grobs = plots, ncol = 1)
  dev.off()
}


# =============================================================================
# SECTION 18 — DIFFERENTIAL GENE TABLES AND HEATMAPS
# =============================================================================
# Differentially expressed genes are classified (up = 1, down = −1,
# unchanged = 0) and summarized in a cross-cell-type table. A hierarchically
# clustered heatmap of log2FC values is generated per comparison.

for (comp in comparaciones) {
  tag       <- comp$tag
  csv_dir   <- file.path(output_dir, "deseq2", tag)
  diff_dir  <- file.path(output_dir, "diff",   tag)
  csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

  dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(csv_files)) { message("No CSV files found for comparison: ", tag); next }

  listas      <- lapply(csv_files, procesar_deseq2_resultado, output_dir = diff_dir)
  tabla_class <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), lapply(listas, `[[`, "class"))
  tabla_logfc <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), lapply(listas, `[[`, "logfc"))

  tabla_class <- tabla_class %>% filter(apply(select(., -gene_id) != 0, 1, any))
  tabla_logfc <- tabla_logfc %>% filter(gene_id %in% tabla_class$gene_id)

  write_tsv(tabla_class, file.path(diff_dir, "tabla_diferenciales.tsv"))
  write_tsv(tabla_logfc, file.path(diff_dir, "tabla_log2FC.tsv"))

  matriz <- as.matrix(column_to_rownames(tabla_logfc, "gene_id"))
  matriz[is.na(matriz)] <- 0

  if (nrow(matriz) > 1) {
    pdf(file.path(diff_dir, paste0("heatmap_", tag, ".pdf")), width = 14, height = 18)
    tryCatch(hacer_heatmap(matriz), error = function(e) message("Heatmap error: ", e$message))
    dev.off()
  }
}


# =============================================================================
# SECTION 19 — GO ENRICHMENT ANALYSIS
# =============================================================================
# Gene Ontology enrichment (clusterProfiler::enrichGO) is run for each cell
# type and pairwise comparison. Both full and semantically simplified (redundancy-
# reduced) result sets are generated and visualized as balloon plots.

universo <- keys(ORGDB, keytype = KEYTYPE)

dir.create(file.path(output_dir, "Enrichment"), showWarnings = FALSE)

for (comp in comparaciones) {
  tag        <- comp$tag
  diff_dir   <- file.path(output_dir, "diff",       tag)
  enr_dir    <- file.path(output_dir, "Enrichment", tag)
  tabla_path <- file.path(diff_dir, "tabla_diferenciales.tsv")

  dir.create(enr_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(tabla_path)) { message("No differential table for: ", tag); next }

  tabla <- read.table(tabla_path, header = TRUE, row.names = 1, sep = "\t")
  tabla <- tabla[, colSums(tabla != 0) > 0, drop = FALSE]

  go_total  <- correr_enriquecimiento_go(tabla, universo, GO_ONT,
                                         orgdb = ORGDB, keytype = KEYTYPE,
                                         simplificar = FALSE, output_dir = enr_dir)
  go_simple <- correr_enriquecimiento_go(tabla, universo, GO_ONT,
                                         orgdb = ORGDB, keytype = KEYTYPE,
                                         simplificar = TRUE,  output_dir = enr_dir)

  go_total_podado  <- podar_go(go_total,  GO_LEVEL, GO_ONT, GO_QVAL,
                               simplificar = FALSE, output_dir = enr_dir)
  go_simple_podado <- podar_go(go_simple, GO_LEVEL, GO_ONT, GO_QVAL,
                               simplificar = TRUE,  output_dir = enr_dir)

  pdf(file.path(enr_dir, paste0("GO_enrichment_", tag, ".pdf")), width = 18, height = 18)
  tryCatch({
    print(graficar_go_balones(go_total))
    print(graficar_go_balones(go_simple))
    print(graficar_go_balones(go_total_podado))
    print(graficar_go_balones(go_simple_podado))
  }, error = function(e) message("GO plot error: ", e$message))
  dev.off()
}


# =============================================================================
# SECTION 20 — GLOBAL PSEUDOBULK AND REPLICATE CORRELATION
# =============================================================================
# A global pseudobulk aggregated by sample (orig.ident) is generated and
# visualized as a Pearson correlation heatmap to verify replicate consistency
# across samples before differential expression.

pseudobulk <- generate_pseudobulk(pbmc_harmony, group_by = "orig.ident")

save_pdf(plot_replicate_correlation(pseudobulk$by_sample),
         "pseudobulk_correlation.pdf", w = 8, h = 8)


# =============================================================================
# END OF PIPELINE
# =============================================================================
