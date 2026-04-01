# =============================================================================
# Single-Cell RNA-seq Analysis Pipeline
# =============================================================================
# Protocol: Methods in Molecular Biology
#
# Description:
#   End-to-end scRNA-seq pipeline covering:
#     Part 1 — QC, integration, clustering, annotation, and export
#     Part 2 — Pseudobulk differential expression and GO enrichment
#
# Helper scripts (fully documented at https://github.com/cliford2001/ScRNASeq-Docker):
#   load_libraries.R           — loads all required R packages
#   ScRNA_Analysis_Functions.R — core functions (QC, clustering, annotation, DE, GO)
#   custom_seurat.R            — custom Seurat plot utilities (cluster bar charts)
#
# Usage:
#   Run sections sequentially. Section 12 (Cell-Type Curation) must be
#   run interactively — do NOT source the entire script with it active.
#
# Organism:
#   Default: Arabidopsis thaliana. Adapting to another organism requires
#   updating only the organelle patterns, GO database, and key type
#   (each marked clearly in the relevant section below).
#
# Reproducibility:
#   Fixed random seed (1807). Intermediate results written to disk.
# =============================================================================

# =============================================================================
# CONFIGURATION — Edit only this block before running the pipeline
# =============================================================================

# Directory containing the pipeline helper scripts.
# Inside the Docker container this is typically /workspace/ScRNASeq-Docker
PIPELINE_DIR <- "/workspace/ScRNASeq-Docker"

# Root directory for your project data and results.
# All result files will be written to DATA_DIR/results/<step>/
DATA_DIR   <- "/workspace/ScRNA"
base_dir   <- file.path(DATA_DIR, "results")

# ── Input format ──────────────────────────────────────────────────────────────
# USE_CELLBENDER = TRUE  → load CellBender-filtered HDF5 files (recommended)
# USE_CELLBENDER = FALSE → load CellRanger filtered_feature_bc_matrix/ directly
#                          (use this if you skipped the CellBender step)
#
# If TRUE  → set samples$file to the .h5 file path relative to DATA_DIR
# If FALSE → set samples$file to the filtered_feature_bc_matrix/ directory
#            path relative to DATA_DIR
USE_CELLBENDER <- TRUE

# ── Sample manifest ───────────────────────────────────────────────────────────
# Add one entry per sample. Each entry needs:
#   file      — path to the input file or directory (relative to DATA_DIR)
#   label     — unique name for this sample (appears in all plots)
#   condition — experimental group this sample belongs to
samples <- list(
  list(file = "cellbender/Sample_0N_cellbender_filtered.h5",      label = "0N",      condition = "0N"),
  list(file = "cellbender/Sample_05N_R1_cellbender_filtered.h5",  label = "0.5N_R1", condition = "0.5N"),
  list(file = "cellbender/Sample_05N_2_cellbender_filtered.h5",   label = "0.5N_R2", condition = "0.5N"),
  list(file = "cellbender/Sample_5N_R1_cellbender_filtered.h5",   label = "5N_R1",   condition = "5N"),
  list(file = "cellbender/Sample_5N_2_cellbender_filtered.h5",    label = "5N_R2",   condition = "5N")
)

# ── Plot colors (one color per sample label) ───────────────────────────────────
colors <- c(
  "0N"      = "#66c2a5",
  "0.5N_R1" = "#fc8d62", "0.5N_R2" = "#fc8d62",
  "5N_R1"   = "#8da0cb", "5N_R2"   = "#8da0cb"
)


# =============================================================================
# INITIALIZATION
# =============================================================================

# ── Load helper scripts ────────────────────────────────────────────────────────
# Each file is fully documented at https://github.com/cliford2001/ScRNASeq-Docker
source(file.path(PIPELINE_DIR, "load_libraries.R"))          # all R packages
source(file.path(PIPELINE_DIR, "custom_seurat.R"))           # plot_integrated_clusters()
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))# analysis functions

set.seed(1807)
options(Seurat.allow.s4 = FALSE)
setwd(DATA_DIR)

# ── Create per-step output directories ────────────────────────────────────────
# mkd() creates the folder and returns its path in one call.
mkd <- function(name) { d <- file.path(base_dir, name); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

dir_00 <- mkd("00_workflow")
dir_01 <- mkd("01_qc")
dir_02 <- mkd("02_filtering")
dir_03 <- mkd("03_integration")
dir_04 <- mkd("04_clustering")
dir_05 <- mkd("05_annotation")
dir_06 <- mkd("06_expression")
dir_07 <- mkd("07_curation")
dir_08 <- mkd("08_export")
dir_09 <- mkd("09_pseudobulk")
dir_10 <- mkd("10_deseq2")
dir_11 <- mkd("11_volcano")
dir_12 <- mkd("12_heatmaps")
dir_13 <- mkd("13_go")

# output_dir is the global variable used by save_pdf / save_qc / save_vln helpers.
# It is reassigned at the start of each section to the appropriate step directory.
output_dir <- base_dir


# =============================================================================
# SECTION 0 — PIPELINE WORKFLOW FIGURE
# =============================================================================
# Generates a visual overview of the full pipeline saved to 00_workflow/.
# Run this section once immediately after initialization.

plot_pipeline_workflow(file.path(dir_00, "pipeline_workflow.pdf"))


# =============================================================================
# ████████████████████████  PART 1 — SINGLE-CELL ANALYSIS  ████████████████████
# =============================================================================


# =============================================================================
# SECTION 1 — DATA LOADING AND PRE-FILTER QC
# =============================================================================
# Each sample is loaded from its input file and mitochondrial / chloroplast
# read fractions are computed per cell to guide filtering thresholds in
# Section 2. Pre-filter violin plots are saved to 01_qc/.
#
# ┌─ CHANGE FOR YOUR ORGANISM ──────────────────────────────────────────────────
#   Arabidopsis : mt_pattern = "^ATMG"  |  cp_pattern = "^ATCG"
#   Human       : mt_pattern = "^MT-"   |  cp_pattern = NULL
#   Mouse       : mt_pattern = "^mt-"   |  cp_pattern = NULL
# └─────────────────────────────────────────────────────────────────────────────
mt_pattern <- "^ATMG"
cp_pattern <- "^ATCG"

output_dir <- dir_01

if (USE_CELLBENDER) {
  # Load CellBender-filtered HDF5 files (recommended)
  seurat_list_raw <- lapply(samples, load_sample,
                            mt_pattern = mt_pattern,
                            cp_pattern = cp_pattern)
} else {
  # Load directly from CellRanger filtered_feature_bc_matrix/ directories.
  # Use this path if you skipped the CellBender step.
  seurat_list_raw <- lapply(samples, function(s) {
    mat <- Read10X(data.dir = file.path(DATA_DIR, s$file))
    obj <- CreateSeuratObject(counts = mat, project = s$label,
                              min.cells = 3, min.features = 200)
    obj$condition             <- s$condition
    obj[["percent.mt"]]       <- PercentageFeatureSet(obj, pattern = mt_pattern)
    if (!is.null(cp_pattern))
      obj[["percent.cp"]]     <- PercentageFeatureSet(obj, pattern = cp_pattern)
    obj
  })
}
names(seurat_list_raw) <- sapply(samples, `[[`, "label")

plots_pre <- imap(seurat_list_raw, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))
save_qc(plots_pre, "qc_prefilter.pdf")


# =============================================================================
# SECTION 2 — CELL FILTERING AND DOUBLET DETECTION
# =============================================================================
# Thresholds are set based on the pre-filter violin plots (01_qc/qc_prefilter.pdf).
# Cells below min_features or above max_mt are removed. DoubletFinder is
# applied per sample to discard putative doublets.
#
# ┌─ ADJUST THRESHOLDS AFTER INSPECTING 01_qc/qc_prefilter.pdf ─────────────────
#   min_features : minimum number of detected genes per cell (default 200)
#   max_mt       : maximum mitochondrial read percentage  (default 5 %)
# └─────────────────────────────────────────────────────────────────────────────
min_features <- 200
max_mt       <- 5

output_dir <- dir_02

seurat_list <- lapply(seurat_list_raw, filter_sample,
                      min_features = min_features, max_mt = max_mt)
names(seurat_list) <- sapply(samples, `[[`, "label")

plots_post <- imap(seurat_list, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))
save_qc(plots_post, "qc_postfilter.pdf")


# =============================================================================
# SECTION 3 — MERGE AND INITIAL PREPROCESSING
# =============================================================================
# Filtered samples are merged and preprocessed: log-normalization, variable
# feature selection (VST, 2,000 features), scaling, PCA (30 PCs), and UMAP.
# The resulting UMAP shows batch effects before integration.

output_dir <- dir_03

pbmc_harmony <- reduce(seurat_list, merge) %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

pbmc_harmony$orig.ident_uni <- pbmc_harmony$condition

message("Cell counts per condition (pre-integration):")
print(table(pbmc_harmony$condition))

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_preharmony.pdf")

# Checkpoint — restore with: pbmc_harmony <- pbmc_harmony.bkp
pbmc_harmony.bkp <- pbmc_harmony


# =============================================================================
# SECTION 4 — HARMONY BATCH CORRECTION
# =============================================================================
# Harmony adjusts cell embeddings to remove sample-level batch effects while
# preserving biological variation. All downstream steps use the "harmony"
# reduction instead of "pca".
#
# ┌─ DIMENSIONALITY PARAMETERS ─────────────────────────────────────────────────
#   dims_use : how many Harmony dimensions to use downstream (default 1:30)
#   k_param  : number of nearest neighbors for the cell graph (default 30)
# └─────────────────────────────────────────────────────────────────────────────
dims_use <- 1:30
k_param  <- 30

pbmc_harmony <- pbmc_harmony %>%
  RunHarmony("orig.ident", plot_convergence = FALSE)


# =============================================================================
# SECTION 5 — RESOLUTION OPTIMIZATION
# =============================================================================
# Two diagnostics guide the choice of clustering resolution:
#   (a) Elbow plot — k-means within-cluster sum of squares across k = 2-40.
#       The inflection point suggests the number of major cell types.
#   (b) Clustree  — tracks cluster stability across Leiden resolutions.
#       Choose the lowest resolution where clusters stop merging.
#
# ┌─ RESOLUTIONS TO TEST ───────────────────────────────────────────────────────
#   Inspect 04_clustering/clustree.pdf and elbow_plot.pdf before choosing
#   cluster_resolution in Section 6.
# └─────────────────────────────────────────────────────────────────────────────
resolutions_test <- c(0.15, 0.30, 0.50, 0.8, 1.0)

output_dir <- dir_04

# ── 5a. Elbow plot ────────────────────────────────────────────────────────────
k_range  <- 2:40
pca_data <- Embeddings(pbmc_harmony, "pca")[, dims_use]
wss      <- sapply(k_range, function(k) {
  kmeans(pca_data, centers = k, nstart = 10)$tot.withinss
})

elbow_plot <- ggplot(data.frame(k = k_range, wss = wss), aes(k, wss)) +
  geom_line() + geom_point() +
  labs(x = "Number of clusters (k)", y = "Within-cluster sum of squares") +
  theme_minimal()

save_pdf(elbow_plot, "elbow_plot.pdf", w = 8, h = 6)

# ── 5b. Clustree ──────────────────────────────────────────────────────────────
clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use,
                k.param = k_param, verbose = FALSE)

for (res in resolutions_test)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(clustree(clu, prefix = "RNA_snn_res."), "clustree.pdf", w = 14, h = 14)


# =============================================================================
# SECTION 6 — FINAL CLUSTERING
# =============================================================================
# Apply the selected resolution for the final cluster assignment.
# After clustering, a UMAP coloured by sample identity (umap_postharmony.pdf)
# and a simple bar chart of cells per sample are saved.
#
# ┌─ SET RESOLUTION AFTER INSPECTING elbow_plot.pdf AND clustree.pdf ──────────
#   cluster_resolution : Leiden resolution for final clustering (default 0.3)
# └─────────────────────────────────────────────────────────────────────────────
cluster_resolution <- 0.3

output_dir <- dir_04

pbmc_harmony <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use,
                k.param = k_param, verbose = FALSE) %>%
  FindClusters(resolution = cluster_resolution, algorithm = 4, verbose = FALSE)

message("Cells per cluster:")
print(table(Idents(pbmc_harmony)))

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_postharmony.pdf")

colors_clusters <- sample(colors(distinct = TRUE),
                          length(unique(pbmc_harmony$seurat_clusters)))
Idents(pbmc_harmony) <- "seurat_clusters"
save_pdf(DimPlot(pbmc_harmony, group.by = "seurat_clusters", cols = colors_clusters),
         "umap_seuratclusters.pdf")

# ── Cell count per sample ─────────────────────────────────────────────────────
cell_count_plot <- ggplot(pbmc_harmony@meta.data,
                          aes(x = orig.ident, fill = orig.ident)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.4, size = 4) +
  scale_fill_manual(values = colors) +
  theme_bw(base_size = 14) +
  labs(title = "Total cells per sample after filtering and integration",
       x = "Sample", y = "Cell count") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_pdf(cell_count_plot, "cell_count_per_sample.pdf", w = 8, h = 6)


# =============================================================================
# SECTION 7 — DOTPLOT: MARKER GENES BY CLUSTER (pre-annotation guide)
# =============================================================================
# Before assigning cell-type labels, this plot helps you identify which
# numbered Seurat cluster corresponds to which cell type by showing the
# expression of bibliography-derived marker genes across all clusters.
# Clusters that strongly express a known marker (e.g., AT5G26000 for Guard
# Cell) should be labelled as that cell type in Section 8.
# Dot size = fraction of expressing cells; color = mean expression level.

output_dir <- dir_05

marcadores <- read.table(file.path(base_dir, "biblio_marks.txt"),
                         header = TRUE, sep = "\t", quote = "")

hacer_dotplot_marcadores(
  pbmc_harmony,
  marcadores,
  annot_col = "seurat_clusters",
  outfile   = file.path(output_dir, "dotplot_marcadores_preannotation.pdf"),
  width = 20, height = 10
)


# =============================================================================
# SECTION 8 — CELL-TYPE ANNOTATION
# =============================================================================
# Two strategies assign cell-type labels to clusters:
#
#   (a) Bibliography-based markers (biblio_marks.txt) — the same marker table
#       read above is crossed with cluster-level differential genes to assign
#       initial labels automatically.
#
#   (b) Reference transfer — labels are projected from a published Seurat
#       object (Arabidopsis leaf atlas, GSE273033) using FindTransferAnchors
#       and TransferData. Result stored in pbmc_harmony$celltype_reference.

output_dir <- dir_05

# ── 8a. Bibliography-based annotation ─────────────────────────────────────────
markers <- find_markers(pbmc_harmony,
                        output_file = file.path(output_dir, "FindAllMarkers.tsv"))

pbmc_harmony <- annotate_by_markers(pbmc_harmony, markers,
                                    reference_file = file.path(base_dir, "biblio_marks.txt"))

# ── 8b. Reference-based annotation ────────────────────────────────────────────
esp          <- readRDS(file.path(base_dir, "GSE273033_seuratObj_for_publication.rds"))
pbmc_harmony <- annotate_by_reference(pbmc_harmony,
                                      reference_obj = esp,
                                      reference_col = "annotation")
# Annotation stored in: pbmc_harmony$celltype_reference


# =============================================================================
# SECTION 9 — ANNOTATED CLUSTREE
# =============================================================================
# Re-runs the resolution sweep with cell-type labels overlaid on each node.
# Confirms that the chosen resolution cleanly separates known cell types.

output_dir <- dir_05

Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use,
                k.param = k_param, verbose = FALSE)

for (res in resolutions_test)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(
  clustree(clu, prefix = "RNA_snn_res.",
           node_label = "celltype_reference", node_label_aggr = "Mode"),
  "clustree_annotated.pdf", w = 14, h = 14
)


# =============================================================================
# SECTION 10 — GENE EXPRESSION VISUALIZATION
# =============================================================================
# Violin and feature plots for individual genes or gene sets, generated both
# across all cell types and within a specific cell type of interest.
# JoinLayers() is required in Seurat 5 before subsetting after merge.
#
# ┌─ SET YOUR GENES AND CELL TYPE OF INTEREST ──────────────────────────────────
#   gene              : single gene to inspect
#   genes_of_interest : gene set to inspect together
#   celltype          : cell type to zoom into (must match celltype_reference)
# └─────────────────────────────────────────────────────────────────────────────
gene              <- "AT5G26000"
genes_of_interest <- c("AT5G26000", "AT5G54250")
celltype          <- "Guard Cell"

output_dir <- dir_06

pbmc_harmony     <- JoinLayers(pbmc_harmony)
Idents(pbmc_harmony) <- "celltype_reference"

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
# SECTION 11 — CELL-TYPE GROUPING
# =============================================================================
# Fine-grained labels are collapsed into broader categories for downstream
# analyses. Cell types NOT listed in 'grouping' keep their original label.
#
# ┌─ EDIT THIS MAP TO MATCH YOUR CELL TYPES ───────────────────────────────────
#   Left side  : original label (must match exactly)
#   Right side : new broader label to assign
# └─────────────────────────────────────────────────────────────────────────────
grouping <- c(
  "Companion Cell"    = "Vascular Cell",
  "Cambium"           = "Vascular Cell",
  "Phloem Parenchyma" = "Vascular Cell",
  "Xylem"             = "Vascular Cell",
  "Sieve Element"     = "Vascular Cell",
  "Meristemoid"       = "Stomatal Line"
)

output_dir <- dir_07

pbmc_harmony$annotation_agrupada <- recode(pbmc_harmony$celltype_reference, !!!grouping)

save_pdf(
  DimPlot(pbmc_harmony, group.by = "annotation_agrupada",
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_annotated.pdf"
)


# =============================================================================
# SECTION 12 — INTERACTIVE CELL-TYPE CURATION
# =============================================================================
# !!! WARNING: Run this section interactively, step by step.
# !!! Do NOT source the entire script with this section active.
#
# Purpose: subcluster populations that appear heterogeneous in the UMAP,
# inspect them, and reassign cells to the correct cell type manually.
#
# Step 1 → subcluster the heterogeneous types
# Step 2 → generate a composite inspection figure and save it to disk
# Step 3 → fill in the reassignment table (reassign) below
# Step 4 → apply corrections to the global object

output_dir <- dir_07
Idents(pbmc_harmony) <- "annotation_agrupada"

# ── Step 1. Subcluster ────────────────────────────────────────────────────────
meristemoid_umap   <- subclustar_tipo(pbmc_harmony, "Stomatal Line")
pavement_cell_umap <- subclustar_tipo(pbmc_harmony, "Pavement Cell")

# ── Step 2. Composite inspection figure (view, then fill in Step 3) ───────────
# All visual outputs for this step are assembled into one large PDF.
# Open 07_curation/subclustering_inspection.pdf, decide on the reassignments,
# then continue to Step 3.

p_meris_dim <- DimPlot(meristemoid_umap, group.by = "cluster_subtipo",
                       label = TRUE, raster = FALSE) +
  ggtitle("Stomatal Line \u2014 subclusters")

p_pave_dim  <- DimPlot(pavement_cell_umap, group.by = "cluster_subtipo",
                       label = TRUE, raster = FALSE) +
  ggtitle("Pavement Cell \u2014 subclusters")

marker_plots <- lapply(seq_len(nrow(marcadores)), function(i) {
  FeaturePlot(pavement_cell_umap, features = marcadores$gene[i]) +
    ggtitle(paste0(marcadores$cell.types[i], "\n", marcadores$gene[i])) +
    theme(plot.title = element_text(size = 8))
})

n_markers    <- length(marker_plots)
ncol_markers <- min(5L, n_markers)
nrow_markers <- ceiling(n_markers / ncol_markers)

composite_inspect <- (p_meris_dim | p_pave_dim) /
  wrap_plots(marker_plots, ncol = ncol_markers)

ggsave(
  file.path(output_dir, "subclustering_inspection.pdf"),
  composite_inspect,
  width     = max(20, ncol_markers * 4),
  height    = 10 + nrow_markers * 4,
  limitsize = FALSE
)
message("Subclustering inspection figure saved to 07_curation/subclustering_inspection.pdf")
message("Open it, decide on subcluster reassignments, fill in Step 3, then continue.")

# ── Step 3. Reassignment table ────────────────────────────────────────────────
# For each subclustered object: map subcluster IDs to final cell-type labels.
# Subcluster IDs come from $cluster_subtipo (values: "0", "1", "2", ...).
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
# SECTION 13 — EXPORT TO H5AD (Scanpy / Python)
# =============================================================================
# Exports the curated object to AnnData h5ad format for Python-based
# trajectory and velocity analyses (Scanpy, scFates, Palantir — all
# pre-installed in the Docker image).

output_dir <- dir_08

exportar_para_scanpy(pbmc_harmony,
                     file.path(output_dir, "pbmc_harmony_curated.h5ad"))

# To export a specific cell type:
# exportar_para_scanpy(
#   subset(pbmc_harmony, subset = celltype_reference_curated == "Guard Cell"),
#   file.path(output_dir, "GuardCell.h5ad")
# )


# =============================================================================
# ████████████████  PART 2 — PSEUDOBULK DE AND GO ENRICHMENT  █████████████████
# =============================================================================


# =============================================================================
# SECTION 14 — GLOBAL PSEUDOBULK AND REPLICATE CORRELATION
# =============================================================================
# Aggregates all cells per sample into a pseudobulk count vector and plots
# a Pearson correlation heatmap. Replicates from the same condition should
# cluster together with r > 0.95, confirming reproducibility before DE.

output_dir <- dir_09

pseudobulk <- generate_pseudobulk(pbmc_harmony, group_by = "orig.ident")

save_pdf(plot_replicate_correlation(pseudobulk$by_sample),
         "pseudobulk_correlation.pdf", w = 8, h = 8)


# =============================================================================
# SECTION 15 — PER-CELL-TYPE SUBSETS
# =============================================================================
# Splits the integrated object into one Seurat object per curated cell type
# for independent pseudobulk differential expression analyses.

output_dir <- dir_09

celular_subsets <- setNames(
  lapply(unique(pbmc_harmony$celltype_reference_curated),
         function(t) subset(pbmc_harmony, subset = celltype_reference_curated == t)),
  gsub("[^[:alnum:]_]", "_", unique(pbmc_harmony$celltype_reference_curated))
)


# =============================================================================
# SECTION 16 — PSEUDO-REPLICATES AND PSEUDOBULK COUNT MATRICES
# =============================================================================
# Because this experiment lacks within-condition biological replicates, cells
# are randomly assigned to n = 3 pseudo-replicates per condition. Counts are
# then summed per pseudo-replicate to build the matrix required by DESeq2.
# Cell types with fewer than 2 conditions are excluded automatically.

output_dir <- dir_09

celular_subsets_replicados <- Filter(Negate(is.null),
                                     lapply(celular_subsets, asignar_pseudoreplicados))
pseudobulk_list            <- lapply(celular_subsets_replicados, hacer_pseudobulk)

dir.create(file.path(output_dir, "pseudobulk_replicas"), recursive = TRUE, showWarnings = FALSE)
for (tipo in names(pseudobulk_list))
  write.csv(pseudobulk_list[[tipo]],
            file.path(output_dir, "pseudobulk_replicas", paste0("Pseudobulk_", tipo, ".csv")),
            row.names = TRUE)


# =============================================================================
# SECTION 17 — DESEQ2 DIFFERENTIAL EXPRESSION
# =============================================================================
# Pairwise contrasts are run for all condition comparisons.
# Results are written as CSV files under 10_deseq2/<tag>/
#
# ┌─ DEFINE YOUR COMPARISONS ───────────────────────────────────────────────────
#   conds : c("reference_condition", "treatment_condition")
#   tag   : short label used for output file names
# └─────────────────────────────────────────────────────────────────────────────
comparaciones <- list(
  list(conds = c("0.5N", "5N"),   tag = "05_5"),
  list(conds = c("0N",   "5N"),   tag = "0_5"),
  list(conds = c("0N",   "0.5N"), tag = "0_05")
)

output_dir <- dir_10

for (tag in sapply(comparaciones, `[[`, "tag"))
  dir.create(file.path(output_dir, tag), recursive = TRUE, showWarnings = FALSE)

for (tipo in names(pseudobulk_list))
  correr_deseq2(as.matrix(pseudobulk_list[[tipo]]), comparaciones,
                output_dir = output_dir, tipo = tipo)


# =============================================================================
# SECTION 18 — VOLCANO PLOTS
# =============================================================================
# One volcano plot per cell type per comparison; two plots per PDF page.
# Points are colored by significance (padj < 0.05, |log2FC| > 1).

output_dir <- dir_11

for (comp in comparaciones) {
  tag       <- comp$tag
  csv_files <- list.files(file.path(dir_10, tag),
                          pattern = "\\.csv$", full.names = TRUE)
  if (!length(csv_files)) { message("No CSV files found for comparison: ", tag); next }

  pdf(file.path(output_dir, paste0("VolcanoPlots_", tag, ".pdf")),
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
# SECTION 19 — DIFFERENTIAL GENE TABLES AND HEATMAPS
# =============================================================================
# Genes are classified (up = 1 / down = -1 / unchanged = 0) and summarized
# in cross-cell-type tables. A hierarchically clustered log2FC heatmap is
# produced per comparison.

output_dir <- dir_12

for (comp in comparaciones) {
  tag       <- comp$tag
  csv_dir   <- file.path(dir_10, tag)
  diff_dir  <- file.path(output_dir, tag)
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
# SECTION 20 — GO ENRICHMENT ANALYSIS
# =============================================================================
# Gene Ontology enrichment (clusterProfiler::enrichGO) is run per cell type
# and comparison. Both full and semantically simplified term sets are
# generated and visualized as balloon plots.
#
# ┌─ CHANGE FOR YOUR ORGANISM ──────────────────────────────────────────────────
#   Arabidopsis thaliana : orgdb = org.At.tair.db, keytype = "TAIR"
#   Homo sapiens         : orgdb = org.Hs.eg.db,   keytype = "ENSEMBL"
#   Mus musculus         : orgdb = org.Mm.eg.db,   keytype = "ENSEMBL"
# └─────────────────────────────────────────────────────────────────────────────
orgdb      <- org.At.tair.db
keytype    <- "TAIR"
espacio    <- "BP"    # "BP" = Biological Process | "MF" = Molecular Function | "CC" = Cellular Component
qval       <- 0.05
nivel_poda <- 6       # Maximum GO hierarchy depth for term pruning

output_dir <- dir_13

universo <- keys(orgdb, keytype = keytype)

for (comp in comparaciones) {
  tag        <- comp$tag
  diff_dir   <- file.path(dir_12, tag)
  enr_dir    <- file.path(output_dir, tag)
  tabla_path <- file.path(diff_dir, "tabla_diferenciales.tsv")

  dir.create(enr_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(tabla_path)) { message("No differential table for: ", tag); next }

  tabla <- read.table(tabla_path, header = TRUE, row.names = 1, sep = "\t")
  tabla <- tabla[, colSums(tabla != 0) > 0, drop = FALSE]

  go_total  <- correr_enriquecimiento_go(tabla, universo, espacio,
                                         orgdb = orgdb, keytype = keytype,
                                         simplificar = FALSE, output_dir = enr_dir)
  go_simple <- correr_enriquecimiento_go(tabla, universo, espacio,
                                         orgdb = orgdb, keytype = keytype,
                                         simplificar = TRUE,  output_dir = enr_dir)

  go_total_podado  <- podar_go(go_total,  nivel_poda, espacio, qval,
                               simplificar = FALSE, output_dir = enr_dir)
  go_simple_podado <- podar_go(go_simple, nivel_poda, espacio, qval,
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
# END OF PIPELINE
# =============================================================================
