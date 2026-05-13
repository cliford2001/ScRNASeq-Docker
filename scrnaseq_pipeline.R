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
PIPELINE_DIR <- "~/projects2/eleo/ScRNA/metodologia/ScRNASeq-Docker/"

# Root directory for your project data and results.
# All result files will be written to DATA_DIR/results/<step>/
DATA_DIR   <- "~/projects2/eleo/ScRNA/"
base_dir   <- file.path(DATA_DIR, "metodologia/resultados")

# ── Input format ──────────────────────────────────────────────────────────────
# USE_CELLBENDER = TRUE  → load CellBender-filtered HDF5 files (recommended)
# USE_CELLBENDER = FALSE → load CellRanger filtered_feature_bc_matrix/ directly
#                          (use this if you skipped the CellBender step)
#
# If TRUE  → set samples$file to the .h5 file path relative to DATA_DIR
# If FALSE → set samples$file to the filtered_feature_bc_matrix/ directory
#            path relative to DATA_DIR
USE_CELLBENDER <- FALSE

# ── Sample manifest ───────────────────────────────────────────────────────────
# Add one entry per sample. Each entry needs:
#   file      — path to the input file or directory (relative to DATA_DIR)
#   label     — unique name for this sample (appears in all plots)
#   condition — experimental group this sample belongs to

# ── OPTION 1: CellBender-filtered HDF5 files (USE_CELLBENDER = TRUE) ─────────
# samples <- list(
#   list(file = "cellbender/Sample_0N_cellbender_filtered.h5",      label = "0N",      condition = "0N"),
#   list(file = "cellbender/Sample_05N_R1_cellbender_filtered.h5",  label = "0.5N_R1", condition = "0.5N"),
#   list(file = "cellbender/Sample_05N_2_cellbender_filtered.h5",   label = "0.5N_R2", condition = "0.5N"),
#   list(file = "cellbender/Sample_5N_R1_cellbender_filtered.h5",   label = "5N_R1",   condition = "5N"),
#   list(file = "cellbender/Sample_5N_2_cellbender_filtered.h5",    label = "5N_R2",   condition = "5N")
# )

# ── OPTION 2: CellRanger filtered_feature_bc_matrix (USE_CELLBENDER = FALSE) ─
samples <- list(
  list(file = "cellranger/Sample_0N/outs/filtered_feature_bc_matrix",      label = "0N",      condition = "0N"),
  list(file = "cellranger/Sample_05N/outs/filtered_feature_bc_matrix",     label = "0.5N_R1", condition = "0.5N"),
  list(file = "cellranger/Sample_05N_2/outs/filtered_feature_bc_matrix",   label = "0.5N_R2", condition = "0.5N"),
  list(file = "cellranger/Sample_5N/outs/filtered_feature_bc_matrix",      label = "5N_R1",   condition = "5N"),
  list(file = "cellranger/Sample_5N_2/outs/filtered_feature_bc_matrix",    label = "5N_R2",   condition = "5N")
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
options(Seurat.allow.s4 = FALSE)  # required for Seurat 5 compatibility
setwd(DATA_DIR)

# ── Create per-step output directories ────────────────────────────────────────
list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)  # creates output folders and loads their paths as variables

# output_dir is the global variable used by save_pdf / save_qc / save_vln helpers.
# It is reassigned at the start of each section to the appropriate step directory.
output_dir <- base_dir


# =============================================================================
# SECTION 0 — PIPELINE WORKFLOW FIGURE
# =============================================================================
# Generates a visual overview of the full pipeline saved to 01_qc/.
# Run this section once immediately after initialization.

plot_pipeline_workflow(file.path(dir_01, "pipeline_workflow.pdf"))


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
output_dir <- dir_01

mt_pattern <- "^ATMG"  # Arabidopsis mitochondrial genes
cp_pattern <- "^ATCG"  # Arabidopsis chloroplast genes


# Load all samples using the helper function
seurat_list_raw <- load_seurat_samples(samples = samples,
                                       DATA_DIR = DATA_DIR,
                                       USE_CELLBENDER = USE_CELLBENDER,
                                       mt_pattern = mt_pattern,
                                       cp_pattern = cp_pattern)

plot_qc_batch(seurat_list_raw, colors, "qc_prefilter.pdf")


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
output_dir <- dir_01   # both pre- and post-filter QC plots go to 01_qc/

seurat_list <- filter_seurat_samples(seurat_list_raw, min_features = 200, max_mt = 5)

plot_qc_batch(seurat_list, colors, "qc_postfilter.pdf")

# Checkpoint — restore with: seurat_list <- readRDS(file.path(dir_objects, "seurat_list_postfilter.rds"))
saveRDS(seurat_list, file.path(dir_objects, "seurat_list_postfilter.rds"))


# =============================================================================
# SECTION 3 — MERGE AND INITIAL PREPROCESSING
# =============================================================================
# Filtered samples are merged and preprocessed: log-normalization, variable
# feature selection (VST, 2,000 features), scaling, PCA (30 PCs), and UMAP.
# The resulting UMAP shows batch effects before integration.

output_dir <- dir_01

pbmc_harmony <- reduce(seurat_list, merge) %>%  # merge all samples into one object
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_preharmony.pdf")

# Checkpoint — restore with: pbmc_harmony <- readRDS(file.path(dir_objects, "pbmc_harmony_preharmony.rds"))
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_preharmony.rds"))


# =============================================================================
# SECTION 4 — HARMONY BATCH CORRECTION
# =============================================================================
# Harmony adjusts cell embeddings to remove sample-level batch effects while
# preserving biological variation. All downstream steps use the "harmony"
# reduction instead of "pca".

output_dir <- dir_01

pbmc_harmony <- pbmc_harmony %>%
  RunHarmony("orig.ident", plot_convergence = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE)

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_postharmony.pdf")

# Checkpoint — restore with: pbmc_harmony <- readRDS(file.path(dir_objects, "pbmc_harmony_postharmony.rds"))
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_postharmony.rds"))


# =============================================================================
# SECTION 5 — RESOLUTION OPTIMIZATION
# =============================================================================
# Two diagnostics guide the choice of clustering resolution:
#   (a) Elbow plot — k-means within-cluster sum of squares across k values.
#       The inflection point suggests the number of major cell types.
#   (b) Clustree  — tracks cluster stability across Leiden resolutions.
#       Choose the lowest resolution where clusters stop merging.
#
# ┌─ PARAMETERS ────────────────────────────────────────────────────────────────
#   k_range          : k values tested in the elbow plot
#   resolutions_test : Leiden resolutions swept by clustree
#   → Inspect elbow_plot.pdf and clustree.pdf before setting cluster_resolution
#     in Section 6.
# └─────────────────────────────────────────────────────────────────────────────
k_range          <- 1:31
resolutions_test <- c(0.15, 0.30, 0.50, 0.8, 1.0)

output_dir <- dir_02

# ── 5a. Elbow plot ────────────────────────────────────────────────────────────
pca_data <- Embeddings(pbmc_harmony, "pca")[, 1:30]
wss      <- sapply(k_range, function(k) kmeans(pca_data, centers = k, nstart = 4)$tot.withinss)

elbow_plot <- ggplot(data.frame(k = k_range, wss = wss), aes(k, wss)) +
  geom_line() + geom_point() +
  labs(x = "Number of clusters (k)", y = "Within-cluster sum of squares") +
  theme_minimal()

save_pdf(elbow_plot, "elbow_plot.pdf", w = 8, h = 6)

# ── 5b. Clustree ──────────────────────────────────────────────────────────────
clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30,
                k.param = 30, verbose = FALSE)

for (res in resolutions_test)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(clustree(clu, prefix = "RNA_snn_res."), "clustree.pdf", w = 14, h = 14)


# =============================================================================
# SECTION 6 — FINAL CLUSTERING
# =============================================================================
# Apply the selected resolution for the final cluster assignment.
# After clustering, a UMAP coloured by Seurat cluster and a bar chart of
# cells per sample are saved to 02_clustering/.
#
# ┌─ SET RESOLUTION AFTER INSPECTING elbow_plot.pdf AND clustree.pdf ──────────
#   cluster_resolution : Leiden resolution for final clustering (default 0.3)
# └─────────────────────────────────────────────────────────────────────────────
cluster_resolution <- 0.3
output_dir <- dir_02

# clu (Section 5) was temporary — re-run on pbmc_harmony to embed the final clusters
pbmc_harmony <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30,
                k.param = 30, verbose = FALSE) %>%
  FindClusters(resolution = cluster_resolution, algorithm = 4, verbose = FALSE)

Idents(pbmc_harmony) <- "seurat_clusters"
save_pdf(DimPlot(pbmc_harmony, group.by = "seurat_clusters", label = TRUE),
         "umap_seuratclusters.pdf")


# =============================================================================
# SECTION 7 — DOTPLOT: MARKER GENES BY CLUSTER (pre-annotation guide)
# =============================================================================
# Before assigning cell-type labels, this plot helps you identify which
# numbered Seurat cluster corresponds to which cell type by showing the
# expression of bibliography-derived marker genes across all clusters.
# Clusters that strongly express a known marker (e.g., AT5G26000 for Guard
# Cell) should be labelled as that cell type in Section 8.
# Dot size = fraction of expressing cells; color = mean expression level.
output_dir <- dir_03

biblio_marks_file <- file.path(DATA_DIR, "metodologia/biblio_marks.txt")
marker_table      <- read.table(biblio_marks_file, header = TRUE, sep = "\t", quote = "")

plot_marker_dotplot(
  pbmc_harmony,
  marker_table,
  annot_col = "seurat_clusters",
  outfile   = file.path(output_dir, "dotplot_marker_table_preannotation.pdf"),
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

output_dir <- dir_03

# ── 8a. Bibliography-based annotation ─────────────────────────────────────────
markers <- find_markers(pbmc_harmony,
                        output_file = file.path(output_dir, "FindAllMarkers.tsv"))

pbmc_harmony <- annotate_by_markers(pbmc_harmony, markers,
                                    reference_file = biblio_marks_file)
# Annotation stored in: pbmc_harmony$celltype

plot_marker_dotplot(
  pbmc_harmony,
  marker_table,
  annot_col = "celltype", # uses the newly assigned annotation column
  outfile   = file.path(output_dir, "dotplot_marker_table_annotation_biblio.pdf"),
  width = 20, height = 10
)

# ── 8b. Reference-based annotation ────────────────────────────────────────────
reference_obj <- readRDS(file.path(DATA_DIR, "metodologia/GSE273033_seuratObj_for_publication.rds"))
pbmc_harmony <- annotate_by_reference(pbmc_harmony,
                                      reference_obj = reference_obj,
                                      reference_col = "annotation")

plot_marker_dotplot(
  pbmc_harmony,
  marker_table,
  annot_col = "celltype_reference", # uses the newly assigned annotation column
  outfile   = file.path(output_dir, "dotplot_marker_table_annotation_reference.pdf"),
  width = 20, height = 10
)
# Annotation stored in: pbmc_harmony$celltype_reference


# =============================================================================
# SECTION 9 — ANNOTATED CLUSTREE
# =============================================================================
# Re-runs the resolution sweep with cell-type labels overlaid on each node.
# Confirms that the chosen resolution cleanly separates known cell types.

output_dir <- dir_03

# Mode: returns the most frequent value in a vector
Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30,
                k.param = 30, verbose = FALSE)

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

output_dir <- dir_04

# JoinLayers is required in Seurat 5 before subsetting after merge
pbmc_harmony <- JoinLayers(pbmc_harmony)
Idents(pbmc_harmony) <- "celltype_reference"

n_genes <- length(genes_of_interest)

# ── 10a. All cell types ───────────────────────────────────────────────────────
save_vln(VlnPlot(pbmc_harmony, features = gene),                  "vln_gene_all.pdf")
save_pdf(FeaturePlot(pbmc_harmony, features = gene),              "feature_gene_all.pdf")
save_vln(VlnPlot(pbmc_harmony, features = genes_of_interest),     "vln_geneset_all.pdf",     n = n_genes)
save_pdf(FeaturePlot(pbmc_harmony, features = genes_of_interest), "feature_geneset_all.pdf", h = 8 * n_genes)

# ── 10b. Cell type of interest ────────────────────────────────────────────────
sub_obj <- subset(pbmc_harmony, idents = celltype)

save_vln(VlnPlot(sub_obj, features = gene),                  "vln_gene_celltype.pdf")
save_pdf(FeaturePlot(sub_obj, features = gene),              "feature_gene_celltype.pdf")
save_vln(VlnPlot(sub_obj, features = genes_of_interest),     "vln_geneset_celltype.pdf",     n = n_genes)
save_pdf(FeaturePlot(sub_obj, features = genes_of_interest), "feature_geneset_celltype.pdf", h = 8 * n_genes)


# =============================================================================
# SECTION 11 — CELL-TYPE GROUPING  [OPTIONAL]
# =============================================================================
# Fine-grained labels are collapsed into broader categories for downstream
# analyses. Cell types NOT listed in 'grouping' keep their original label.
# Skip this section if you do not need to merge cell types.
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

output_dir <- dir_05

# !!! unpacks the grouping vector as named arguments to recode()
pbmc_harmony$celltype_grouped <- recode(pbmc_harmony$celltype_reference, !!!grouping)

save_pdf(
  DimPlot(pbmc_harmony, group.by = "celltype_grouped",  # grouped broad cell types
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_grouped.pdf"
)


# =============================================================================
# SECTION 12 — INTERACTIVE CELL-TYPE CURATION  [OPTIONAL]
# =============================================================================
# !!! WARNING: Run this section interactively, step by step.
# !!! Do NOT source the entire script with this section active.
# Skip this section if you are satisfied with the annotation from Section 8.
#
# Purpose: subcluster populations that appear heterogeneous in the UMAP,
# inspect them, and reassign cells to the correct cell type manually.
#
# Step 1 → subcluster the heterogeneous types
# Step 2 → save inspection figures; open them and decide on reassignments
# Step 3 → fill in the reassignment table below
# Step 4 → apply corrections to the global object

output_dir   <- dir_05
curation_col <- "celltype_grouped"   # starting annotation column for curation
Idents(pbmc_harmony) <- curation_col
table(pbmc_harmony[[curation_col]])

# ── Step 1. Subcluster ────────────────────────────────────────────
mesophyll_umap     <- subcluster_cell_type(pbmc_harmony, "Mesophyll",     annot_col = curation_col)
pavement_cell_umap <- subcluster_cell_type(pbmc_harmony, "Pavement Cell", annot_col = curation_col)

# Check how many subclusters each type produced:
# table(mesophyll_umap$cluster_subtipo)
# table(pavement_cell_umap$cluster_subtipo)

# ── Step 2. Inspection figures ────────────────────────────────────
# Each call creates the DimPlot, saves it as PDF, and returns it for the composite
p_meso_dim <- plot_subcluster_umap(mesophyll_umap,     "Mesophyll",     output_dir)
p_pave_dim <- plot_subcluster_umap(pavement_cell_umap, "Pavement Cell", output_dir)

# Composite: each row = [ UMAP | marker genes ] for one cell type
save_subcluster_composite(
  subcluster_list = list(
    list(umap_plot = p_meso_dim, obj = mesophyll_umap),
    list(umap_plot = p_pave_dim, obj = pavement_cell_umap)
  ),
  marker_table = marker_table,
  output_dir   = output_dir
)

# ── Step 3. Reassignment table ────────────────────────────────────
# Map subcluster IDs \u2192 final cell-type labels.
# "others" is a catch-all for any subcluster ID not listed.
# Names must match the variable names used in Step 1 exactly.
reassign <- list(
  mesophyll_umap = c(
    "0"      = "Mesophyll",
    "1"      = "Mesophyll",
    "2"      = "Mesophyll",
    "others" = "Mesophyll"
  ),
  pavement_cell_umap = c(
    "0"      = "Pavement Cell",
    "1"      = "Pavement Cell",
    "2"      = "Pavement Cell",
    "3"      = "Pavement Cell",
    "4"      = "Pavement Cell",
    "others" = "Pavement Cell"
  )
)

# ── Step 4. Apply corrections ───────────────────────────────────────────────
subcluster_list <- list(
  mesophyll_umap     = mesophyll_umap,
  pavement_cell_umap = pavement_cell_umap
)

pbmc_harmony <- apply_subcluster_reassignment(
  obj             = pbmc_harmony,
  subcluster_list = subcluster_list,
  reassign        = reassign,
  source_col      = curation_col,
  dest_col        = "celltype_curated"
)

save_pdf(
  DimPlot(pbmc_harmony, group.by = "celltype_curated",
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_curated.pdf"
)

# Checkpoint — restore with: pbmc_harmony <- readRDS(file.path(dir_objects, "pbmc_harmony_curated.rds"))
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_curated.rds"))

# =============================================================================
# SECTION 13 — EXPORT TO H5AD (Scanpy / Python)
# =============================================================================
# Exports the curated object to AnnData h5ad format for Python-based
# trajectory and velocity analyses (Scanpy, scFates, Palantir — all
# pre-installed in the Docker image).

export_to_scanpy(pbmc_harmony,
                 file.path(dir_objects, "pbmc_harmony_curated.h5ad"))

# To export a specific cell type:
# export_to_scanpy(
#   subset(pbmc_harmony, subset = celltype_curated == "Guard Cell"),
#   file.path(dir_objects, "GuardCell.h5ad")
# )


# =============================================================================
# ████████████████████████  PART 2 — PSEUDOBULK DE & GO  ██████████████████████
# =============================================================================


# =============================================================================
# SECTION 14 — CELL-TYPE SUBSETS
# =============================================================================
# Creates one Seurat subset per curated cell type for downstream pseudobulk
# analysis. Object names are sanitised so they can be used safely as list names
# and output filenames.
#
# ┌─ SET THE ANNOTATION COLUMN TO USE FOR PART 2 ───────────────────────────────
#   pseudobulk_annot_col : metadata column containing the final cell-type labels
# └─────────────────────────────────────────────────────────────────────────────
pseudobulk_annot_col <- "celltype_curated"

# Create cell-type subsets for pseudobulk analysis
cell_type_subsets <- create_cell_type_subsets(pbmc_harmony, annot_col = pseudobulk_annot_col)


# =============================================================================
# SECTION 15 — PSEUDO-REPLICATE ASSIGNMENT
# =============================================================================
# Assigns random pseudo-replicates within each condition for every cell type.
# Only subsets containing at least two conditions are kept for Part 2.
#
# ┌─ PSEUDOBULK PARAMETERS ──────────────────────────────────────────────────────
#   pseudobulk_conditions : optional condition subset to retain (NULL = all)
#                           Examples:
#                             NULL              → use all conditions (0N, 0.5N, 5N)
#                             c("0N", "0.5N")   → use only 0N and 0.5N
#                             "5N"              → use only 5N
#   n_pseudoreps          : number of pseudo-replicates per condition (per cell type)
# └─────────────────────────────────────────────────────────────────────────────

# Use all conditions
#pseudobulk_conditions <- NULL
# Uncomment below to use only specific conditions:
# pseudobulk_conditions <- c("0N", "0.5N")      # Compare control vs low nitrogen
#pseudobulk_conditions <- c("0.5N", "5N")      # Compare nitrogen treatments
# pseudobulk_conditions <- "5N"                  # Single condition (rarely useful)

n_pseudoreps <- 3

# Assign pseudo-replicates (uses global random seed set in INITIALIZATION)
cell_type_subsets_replicates <- assign_pseudoreplicates_batch(cell_type_subsets,
                                                             pseudobulk_conditions = pseudobulk_conditions,
                                                             n_pseudoreps = n_pseudoreps)


# Example QC checks for one subset:
table(cell_type_subsets_replicates$Pavement_Cell$replicate)
table(cell_type_subsets_replicates$Pavement_Cell$orig.ident)


# =============================================================================
# SECTION 16 — PSEUDOBULK TABLES AND DESEQ2
# =============================================================================
# Aggregates counts by pseudo-replicate for each cell type, saves the resulting
# count tables, and runs DESeq2 for the user-defined pairwise contrasts.
#
# ┌─ DEFINE YOUR CONDITION CONTRASTS HERE ──────────────────────────────────────
#   comparaciones : each entry must contain
#                   conds = c("reference", "treatment")
#                   tag   = output folder / file label
# └─────────────────────────────────────────────────────────────────────────────
comparaciones <- list(
  list(conds = c("0.5N", "5N"), tag = "0.5N_vs_5N"),
  list(conds = c("0N",   "5N"), tag = "0N_vs_5N")
)

output_dir <- dir_06

# ┌─ SELECT WHICH CELL TYPES TO ANALYZE ────────────────────────────────────────
#   NULL = analyze all cell types
#   c("Epidermis", "Cortex") = analyze only these
# └─────────────────────────────────────────────────────────────────────────────
cell_types_to_analyze <- NULL  # Change to c("Epidermis", "Cortex") to filter

# Run pseudobulk aggregation and DESeq2 analysis
deseq2_results <- run_pseudobulk_deseq2_analysis(
  cell_type_subsets_replicates = cell_type_subsets_replicates,
  comparisons = comparaciones,
  output_dir = output_dir,
  cell_types = cell_types_to_analyze,
  pseudobulk_dir = file.path(dir_objects, "pseudobulk_replicas")
)


# =============================================================================
# SECTION 17 — VOLCANO PLOTS
# =============================================================================
# Renders one PNG volcano plot per cell type for a selected contrast and also
# combines them into a single PDF.
#
# ┌─ VOLCANO PARAMETERS ─────────────────────────────────────────────────────────
#   volcano_tag : tag of the contrast to visualize (must match comparaciones)
#   padj_cut    : adjusted p-value threshold
#   lfc_cut     : absolute log2 fold-change threshold
# └─────────────────────────────────────────────────────────────────────────────
volcano_tag <- "0.5N_vs_5N"
padj_cut    <- 0.05
lfc_cut     <- 1

render_volcano_plots(
  results_dir = file.path(dir_06, volcano_tag),
  output_dir  = file.path(dir_06, volcano_tag, "volcano"),
  pdf_name    = paste0("VolcanoPlots_", volcano_tag, ".pdf"),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut
)


# =============================================================================
# SECTION 18 — DIFFERENTIAL GENE TABLES
# =============================================================================
# Builds combined differential-expression summary tables for one selected
# contrast across all cell types.
#
# ┌─ DIFFERENTIAL TABLE PARAMETERS ──────────────────────────────────────────────
#   diff_tag    : tag of the contrast to summarize
#   diff_prefix : output filename prefix for the discrete matrix
# └─────────────────────────────────────────────────────────────────────────────
diff_tag    <- volcano_tag
diff_prefix <- paste0("tabla_diferenciales_", diff_tag)

diff_tables <- build_differential_tables(
  results_dir = file.path(dir_06, diff_tag),
  output_dir  = file.path(dir_06, diff_tag),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut,
  prefix      = diff_prefix
)


# =============================================================================
# SECTION 19 — GO ENRICHMENT (SIMPLE)
# =============================================================================
# Gene Ontology enrichment per cell type for the selected contrast.

go_space    <- "BP"          # Change to "MF" or "CC" if desired
padj_cutoff <- 0.05

deseq2_files <- list.files(file.path(dir_06, diff_tag),
                           pattern = "^DESeq2_.*\\.csv$",
                           full.names = TRUE)

for (deseq2_file in deseq2_files) {
  cell_type <- gsub("^DESeq2_|\\.csv$", "", basename(deseq2_file))

  deseq2_results <- read.csv(deseq2_file, row.names = 1)
  sig_genes <- rownames(deseq2_results)[deseq2_results$padj < padj_cutoff]

  if (length(sig_genes) > 0) {
    run_simple_go_enrichment(
      diff_table = data.frame(gene_id = sig_genes),
      output_dir = file.path(dir_07, diff_tag),
      orgdb = org.At.tair.db,
      keytype = "TAIR",
      go_space = go_space,
      padj_cutoff = padj_cutoff,
      cell_type = cell_type,
      contrast_tag = diff_tag
    )
  }
}

# =============================================================================
# SECTION 20 — LOG2FC HEATMAP + CLUSTERING
# =============================================================================
# Heatmap of log2FC values across all cell types for the selected contrast.
#
# ┌─ PARAMETERS ─────────────────────────────────────────────────────────────────
#   CLUSTER_METHOD : "hclust" — hierarchical (euclidean + complete + cutreeDynamic)
#                    "wgcna"  — coexpression network (TOM + mergeCloseModules)
#   heatmap_limits : color scale range
#   wgcna_merge_cut: merge WGCNA modules with correlation > (1 - wgcna_merge_cut)
# └─────────────────────────────────────────────────────────────────────────────
CLUSTER_METHOD  <- "wgcna"   # "hclust" or "wgcna"
heatmap_limits  <- c(-5, 5)
wgcna_merge_cut <- 0.25

heatmap_results <- build_logfc_heatmap(
  logfc_table  = diff_tables$logfc,
  contrast_tag = diff_tag,
  output_dir   = file.path(dir_06, diff_tag),
  method       = CLUSTER_METHOD,
  limits       = heatmap_limits,
  merge_cut    = wgcna_merge_cut
)


# =============================================================================
# SECTION 21 — GO ENRICHMENT PER CLUSTER
# =============================================================================
# Runs GO enrichment for each cluster identified in Section 20.
# Uses the cluster assignments from heatmap_results.

go_clusters_padj <- 0.05

for (clust_id in unique(heatmap_results$cluster)) {
  genes <- heatmap_results$gene_id[heatmap_results$cluster == clust_id]

  run_simple_go_enrichment(
    diff_table   = data.frame(gene_id = genes),
    output_dir   = file.path(dir_06, diff_tag, paste0("GO_clusters_", CLUSTER_METHOD)),
    orgdb        = org.At.tair.db,
    keytype      = "TAIR",
    go_space     = "BP",
    padj_cutoff  = go_clusters_padj,
    cell_type    = as.character(clust_id),
    contrast_tag = diff_tag
  )
}


# =============================================================================
# SECTION 22 — NETWORK INFERENCE PER CLUSTER (3-METHOD MIX STRATEGY)
# =============================================================================
# THREE COMPLEMENTARY network inference analyses per cluster from Section 20:
# (Analogous to Sección 20's dual-clustering strategy: hclust vs WGCNA)
#
# ┌─ THE MIX STRATEGY ────────────────────────────────────────────────────────────
#   We do NOT choose between GENIE3, WGCNA, and SYNERGY. Instead, we run all 3
#   because they answer different questions:
#
#   • GENIE3 (directed, like hclust geometry):
#       Finds TF → target edges with predictive power.
#       Filter: Pearson |r| ≥ 0.90 (avoids noise via correlation).
#       ✓ Strength: directionality, interpretable as causality.
#       ⚠ Weakness: only TFs regulate; less robust to outliers.
#
#   • WGCNA (undirected, like WGCNA coexpression):
#       Finds coexpressed gene pairs via TOM (Topological Overlap).
#       Filter: TOM ≥ 0.15 (genes that share ≥15% of neighbors).
#       ✓ Strength: coexpression robustness; ANY gene can regulate ANY gene.
#       ⚠ Weakness: no directionality; local structure matters, not global.
#
#   • SYNERGY (mix, high-confidence):
#       Combines GENIE3 directionality + WGCNA coexpression validation.
#       Filter: Pearson ≥0.90 AND TOM ≥0.15 (both layers must pass).
#       Score: geometric mean of rank-normalized GENIE3 × TOM.
#       ✓ Strength: high-confidence TF→target edges backed by coexpression.
#       ⚠ Weakness: very restrictive; fewer edges (top tier only).
#
# ┌─ WHAT IS TOM? (Topological Overlap Matrix) ──────────────────────────────────
#   TOM measures "robustness" of a gene pair's connection by counting shared
#   neighbors: if gene A and gene B both correlate with genes {C, D, E, ...},
#   they are truly connected, not by chance.
#
#   Formula:   TOM(A,B) = (# shared neighbors) / min(neighbors(A), neighbors(B))
#   Range:     0 to 1 (0 = no shared neighbors, 1 = identical neighborhoods)
#   Severity:  TOM ≥ 0.25 is "strict", 0.15 is "moderate", 0.08 is "loose"
#   Your threshold (0.15) = top 5% of edge confidences = MODERATE
#
# Both functions read pseudobulk replicate counts from Section 16, normalize
# (CPM + log2) and run independently. Outputs are saved in dir_08/<contrast>/.
#
# ┌─ COMMON PARAMETERS ──────────────────────────────────────────────────────────
#   n_top_clusters : how many largest clusters to analyze
#   min_var_filter : drop genes with variance below this across samples
# ├─ GENIE3 PARAMETERS ──────────────────────────────────────────────────────────
#   net_orgdb      : Bioconductor OrgDb (e.g. org.At.tair.db, org.Hs.eg.db)
#   net_keytype    : key type matching gene IDs (e.g. "TAIR", "ENSEMBL")
#   custom_tfs     : optional vector of TF IDs to override GO-based detection
#   cor_min        : Pearson |r| correlation threshold (≥0.90 = MODERATE)
#   genie3_ntrees  : Random Forest trees (more = stabler, slower)
#   n_cores        : parallel cores for GENIE3
# ├─ WGCNA PARAMETERS ───────────────────────────────────────────────────────────
#   soft_power     : power for adjacency (default 6; higher = fewer edges)
#   network_type   : "signed" (correlation direction matters) or "unsigned"
#   tom_threshold  : TOM threshold (≥0.15 = MODERATE, equivalent to Pearson 0.90)
# └─────────────────────────────────────────────────────────────────────────────
# ┌─ CHOOSE WHICH METHODS TO RUN ────────────────────────────────────────────────
#   Edit the line below to run only desired methods. Options:
#   c("GENIE3", "WGCNA", "SYNERGY")  — run all 3
#   c("GENIE3", "SYNERGY")            — skip WGCNA
#   c("SYNERGY")                      — only high-confidence
#   c("GENIE3", "WGCNA")             — skip SYNERGY
# └─────────────────────────────────────────────────────────────────────────────
network_methods <- c("GENIE3", "WGCNA", "SYNERGY")  # CHANGE AS NEEDED

# ── Common parameters ────────────────────────────────────────────────────────
n_top_clusters <- 3
min_var_filter <- 0.01

# GENIE3 parameters
net_orgdb     <- org.At.tair.db
net_keytype   <- "TAIR"
custom_tfs    <- NULL
cor_min       <- 0.75        # EXPLORATORY — top 15% (Pearson |r| >= 0.75)
genie3_ntrees <- 100
n_cores       <- 4

# WGCNA parameters
soft_power    <- 6
network_type  <- "signed"
tom_threshold <- 0.05        # EXPLORATORY — top 15% (TOM >= 0.05)

# ── Run all selected methods in one call ─────────────────────────────────────
net_pipeline <- run_network_inference_pipeline(
  heatmap_results      = heatmap_results,
  pseudobulk_dir       = file.path(dir_objects, "pseudobulk_replicas"),
  output_base_dir      = file.path(dir_08, diff_tag),
  methods              = network_methods,
  orgdb                = net_orgdb,
  keytype              = net_keytype,
  custom_tfs           = custom_tfs,
  cor_min              = cor_min,
  genie3_ntrees        = genie3_ntrees,
  n_cores              = n_cores,
  soft_power           = soft_power,
  network_type         = network_type,
  tom_threshold        = tom_threshold,
  n_top_clusters       = n_top_clusters,
  min_var_filter       = min_var_filter
)

# ── Extract results for downstream sections ─────────────────────────────────
genie3_results  <- net_pipeline$results$GENIE3
wgcna_results   <- net_pipeline$results$WGCNA
synergy_results <- net_pipeline$results$SYNERGY


# =============================================================================
# SECTION 22B — THRESHOLD TESTING (OPTIONAL)
# =============================================================================
# Test 5 different threshold combinations to find optimal settings.
# Generates PDF comparing edge counts and cluster coverage.
#
# ┌─ UNCOMMENT TO RUN ────────────────────────────────────────────────────────
#   Set RUN_THRESHOLD_TEST <- TRUE below to execute
# └─────────────────────────────────────────────────────────────────────────────

RUN_THRESHOLD_TEST <- FALSE  # Set to TRUE to run threshold exploration

if (RUN_THRESHOLD_TEST) {
  message("\n🔬 TESTING NETWORK THRESHOLDS (5 combinations)...")

  threshold_test <- test_network_thresholds(
    heatmap_results = heatmap_results,
    pseudobulk_dir  = file.path(dir_objects, "pseudobulk_replicas"),
    output_dir      = file.path(dir_08, diff_tag, "THRESHOLD_TEST"),
    method          = "SYNERGY",  # Change to "GENIE3" or "WGCNA" if desired
    orgdb           = net_orgdb,
    keytype         = net_keytype,
    custom_tfs      = custom_tfs,
    genie3_ntrees   = genie3_ntrees,
    n_cores         = n_cores,
    soft_power      = soft_power,
    network_type    = network_type,
    n_top_clusters  = n_top_clusters,
    min_var_filter  = min_var_filter
  )

  message("\n✓ Threshold test complete.")
  message("  PDF: ", threshold_test$pdf)
  message("  Recommendation: ", threshold_test$recommendation)
}


# =============================================================================
# SECTION 23 — NETWORK VISUALIZATION (FORCE-DIRECTED LAYOUT)
# =============================================================================
# Clean network visualization using force-directed (Fruchterman-Reingold) layout.
# Choose ONE method below; visualizes all filtered edges with igraph.
#
# ┌─ CHOOSE METHOD ─────────────────────────────────────────────────────────────
#   Selected method will be visualized in detail (force-directed, node sizes by
#   degree, edge widths by weight). This complements SEC 22's PDF summaries.
#
#   Recommendations:
#   • GENIE3  → see TF directionality visually
#   • WGCNA   → see coexpression structure and modules
#   • SYNERGY → highest confidence edges (TF→target validated by coexpression)
# └─────────────────────────────────────────────────────────────────────────────

viz_method <- "SYNERGY"   # "GENIE3", "WGCNA", or "SYNERGY"

viz_results <- switch(viz_method,
  "GENIE3"  = genie3_results,
  "WGCNA"   = wgcna_results,
  "SYNERGY" = synergy_results,
  synergy_results  # default to SYNERGY
)

viz_weight_col <- switch(viz_method,
  "GENIE3"  = "weight",
  "WGCNA"   = "TOM",
  "SYNERGY" = "score_synergy",
  "score_synergy"
)

viz_directed <- switch(viz_method,
  "GENIE3"  = TRUE,
  "WGCNA"   = FALSE,
  "SYNERGY" = TRUE,
  TRUE
)

viz_edge_color <- switch(viz_method,
  "GENIE3"  = "#2ca02c",   # green
  "WGCNA"   = "#1f77b4",   # blue
  "SYNERGY" = "#d62728",   # red
  "#1f77b4"
)

visualize_network_per_cluster(
  network_results     = viz_results,
  cluster_assignments = heatmap_results,
  output_dir          = file.path(dir_08, diff_tag, "VISUALIZATION"),
  method_name         = viz_method,
  weight_col          = viz_weight_col,
  directed            = viz_directed,
  edge_color          = viz_edge_color
)

message("\n✓ SECTION 23 COMPLETE: Network visualization saved")


# =============================================================================
# SECTION 24 — CLUSTER PROFILE REPORTS
# =============================================================================
# Per-cluster profiles: heatmaps, expression statistics, functional annotation.
#
# For each cluster from SEC 20:
#   • Expression heatmap (pseudobulk × genes, with row clustering)
#   • Expression statistics (mean, SD, range)
#   • Gene count and composition
#
# Links to GO enrichment results from SEC 21 for functional context.

exprMatr_pseudobulk <- load_pseudobulk_matrix(
  file.path(dir_objects, "pseudobulk_replicas"),
  normalize = TRUE
)

generate_cluster_profile_report(
  cluster_assignments = heatmap_results,
  pseudobulk_matrix   = exprMatr_pseudobulk,
  output_dir          = file.path(dir_06, diff_tag, "CLUSTER_PROFILES"),
  method_name         = "WGCNA"  # clustering method used in SEC 20
)

message("\n✓ SECTION 24 COMPLETE: Cluster profiles saved")


