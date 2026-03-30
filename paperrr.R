# =============================================================================
# scRNA-seq Analysis — Arabidopsis thaliana Nitrogen Conditions
# =============================================================================
# Description : Full pipeline from CellBender HDF5 to annotated Seurat object.
#               Conditions: 0N, 0.5N (R1/R2), 5N (R1/R2).
# Requires    : load_libraries.R, ScRNA_Analysis_Functions.R, custom_seurat.R
# =============================================================================

# ── Libraries and functions ───────────────────────────────────────────────────
source("load_libraries.R")
source("metodologia/custom_seurat.R")
source("metodologia/ScRNA_Analysis_Functions.R")

sessionInfo()


# ── Global settings ────────────────────────────────────────────────────────────
options(Seurat.allow.s4 = FALSE)
knitr::opts_chunk$set(cache = TRUE)
setwd("/home/mvergara/projects2/eleo/ScRNA/")


# =============================================================================
# CONFIGURATION
# =============================================================================

samples <- list(
  list(file = "cellbender/Sample_0N_cellbender_filtered.h5",      label = "0N",      condition = "0N"),
  list(file = "cellbender/Sample_05N_R1_cellbender_filtered.h5",  label = "0.5N_R1", condition = "0.5N"),
  list(file = "cellbender/Sample_05N_2_cellbender_filtered.h5",   label = "0.5N_R2", condition = "0.5N"),
  list(file = "cellbender/Sample_5N_R1_cellbender_filtered.h5",   label = "5N_R1",   condition = "5N"),
  list(file = "cellbender/Sample_5N_2_cellbender_filtered.h5",    label = "5N_R2",   condition = "5N")
)

colors <- c(
  "0N"      = "#66c2a5",
  "0.5N_R1" = "#fc8d62",
  "0.5N_R2" = "#fc8d62",
  "5N_R1"   = "#8da0cb",
  "5N_R2"   = "#8da0cb"
)

output_dir <- "metodologia"


# =============================================================================
# QC — PRE-FILTER
# =============================================================================

seurat_list_raw <- lapply(samples, process_sample,
                          mt_pattern = "^ATMG",   # Arabidopsis mitochondria
                          cp_pattern = "^ATCG")   # Arabidopsis chloroplast
                          # Human : mt_pattern = "^MT-",  cp_pattern = NULL
                          # Mouse : mt_pattern = "^mt-",  cp_pattern = NULL
names(seurat_list_raw) <- sapply(samples, `[[`, "label")

plots_pre <- imap(seurat_list_raw, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))

ggsave(
  filename = file.path(output_dir, "qc_violin_all_samples_prefilter.png"),
  plot     = wrap_plots(plots_pre, ncol = 1),
  width    = 14,
  height   = 6 * length(plots_pre),
  dpi      = 300,
  bg       = "white"
)


# =============================================================================
# QC — FILTERING AND POST-FILTER CHECK
# =============================================================================

seurat_list <- lapply(samples, process_sample,
                      min_features = 200,
                      max_mt       = 5,
                      mt_pattern   = "^ATMG",
                      cp_pattern   = "^ATCG")
names(seurat_list) <- sapply(samples, `[[`, "label")

plots_post <- imap(seurat_list, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))

ggsave(
  filename = file.path(output_dir, "qc_violin_all_samples_postfilter.png"),
  plot     = wrap_plots(plots_post, ncol = 1),
  width    = 14,
  height   = 6 * length(plots_post),
  dpi      = 300,
  bg       = "white"
)


# =============================================================================
# MERGE AND PREPROCESSING
# =============================================================================

pbmc_harmony <- reduce(seurat_list, merge)

pbmc_harmony <- pbmc_harmony %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

pbmc_harmony$orig.ident_uni <- pbmc_harmony$condition
pbmc_harmony.bkp            <- pbmc_harmony

cat("\n=== SUMMARY ===\n")
cat("Samples processed:", length(seurat_list), "\n")
cat("Total cells:      ", ncol(pbmc_harmony), "\n")
cat("Conditions:       ", paste(unique(pbmc_harmony$condition), collapse = ", "), "\n")
table(pbmc_harmony$condition)
table(pbmc_harmony$orig.ident)

ggsave(
  filename  = file.path(output_dir, "Umap_preharmony.pdf"),
  plot      = DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
  width     = 18, height = 18, dpi = 500, limitsize = FALSE
)


# =============================================================================
# HARMONY BATCH CORRECTION
# =============================================================================

harmony_var      <- "orig.ident"
dims_use         <- 1:30
k_param          <- 30
resolutions_test <- c(0.05, 0.15, 0.30, 0.40, 0.50, 0.60, 0.7, 0.8, 0.9, 1.0)
k_range          <- 2:40

pbmc_harmony <- pbmc_harmony %>%
  RunHarmony(harmony_var, plot_convergence = FALSE)


# =============================================================================
# ELBOW PLOT (K-MEANS WSS)
# =============================================================================

pca_data <- Embeddings(pbmc_harmony, "pca")[, dims_use]
wss      <- sapply(k_range, function(k) kmeans(pca_data, centers = k, nstart = 10)$tot.withinss)

elbow_plot <- ggplot(data.frame(k = k_range, wss = wss), aes(k, wss)) +
  geom_line() + geom_point() +
  labs(x = "Number of clusters (k)", y = "Within-cluster sum of squares") +
  theme_minimal()

ggsave(file.path(output_dir, "elbow_plot.pdf"), elbow_plot, width = 8, height = 6)


# =============================================================================
# CLUSTREE — RESOLUTION SWEEP
# =============================================================================

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use, k.param = k_param, verbose = FALSE)

for (res in resolutions_test) {
  clu        <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)
  n_clusters <- length(unique(clu[[paste0("RNA_snn_res.", res)]][, 1]))
  message("Resolution ", res, " -> ", n_clusters, " clusters")
}

ggsave(
  file.path(output_dir, "clustree.pdf"),
  clustree(clu, prefix = "RNA_snn_res."),
  width = 14, height = 14
)


# =============================================================================
# FINAL CLUSTERING
# =============================================================================

cluster_resolution <- 0.3

pbmc_harmony <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use, k.param = k_param, verbose = FALSE) %>%
  FindClusters(resolution = cluster_resolution, algorithm = 4, verbose = FALSE)

cat("\n=== CLUSTERING ===\n")
cat("Resolution:", cluster_resolution, "\n")
cat("Clusters:  ", length(unique(Idents(pbmc_harmony))), "\n")
table(Idents(pbmc_harmony))

ggsave(
  filename  = file.path(output_dir, "Umap_postharmony.pdf"),
  plot      = DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
  width     = 18, height = 18, dpi = 500, limitsize = FALSE
)


# =============================================================================
# CELL COUNT PER CLUSTER
# =============================================================================

Idents(pbmc_harmony) <- "orig.ident"
ggsave(
  filename  = file.path(output_dir, "conteocelulas.pdf"),
  plot      = plot_integrated_clusters(pbmc_harmony),
  width     = 20, height = 20, dpi = 500, limitsize = FALSE
)

n_groups        <- length(unique(pbmc_harmony$seurat_clusters))
colors_clusters <- sample(colors(distinct = TRUE), n_groups)

Idents(pbmc_harmony) <- "seurat_clusters"
ggsave(
  file.path(output_dir, "Umap_seuratclusters.pdf"),
  DimPlot(pbmc_harmony, group.by = "seurat_clusters", cols = colors_clusters),
  width = 18, height = 18, dpi = 500, limitsize = FALSE
)
ggsave(
  file.path(output_dir, "conteocelulas_seurat.pdf"),
  plot_integrated_clusters(pbmc_harmony),
  width = 18, height = 18, dpi = 500, limitsize = FALSE
)


# =============================================================================
# PSEUDOBULK
# =============================================================================

pseudobulk <- generate_pseudobulk(pbmc_harmony, group_by = "orig.ident")
pseudobulk$by_sample     # matrix: one column per sample
pseudobulk$by_condition  # matrix: one column per condition

plot_replicate_correlation(pseudobulk$by_sample)
ggsave(file.path(output_dir, "pseudobulk_correlation.pdf"), width = 10, height = 10)


# =============================================================================
# MARKER GENES AND ANNOTATION
# =============================================================================

markers      <- find_markers(pbmc_harmony,
                             output_file = file.path(output_dir, "FindAllMarkers_0.3C.tsv"))

pbmc_harmony <- annotate_by_markers(pbmc_harmony, markers,
                                    reference_file = file.path(output_dir, "biblio_marks.txt"))

esp          <- readRDS(file.path(output_dir, "GSE273033_seuratObj_for_publication.rds"))
pbmc_harmony <- annotate_by_reference(pbmc_harmony,
                                      reference_obj = esp,
                                      reference_col = "annotation")


# =============================================================================
# CLUSTREE — ANNOTATED
# =============================================================================

Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use, k.param = k_param, verbose = FALSE)

for (res in resolutions_test) {
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)
}

ggsave(
  file.path(output_dir, "clustree_annotated.pdf"),
  clustree(clu, prefix = "RNA_snn_res.",
           node_label      = "celltype_reference",
           node_label_aggr = "Mode"),
  width = 18, height = 18
)


# =============================================================================
# GENE EXPRESSION VISUALISATION
# =============================================================================

# Join layers before subsetting (required in Seurat 5 after merge)
pbmc_harmony <- JoinLayers(pbmc_harmony)

Idents(pbmc_harmony) <- "celltype_reference"

genes_of_interest    <- c("AT5G26000", "AT5G54250")
celltype_of_interest <- "Guard Cell"   # copy exact label from annotation table
gene                 <- "AT5G26000"

# ── Single gene, all clusters ─────────────────────────────────────────────────
VlnPlot(pbmc_harmony, features = gene)
FeaturePlot(pbmc_harmony, features = gene)

# ── Gene set, all clusters ────────────────────────────────────────────────────
VlnPlot(pbmc_harmony, features = genes_of_interest)
FeaturePlot(pbmc_harmony, features = genes_of_interest)

# ── Single gene, one cell type ────────────────────────────────────────────────
sub_obj <- subset(pbmc_harmony, idents = celltype_of_interest)
VlnPlot(sub_obj, features = gene)
FeaturePlot(sub_obj, features = gene)

# ── Gene set, one cell type ───────────────────────────────────────────────────
VlnPlot(sub_obj, features = genes_of_interest)
FeaturePlot(sub_obj, features = genes_of_interest)
