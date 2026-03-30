# =============================================================================
# scRNA-seq Analysis — Arabidopsis thaliana Nitrogen Conditions
# =============================================================================
# Description : Full pipeline from CellBender HDF5 to annotated Seurat object.
#               Conditions: 0N, 0.5N (R1/R2), 5N (R1/R2).
# Requires    : load_libraries.R, ScRNA_Analysis_Functions.R, custom_seurat.R
#
# Organism-agnostic: Change mt_pattern / cp_pattern and orgdb / keytype below
#   to adapt this pipeline to any species.
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

# ── Samples ───────────────────────────────────────────────────────────────────
# Adjust file paths, labels, and condition tags to match your experiment.

samples <- list(
  list(file = "cellbender/Sample_0N_cellbender_filtered.h5",      label = "0N",      condition = "0N"),
  list(file = "cellbender/Sample_05N_R1_cellbender_filtered.h5",  label = "0.5N_R1", condition = "0.5N"),
  list(file = "cellbender/Sample_05N_2_cellbender_filtered.h5",   label = "0.5N_R2", condition = "0.5N"),
  list(file = "cellbender/Sample_5N_R1_cellbender_filtered.h5",   label = "5N_R1",   condition = "5N"),
  list(file = "cellbender/Sample_5N_2_cellbender_filtered.h5",    label = "5N_R2",   condition = "5N")
)

# ── Colours (one per label) ───────────────────────────────────────────────────

colors <- c(
  "0N"      = "#66c2a5",
  "0.5N_R1" = "#fc8d62",
  "0.5N_R2" = "#fc8d62",
  "5N_R1"   = "#8da0cb",
  "5N_R2"   = "#8da0cb"
)

# ── Output directory ──────────────────────────────────────────────────────────

output_dir <- "metodologia"


# =============================================================================
# QC — PRE-FILTER
# =============================================================================
# Organism patterns:
#   Arabidopsis : mt_pattern = "^ATMG",  cp_pattern = "^ATCG"
#   Human       : mt_pattern = "^MT-",   cp_pattern = NULL
#   Mouse       : mt_pattern = "^mt-",   cp_pattern = NULL

seurat_list_raw <- lapply(samples, process_sample,
                          mt_pattern = "^ATMG",
                          cp_pattern = "^ATCG")
names(seurat_list_raw) <- sapply(samples, `[[`, "label")

plots_pre <- imap(seurat_list_raw, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))

ggsave(
  filename = file.path(output_dir, "qc_violin_all_samples_prefilter.pdf"),
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
  filename = file.path(output_dir, "qc_violin_all_samples_postfilter.pdf"),
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
  filename = file.path(output_dir, "Umap_preharmony.pdf"),
  plot     = DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
  width    = 18, height = 18, dpi = 500, limitsize = FALSE
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
  filename = file.path(output_dir, "Umap_postharmony.pdf"),
  plot     = DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
  width    = 18, height = 18, dpi = 500, limitsize = FALSE
)


# =============================================================================
# CELL COUNT PER CLUSTER
# =============================================================================

Idents(pbmc_harmony) <- "orig.ident"

ggsave(
  filename = file.path(output_dir, "conteocelulas.pdf"),
  plot     = plot_integrated_clusters(pbmc_harmony),
  width    = 20, height = 20, dpi = 500, limitsize = FALSE
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
# PSEUDOBULK (GLOBAL)
# =============================================================================

pseudobulk <- generate_pseudobulk(pbmc_harmony, group_by = "orig.ident")
pseudobulk$by_sample     # matrix: one column per sample
pseudobulk$by_condition  # matrix: one column per condition

plot_replicate_correlation(pseudobulk$by_sample)
ggsave(file.path(output_dir, "pseudobulk_correlation.pdf"), width = 10, height = 10)


# =============================================================================
# MARKER GENES AND ANNOTATION
# =============================================================================

markers <- find_markers(pbmc_harmony,
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
  clustree(clu, prefix       = "RNA_snn_res.",
               node_label      = "celltype_reference",
               node_label_aggr = "Mode"),
  width = 18, height = 18
)


# =============================================================================
# GENE EXPRESSION VISUALISATION
# =============================================================================
# Join layers before subsetting (required in Seurat 5 after merge).

pbmc_harmony <- JoinLayers(pbmc_harmony)

Idents(pbmc_harmony) <- "celltype_reference"

genes_of_interest    <- c("AT5G26000", "AT5G54250")
celltype_of_interest <- "Guard Cell"   # copy exact label from annotation table
gene                 <- "AT5G26000"

# ── Single gene, all clusters ─────────────────────────────────────────────────

ggsave(
  file.path(output_dir, "VlnPlot_single_gene_all.pdf"),
  VlnPlot(pbmc_harmony, features = gene),
  width = 14, height = 6
)

ggsave(
  file.path(output_dir, "FeaturePlot_single_gene_all.pdf"),
  FeaturePlot(pbmc_harmony, features = gene),
  width = 10, height = 8
)

# ── Gene set, all clusters ────────────────────────────────────────────────────

ggsave(
  file.path(output_dir, "VlnPlot_geneset_all.pdf"),
  VlnPlot(pbmc_harmony, features = genes_of_interest),
  width = 14, height = 6 * length(genes_of_interest)
)

ggsave(
  file.path(output_dir, "FeaturePlot_geneset_all.pdf"),
  FeaturePlot(pbmc_harmony, features = genes_of_interest),
  width = 10, height = 8 * length(genes_of_interest)
)

# ── Single gene, one cell type ────────────────────────────────────────────────

sub_obj <- subset(pbmc_harmony, idents = celltype_of_interest)

ggsave(
  file.path(output_dir, paste0("VlnPlot_", gsub(" ", "_", celltype_of_interest), "_single.pdf")),
  VlnPlot(sub_obj, features = gene),
  width = 10, height = 6
)

ggsave(
  file.path(output_dir, paste0("FeaturePlot_", gsub(" ", "_", celltype_of_interest), "_single.pdf")),
  FeaturePlot(sub_obj, features = gene),
  width = 10, height = 8
)

# ── Gene set, one cell type ───────────────────────────────────────────────────

ggsave(
  file.path(output_dir, paste0("VlnPlot_", gsub(" ", "_", celltype_of_interest), "_geneset.pdf")),
  VlnPlot(sub_obj, features = genes_of_interest),
  width = 10, height = 6 * length(genes_of_interest)
)

ggsave(
  file.path(output_dir, paste0("FeaturePlot_", gsub(" ", "_", celltype_of_interest), "_geneset.pdf")),
  FeaturePlot(sub_obj, features = genes_of_interest),
  width = 10, height = 8 * length(genes_of_interest)
)


# =============================================================================
# CELL TYPE GROUPING (OPTIONAL)
# =============================================================================
# Map fine annotations to broader groups. Cell types not listed here keep
# their original name.

grouping <- c(
  "Companion Cell"    = "Vascular Cell",
  "Cambium"           = "Vascular Cell",
  "Phloem Parenchyma" = "Vascular Cell",
  "Xylem"             = "Vascular Cell",
  "Sieve Element"     = "Vascular Cell",
  "Meristemoid"       = "Stomatal Line"
)

pbmc_harmony$annotation_agrupada <- recode(pbmc_harmony$celltype_reference, !!!grouping)

# ── UMAP with grouped annotations ────────────────────────────────────────────

figure_grouped <- DimPlot(pbmc_harmony,
                          group.by = "annotation_agrupada",
                          label = TRUE, repel = TRUE, raster = FALSE)

ggsave(
  file.path(output_dir, "umap_annotated.pdf"),
  figure_grouped,
  width = 18, height = 18, dpi = 500, limitsize = FALSE
)


# =============================================================================
# CELL TYPE CURATION — INTERACTIVE (DO NOT AUTOMATE)
# =============================================================================
# Run interactively: inspect subclusters, then fill in the corrections block.

Idents(pbmc_harmony) <- "annotation_agrupada"

marcadores <- read.table("recursos/biblio_marks.txt", header = TRUE, sep = "\t", quote = "")

# ── Step 1: Subcluster cell types of interest ─────────────────────────────────

meristemoid_umap   <- subclustar_tipo(pbmc_harmony, "Stomatal Line")
pavement_cell_umap <- subclustar_tipo(pbmc_harmony, "Pavement Cell")

# ── Step 2: Inspect subclusters ───────────────────────────────────────────────

DimPlot(meristemoid_umap,   group.by = "cluster_subtipo", label = TRUE, raster = FALSE)
DimPlot(pavement_cell_umap, group.by = "cluster_subtipo", label = TRUE, raster = FALSE)

for (i in 1:nrow(marcadores)) {
  print(FeaturePlot(pavement_cell_umap, features = marcadores$gene[i]) +
          ggtitle(paste(marcadores$cell.type[i], "-", marcadores$gene[i])))
}

# ── Step 3: View cluster-level markers ───────────────────────────────────────

markers_meristemoid   <- FindAllMarkers(meristemoid_umap,   only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_pavement_cell <- FindAllMarkers(pavement_cell_umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers_meristemoid   %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)
markers_pavement_cell %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)

# ── Step 4: Define corrections (fill after inspecting) ───────────────────────

correcciones <- list(

  # Meristemoid = list(
  #   obj = meristemoid_umap,
  #   mapa = c(
  #     "0" = "Stomatal Line",
  #     "1" = "Stomatal Line",
  #     "2" = "Pavement Cell",
  #     "3" = "Stomatal Line",
  #     "4" = "Stomatal Line"
  #   )
  # ),

  PavementCell = list(
    obj = pavement_cell_umap,
    mapa = c(
      "0" = "Pavement Cell",
      "1" = "Pavement Cell",
      "2" = "Pavement Cell",
      "3" = "Mesophyll",
      "4" = "Pavement Cell"
    )
  )
)

# ── Step 5: Apply corrections to the global object ───────────────────────────

pbmc_harmony$celltype_reference_curated <- pbmc_harmony$annotation_agrupada

for (tipo in names(correcciones)) {
  obj_local  <- correcciones[[tipo]]$obj
  mapa_local <- correcciones[[tipo]]$mapa
  celdas     <- colnames(obj_local)
  nuevos     <- mapa_local[obj_local$cluster_subtipo]
  pbmc_harmony$celltype_reference_curated[celdas] <- nuevos
}

# ── UMAP with curated annotations ────────────────────────────────────────────

figure_curated <- DimPlot(pbmc_harmony, group.by = "celltype_reference_curated",
                          label = TRUE, repel = TRUE, raster = FALSE)

ggsave(
  file.path(output_dir, "umap_curada.pdf"),
  figure_curated,
  width = 18, height = 18, dpi = 500, limitsize = FALSE
)


# =============================================================================
# CELL TYPE SUBSETS
# =============================================================================

celular_subsets <- setNames(
  lapply(unique(pbmc_harmony$celltype_reference_curated), function(tipo) {
    subset(pbmc_harmony, subset = celltype_reference_curated == tipo)
  }),
  gsub("[^[:alnum:]_]", "_", unique(pbmc_harmony$celltype_reference_curated))
)


# =============================================================================
# PSEUDOREPLICATES AND PSEUDOBULK PER CELL TYPE
# =============================================================================

celular_subsets_replicados <- Filter(Negate(is.null),
                                     lapply(celular_subsets, asignar_pseudoreplicados))

pseudobulk_list <- lapply(celular_subsets_replicados, hacer_pseudobulk)

# ── Save pseudobulk matrices ──────────────────────────────────────────────────

dir.create("results/pseudobulk_replicas", recursive = TRUE, showWarnings = FALSE)

for (tipo in names(pseudobulk_list)) {
  write.csv(pseudobulk_list[[tipo]],
            file.path("results/pseudobulk_replicas", paste0("Pseudobulk_", tipo, ".csv")),
            row.names = TRUE)
}


# =============================================================================
# DESEQ2 — ALL COMPARISONS
# =============================================================================
# Add or remove comparisons as needed.

comparaciones <- list(
  list(conds = c("0.5N", "5N"), tag = "05_5"),
  list(conds = c("0N",   "5N"), tag = "0_5"),
  list(conds = c("0N", "0.5N"), tag = "0_05")
)

for (tag in sapply(comparaciones, `[[`, "tag")) {
  dir.create(file.path("results/deseq2", tag), recursive = TRUE, showWarnings = FALSE)
}

for (tipo in names(pseudobulk_list)) {
  cat("DESeq2:", tipo, "\n")
  correr_deseq2(as.matrix(pseudobulk_list[[tipo]]),
                comparaciones = comparaciones,
                output_dir    = "results/deseq2")
}


# =============================================================================
# VOLCANO PLOTS — ALL COMPARISONS
# =============================================================================

dir.create("results/volcano_plots", showWarnings = FALSE)

for (comp in comparaciones) {

  tag      <- comp$tag
  csv_dir  <- file.path("results/deseq2", tag)
  csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(csv_files) == 0) {
    message("No CSV files found for comparison: ", tag)
    next
  }

  pdf_path <- file.path("results/volcano_plots", paste0("VolcanoPlots_", tag, ".pdf"))
  pdf(pdf_path, width = 12, height = 6)

  plots <- list()
  for (file in csv_files) {
    plots <- append(plots, list(hacer_volcano(file, file.path("results/volcano_plots", tag))))
    if (length(plots) == 2) {
      grid.arrange(grobs = plots, ncol = 2)
      plots <- list()
    }
  }
  if (length(plots) == 1) grid.arrange(grobs = plots, ncol = 1)

  dev.off()
  message("Volcano plots saved: ", pdf_path)
}


# =============================================================================
# DIFFERENTIAL EXPRESSION TABLE AND HEATMAP — ALL COMPARISONS
# =============================================================================

for (comp in comparaciones) {

  tag       <- comp$tag
  csv_dir   <- file.path("results/deseq2", tag)
  diff_dir  <- file.path("results/diff", tag)
  csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

  dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(csv_files) == 0) {
    message("No CSV files for comparison: ", tag)
    next
  }

  listas     <- lapply(csv_files, procesar_deseq2_resultado, output_dir = diff_dir)
  tabla_class <- Reduce(function(x, y) full_join(x, y, by = "AGI"), lapply(listas, `[[`, "class"))
  tabla_logfc <- Reduce(function(x, y) full_join(x, y, by = "AGI"), lapply(listas, `[[`, "logfc"))

  tabla_class <- tabla_class %>% filter(apply(select(., -AGI) != 0, 1, any))
  tabla_logfc <- tabla_logfc %>% filter(AGI %in% tabla_class$AGI)

  write_tsv(tabla_class, file.path(diff_dir, "tabla_diferenciales.tsv"))
  write_tsv(tabla_logfc, file.path(diff_dir, "tabla_log2FC.tsv"))

  # ── Heatmap ──────────────────────────────────────────────────────────────────

  matriz <- as.matrix(column_to_rownames(tabla_logfc, "AGI"))
  matriz[is.na(matriz)] <- 0

  if (nrow(matriz) > 1) {
    pdf(file.path(diff_dir, paste0("heatmap_", tag, ".pdf")), width = 14, height = 18)
    tryCatch(hacer_heatmap(matriz), error = function(e) message("Heatmap error: ", e$message))
    dev.off()
    message("Heatmap saved for: ", tag)
  } else {
    message("Not enough genes for heatmap in comparison: ", tag)
  }
}


# =============================================================================
# GO ENRICHMENT — ALL COMPARISONS
# =============================================================================
# Organism-specific parameters:
#   Arabidopsis : orgdb = org.At.tair.db, keytype = "TAIR"
#   Human       : orgdb = org.Hs.eg.db,   keytype = "ENSEMBL" or "SYMBOL"
#   Mouse       : orgdb = org.Mm.eg.db,   keytype = "ENSEMBL" or "SYMBOL"

orgdb   <- org.At.tair.db
keytype <- "TAIR"

# Background gene universe (TAIR format for Arabidopsis)
universo <- read.table("recursos/gene_association.tair.tsv",
                       sep = "\t", quote = "", fill = TRUE, comment.char = "!")[, -1]
universo <- as.matrix(universo)[, 1]

# GO parameters
espacio    <- "BP"
qval       <- 0.05
nivel_poda <- 6

dir.create("results/Enrichment", showWarnings = FALSE)

for (comp in comparaciones) {

  tag      <- comp$tag
  diff_dir <- file.path("results/diff", tag)
  enr_dir  <- file.path("results/Enrichment", tag)
  dir.create(enr_dir, showWarnings = FALSE, recursive = TRUE)

  tabla_path <- file.path(diff_dir, "tabla_diferenciales.tsv")
  if (!file.exists(tabla_path)) {
    message("No differential table for comparison: ", tag)
    next
  }

  tabla <- read.table(tabla_path, header = TRUE, row.names = 1, sep = "\t")
  tabla <- tabla[, colSums(tabla != 0) > 0, drop = FALSE]

  # ── Run enrichment ───────────────────────────────────────────────────────────

  go_total  <- correr_enriquecimiento_go(tabla, universo, espacio,
                                         orgdb = orgdb, keytype = keytype,
                                         simplificar = FALSE,
                                         output_dir  = enr_dir)

  go_simple <- correr_enriquecimiento_go(tabla, universo, espacio,
                                         orgdb = orgdb, keytype = keytype,
                                         simplificar = TRUE,
                                         output_dir  = enr_dir)

  # ── Prune and plot ───────────────────────────────────────────────────────────

  go_total_podado  <- podar_go(go_total,  nivel_poda, espacio, qval, simplificar = FALSE)
  go_simple_podado <- podar_go(go_simple, nivel_poda, espacio, qval, simplificar = TRUE)

  p1 <- graficar_go_balones(go_total)
  p2 <- graficar_go_balones(go_simple)
  p3 <- graficar_go_balones(go_total_podado)
  p4 <- graficar_go_balones(go_simple_podado)

  pdf(file.path(enr_dir, paste0("GO_enrichment_", tag, ".pdf")), width = 18, height = 18)
  tryCatch({ print(p1); print(p2); print(p3); print(p4) }, silent = TRUE)
  dev.off()

  message("GO enrichment saved for: ", tag)
}
