# =============================================================================
# scRNA-seq Analysis Pipeline
# =============================================================================
# Requires: load_libraries.R  |  ScRNA_Analysis_Functions.R  |  custom_seurat.R
# Organism-agnostic: adjust mt_pattern, cp_pattern, orgdb, keytype below.
# =============================================================================

source("metodologia/ScRNASeq-Docker/load_libraries.R")
source("metodologia/ScRNASeq-Docker/custom_seurat.R")
source("metodologia/ScRNASeq-Docker/ScRNA_Analysis_Functions.R")

options(Seurat.allow.s4 = FALSE)
setwd("/home/mvergara/projects2/eleo/ScRNA/")


# =============================================================================
# CONFIGURATION — edit this block for your experiment
# =============================================================================

samples <- list(
  list(file = "cellbender/Sample_0N_cellbender_filtered.h5",      label = "0N",      condition = "0N"),
  list(file = "cellbender/Sample_05N_R1_cellbender_filtered.h5",  label = "0.5N_R1", condition = "0.5N"),
  list(file = "cellbender/Sample_05N_2_cellbender_filtered.h5",   label = "0.5N_R2", condition = "0.5N"),
  list(file = "cellbender/Sample_5N_R1_cellbender_filtered.h5",   label = "5N_R1",   condition = "5N"),
  list(file = "cellbender/Sample_5N_2_cellbender_filtered.h5",    label = "5N_R2",   condition = "5N")
)

colors <- c("0N" = "#66c2a5", "0.5N_R1" = "#fc8d62", "0.5N_R2" = "#fc8d62",
            "5N_R1" = "#8da0cb", "5N_R2" = "#8da0cb")

output_dir       <- "metodologia"
resolutions_test <- c(0.15, 0.30, 0.50, 0.8, 1.0)
marcadores       <- read.table("recursos/biblio_marks.txt", header = TRUE, sep = "\t", quote = "")


# ── Plot-saving helpers ────────────────────────────────────────────────────────
# save_pdf(plot, "name.pdf")             — UMAP / FeaturePlot  (10 × 8)
# save_vln(plot, "name.pdf")             — VlnPlot single gene  (14 × 6)
# save_vln(plot, "name.pdf", n = k)      — VlnPlot k genes      (14 × 6k)
# save_qc(plot_list, "name.pdf")         — stacked QC grid

save_pdf <- function(plot, file, w = 10, h = 8)
  ggsave(file.path(output_dir, file), plot, width = w, height = h,
         dpi = 300, limitsize = FALSE)

save_vln <- function(plot, file, n = 1)
  save_pdf(plot, file, w = 14, h = 6 * n)

save_qc <- function(plot_list, file)
  ggsave(file.path(output_dir, file), wrap_plots(plot_list, ncol = 1),
         width = 14, height = 6 * length(plot_list), dpi = 300, bg = "white")


# =============================================================================
# QC — LOAD AND ANNOTATE
# =============================================================================
# mt_pattern / cp_pattern:  Arabidopsis = "^ATMG" / "^ATCG"
#                           Human = "^MT-" / NULL   |   Mouse = "^mt-" / NULL

seurat_list_raw <- lapply(samples, load_sample, mt_pattern = "^ATMG", cp_pattern = "^ATCG")
names(seurat_list_raw) <- sapply(samples, `[[`, "label")

plots_pre <- imap(seurat_list_raw, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))
save_qc(plots_pre, "qc_prefilter.pdf")


# =============================================================================
# QC — FILTER AND DOUBLETFINDER
# =============================================================================
# Adjust min_features and max_mt based on the pre-filter plots above.

seurat_list <- lapply(seurat_list_raw, filter_sample, min_features = 200, max_mt = 5)
names(seurat_list) <- sapply(samples, `[[`, "label")

plots_post <- imap(seurat_list, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))
save_qc(plots_post, "qc_postfilter.pdf")


# =============================================================================
# MERGE AND PREPROCESSING
# =============================================================================

pbmc_harmony <- reduce(seurat_list, merge) %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

pbmc_harmony$orig.ident_uni <- pbmc_harmony$condition

table(pbmc_harmony$condition)
table(pbmc_harmony$orig.ident)

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors), "umap_preharmony.pdf")  # 10×8


# =============================================================================
# HARMONY BATCH CORRECTION
# =============================================================================

dims_use <- 1:30
k_param  <- 30

pbmc_harmony <- pbmc_harmony %>% RunHarmony("orig.ident", plot_convergence = FALSE)


# =============================================================================
# ELBOW PLOT (K-MEANS WSS)
# =============================================================================

k_range  <- 2:40
pca_data <- Embeddings(pbmc_harmony, "pca")[, dims_use]
wss      <- sapply(k_range, function(k) kmeans(pca_data, centers = k, nstart = 10)$tot.withinss)

elbow_plot <- ggplot(data.frame(k = k_range, wss = wss), aes(k, wss)) +
  geom_line() + geom_point() +
  labs(x = "Number of clusters (k)", y = "Within-cluster sum of squares") +
  theme_minimal()

save_pdf(elbow_plot, "elbow_plot.pdf", w = 8, h = 6)


# =============================================================================
# CLUSTREE — RESOLUTION SWEEP
# =============================================================================

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use, k.param = k_param, verbose = FALSE)

for (res in resolutions_test)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(clustree(clu, prefix = "RNA_snn_res."), "clustree.pdf", w = 14, h = 14)


# =============================================================================
# FINAL CLUSTERING
# =============================================================================

cluster_resolution <- 0.3

pbmc_harmony <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use, k.param = k_param, verbose = FALSE) %>%
  FindClusters(resolution = cluster_resolution, algorithm = 4, verbose = FALSE)

table(Idents(pbmc_harmony))

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors), "umap_postharmony.pdf")


# =============================================================================
# CELL COUNT PER CLUSTER
# =============================================================================

Idents(pbmc_harmony) <- "orig.ident"
save_pdf(plot_integrated_clusters(pbmc_harmony), "conteocelulas.pdf", w = 14, h = 10)

colors_clusters <- sample(colors(distinct = TRUE), length(unique(pbmc_harmony$seurat_clusters)))
Idents(pbmc_harmony) <- "seurat_clusters"

save_pdf(DimPlot(pbmc_harmony, group.by = "seurat_clusters", cols = colors_clusters), "umap_seuratclusters.pdf")
save_pdf(plot_integrated_clusters(pbmc_harmony), "conteocelulas_seurat.pdf", w = 14, h = 10)


# =============================================================================
# PSEUDOBULK (GLOBAL)
# =============================================================================

pseudobulk <- generate_pseudobulk(pbmc_harmony, group_by = "orig.ident")

save_pdf(plot_replicate_correlation(pseudobulk$by_sample), "pseudobulk_correlation.pdf", w = 8, h = 8)


# =============================================================================
# MARKER GENES AND ANNOTATION
# =============================================================================

markers      <- find_markers(pbmc_harmony,
                             output_file = file.path(output_dir, "FindAllMarkers.tsv"))

pbmc_harmony <- annotate_by_markers(pbmc_harmony, markers,
                                    reference_file = file.path(output_dir, "biblio_marks.txt"))

esp          <- readRDS(file.path(output_dir, "GSE273033_seuratObj_for_publication.rds"))
pbmc_harmony <- annotate_by_reference(pbmc_harmony, reference_obj = esp, reference_col = "annotation")


# =============================================================================
# CLUSTREE — ANNOTATED
# =============================================================================

Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = dims_use, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = dims_use, k.param = k_param, verbose = FALSE)

for (res in resolutions_test)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(
  clustree(clu, prefix = "RNA_snn_res.", node_label = "celltype_reference", node_label_aggr = "Mode"),
  "clustree_annotated.pdf", w = 14, h = 14
)


# =============================================================================
# GENE EXPRESSION VISUALISATION
# =============================================================================
# Join layers before subsetting (required in Seurat 5 after merge).

pbmc_harmony     <- JoinLayers(pbmc_harmony)
Idents(pbmc_harmony) <- "celltype_reference"

gene             <- "AT5G26000"
genes_of_interest <- c("AT5G26000", "AT5G54250")
celltype         <- "Guard Cell"   # exact label from annotation table
sub_obj          <- subset(pbmc_harmony, idents = celltype)

save_vln(VlnPlot(pbmc_harmony, features = gene),                   "vln_gene_all.pdf")
save_pdf(FeaturePlot(pbmc_harmony, features = gene),               "feature_gene_all.pdf")
save_vln(VlnPlot(pbmc_harmony, features = genes_of_interest),      "vln_geneset_all.pdf",      n = length(genes_of_interest))
save_pdf(FeaturePlot(pbmc_harmony, features = genes_of_interest),  "feature_geneset_all.pdf",  h = 8 * length(genes_of_interest))

save_vln(VlnPlot(sub_obj, features = gene),                        "vln_gene_celltype.pdf")
save_pdf(FeaturePlot(sub_obj, features = gene),                    "feature_gene_celltype.pdf")
save_vln(VlnPlot(sub_obj, features = genes_of_interest),           "vln_geneset_celltype.pdf", n = length(genes_of_interest))
save_pdf(FeaturePlot(sub_obj, features = genes_of_interest),       "feature_geneset_celltype.pdf", h = 8 * length(genes_of_interest))


# =============================================================================
# CELL TYPE GROUPING (OPTIONAL)
# =============================================================================
# Map fine annotations to broader groups; unlisted types keep their name.

grouping <- c(
  "Companion Cell" = "Vascular Cell", "Cambium"  = "Vascular Cell",
  "Phloem Parenchyma" = "Vascular Cell", "Xylem" = "Vascular Cell",
  "Sieve Element"  = "Vascular Cell", "Meristemoid" = "Stomatal Line"
)

pbmc_harmony$annotation_agrupada <- recode(pbmc_harmony$celltype_reference, !!!grouping)

save_pdf(DimPlot(pbmc_harmony, group.by = "annotation_agrupada", label = TRUE, repel = TRUE, raster = FALSE),
         "umap_annotated.pdf")


# =============================================================================
# CELL TYPE CURATION — INTERACTIVE (run step by step, do not source all at once)
# =============================================================================
# This section is intentionally manual: cell types differ across experiments.
# Workflow:
#   1. Subcluster types that look heterogeneous in the UMAP
#   2. Inspect the subclusters with DimPlot + FeaturePlot
#   3. Fill in the reassignment table below
#   4. Apply corrections to the global object

Idents(pbmc_harmony) <- "annotation_agrupada"

# ── 1. Subcluster ─────────────────────────────────────────────────────────────
# Add or remove cell types as needed for your dataset.

meristemoid_umap   <- subclustar_tipo(pbmc_harmony, "Stomatal Line")
pavement_cell_umap <- subclustar_tipo(pbmc_harmony, "Pavement Cell")

# ── 2. Inspect ────────────────────────────────────────────────────────────────

DimPlot(meristemoid_umap,   group.by = "cluster_subtipo", label = TRUE, raster = FALSE)
DimPlot(pavement_cell_umap, group.by = "cluster_subtipo", label = TRUE, raster = FALSE)

for (i in seq_len(nrow(marcadores)))
  print(FeaturePlot(pavement_cell_umap, features = marcadores$gene[i]) +
        ggtitle(paste(marcadores$cell.type[i], "-", marcadores$gene[i])))

FindAllMarkers(meristemoid_umap,   only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)

FindAllMarkers(pavement_cell_umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)

# ── 3. Reassignment table ─────────────────────────────────────────────────────
# For each subclustered object, map cluster IDs to cell type names.
# Cluster IDs come from $cluster_subtipo (values like "0", "1", "2", ...).
# Add one entry per cell type you subclustered.

reassign <- list(
  meristemoid_umap = c(
    "0" = "Stomatal Line",
    "1" = "Stomatal Line",
    "2" = "Pavement Cell",
    "3" = "Stomatal Line",
    "4" = "Stomatal Line"
  ),
  pavement_cell_umap = c(
    "0" = "Pavement Cell",
    "1" = "Pavement Cell",
    "2" = "Pavement Cell",
    "3" = "Mesophyll",
    "4" = "Pavement Cell"
  )
)

# ── 4. Apply ──────────────────────────────────────────────────────────────────

pbmc_harmony$celltype_reference_curated <- pbmc_harmony$annotation_agrupada

for (obj_name in names(reassign)) {
  obj   <- get(obj_name)
  pbmc_harmony$celltype_reference_curated[colnames(obj)] <- reassign[[obj_name]][obj$cluster_subtipo]
}

save_pdf(DimPlot(pbmc_harmony, group.by = "celltype_reference_curated", label = TRUE, repel = TRUE, raster = FALSE),
         "umap_curada.pdf")


# =============================================================================
# EXPORT TO H5AD (SCANPY)
# =============================================================================

exportar_para_scanpy(pbmc_harmony, "results/objs/pbmc_harmony_curated.h5ad")

# To export a single cell type:
# exportar_para_scanpy(subset(pbmc_harmony, subset = celltype_reference_curated == "Guard Cell"),
#                      "results/objs/GuardCell.h5ad")


# =============================================================================
# DOTPLOT — MARKER GENES BY CELL TYPE
# =============================================================================

cell_order_dotplot <- c(
  "Pavement Cell", "Stomatal lineage", "Guard Cell", "Mesophyll", "Bundle Sheath",
  "Phloem Parenchyma", "Cambium", "Xylem", "Companion Cell",
  "Hydathode", "Cell Cycle: G1-S", "Cell Cycle: G2-M"
)

hacer_dotplot_marcadores(
  pbmc_harmony, marcadores,
  annot_col       = "celltype_reference_curated",
  cell_order      = cell_order_dotplot,
  clusters_remove = c("Sieve Element", "Myrosin Idioblast"),
  outfile         = file.path(output_dir, "dotplot_marcadores.pdf"),
  width = 20, height = 10
)


# =============================================================================
# CELL TYPE SUBSETS
# =============================================================================

celular_subsets <- setNames(
  lapply(unique(pbmc_harmony$celltype_reference_curated),
         function(t) subset(pbmc_harmony, subset = celltype_reference_curated == t)),
  gsub("[^[:alnum:]_]", "_", unique(pbmc_harmony$celltype_reference_curated))
)


# =============================================================================
# PSEUDOREPLICATES AND PSEUDOBULK PER CELL TYPE
# =============================================================================

celular_subsets_replicados <- Filter(Negate(is.null), lapply(celular_subsets, asignar_pseudoreplicados))
pseudobulk_list            <- lapply(celular_subsets_replicados, hacer_pseudobulk)

dir.create("results/pseudobulk_replicas", recursive = TRUE, showWarnings = FALSE)
for (tipo in names(pseudobulk_list))
  write.csv(pseudobulk_list[[tipo]],
            file.path("results/pseudobulk_replicas", paste0("Pseudobulk_", tipo, ".csv")),
            row.names = TRUE)


# =============================================================================
# DESEQ2 — ALL COMPARISONS
# =============================================================================

comparaciones <- list(
  list(conds = c("0.5N", "5N"), tag = "05_5"),
  list(conds = c("0N",   "5N"), tag = "0_5"),
  list(conds = c("0N", "0.5N"), tag = "0_05")
)

for (tag in sapply(comparaciones, `[[`, "tag"))
  dir.create(file.path("results/deseq2", tag), recursive = TRUE, showWarnings = FALSE)

for (tipo in names(pseudobulk_list))
  correr_deseq2(as.matrix(pseudobulk_list[[tipo]]), comparaciones,
                output_dir = "results/deseq2", tipo = tipo)


# =============================================================================
# VOLCANO PLOTS — ALL COMPARISONS
# =============================================================================

dir.create("results/volcano_plots", showWarnings = FALSE)

for (comp in comparaciones) {
  tag       <- comp$tag
  csv_files <- list.files(file.path("results/deseq2", tag), pattern = "\\.csv$", full.names = TRUE)
  if (!length(csv_files)) { message("No CSVs: ", tag); next }

  pdf(file.path("results/volcano_plots", paste0("VolcanoPlots_", tag, ".pdf")), width = 12, height = 6)
  plots <- list()
  for (f in csv_files) {
    plots <- c(plots, list(hacer_volcano(f)))
    if (length(plots) == 2) { grid.arrange(grobs = plots, ncol = 2); plots <- list() }
  }
  if (length(plots)) grid.arrange(grobs = plots, ncol = 1)
  dev.off()
}


# =============================================================================
# DIFFERENTIAL TABLE AND HEATMAP — ALL COMPARISONS
# =============================================================================

for (comp in comparaciones) {
  tag       <- comp$tag
  csv_dir   <- file.path("results/deseq2", tag)
  diff_dir  <- file.path("results/diff",   tag)
  csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

  dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(csv_files)) { message("No CSVs: ", tag); next }

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
# GO ENRICHMENT — ALL COMPARISONS
# =============================================================================
# Organism:  Arabidopsis = org.At.tair.db / "TAIR"
#            Human       = org.Hs.eg.db   / "ENSEMBL"
#            Mouse       = org.Mm.eg.db   / "ENSEMBL"

orgdb      <- org.At.tair.db
keytype    <- "TAIR"
universo   <- keys(orgdb, keytype = keytype)
espacio    <- "BP"   # "BP" | "MF" | "CC"
qval       <- 0.05
nivel_poda <- 6

dir.create("results/Enrichment", showWarnings = FALSE)

for (comp in comparaciones) {
  tag        <- comp$tag
  diff_dir   <- file.path("results/diff",       tag)
  enr_dir    <- file.path("results/Enrichment", tag)
  tabla_path <- file.path(diff_dir, "tabla_diferenciales.tsv")

  dir.create(enr_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(tabla_path)) { message("No table: ", tag); next }

  tabla <- read.table(tabla_path, header = TRUE, row.names = 1, sep = "\t")
  tabla <- tabla[, colSums(tabla != 0) > 0, drop = FALSE]

  go_total  <- correr_enriquecimiento_go(tabla, universo, espacio, orgdb = orgdb, keytype = keytype,
                                         simplificar = FALSE, output_dir = enr_dir)
  go_simple <- correr_enriquecimiento_go(tabla, universo, espacio, orgdb = orgdb, keytype = keytype,
                                         simplificar = TRUE,  output_dir = enr_dir)

  go_total_podado  <- podar_go(go_total,  nivel_poda, espacio, qval, simplificar = FALSE)
  go_simple_podado <- podar_go(go_simple, nivel_poda, espacio, qval, simplificar = TRUE)

  pdf(file.path(enr_dir, paste0("GO_enrichment_", tag, ".pdf")), width = 18, height = 18)
  tryCatch({ print(graficar_go_balones(go_total));        print(graficar_go_balones(go_simple))
             print(graficar_go_balones(go_total_podado)); print(graficar_go_balones(go_simple_podado)) },
           error = function(e) message("GO plot error: ", e$message))
  dev.off()
}
