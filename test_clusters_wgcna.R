# =============================================================================
# TEST: Hierarchical clusters + WGCNA modules + double-annotation heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

# ── INPUT ─────────────────────────────────────────────────────────────────────
logfc_file <- "~/projects2/eleo/ScRNA/metodologia/resultados/12_heatmaps/0.5N_vs_5N/tabla_log2FC_fc1_padj_005.tsv"
output_dir <- "~/projects2/eleo/ScRNA/metodologia/resultados/12_heatmaps/0.5N_vs_5N/coexpression"
contrast   <- "0.5N_vs_5N"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── LOAD MATRIX ───────────────────────────────────────────────────────────────
df  <- read.table(logfc_file, header = TRUE, sep = "\t", check.names = FALSE)
mat <- as.matrix(df[, -1])
rownames(mat) <- df[[1]]
colnames(mat) <- gsub(paste0("_", contrast, "$"), "", colnames(mat))
mat[is.na(mat)] <- 0

cat("Matrix dimensions:", nrow(mat), "genes x", ncol(mat), "cell types\n")

# =============================================================================
# PART 1 — HIERARCHICAL CLUSTERS (hclust + cutreeDynamic)
# =============================================================================
cat("\n── Part 1: hierarchical clustering...\n")

dist_rows  <- dist(mat, method = "euclidean")
hc_rows    <- hclust(dist_rows, method = "complete")

hclust_clusters <- cutreeDynamic(
  dendro            = hc_rows,
  distM             = as.matrix(dist_rows),
  deepSplit         = 1,
  minClusterSize    = 10,
  pamRespectsDendro = FALSE
)
names(hclust_clusters) <- rownames(mat)

cat("Hierarchical clusters found:", length(unique(hclust_clusters[hclust_clusters > 0])), "\n")
print(table(hclust_clusters))

# =============================================================================
# PART 2 — WGCNA MODULES (TOM-based)
# =============================================================================
cat("\n── Part 2: WGCNA modules...\n")

enableWGCNAThreads(2)

# Soft threshold power
powers    <- c(1:10, seq(12, 20, 2))
sft       <- pickSoftThreshold(t(mat), powerVector = powers, verbose = 0,
                                networkType = "signed")
power_val <- sft$powerEstimate
if (is.na(power_val)) power_val <- 6
cat("Soft threshold power selected:", power_val, "\n")

# Adjacency + TOM
adj <- adjacency(t(mat), power = power_val, type = "signed")
TOM <- TOMsimilarity(adj, TOMType = "signed")
dissTOM <- 1 - TOM
rownames(dissTOM) <- colnames(dissTOM) <- rownames(mat)

# Cluster on TOM dissimilarity
gene_tree <- hclust(as.dist(dissTOM), method = "average")

wgcna_modules <- cutreeDynamic(
  dendro            = gene_tree,
  distM             = dissTOM,
  deepSplit         = 1,
  minClusterSize    = 10,
  pamRespectsDendro = FALSE
)
names(wgcna_modules) <- rownames(mat)

# Convert module numbers to colors (WGCNA standard)
module_colors <- labels2colors(wgcna_modules)
names(module_colors) <- rownames(mat)

cat("WGCNA modules found:", length(unique(module_colors[module_colors != "grey"])), "\n")
print(table(module_colors))

# =============================================================================
# PART 3 — DOUBLE ANNOTATION HEATMAP
# =============================================================================
cat("\n── Part 3: double annotation heatmap...\n")

# Palettes
n_hclust <- length(unique(hclust_clusters[hclust_clusters > 0]))
pal_hclust <- setNames(
  colorRampPalette(brewer.pal(min(n_hclust, 12), "Dark2"))(n_hclust),
  sort(unique(hclust_clusters[hclust_clusters > 0]))
)
pal_hclust["0"] <- "grey85"   # unassigned

unique_mods <- unique(module_colors)
pal_wgcna   <- setNames(unique_mods, unique_mods)  # WGCNA colors are already named

# Left annotations
left_ha <- rowAnnotation(
  `Hierarchical` = as.character(hclust_clusters),
  `WGCNA`        = module_colors,
  col = list(
    `Hierarchical` = pal_hclust,
    `WGCNA`        = pal_wgcna
  ),
  annotation_name_gp  = gpar(fontsize = 11, fontface = "bold"),
  annotation_width    = unit(c(0.4, 0.4), "cm")
)

# Order rows by WGCNA module, then by gene_tree order within each module
wgcna_order <- gene_tree$order  # WGCNA dendrogram order
genes_ordered <- rownames(mat)[wgcna_order]
module_order  <- module_colors[wgcna_order]

# Sort by module name so same-module genes are grouped, preserving intra-module order
row_order <- genes_ordered[order(module_order, method = "radix")]

# Heatmap
ht <- Heatmap(
  mat,
  name = "log2FC",
  col  = colorRamp2(c(-5, 0, 5), c("blue", "black", "yellow")),

  cluster_rows    = FALSE,
  row_order       = row_order,
  cluster_columns = FALSE,

  left_annotation   = left_ha,
  show_row_names    = FALSE,
  show_column_names = TRUE,
  column_names_gp   = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot  = 45,

  column_title    = sprintf("log2FC — %s  (%d genes)", contrast, nrow(mat)),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),

  heatmap_legend_param = list(
    title         = "log2FC",
    title_gp      = gpar(fontsize = 11, fontface = "bold"),
    legend_height = unit(4, "cm")
  )
)

pdf(file.path(output_dir, paste0("heatmap_double_annot_", contrast, ".pdf")),
    width = 14, height = 18)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 10), "mm"))
dev.off()

# Save cluster assignments
assignments <- data.frame(
  gene_id     = rownames(mat),
  hclust_cluster = hclust_clusters,
  wgcna_module   = module_colors
)
write.table(assignments,
            file.path(output_dir, paste0("cluster_assignments_", contrast, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n✓ Done. Output in:", output_dir, "\n")
