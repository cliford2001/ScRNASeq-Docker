# =============================================================================
# TEST: 5 clustering methods comparison on log2FC matrix
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(cluster)    # PAM
})

# ── INPUT ─────────────────────────────────────────────────────────────────────
logfc_file <- "~/projects2/eleo/ScRNA/metodologia/resultados/12_heatmaps/0.5N_vs_5N/tabla_log2FC_fc1_padj_005.tsv"
output_dir <- "~/projects2/eleo/ScRNA/metodologia/resultados/12_heatmaps/0.5N_vs_5N/coexpression"
contrast   <- "0.5N_vs_5N"
k_fixed    <- 8     # number of clusters for methods that need k (ward + kmeans + PAM)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── LOAD MATRIX ───────────────────────────────────────────────────────────────
df  <- read.table(logfc_file, header = TRUE, sep = "\t", check.names = FALSE)
mat <- as.matrix(df[, -1])
rownames(mat) <- df[[1]]
colnames(mat) <- gsub(paste0("_", contrast, "$"), "", colnames(mat))
mat[is.na(mat)] <- 0

cat("Matrix:", nrow(mat), "genes x", ncol(mat), "cell types\n")

# =============================================================================
# METHOD 1 — hclust complete + cutreeDynamic
# =============================================================================
cat("\n── Method 1: hclust complete + cutreeDynamic...\n")

dist_mat     <- dist(mat, method = "euclidean")
hc_complete  <- hclust(dist_mat, method = "complete")

clust_m1 <- cutreeDynamic(
  dendro = hc_complete, distM = as.matrix(dist_mat),
  deepSplit = 1, minClusterSize = 10, pamRespectsDendro = FALSE
)
names(clust_m1) <- rownames(mat)
cat("Clusters:", length(unique(clust_m1[clust_m1 > 0])), "\n")

# =============================================================================
# METHOD 2 — hclust ward.D2 + cutree(k)
# =============================================================================
cat("\n── Method 2: hclust ward.D2 + cutree(k =", k_fixed, ")...\n")

hc_ward  <- hclust(dist_mat, method = "ward.D2")
clust_m2 <- cutree(hc_ward, k = k_fixed)
names(clust_m2) <- rownames(mat)
cat("Clusters:", length(unique(clust_m2)), "\n")

# =============================================================================
# METHOD 3 — WGCNA (TOM)
# =============================================================================
cat("\n── Method 3: WGCNA (TOM)...\n")

enableWGCNAThreads(2)

sft       <- pickSoftThreshold(t(mat), powerVector = c(1:10, seq(12,20,2)),
                                verbose = 0, networkType = "signed")
power_val <- sft$powerEstimate
if (is.na(power_val)) power_val <- 6
cat("Soft power:", power_val, "\n")

adj      <- adjacency(t(mat), power = power_val, type = "signed")
TOM      <- TOMsimilarity(adj, TOMType = "signed")
dissTOM  <- 1 - TOM
rownames(dissTOM) <- colnames(dissTOM) <- rownames(mat)

gene_tree   <- hclust(as.dist(dissTOM), method = "average")
clust_m3_num <- cutreeDynamic(
  dendro = gene_tree, distM = dissTOM,
  deepSplit = 1, minClusterSize = 10, pamRespectsDendro = FALSE
)
clust_m3 <- labels2colors(clust_m3_num)
names(clust_m3) <- rownames(mat)
cat("WGCNA modules:", length(unique(clust_m3[clust_m3 != "grey"])), "\n")

# =============================================================================
# METHOD 4 — K-means
# =============================================================================
cat("\n── Method 4: K-means (k =", k_fixed, ")...\n")

set.seed(1807)
km       <- kmeans(mat, centers = k_fixed, nstart = 25, iter.max = 100)
clust_m4 <- km$cluster
names(clust_m4) <- rownames(mat)
cat("Clusters:", length(unique(clust_m4)), "\n")

# =============================================================================
# METHOD 5 — PAM (Partitioning Around Medoids)
# =============================================================================
cat("\n── Method 5: PAM (k =", k_fixed, ")...\n")

pam_res  <- pam(mat, k = k_fixed, metric = "euclidean")
clust_m5 <- pam_res$clustering
names(clust_m5) <- rownames(mat)
cat("Clusters:", length(unique(clust_m5)), "\n")

# =============================================================================
# HEATMAP WITH 5 ANNOTATIONS
# =============================================================================
cat("\n── Building comparison heatmap...\n")

make_pal <- function(clusters, palette = "Set2") {
  u <- sort(unique(clusters))
  setNames(colorRampPalette(brewer.pal(min(length(u), 8), palette))(length(u)), u)
}

# WGCNA colors are self-named
pal_m3 <- setNames(unique(clust_m3), unique(clust_m3))

left_ha <- rowAnnotation(
  `1-hclust`  = as.character(clust_m1),
  `2-ward`    = as.character(clust_m2),
  `3-WGCNA`   = clust_m3,
  `4-kmeans`  = as.character(clust_m4),
  `5-PAM`     = as.character(clust_m5),
  col = list(
    `1-hclust` = make_pal(as.character(clust_m1), "Dark2"),
    `2-ward`   = make_pal(as.character(clust_m2), "Set1"),
    `3-WGCNA`  = pal_m3,
    `4-kmeans` = make_pal(as.character(clust_m4), "Paired"),
    `5-PAM`    = make_pal(as.character(clust_m5), "Set3")
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_width   = unit(rep(0.4, 5), "cm")
)

# Order rows by WGCNA module (same as previous script)
wgcna_order <- gene_tree$order
genes_ordered <- rownames(mat)[wgcna_order]
module_order  <- clust_m3[wgcna_order]
row_order     <- genes_ordered[order(module_order)]

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

  column_title    = sprintf("5 clustering methods — %s  (%d genes)", contrast, nrow(mat)),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),

  heatmap_legend_param = list(
    title         = "log2FC",
    title_gp      = gpar(fontsize = 11, fontface = "bold"),
    legend_height = unit(4, "cm")
  )
)

pdf(file.path(output_dir, paste0("heatmap_5methods_", contrast, ".pdf")),
    width = 16, height = 18)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 10), "mm"))
dev.off()

# Save all assignments
assignments <- data.frame(
  gene_id       = rownames(mat),
  hclust_dyn    = clust_m1,
  ward_cutree   = clust_m2,
  wgcna_module  = clust_m3,
  kmeans        = clust_m4,
  pam           = clust_m5
)
write.table(assignments,
            file.path(output_dir, paste0("cluster_assignments_5methods_", contrast, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n✓ Done. Output in:", output_dir, "\n")
