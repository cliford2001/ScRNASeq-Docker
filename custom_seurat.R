# =============================================================================
# custom_seurat.R — Custom Seurat plot utilities
# =============================================================================

#' Cell Count Bar Chart per Cluster and Sample
#'
#' Produces a two-panel figure: (left) fraction of cells per dataset within
#' each identity/cluster; (right) total cell count per cluster on a log10 scale
#' with count labels.
#'
#' @param srat Seurat object. Uses active Idents() as the grouping variable.
#' @return A patchwork ggplot object.
#' @export
plot_integrated_clusters <- function(srat) {

  count_table <- table(Idents(srat), srat@meta.data$orig.ident)
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$identity <- rownames(count_mtx)
  melt_mtx    <- reshape2::melt(count_mtx)
  melt_mtx$identity <- as.factor(melt_mtx$identity)

  cluster_size <- aggregate(value ~ identity, data = melt_mtx, FUN = sum)

  sorted_labels        <- cluster_size$identity[order(cluster_size$value, decreasing = TRUE)]
  cluster_size$identity <- factor(cluster_size$identity, levels = sorted_labels)
  melt_mtx$identity    <- factor(melt_mtx$identity,     levels = sorted_labels)

  colnames(melt_mtx)[2] <- "dataset"

  p_fraction <- ggplot(melt_mtx, aes(x = identity, y = value, fill = dataset)) +
    geom_bar(position = "fill", stat = "identity") +
    theme_bw() + coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    ylab("Fraction of cells per dataset") +
    xlab("Identity") +
    theme(legend.position = "top")

  p_count <- ggplot(cluster_size, aes(y = identity, x = value)) +
    geom_bar(stat = "identity", fill = "grey60") +
    geom_text(aes(label = value),
              position = position_identity(),
              hjust = 1.1, size = 8, color = "blue") +
    scale_x_log10() +
    xlab("Cells per identity (log10 scale)") +
    ylab("") +
    theme_bw()

  p_fraction + p_count + plot_layout(widths = c(3, 1))
}
