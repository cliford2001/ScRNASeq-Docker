# =============================================================================
# ScRNA-Seq Single-Cell Analysis - Custom Functions Library
# =============================================================================
# Author: Ellie Malcolm et al.
# Date: 2025-03
# Description: Reusable functions for QC, preprocessing, clustering, and DE analysis
# =============================================================================

# =============================================================================
# TABLE OF CONTENTS
# =============================================================================
#
#  1. QC AND VISUALIZATION FUNCTIONS
#     - load_cellbender_filtered_h5
#     - plot_qc_violin_grid
#     - summarize_nfeature_plot
#
#  2. PREPROCESSING AND DOUBLET DETECTION
#     - preprocesar_y_doubletfinder
#     - doubletfinder_pipeline
#     - load_sample          (load + annotate only, no filtering)
#     - filter_sample        (filter + DoubletFinder on annotated object)
#     - process_sample       (shortcut: load_sample + filter_sample)
#
#  3. BULK / PSEUDOBULK UTILITIES
#     - normalize_bulk_pseudobulk
#     - clasificar_residuos
#     - generate_pseudobulk
#     - plot_replicate_correlation
#
#  4. SEURAT UTILITIES
#     - unificar_nombres
#     - mostrar_tabla
#     - export_to_scanpy
#     - safe_vln
#     - unir_layers_counts
#
#  5. ANNOTATION
#     - find_markers
#     - annotate_by_markers
#     - annotate_by_reference
#     - subcluster_cell_type
#
#  6. PSEUDOBULK, DESEQ2, VOLCANO, HEATMAP
#     - assign_pseudo_replicates
#     - run_pseudobulk
#     - correr_deseq2
#     - plot_volcano
#     - procesar_deseq2_resultado
#     - plot_heatmap
#     - plot_marker_dotplot
#
#  7. GO ENRICHMENT
#     - correr_enriquecimiento_go
#     - podar_go
#     - graficar_go_balones
#
# =============================================================================


# =============================================================================
# 1. QC AND VISUALIZATION FUNCTIONS
# =============================================================================

#' Load CellBender Filtered HDF5 Data
#'
#' Reads filtered expression matrix from CellBender HDF5 output.
#'
#' @param h5_path Path to filtered HDF5 file.
#' @param project  Project name for Seurat object metadata.
#' @return A Seurat object containing raw counts.
#' @export
load_cellbender_filtered_h5 <- function(h5_path, project = "Sample") {

  f <- H5File$new(h5_path, mode = "r")

  message("Reading CSR components...")
  data    <- f[["matrix/data"]]$read()
  indices <- f[["matrix/indices"]]$read()
  indptr  <- f[["matrix/indptr"]]$read()
  shape   <- f[["matrix/shape"]]$read()

  message("Reading gene IDs and barcodes...")
  gene_ids <- f[["matrix/features/id"]]$read()
  barcodes <- f[["matrix/barcodes"]]$read()

  message("Creating sparse gene x cell matrix...")
  mat <- new("dgCMatrix",
             x         = as.numeric(data),
             i         = indices,
             p         = indptr,
             Dim       = shape,
             Dimnames  = list(gene_ids, barcodes))

  seu <- CreateSeuratObject(counts = mat, project = project)

  return(seu)
}


#' QC Violin Plot Grid
#'
#' Visualizes nFeature_RNA, nCount_RNA, percent.mt, and percent.cp by condition.
#'
#' @param obj1  Seurat object.
#' @param label Condition label.
#' @param color Color for plotting.
#' @return A ggplot object.
#' @export
plot_qc_batch <- function(seurat_list, colors, file) {
  plots <- imap(seurat_list, ~ plot_qc_violin_grid(.x, .y, colors[[.y]]))
  save_qc(plots, file)
  invisible(plots)
}


plot_qc_violin_grid <- function(obj1, label, color) {

  n1        <- ncol(obj1)
  obj1$cond <- label

  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  if ("percent.cp" %in% colnames(obj1@meta.data)) {
    features <- c(features, "percent.cp")
  }

  p1 <- VlnPlot(obj1,
                features = features,
                pt.size  = 0.1,
                ncol     = length(features),
                group.by = "cond",
                cols     = color) +
    ggtitle(paste0(label, " (", n1, " cells)")) +
    theme_minimal(base_size = 12)

  return(p1)
}


#' Summary of nFeature_RNA Distribution
#'
#' Creates a boxplot alongside quartile and quintile summary tables.
#'
#' @param obj_list  List of Seurat objects.
#' @param labels Labels for each object.
#' @param colores   Color vector (named or positional).
#' @export
summarize_nfeature_plot <- function(obj_list, labels = NULL, colores = NULL) {

  if (is.null(labels)) labels <- paste0("Group", seq_along(obj_list))
  if (length(labels) != length(obj_list)) stop("Labels must match objects.")

  if (is.null(colores)) {
    colores       <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")[1:length(obj_list)]
    names(colores) <- labels
  }

  lista_df <- lapply(seq_along(obj_list), function(i) {
    obj <- obj_list[[i]]
    data.frame(nFeature_RNA = obj@meta.data$nFeature_RNA, grupo = labels[i])
  })

  meta_comb       <- bind_rows(lista_df)
  meta_comb$grupo <- factor(meta_comb$grupo, levels = labels)

  p_box <- ggplot(meta_comb, aes(x = grupo, y = nFeature_RNA, fill = grupo)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
    scale_fill_manual(values = colores) +
    labs(title = "nFeature_RNA Distribution", x = "Condition", y = "nFeature_RNA") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  cuartiles <- meta_comb %>%
    group_by(grupo) %>%
    summarise(
      Min    = quantile(nFeature_RNA, 0),
      Q1     = quantile(nFeature_RNA, 0.25),
      Median = quantile(nFeature_RNA, 0.5),
      Q3     = quantile(nFeature_RNA, 0.75),
      Max    = quantile(nFeature_RNA, 1),
      .groups = "drop"
    ) %>%
    arrange(factor(grupo, levels = labels))

  quintiles <- meta_comb %>%
    group_by(grupo) %>%
    summarise(
      `0%`   = quantile(nFeature_RNA, 0.0),
      `20%`  = quantile(nFeature_RNA, 0.2),
      `40%`  = quantile(nFeature_RNA, 0.4),
      `60%`  = quantile(nFeature_RNA, 0.6),
      `80%`  = quantile(nFeature_RNA, 0.8),
      `100%` = quantile(nFeature_RNA, 1.0),
      .groups = "drop"
    ) %>%
    arrange(factor(grupo, levels = labels))

  tabla_cuartiles <- tableGrob(cuartiles)
  tabla_quintiles <- tableGrob(quintiles)

  panel_tablas <- plot_grid(
    ggdraw() + draw_label("Quartiles", fontface = "bold", size = 13),
    ggdraw() + draw_grob(tabla_cuartiles),
    ggdraw() + draw_label("Quintiles", fontface = "bold", size = 13),
    ggdraw() + draw_grob(tabla_quintiles),
    ncol        = 1,
    rel_heights = c(0.15, 1, 0.15, 1)
  )

  final_plot <- plot_grid(p_box, panel_tablas, ncol = 2, rel_widths = c(1.5, 1))
  print(final_plot)
}


# =============================================================================
# 2. PREPROCESSING AND DOUBLET DETECTION
# =============================================================================

#' Preprocessing + DoubletFinder Pipeline
#'
#' Normalizes, scales, runs PCA, and performs doublet detection via DoubletFinder.
#'
#' @param seurat_obj             Seurat object.
#' @param pcs                    PCs to use (e.g. 1:20).
#' @param expected_doublet_rate  Expected doublet rate (default 0.075).
#' @param project_id             Project ID label.
#' @return Seurat object with doublet classifications in metadata.
#' @export
preprocesar_y_doubletfinder <- function(seurat_obj,
                                        pcs                   = 1:20,
                                        expected_doublet_rate = 0.075,
                                        project_id            = "sample") {

  message("Normalizing...")
  seurat_obj <- NormalizeData(seurat_obj)

  message("Finding variable features...")
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

  message("Scaling...")
  seurat_obj <- ScaleData(seurat_obj)

  message("PCA...")
  seurat_obj <- RunPCA(seurat_obj, npcs = max(pcs))

  message("Running DoubletFinder...")
  sweep.res   <- paramSweep(seurat_obj, PCs = pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn       <- find.pK(sweep.stats)
  best.pK     <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  nExp        <- round(expected_doublet_rate * ncol(seurat_obj))

  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs       = pcs,
    pN        = 0.25,
    pK        = best.pK,
    nExp      = nExp,
    reuse.pANN = FALSE,
    sct       = FALSE
  )

  return(seurat_obj)
}


#' Full DoubletFinder Pipeline with Clustering
#'
#' Comprehensive doublet detection pipeline including normalization, PCA,
#' clustering, parameter sweep, and optional singlet filtering.
#'
#' @param obj              Seurat object.
#' @param etiqueta         Sample label for messages.
#' @param PCs              PCs for analysis (e.g. 1:20).
#' @param resolution       Leiden/Louvain resolution.
#' @param return_singlets  If TRUE, return only singlets.
#' @param sct              Whether to use SCT normalization.
#' @return Seurat object, optionally filtered to singlets.
#' @export
doubletfinder_pipeline <- function(obj,
                                   etiqueta        = "Sample",
                                   PCs             = 1:20,
                                   resolution      = 0.5,
                                   return_singlets = TRUE,
                                   sct             = FALSE) {

  message("Processing: ", etiqueta)

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = max(PCs))
  obj <- FindNeighbors(obj, dims = PCs)
  obj <- FindClusters(obj, resolution = resolution)

  sweep.res          <- paramSweep(obj, PCs = PCs, sct = sct)
  sweep.stats        <- summarizeSweep(sweep.res, GT = FALSE)
  sweep.stats$pK     <- as.numeric(as.character(sweep.stats$pK))
  sweep.stats$pN     <- as.numeric(as.character(sweep.stats$pN))

  best_row <- sweep.stats[which.max(sweep.stats$BCreal), ]
  best.pK  <- best_row$pK
  best.pN  <- best_row$pN
  nExp     <- round(best.pN * ncol(obj))

  message("Best pK: ", best.pK, ", pN: ", best.pN, ", nExp: ", nExp)

  obj <- doubletFinder(
    obj,
    PCs        = PCs,
    pN         = best.pN,
    pK         = best.pK,
    nExp       = nExp,
    reuse.pANN = NULL,
    sct        = sct
  )

  # Fix any data.frame columns returned by doubletFinder
  for (col in colnames(obj@meta.data)) {
    if (is.data.frame(obj@meta.data[[col]])) {
      message("Fixing column: ", col)
      obj@meta.data[[col]] <- obj@meta.data[[col]][, 1]
    }
  }

  df_col        <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  obj$doublet_class <- obj[[df_col]]

  tab <- table(obj$doublet_class)
  message("Summary for ", etiqueta, ":")
  print(tab)

  if (return_singlets) {
    obj <- subset(obj, subset = doublet_class == "Singlet")
    message("Retained singlets: ", ncol(obj))
  }

  return(obj)
}


#' Load and Annotate a Single Sample
#'
#' Loads a CellBender h5 file and computes mitochondrial / chloroplast
#' percentages. No filtering or doublet detection — use this to inspect raw
#' QC metrics before deciding thresholds.
#'
#' @param sample_info Named list with fields: file, label, condition.
#' @param mt_pattern  Regex for mitochondrial genes (e.g. "^MT-", "^ATMG").
#' @param cp_pattern  Regex for chloroplast genes (e.g. "^ATCG"); NULL to skip.
#' @return Seurat object with percent.mt (and percent.cp) in metadata.
#' @export
load_sample <- function(sample_info,
                        mt_pattern = "^ATMG",
                        cp_pattern = "^ATCG") {

  obj <- load_cellbender_filtered_h5(sample_info$file, sample_info$label)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  if (!is.null(cp_pattern))
    obj[["percent.cp"]] <- PercentageFeatureSet(obj, pattern = cp_pattern)

  obj <- RenameCells(obj, add.cell.id = sample_info$condition)
  obj$condition <- sample_info$condition

  return(obj)
}


#' Filter and Run DoubletFinder on an Annotated Sample
#'
#' Applies QC thresholds to an already-annotated Seurat object (output of
#' load_sample) and optionally runs DoubletFinder.
#'
#' @param obj             Seurat object with percent.mt (and percent.cp).
#' @param min_features    Minimum nFeature_RNA.
#' @param max_features    Maximum nFeature_RNA.
#' @param min_counts      Minimum nCount_RNA.
#' @param max_counts      Maximum nCount_RNA.
#' @param max_mt          Maximum mitochondrial percent.
#' @param max_cp          Maximum chloroplast percent (ignored if percent.cp absent).
#' @param run_doubletfinder Whether to run DoubletFinder (default TRUE).
#' @return Filtered Seurat object.
#' @export
filter_seurat_samples <- function(seurat_list, ...) {
  result        <- lapply(seurat_list, filter_sample, ...)
  names(result) <- names(seurat_list)
  result
}


filter_sample <- function(obj,
                          min_features      = 200,
                          max_features      = Inf,
                          min_counts        = 0,
                          max_counts        = Inf,
                          max_mt            = 5,
                          max_cp            = 100,
                          run_doubletfinder = TRUE) {

  has_cp <- "percent.cp" %in% colnames(obj@meta.data)

  if (has_cp) {
    obj <- subset(obj, subset =
                    nFeature_RNA > min_features &
                    nFeature_RNA < max_features &
                    nCount_RNA   > min_counts   &
                    nCount_RNA   < max_counts   &
                    percent.mt   < max_mt       &
                    percent.cp   < max_cp)
  } else {
    obj <- subset(obj, subset =
                    nFeature_RNA > min_features &
                    nFeature_RNA < max_features &
                    nCount_RNA   > min_counts   &
                    nCount_RNA   < max_counts   &
                    percent.mt   < max_mt)
  }

  if (run_doubletfinder)
    obj <- doubletfinder_pipeline(obj, etiqueta = Project(obj))

  return(obj)
}


#' Load, Annotate, Filter and Run DoubletFinder (full pipeline shortcut)
#'
#' Convenience wrapper that calls load_sample() then filter_sample().
#' Useful when you do not need to inspect raw QC plots before filtering.
#'
#' @inheritParams load_sample
#' @inheritParams filter_sample
#' @return Filtered Seurat object.
#' @export
process_sample <- function(sample_info,
                           mt_pattern        = "^ATMG",
                           cp_pattern        = "^ATCG",
                           min_features      = 200,
                           max_features      = Inf,
                           min_counts        = 0,
                           max_counts        = Inf,
                           max_mt            = 5,
                           max_cp            = 100,
                           run_doubletfinder = TRUE) {

  obj <- load_sample(sample_info, mt_pattern = mt_pattern, cp_pattern = cp_pattern)
  obj <- filter_sample(obj,
                       min_features      = min_features,
                       max_features      = max_features,
                       min_counts        = min_counts,
                       max_counts        = max_counts,
                       max_mt            = max_mt,
                       max_cp            = max_cp,
                       run_doubletfinder = run_doubletfinder)
  return(obj)
}


# =============================================================================
# 3. BULK / PSEUDOBULK UTILITIES
# =============================================================================

#' Normalize Pseudobulk vs Bulk Counts
#'
#' DESeq2-based normalization followed by log2 transformation of pseudobulk
#' and bulk count vectors, restricted to their common gene set.
#'
#' @param pseudobulk_counts Named numeric vector of pseudobulk counts.
#' @param bulk_counts        Named numeric vector of bulk counts.
#' @return Data frame with columns gene, pseudobulk, bulk (log2-normalized).
#' @export
normalize_bulk_pseudobulk <- function(pseudobulk_counts, bulk_counts) {

  common_genes <- intersect(names(pseudobulk_counts), names(bulk_counts))

  if (length(common_genes) < 10) {
    stop("Too few common genes.")
  }

  counts_matrix <- data.frame(
    pseudobulk = round(pseudobulk_counts[common_genes]),
    bulk       = round(bulk_counts[common_genes])
  )
  rownames(counts_matrix) <- common_genes

  condition <- factor(c("pseudobulk", "bulk"))
  col_data  <- data.frame(condition = condition)

  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData   = col_data,
                                design    = ~ condition)
  dds        <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)

  log_norm_counts <- log2(norm_counts + 1)

  df <- data.frame(
    gene       = rownames(log_norm_counts),
    pseudobulk = log_norm_counts[, "pseudobulk"],
    bulk       = log_norm_counts[, "bulk"]
  )

  return(df)
}


#' Classify Genes by Residuals
#'
#' Fits a linear model (bulk ~ pseudobulk) and classifies genes as
#' Upregulated, Downregulated, or Consistent based on residual magnitude.
#'
#' @param df     Data frame with columns pseudobulk and bulk.
#' @param umbral Residual threshold for classification.
#' @return The input data frame augmented with residuals and status columns.
#' @export
clasificar_residuos <- function(df, umbral = 5) {

  modelo      <- lm(bulk ~ pseudobulk, data = df)
  df$residuals <- resid(modelo)

  df$status <- case_when(
    df$residuals >  umbral ~ "Upregulated",
    df$residuals < -umbral ~ "Downregulated",
    TRUE                   ~ "Consistent"
  )

  return(df)
}


#' Generate Pseudobulk Counts Matrix
#'
#' Aggregates single-cell counts by a grouping variable. Optionally merges
#' replicates belonging to the same condition.
#'
#' @param seurat_obj        Seurat object.
#' @param group_by          Metadata column to group cells by.
#' @param merge_replicates  If TRUE, sum columns matching each unique condition.
#' @return Matrix (genes x samples), or a named list with by_sample and
#'         by_condition matrices when merge_replicates = TRUE.
#' @export
generate_pseudobulk <- function(seurat_obj,
                                group_by          = "orig.ident",
                                merge_replicates  = TRUE) {

  groups <- unique(seurat_obj@meta.data[[group_by]])
  cat("Generating pseudobulk for:", paste(groups, collapse = ", "), "\n")

  process_group <- function(group_name) {
    cells  <- subset(seurat_obj,
                     cells = colnames(seurat_obj)[seurat_obj@meta.data[[group_by]] == group_name])
    layers <- grep("^counts", Layers(cells[["RNA"]]), value = TRUE)

    if (length(layers) == 0) {
      counts <- GetAssayData(cells[["RNA"]], layer = "data")
    } else if (length(layers) == 1) {
      counts <- GetAssayData(cells[["RNA"]], layer = layers)
    } else {
      mats   <- lapply(layers, function(x) GetAssayData(cells[["RNA"]], layer = x))
      counts <- Reduce(RowMergeSparseMatrices, mats)
    }

    gene_sums <- Matrix::rowSums(counts)
    cat(" ", group_name, "->", ncol(cells), "celulas,", length(gene_sums), "genes\n")
    return(gene_sums)
  }

  pseudobulk_list       <- lapply(groups, process_group)
  names(pseudobulk_list) <- groups

  all_genes <- unique(unlist(lapply(pseudobulk_list, names)))

  pseudobulk_matrix <- sapply(pseudobulk_list, function(x) {
    v       <- x[all_genes]
    v[is.na(v)] <- 0
    return(v)
  })
  rownames(pseudobulk_matrix) <- all_genes

  if (merge_replicates) {
    conditions <- unique(seurat_obj$condition)

    merged_matrix <- sapply(conditions, function(cond) {
      cols <- grep(paste0("^", cond), colnames(pseudobulk_matrix), value = TRUE)
      if (length(cols) == 0) {
        cols <- colnames(pseudobulk_matrix)[grepl(cond, colnames(pseudobulk_matrix))]
      }
      if (length(cols) == 1) return(pseudobulk_matrix[, cols])
      return(rowSums(pseudobulk_matrix[, cols, drop = FALSE]))
    })
    colnames(merged_matrix) <- conditions

    cat("\nReplicas fusionadas:\n")
    print(colnames(merged_matrix))

    return(list(
      by_sample    = pseudobulk_matrix,
      by_condition = merged_matrix
    ))
  }

  return(pseudobulk_matrix)
}


#' Plot Replicate Correlation Heatmap
#'
#' Computes pairwise Pearson correlations across columns (samples/replicates)
#' of a pseudobulk count matrix and displays a pheatmap of the correlation
#' matrix.
#'
#' @param pseudobulk_mat A numeric matrix with genes as rows and samples as
#'   columns. Typically the output of generate_pseudobulk() or
#'   run_pseudobulk().
#' @param main           Title for the heatmap (default: "Replicate Correlation").
#' @return Invisible: the correlation matrix.
#' @export
plot_replicate_correlation <- function(pseudobulk_mat,
                                       main = "Replicate Correlation") {

  cor_mat <- cor(pseudobulk_mat, method = "pearson", use = "pairwise.complete.obs")

  p <- pheatmap(cor_mat,
                display_numbers = TRUE,
                number_format   = "%.2f",
                color           = colorRampPalette(c("white", "steelblue"))(50),
                main            = main,
                border_color    = NA)

  invisible(p)
}


# =============================================================================
# 4. SEURAT UTILITIES
# =============================================================================

#' Unify Ident Names
#'
#' Removes numeric suffixes (e.g. ".1", "_2") from cluster identity names.
#'
#' @param obj Seurat object.
#' @return Seurat object with updated Idents.
#' @export
unificar_nombres <- function(obj) {

  old_levels <- levels(obj)
  new_levels <- gsub("[._][0-9]+$", "", old_levels)
  new_ids    <- setNames(new_levels, old_levels)
  obj        <- RenameIdents(obj, new_ids)

  return(obj)
}


#' Display Annotation Table as Grid
#'
#' Creates and displays a cell type count comparison table using grid graphics.
#'
#' @param filtered_vec   Filtered annotation vector.
#' @param cellbender_vec CellBender annotation vector.
#' @param titulo         Table title.
#' @export
mostrar_tabla <- function(filtered_vec, cellbender_vec, titulo = "Annotations") {

  t1       <- table(filtered_vec)
  t2       <- table(cellbender_vec)
  all_types <- union(names(t1), names(t2))

  df <- data.frame(
    celltype   = all_types,
    filtered   = as.integer(t1[all_types]),
    cellbender = as.integer(t2[all_types]),
    stringsAsFactors = FALSE
  )
  df[is.na(df)] <- 0

  total_row <- data.frame(
    celltype   = "Total",
    filtered   = sum(df$filtered),
    cellbender = sum(df$cellbender),
    stringsAsFactors = FALSE
  )
  df <- rbind(df, total_row)

  grid.newpage()
  grid.draw(tableGrob(df, rows = NULL, theme = ttheme_minimal()))
}


#' Export Seurat to Scanpy h5ad Format
#'
#' Converts a Seurat object to SingleCellExperiment and writes it as an h5ad
#' file compatible with Scanpy. Uses zellkonverter if available, falling back
#' to SeuratDisk.
#'
#' @param seurat_obj Seurat object.
#' @param outfile    Output file path (should end in .h5ad).
#' @param assay_name Assay to export (default "RNA").
#' @param use_reduc  Reductions to include (default c("pca","umap","harmony")).
#' @param X_name     Layer name to store as .X in Scanpy (default "logcounts").
#' @param overwrite  Overwrite existing file (default TRUE).
#' @return Invisible SingleCellExperiment.
#' @export
export_to_scanpy <- function(seurat_obj,
                                 outfile,
                                 assay_name = "RNA",
                                 use_reduc  = c("pca", "umap", "harmony"),
                                 X_name     = "logcounts",
                                 overwrite  = TRUE) {

  stopifnot(inherits(seurat_obj, "Seurat"))
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)

  # ── Deep clean before conversion ─────────────────────────────────────────────
  # These slots often contain non-serialisable objects that break h5ad export.
  message("Cleaning Seurat object before conversion...")
  seurat_obj@misc  <- list()
  seurat_obj@tools <- list()

  cols_df   <- grep("DF.classifications|^pANN_|^doublet_class",
                    colnames(seurat_obj@meta.data), value = TRUE)
  if (length(cols_df) > 0) seurat_obj@meta.data[, cols_df] <- NULL

  # Deduplicate column names (can appear after multi-sample merge)
  seurat_obj@meta.data <- seurat_obj@meta.data[
    , !duplicated(names(seurat_obj@meta.data)), drop = FALSE]

  seurat_obj[[assay_name]]@meta.data <- data.frame(row.names = rownames(seurat_obj))

  for (rd in names(seurat_obj@reductions)) {
    seurat_obj@reductions[[rd]]@misc <- list()
  }

  # ── Convert to SCE ────────────────────────────────────────────────────────────
  message("Converting Seurat -> SCE...")
  sce <- as.SingleCellExperiment(seurat_obj, assay = assay_name)

  if (is.null(rownames(sce))) stop("Object missing rownames.")
  rownames(sce) <- make.unique(rownames(sce))

  if (is.null(colnames(sce)) || any(!nzchar(colnames(sce)))) {
    colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
  }
  stopifnot(!anyDuplicated(colnames(sce)))

  # ── Assays ────────────────────────────────────────────────────────────────────
  message("Extracting counts and logcounts...")
  assay(sce, "counts")    <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
  assay(sce, "logcounts") <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = "data")

  # ── Cell metadata ─────────────────────────────────────────────────────────────
  # as.SingleCellExperiment already transfers meta.data → colData.
  # Deduplicate any repeated columns (can arise after multi-sample merge).
  cd <- as.data.frame(colData(sce))
  colData(sce) <- S4Vectors::DataFrame(cd[, !duplicated(names(cd)), drop = FALSE])

  # ── Reductions ────────────────────────────────────────────────────────────────
  message("Exporting reductions...")
  reds <- Seurat::Reductions(seurat_obj)
  for (red in use_reduc) {
    if (red %in% reds) {
      message("   Including: ", red)
      reducedDims(sce)[[toupper(red)]] <- seurat_obj@reductions[[red]]@cell.embeddings
    }
  }

  if (file.exists(outfile) && overwrite) file.remove(outfile)

  ok <- FALSE
  if (requireNamespace("anndataR", quietly = TRUE)) {
    message("Writing h5ad with anndataR (native R)...")
    adata <- anndataR::as_AnnData(sce)
    anndataR::write_h5ad(adata, path = outfile)
    ok <- TRUE
  } else if (requireNamespace("zellkonverter", quietly = TRUE)) {
    message("Writing h5ad with zellkonverter...")
    zellkonverter::writeH5AD(sce, file = outfile, X_name = X_name)
    ok <- TRUE
  } else if (requireNamespace("SeuratDisk", quietly = TRUE)) {
    message("Using SeuratDisk fallback...")
    tmp_h5seu <- file.path(tempdir(), paste0(basename(outfile), ".h5seurat"))
    if (file.exists(tmp_h5seu)) file.remove(tmp_h5seu)
    SeuratDisk::SaveH5Seurat(seurat_obj, filename = tmp_h5seu, overwrite = TRUE)
    SeuratDisk::Convert(source = tmp_h5seu, dest = "h5ad", overwrite = TRUE)
    gen_h5ad <- sub("\\.h5seurat$", ".h5ad", tmp_h5seu)
    if (!file.exists(gen_h5ad)) stop("Failed to generate h5ad.")
    file.rename(gen_h5ad, outfile)
    ok <- TRUE
  } else {
    stop("Install 'anndataR' to export h5ad: BiocManager::install('anndataR')")
  }

  if (!ok || !file.exists(outfile)) stop("Export failed: ", outfile)

  message("Export complete: ", outfile)
  invisible(sce)
}


#' Safe VlnPlot for RMarkdown
#'
#' Wrapper around VlnPlot with custom fill colors, grouped by orig.ident.
#'
#' @param obj     Seurat object.
#' @param feature Gene or metadata feature to plot.
#' @param colors  Named or positional color palette.
#' @return A ggplot object.
#' @export
safe_vln <- function(obj, feature, colors) {

  p <- VlnPlot(
    obj,
    features = feature,
    group.by = "orig.ident",
    pt.size  = 0
  )
  p <- p + scale_fill_manual(values = colors)

  return(p)
}


#' Unify and Merge Seurat Layers
#'
#' Combines multiple sparse count matrices from different layers of the RNA
#' assay into a single merged matrix.
#'
#' @param obj   Seurat object.
#' @param capas Character vector of layer names to merge.
#' @return A merged sparse matrix.
#' @export
unir_layers_counts <- function(obj, capas) {

  if (length(capas) == 1) {
    return(GetAssayData(obj[["RNA"]], layer = capas))
  }

  message("Merging ", length(capas), " layers...")
  mats   <- lapply(capas, function(x) GetAssayData(obj[["RNA"]], layer = x))
  merged <- Reduce(RowMergeSparseMatrices, mats)

  return(merged)
}


# =============================================================================
# 5. ANNOTATION
# =============================================================================

#' Find Cluster Markers
#'
#' Runs FindAllMarkers on Seurat clusters, caching results to disk.
#'
#' @param seurat_obj      Seurat object.
#' @param output_file     Path to cache markers as TSV.
#' @param only_pos        Only return positive markers.
#' @param min_pct         Minimum cell fraction expressing the gene.
#' @param logfc_threshold Log fold-change threshold.
#' @param force           Recompute even if cache exists.
#' @return Data frame of markers.
#' @export
find_markers <- function(seurat_obj,
                         output_file     = "results/FindAllMarkers.tsv",
                         only_pos        = TRUE,
                         min_pct         = 0.25,
                         logfc_threshold = 0.25,
                         force           = FALSE) {

  seurat_obj        <- JoinLayers(seurat_obj)
  Idents(seurat_obj) <- "seurat_clusters"

  if (file.exists(output_file) && !force) {
    cat("Loading existing markers:", output_file, "\n")
    markers <- read.table(output_file, header = TRUE, sep = "\t", quote = "")
  } else {
    cat("Computing markers...\n")
    markers <- FindAllMarkers(
      seurat_obj,
      only.pos        = only_pos,
      min.pct         = min_pct,
      logfc.threshold = logfc_threshold
    )
    write.table(markers, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
    cat("Guardado en:", output_file, "\n")
  }

  return(markers)
}


#' Annotate Clusters by Marker List
#'
#' Crosses FindAllMarkers output with a reference marker table to assign cell
#' type labels to clusters.
#'
#' @param seurat_obj     Seurat object.
#' @param markers        Data frame from find_markers().
#' @param reference_file Path to reference table (gene | cell.types).
#'   If NULL, a file chooser dialog is shown.
#' @return Seurat object with celltype metadata and updated Idents.
#' @export
annotate_by_markers <- function(seurat_obj,
                                markers,
                                reference_file = NULL) {

  if (is.null(reference_file)) {
    reference_file <- file.choose(caption = "Select reference file (gene | cell.types)")
  }

  cat("Using reference:", reference_file, "\n")

  reference <- read.table(reference_file, header = TRUE, sep = "\t", quote = "")

  merged <- merge(markers, reference, by.x = "gene", by.y = "gene")
  merged <- merged[order(merged$cluster, merged$p_val_adj), ]
  merged <- merged[!duplicated(merged$cluster), ]

  cat("\nCoincidencias encontradas:\n")
  print(merged[, c("cluster", "gene", "cell.types")])

  Idents(seurat_obj)  <- "seurat_clusters"
  new_ids             <- merged$cell.types
  names(new_ids)      <- merged$cluster
  seurat_obj          <- RenameIdents(seurat_obj, new_ids)
  seurat_obj$celltype <- Idents(seurat_obj)

  cat("\nAnotacion final:\n")
  print(table(seurat_obj$celltype))

  return(seurat_obj)
}


#' Annotate Clusters by Reference Transfer
#'
#' Uses Seurat label transfer (FindTransferAnchors + TransferData) to project
#' cell type annotations from a reference Seurat object onto the query.
#'
#' @param seurat_obj   Query Seurat object.
#' @param reference_obj Reference Seurat object. If NULL, a file chooser is shown.
#' @param reference_col Metadata column in reference to transfer. If NULL,
#'   an interactive selection prompt is shown.
#' @param dims         Dimensions to use for anchor finding.
#' @return Seurat object with celltype_reference metadata column.
#' @export
annotate_by_reference <- function(seurat_obj,
                                  reference_obj = NULL,
                                  reference_col = NULL,
                                  dims          = 1:30) {

  if (is.null(reference_obj)) {
    ref_file      <- file.choose(caption = "Select reference Seurat object (.rds) (.rds)")
    cat("Loading reference:", ref_file, "\n")
    reference_obj <- readRDS(ref_file)
  }

  if (is.null(reference_col)) {
    cat("\nAvailable columns in reference:\n")
    cols <- colnames(reference_obj@meta.data)
    for (i in seq_along(cols)) {
      cat(" ", i, "->", cols[i], "\n")
    }
    selection     <- as.integer(readline("Select column number: "))
    reference_col <- cols[selection]
  }

  cat("Using column:", reference_col, "\n")

  reference_obj <- UpdateSeuratObject(reference_obj)
  reference_obj <- DietSeurat(reference_obj, layers = c("counts", "data"), dimreducs = NULL)
  reference_obj <- NormalizeData(reference_obj, verbose = FALSE)
  reference_obj <- FindVariableFeatures(reference_obj, verbose = FALSE)
  shared_var_features <- intersect(VariableFeatures(reference_obj), rownames(seurat_obj))
  cat("Shared variable features for transfer:", length(shared_var_features), "\n")
  VariableFeatures(reference_obj) <- shared_var_features
  reference_obj <- ScaleData(reference_obj, features = shared_var_features, verbose = FALSE)

  anchors <- FindTransferAnchors(
    reference = reference_obj,
    query     = seurat_obj,
    dims      = dims,
    reduction = "cca",
    features  = shared_var_features
  )

  predictions <- TransferData(
    anchorset        = anchors,
    refdata          = reference_obj@meta.data[[reference_col]],
    dims             = dims,
    weight.reduction = "cca"
  )

  seurat_obj$celltype_reference <- predictions$predicted.id

  cat("\nReference annotation:\n")
  print(table(seurat_obj$celltype_reference))

  return(seurat_obj)
}


#' Subcluster a Cell Type
#'
#' Subsets to a specific annotation, re-runs PCA/UMAP/clustering at the
#' given resolution.
#'
#' @param obj        Seurat object.
#' @param tipo       Cell type(s) to subset (must match values in annot_col).
#' @param annot_col  Metadata column holding cell-type labels.
#' @param resolution Clustering resolution.
#' @param dims       Dimensions for UMAP and neighbor finding.
#' @return Seurat object with cluster_subtipo metadata.
#' @export
#' Plot and Save a Subcluster UMAP
#'
#' Creates a labeled DimPlot for a subclustered object, saves it to disk,
#' and returns the plot invisibly for use in composite figures.
#'
#' @param obj        Subclustered Seurat object (output of subcluster_cell_type).
#' @param label      Cell-type label used as the plot title and to derive
#'   the output filename (spaces replaced with underscores).
#' @param output_dir Directory where the PDF will be saved.
#'
#' @return ggplot object (invisibly).
#' @export
plot_subcluster_umap <- function(obj, label, output_dir) {
  p <- DimPlot(obj, group.by = "cluster_subtipo", label = TRUE, raster = FALSE) +
    ggtitle(paste0(label, " — subclusters"))
  filename <- paste0("subcluster_", tolower(gsub(" ", "_", label)), ".pdf")
  ggsave(file.path(output_dir, filename), p, width = 10, height = 8, dpi = 300)
  invisible(p)
}


subcluster_cell_type <- function(obj, tipo, annot_col = "celltype_grouped",
                            resolution = 0.3, dims = 1:20) {

  sub <- subset(obj, cells = colnames(obj)[obj@meta.data[[annot_col]] %in% tipo])

  # For small subsets, recompute variable features
  ncells <- ncol(sub)
  npcs <- min(ncells - 1, 30)

  sub <- sub %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE) %>%
    RunUMAP(dims = dims, verbose = FALSE) %>%
    FindNeighbors(dims = dims, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE)

  sub$cluster_subtipo <- as.character(sub$seurat_clusters)

  return(sub)
}


#' Generate Marker Gene Feature Plots for Subset
#'
#' Creates a list of FeaturePlots for genes in a marker table on a subset object.
#'
#' @param subset_obj Seurat object (e.g., subclustered population)
#' @param marker_table Data frame with columns: gene, cell.types
#'
#' @return List of ggplot objects
#' @export
plot_markers_for_subset <- function(subset_obj, marker_table) {
  lapply(seq_len(nrow(marker_table)), function(i) {
    FeaturePlot(subset_obj, features = marker_table$gene[i]) +
      ggtitle(paste0(marker_table$cell.types[i], "\n", marker_table$gene[i])) +
      theme(plot.title = element_text(size = 8))
  })
}


#' Save Subcluster Composite Inspection Figure (multi-page PDF)
#'
#' Saves a multi-page PDF:
#'   Page 1 — all UMAP subcluster panels side by side
#'   Page 2+ — marker gene FeaturePlots for each cell type (one page per type)
#'
#' @param subcluster_list List of entries, each with:
#'   $umap_plot — DimPlot ggplot for that cell type.
#'   $obj       — Subclustered Seurat object for FeaturePlots.
#' @param marker_table  Data frame with columns: gene, cell.types.
#' @param output_dir    Directory where the PDF will be saved.
#' @param filename      Output filename (default "subclustering_inspection.pdf").
#' @param n_marker_cols Number of marker columns per page (default 4).
#' @export
save_subcluster_composite <- function(subcluster_list, marker_table, output_dir,
                                       filename      = "subclustering_inspection.pdf",
                                       n_marker_cols = 4) {
  n_marker_rows <- ceiling(nrow(marker_table) / n_marker_cols)
  path <- file.path(output_dir, filename)

  pdf(path, width = max(20, n_marker_cols * 4), height = 10 + n_marker_rows * 4)

  # For each cell type: UMAP page → markers page
  for (x in subcluster_list) {
    # Center the UMAP on the wide page so it doesn't stretch
    print(plot_spacer() | x$umap_plot | plot_spacer())
    markers <- plot_markers_for_subset(x$obj, marker_table)
    print(wrap_plots(markers, ncol = n_marker_cols))
  }

  dev.off()
  message("Saved → ", path)
}


#' Apply Subcluster Reassignment to Global Object
#'
#' Maps subcluster IDs from subclustered objects back to the global Seurat
#' object, updating cell-type labels based on a user-defined reassignment table.
#'
#' @param obj           Global Seurat object to update.
#' @param subcluster_list Named list of subclustered Seurat objects (names must
#'   match keys in reassign).
#' @param reassign      Named list of named character vectors. Each key is an
#'   object name in subcluster_list; each value maps subcluster IDs to cell-type
#'   labels. Use "others" as a catch-all for unlisted subcluster IDs.
#' @param source_col    Metadata column to copy as baseline before reassigning.
#' @param dest_col      Metadata column to write final labels into.
#'
#' @return Updated global Seurat object with dest_col populated.
#' @export
apply_subcluster_reassignment <- function(obj, subcluster_list, reassign,
                                           source_col = "celltype_grouped",
                                           dest_col   = "celltype_curated") {
  obj[[dest_col]] <- as.character(obj@meta.data[[source_col]])

  for (obj_name in names(reassign)) {
    sub    <- subcluster_list[[obj_name]]
    ids    <- as.character(sub$cluster_subtipo)
    labels <- reassign[[obj_name]][ids]

    if ("others" %in% names(reassign[[obj_name]])) {
      labels[is.na(labels)] <- reassign[[obj_name]]["others"]
    }

    obj@meta.data[colnames(sub), dest_col] <- unname(labels)
  }

  obj
}


# =============================================================================
# 6. PSEUDOBULK, DESEQ2, VOLCANO, HEATMAP
# =============================================================================

#' Assign Pseudo-replicates
#'
#' Randomly assigns cells within each condition to pseudo-replicate groups.
#' Conditions are auto-detected from the 'condition' column unless explicitly provided.
#'
#' @param obj         Seurat object with a 'condition' metadata column.
#' @param conditions Character vector of condition names to include. NULL
#'   (default) uses all conditions present in the data.
#' @param n_reps      Number of pseudo-replicates per condition.
#' @param seed        Random seed for reproducibility.
#' @return Seurat object with a replicate metadata column, or NULL if fewer
#'   than 2 conditions are present.
#' @export
assign_pseudo_replicates <- function(obj,
                                     conditions = NULL,
                                     n_reps      = 3,
                                     seed        = 1807) {

  set.seed(seed)

  # Auto-detect conditions from data if not provided
  all_conds <- unique(obj$condition)
  present_conditions <- if (!is.null(conditions)) intersect(all_conds, conditions) else all_conds

  if (length(present_conditions) < 2) return(NULL)

  obj$replicate <- NA
  for (cond in present_conditions) {
    idx              <- obj$condition == cond
    obj$replicate[idx] <- sample(paste0(cond, "_rep", 1:n_reps), sum(idx), replace = TRUE)
  }

  return(obj)
}


#' Create Pseudobulk Count Matrix from Seurat Object
#'
#' Aggregates counts by the replicate column using AggregateExpression.
#'
#' @param obj Seurat object with a replicate metadata column.
#' @return Data frame with genes as rows and replicate groups as columns.
#' @export
run_pseudobulk <- function(obj) {

  if (!"replicate" %in% colnames(obj@meta.data)) {
    stop("Object lacks 'replicate' metadata.")
  }

  keep <- !is.na(obj$replicate)
  if (!any(keep)) {
    stop("Object has no cells with assigned pseudo-replicates.")
  }

  obj <- subset(obj, cells = colnames(obj)[keep])

  assay_obj <- obj[["RNA"]]
  assay_layers <- Layers(assay_obj)
  count_layers <- grep("^counts", assay_layers, value = TRUE)

  if ("counts" %in% assay_layers) {
    counts <- GetAssayData(assay_obj, layer = "counts")
  } else if (length(count_layers) == 1) {
    counts <- GetAssayData(assay_obj, layer = count_layers)
  } else if (length(count_layers) > 1) {
    mats <- lapply(count_layers, function(x) GetAssayData(assay_obj, layer = x))
    counts <- Reduce(RowMergeSparseMatrices, mats)
  } else {
    stop("No counts layer found in RNA assay.")
  }

  rep_ids <- as.character(obj$replicate)

  if (ncol(counts) != length(rep_ids)) {
    stop("Counts matrix and replicate metadata have incompatible dimensions.")
  }

  rep_levels <- sort(unique(rep_ids))
  pseudobulk_mat <- vapply(rep_levels, function(rep_id) {
    idx <- rep_ids == rep_id
    Matrix::rowSums(counts[, idx, drop = FALSE])
  }, numeric(nrow(counts)))

  if (is.null(dim(pseudobulk_mat))) {
    pseudobulk_mat <- matrix(pseudobulk_mat, ncol = 1)
    colnames(pseudobulk_mat) <- rep_levels
  } else {
    colnames(pseudobulk_mat) <- rep_levels
  }

  rownames(pseudobulk_mat) <- rownames(counts)

  counts_df <- as.data.frame(pseudobulk_mat, check.names = FALSE)
  colnames(counts_df) <- sub("^g", "", colnames(counts_df))
  counts_df[, sort(colnames(counts_df)), drop = FALSE]
}


#' Save Pseudobulk Replicate Tables
#'
#' Builds one pseudobulk count table per Seurat object in a named list and
#' writes each table to disk as a CSV file. Intended for per-cell-type
#' pseudobulk objects carrying a `replicate` metadata column.
#'
#' @param obj_list    Named list of Seurat objects.
#' @param output_dir  Directory where CSV files will be written.
#' @param prefix      Filename prefix (default "Pseudobulk_Reps_").
#' @return Named list of pseudobulk count data frames.
#' @export
guardar_tablas_pseudobulk <- function(obj_list,
                                      output_dir,
                                      prefix = "Pseudobulk_Reps_",
                                      min_cells = 10,
                                      min_replicates = 2) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pseudobulk_list <- list()

  for (tipo in names(obj_list)) {
    obj <- obj_list[[tipo]]

    if (!"replicate" %in% colnames(obj@meta.data)) {
      warning("Skipping ", tipo, " because it lacks 'replicate' metadata.")
      next
    }

    obj <- subset(obj, cells = colnames(obj)[!is.na(obj$replicate)])
    if (ncol(obj) < min_cells) {
      warning("Skipping ", tipo, " because it has fewer than ", min_cells, " cells with pseudo-replicates.")
      next
    }

    rep_tab <- table(obj$replicate)
    if (length(rep_tab) < min_replicates) {
      warning("Skipping ", tipo, " because it has fewer than ", min_replicates, " pseudo-replicates.")
      next
    }

    counts_reps_df <- tryCatch(
      run_pseudobulk(obj),
      error = function(e) {
        warning("Skipping ", tipo, " due to pseudobulk error: ", e$message)
        NULL
      }
    )
    if (is.null(counts_reps_df)) next

    pseudobulk_list[[tipo]] <- counts_reps_df

    tipo_clean <- gsub("[^[:alnum:]_]", "_", tipo)
    file_name  <- paste0(prefix, tipo_clean, ".csv")

    write.csv(counts_reps_df,
              file = file.path(output_dir, file_name),
              row.names = TRUE)

    message("Saved pseudobulk table for ", tipo, " (", ncol(counts_reps_df), " replicates).")
  }

  pseudobulk_list
}


#' Run DESeq2 Differential Expression
#'
#' Builds a DESeqDataSet from a pseudobulk count matrix, detects condition
#' levels automatically from the column names, and writes per-comparison
#' result CSV files.
#'
#' @param counts_mat   Genes x samples count matrix (integer).
#' @param comparaciones List of lists, each with fields:
#'   \describe{
#'     \item{conds}{Character vector of length 2: c(reference, treatment).}
#'     \item{tag}{String label used for output file naming.}
#'   }
#' @param output_dir   Base directory; results are written to
#'   output_dir/tag/DESeq2_tag.csv.
#' @return Invisible NULL (side effect: writes CSV files).
#' @export
correr_deseq2 <- function(counts_mat, comparaciones, output_dir, tipo = NULL) {

  rep_names <- colnames(counts_mat)
  condition <- gsub("[_-]rep[0-9]+$", "", sub("^g", "", rep_names))

  if (length(unique(condition)) < 2) return(invisible(NULL))

  for (comp in comparaciones) {
    conds <- comp$conds
    tag   <- comp$tag
    keep   <- condition %in% conds
    if (!any(keep)) next

    counts_sub <- counts_mat[, keep, drop = FALSE]
    cond_sub   <- condition[keep]
    rep_tab    <- table(cond_sub)

    if (!all(conds %in% names(rep_tab))) {
      message("Skipping ", tag, " for ", ifelse(is.null(tipo), "global", tipo),
              ": missing condition(s).")
      next
    }

    if (any(rep_tab[conds] < 2) || ncol(counts_sub) <= length(conds)) {
      message("Skipping ", tag, " for ", ifelse(is.null(tipo), "global", tipo),
              ": insufficient pseudo-replicates per condition.")
      next
    }

    colData <- data.frame(
      row.names = colnames(counts_sub),
      condition = factor(cond_sub, levels = conds)
    )

    dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                                  colData   = colData,
                                  design    = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", conds[2], conds[1]))

    prefix <- if (!is.null(tipo)) paste0("DESeq2_", tipo, "_") else "DESeq2_"
    write.csv(as.data.frame(res),
              file = file.path(output_dir, tag, paste0(prefix, tag, ".csv")))
  }
}


#' Make Volcano Plot from DESeq2 Results CSV
#'
#' Reads a DESeq2 CSV output file and produces a volcano plot colored by
#' significance category.
#'
#' @param file       Path to the DESeq2 CSV file.
#' @param output_dir Output directory (currently unused; plot is returned).
#' @param padj_cut   Adjusted p-value cutoff.
#' @param lfc_cut    Log2 fold-change cutoff.
#' @return A ggplot object.
#' @export
plot_volcano <- function(file, padj_cut = 0.05, lfc_cut = 1) {

  nombre_base <- tools::file_path_sans_ext(basename(file))
  titulo      <- gsub("DESeq2_", "", nombre_base)

  df <- read.csv(file) %>%
    rownames_to_column("gene") %>%
    mutate(
      neg_log10_padj = -log10(padj),
      sig = case_when(
        padj <= padj_cut & log2FoldChange >=  lfc_cut ~ "Upregulated",
        padj <= padj_cut & log2FoldChange <= -lfc_cut ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    ) %>%
    filter(!is.na(log2FoldChange), is.finite(neg_log10_padj))

  ggplot(df, aes(log2FoldChange, neg_log10_padj, color = sig)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c(
      "Upregulated"     = "red",
      "Downregulated"   = "blue",
      "Not significant" = "gray"
    )) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut),      linetype = "dashed") +
    labs(title  = titulo,
         x      = "Log2 Fold Change",
         y      = "-Log10 adj p-value",
         color  = "Significance") +
    theme_minimal()
}


#' Process DESeq2 Result File
#'
#' Reads a DESeq2 CSV, classifies genes as up/down/unchanged, extracts log
#' fold-change values, and writes a filtered significant-gene CSV.
#'
#' @param file_path  Path to DESeq2 CSV file.
#' @param output_dir Directory for the filtered output CSV.
#' @param padj_cut   Adjusted p-value cutoff.
#' @param lfc_cut    Log2 fold-change cutoff.
#' @return Named list with elements class and logfc (data frames).
#' @export
procesar_deseq2_resultado <- function(file_path,
                                      output_dir,
                                      padj_cut = 0.05,
                                      lfc_cut  = 1) {

  df          <- read_csv(file_path, show_col_types = FALSE)
  comparacion <- gsub("^DESeq2_(.*)\\.csv$", "\\1", basename(file_path))

  # First column is always the gene ID (written as rownames by write.csv)
  gene_col <- colnames(df)[1]

  df_class <- df %>%
    mutate(
      gene_id       = .data[[gene_col]],
      clasificacion = case_when(
        padj <= padj_cut & log2FoldChange >  lfc_cut ~  1,
        padj <= padj_cut & log2FoldChange < -lfc_cut ~ -1,
        TRUE ~ 0
      )
    ) %>%
    dplyr::select(gene_id, clasificacion) %>%
    setNames(c("gene_id", comparacion))

  df_logfc <- df %>%
    mutate(
      gene_id = .data[[gene_col]],
      logfc   = ifelse(padj <= padj_cut & abs(log2FoldChange) > lfc_cut,
                       log2FoldChange, NA_real_)
    ) %>%
    dplyr::select(gene_id, logfc) %>%
    setNames(c("gene_id", comparacion))

  df_filt <- df %>% filter(padj <= padj_cut, abs(log2FoldChange) > lfc_cut)
  write_csv(df_filt, file.path(output_dir, paste0(comparacion, "_filtrado.csv")))

  list(class = df_class, logfc = df_logfc)
}


#' Render Volcano Plots for One DESeq2 Contrast Directory
#'
#' Reads all DESeq2 CSV files from a contrast-specific directory, writes one PNG
#' per result, and combines them into a single PDF.
#'
#' @param results_dir Directory containing DESeq2 CSV files for one contrast.
#' @param output_dir  Directory where volcano plot files will be written.
#' @param pdf_name    Name of the combined PDF file.
#' @param padj_cut    Adjusted p-value cutoff.
#' @param lfc_cut     Absolute log2 fold-change cutoff.
#' @return Invisible list of ggplot objects.
#' @export
render_volcano_plots <- function(results_dir,
                                 output_dir,
                                 pdf_name = "VolcanoPlots.pdf",
                                 padj_cut = 0.05,
                                 lfc_cut  = 1) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  csv_files <- list.files(results_dir, pattern = "^DESeq2_.*\\.csv$", full.names = TRUE)
  if (!length(csv_files)) {
    message("No DESeq2 CSV files found in: ", results_dir)
    return(invisible(list()))
  }

  pdf(file.path(output_dir, pdf_name), width = 12, height = 6)
  on.exit(dev.off(), add = TRUE)

  plots <- list()

  for (file in csv_files) {
    p <- plot_volcano(file, padj_cut = padj_cut, lfc_cut = lfc_cut) +
      labs(title = paste("Volcano Plot:", gsub("DESeq2_", "", tools::file_path_sans_ext(basename(file))))) +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(
      filename = file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file)), ".png")),
      plot     = p,
      width    = 8,
      height   = 6,
      dpi      = 300
    )

    plots <- c(plots, list(p))
    if (length(plots) == 2) {
      grid.arrange(grobs = plots, ncol = 2)
      plots <- list()
    }
  }

  if (length(plots) == 1) {
    grid.arrange(grobs = plots, ncol = 1)
  }

  invisible(plots)
}


#' Build Combined Differential Expression Tables for One Contrast
#'
#' Processes all DESeq2 CSV files in a contrast directory, writes filtered
#' per-cell-type CSV files, and assembles combined classification and log2FC
#' matrices across cell types.
#'
#' @param results_dir Directory containing DESeq2 CSV files for one contrast.
#' @param output_dir  Directory where combined tables will be written.
#' @param padj_cut    Adjusted p-value cutoff.
#' @param lfc_cut     Absolute log2 fold-change cutoff.
#' @param prefix      Filename prefix for combined output tables.
#' @return Named list with `class` and `logfc` tibbles.
#' @export
build_differential_tables <- function(results_dir,
                                      output_dir,
                                      padj_cut = 0.05,
                                      lfc_cut  = 1,
                                      prefix   = "tabla_diferenciales") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(results_dir, pattern = "^DESeq2_.*\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No DESeq2 CSV files found in: ", results_dir)

  listas <- lapply(files,
                   procesar_deseq2_resultado,
                   output_dir = output_dir,
                   padj_cut   = padj_cut,
                   lfc_cut    = lfc_cut)

  tabla_class <- Reduce(function(x, y) full_join(x, y, by = "gene_id"),
                        lapply(listas, `[[`, "class")) %>%
    arrange(gene_id)

  tabla_logfc <- Reduce(function(x, y) full_join(x, y, by = "gene_id"),
                        lapply(listas, `[[`, "logfc")) %>%
    arrange(gene_id)

  tabla_filtrada <- tabla_class %>%
    filter(apply(dplyr::select(., -gene_id) != 0, 1, any))

  tabla_logfc_filtrada <- tabla_logfc %>%
    filter(gene_id %in% tabla_filtrada$gene_id)

  suffix <- paste0("fc", lfc_cut, "_padj_", gsub("\\.", "", as.character(padj_cut)))

  write_tsv(tabla_filtrada,
            file.path(output_dir, paste0(prefix, "_", suffix, ".tsv")))
  write_tsv(tabla_logfc_filtrada,
            file.path(output_dir, paste0("tabla_log2FC_", suffix, ".tsv")))

  list(class = tabla_filtrada, logfc = tabla_logfc_filtrada)
}


#' Hierarchically Clustered Heatmap of DE Results
#'
#' Clusters rows (genes) by Euclidean distance + dynamic tree cut, clusters
#' columns (conditions) by PCA-based distance, and renders a pheatmap.
#'
#' @param matriz       Numeric matrix (genes x conditions), e.g. log2FC values.
#' @param min_genes    Minimum cluster size for dynamic tree cut.
#' @param deepSplit_val deepSplit parameter for cutreeDynamic.
#' @param breaks       Two-element vector c(min, max) for the color scale.
#' @export
plot_heatmap <- function(matriz,
                          min_genes    = 1,
                          deepSplit_val = 0,
                          breaks       = c(-5, 5)) {

  dist_rows <- dist(matriz, method = "euclidean")
  hc_rows   <- hclust(dist_rows, method = "complete")

  clust <- cutreeDynamic(
    dendro            = hc_rows,
    distM             = as.matrix(dist_rows),
    deepSplit         = deepSplit_val,
    minClusterSize    = min_genes,
    pamRespectsDendro = FALSE
  )

  pca_res <- prcomp(t(matriz), scale. = FALSE)
  var_exp <- summary(pca_res)$importance[3, ]
  n_pcs   <- which(var_exp >= 0.90)[1]
  hc_cols <- hclust(dist(pca_res$x[, 1:n_pcs]), method = "complete")

  paleta         <- colorRampPalette(brewer.pal(12, "Dark2"))(length(unique(clust[clust > 0])))
  annotation_row <- data.frame(Cluster = as.factor(clust))
  rownames(annotation_row) <- rownames(matriz)

  breaks_seq  <- seq(breaks[1], breaks[2], length.out = 80)
  color_scale <- colorRampPalette(c("blue", "black", "yellow"))(length(breaks_seq) - 1)

  pheatmap(matriz,
           cluster_rows    = hc_rows,
           cluster_cols    = hc_cols,
           annotation_row  = annotation_row,
           annotation_colors = list(
             Cluster = setNames(paleta, sort(unique(clust[clust > 0])))
           ),
           color           = color_scale,
           breaks          = breaks_seq,
           show_rownames   = TRUE,
           border_color    = NA,
           fontsize_row    = 1,
           fontsize_col    = 20,
           fontsize        = 22,
           main            = sprintf("Heatmap (%d genes)", nrow(matriz)))
}


#' Marker Gene DotPlot with Custom Cell-Type Order
#'
#' Builds a DotPlot where cell types (Y-axis) and marker genes (X-axis, coord-
#' flipped) follow user-defined orders, producing a near-diagonal expression
#' pattern useful for cell-type validation figures.
#'
#' @param seurat_obj        Seurat object with annotations in `annot_col`.
#' @param marks             Data frame with columns `gene` and `cell.types`.
#' @param annot_col         Metadata column holding cell-type labels.
#' @param cell_order        Character vector: desired order of cell types
#'                          (top to bottom). Types not listed appear at the end.
#' @param clusters_remove   Cell-type labels to exclude (default NULL).
#' @param rename_map        Named character vector for renaming cell types before
#'                          plotting, e.g. c("Meristemoid" = "Stomatal lineage").
#' @param outfile           PDF output path (NULL = no save).
#' @param width             PDF width in inches.
#' @param height            PDF height in inches.
#' @param dot_scale         Dot size scaling factor.
#' @param base_size         Base font size.
#' @return A ggplot object.
#' @export
plot_marker_dotplot <- function(seurat_obj,
                                     marks,
                                     annot_col       = "celltype_grouped",
                                     cell_order      = NULL,
                                     clusters_remove = NULL,
                                     rename_map      = NULL,
                                     outfile         = NULL,
                                     width           = 20,
                                     height          = 10,
                                     dot_scale       = 12,
                                     base_size       = 18) {

  obj <- seurat_obj

  # ── Rename cell types if requested ───────────────────────────────────────────
  if (!is.null(rename_map)) {
    for (old_name in names(rename_map)) {
      obj@meta.data[[annot_col]][obj@meta.data[[annot_col]] == old_name] <- rename_map[[old_name]]
    }
    if (!is.null(marks) && "cell.types" %in% colnames(marks)) {
      for (old_name in names(rename_map)) {
        marks$cell.types[marks$cell.types == old_name] <- rename_map[[old_name]]
      }
    }
  }

  # ── Remove unwanted clusters ──────────────────────────────────────────────────
  if (!is.null(clusters_remove)) {
    obj <- subset(obj, subset = !!sym(annot_col) %in% clusters_remove, invert = TRUE)
  }

  # ── Build ordered factor ──────────────────────────────────────────────────────
  all_types <- unique(obj@meta.data[[annot_col]])
  if (is.null(cell_order)) {
    ordered_levels <- all_types
  } else {
    remaining      <- setdiff(all_types, cell_order)
    ordered_levels <- c(cell_order, remaining)
  }

  obj@meta.data[["annotation_orden"]] <- factor(
    obj@meta.data[[annot_col]],
    levels = ordered_levels
  )
  Idents(obj) <- "annotation_orden"

  # ── Filter genes present in the object ───────────────────────────────────────
  genes_use <- unique(intersect(marks$gene, rownames(obj)))
  if (length(genes_use) == 0) stop("No marker genes found in the Seurat object.")

  # ── Build DotPlot ─────────────────────────────────────────────────────────────
  figure <- DotPlot(
    obj,
    features = genes_use,
    dot.scale = dot_scale,
    cols      = c("yellow", "darkblue")
  ) +
    coord_flip() +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = base_size - 4),
      axis.text.y  = element_text(size = base_size - 4, face = "italic"),
      axis.title   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "right"
    )

  figure$layers[[1]]$aes_params$alpha <- 1   # solid dots

  # ── Save ──────────────────────────────────────────────────────────────────────
  if (!is.null(outfile)) {
    if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)
    ggsave(outfile, figure, width = width, height = height, dpi = 500)
    message("DotPlot saved: ", outfile)
  }

  invisible(figure)
}


# =============================================================================
# 7. GO ENRICHMENT
# =============================================================================

#' Run GO Enrichment Analysis
#'
#' Iterates over columns of a binary classification matrix and runs clusterProfiler
#' enrichGO for each set of upregulated genes. Writes raw and gene-symbol-readable
#' result tables, optionally simplifying redundant terms.
#'
#' Common OrgDb / keytype combinations:
#'   - Arabidopsis thaliana : OrgDb = org.At.tair.db, keytype = "TAIR"
#'   - Homo sapiens         : OrgDb = org.Hs.eg.db,   keytype = "ENSEMBL" or "ENTREZID"
#'   - Mus musculus         : OrgDb = org.Mm.eg.db,   keytype = "ENSEMBL" or "ENTREZID"
#'   - Oryza sativa         : OrgDb = org.Os.eg.db,   keytype = "GID"
#'
#' @param tabla          Binary matrix (genes x comparisons); genes with value 1
#'   are tested for enrichment.
#' @param universo       Character vector of background gene IDs.
#' @param espacio        GO namespace: "BP", "MF", or "CC".
#' @param orgdb          OrgDb annotation object (default org.At.tair.db).
#' @param keytype        Key type matching rownames of tabla (default "TAIR").
#' @param qvalueCutoff   Q-value cutoff for enrichment (default 0.05).
#' @param pvalueCutoff   P-value cutoff for enrichment (default 0.05).
#' @param simplificar    If TRUE, simplify redundant GO terms before saving.
#' @param umbral_simply  Similarity cutoff for simplify() (default 0.7).
#' @param output_dir     Directory for output text files.
#' @return Named list of enrichResult objects (one per column of tabla).
#' @export
correr_enriquecimiento_go <- function(tabla,
                                       universo,
                                       espacio,
                                       orgdb          = org.At.tair.db,
                                       keytype        = "TAIR",
                                       qvalueCutoff   = 0.05,
                                       pvalueCutoff   = 0.05,
                                       simplificar    = FALSE,
                                       umbral_simply  = 0.7,
                                       output_dir     = "results/Enrichment") {

  salida        <- vector("list", ncol(tabla))
  names(salida) <- colnames(tabla)

  for (n in seq_len(ncol(tabla))) {

    gene <- unique(trimws(gsub("\\..*", "", rownames(tabla)[tabla[, n] == 1])))
    if (length(gene) == 0) {
      message("Sin genes: ", colnames(tabla)[n])
      next
    }

    enri <- tryCatch(
      enrichGO(gene          = gene,
               universe      = universo,
               OrgDb         = orgdb,
               keyType       = keytype,
               ont           = espacio,
               pAdjustMethod = "BH",
               pvalueCutoff  = pvalueCutoff,
               qvalueCutoff  = qvalueCutoff,
               readable      = FALSE),
      error = function(e) NULL
    )

    if (is.null(enri) || nrow(enri@result) == 0) {
      message("Sin GO: ", colnames(tabla)[n])
      next
    }

    # Save raw and gene-symbol-readable results
    sufijo <- paste(colnames(tabla)[n], espacio, qvalueCutoff, sep = ".")

    write.table(as.data.frame(enri),
                file.path(output_dir, paste0(sufijo, ".txt")),
                sep = "\t", col.names = NA, quote = FALSE)

    write.table(as.data.frame(setReadable(enri, OrgDb = orgdb)),
                file.path(output_dir, paste0(sufijo, ".symbol.txt")),
                sep = "\t", col.names = NA, quote = FALSE)

    if (simplificar) {
      enri_s <- simplify(enri, cutoff = umbral_simply, by = "p.adjust", select_fun = min)
      if (!is.null(enri_s) && nrow(enri_s@result) > 0) {
        write.table(as.data.frame(enri_s),
                    file.path(output_dir, paste0(sufijo, ".simply.", umbral_simply, ".txt")),
                    sep = "\t", col.names = NA, quote = FALSE)
        salida[[n]] <- enri_s
      }
    } else {
      salida[[n]] <- enri
    }
  }

  return(salida)
}


#' Filter GO Results by Ontology Level
#'
#' Applies gofilter to each enrichResult in a list, keeping only terms at or
#' below the specified GO level, and writes filtered tables to disk.
#'
#' @param resuGO    Named list of enrichResult objects.
#' @param nivel     Maximum GO level to retain.
#' @param espacio   GO namespace string (used in output filenames).
#' @param qvalueCutoff Q-value cutoff (used in output filenames).
#' @param simplificar Logical; affects output filename suffix.
#' @param output_dir Directory for output files.
#' @return Named list of filtered enrichResult objects.
#' @export
podar_go <- function(resuGO,
                     nivel,
                     espacio,
                     qvalueCutoff,
                     simplificar  = FALSE,
                     output_dir   = "results/Enrichment") {

  salida        <- vector("list", length(resuGO))
  names(salida) <- names(resuGO)

  for (k in seq_along(resuGO)) {

    if (is.null(resuGO[[k]])) next

    res <- tryCatch(gofilter(resuGO[[k]], nivel), error = function(e) NULL)
    if (is.null(res) || nrow(res@result) == 0) next

    salida[[k]] <- res

    sufijo <- paste(
      names(resuGO)[k], espacio, qvalueCutoff,
      if (simplificar) "simply" else "total",
      paste0("nivel_", nivel), "txt",
      sep = "."
    )
    write.table(as.data.frame(res),
                file.path(output_dir, sufijo),
                sep = "\t", col.names = NA, quote = FALSE)
  }

  return(salida)
}


#' Balloon Plot of GO Enrichment Results
#'
#' Visualizes enrichment results as a balloon/bubble chart where bubble size
#' encodes fold enrichment and fill color encodes -log10(q-value).
#'
#' @param resuGO Named list of enrichResult objects (one per comparison).
#' @return A ggplot object.
#' @export
graficar_go_balones <- function(resuGO) {

  nombres <- names(resuGO)
  if (is.null(nombres)) nombres <- as.character(seq_along(resuGO))

  bloques <- lapply(seq_along(resuGO), function(k) {
    if (is.null(resuGO[[k]])) return(NULL)
    df <- as.data.frame(resuGO[[k]])
    if (!nrow(df)) return(NULL)

    gr <- as.numeric(unlist(strsplit(df$GeneRatio, "/")))
    br <- as.numeric(unlist(strsplit(df$BgRatio,   "/")))
    gr <- gr[seq(1, length(gr), 2)] / gr[seq(2, length(gr), 2)]
    br <- br[seq(1, length(br), 2)] / br[seq(2, length(br), 2)]

    data.frame(
      Exp          = nombres[k],
      GOid         = df$ID,
      GODesc       = df$Description,
      Log10Qvalue  = -log10(df$qvalue),
      Enrichment   = gr / br
    )
  })

  bloques <- Filter(Negate(is.null), bloques)
  if (!length(bloques)) stop("No GO enrichment results available to plot.")

  dat     <- na.omit(do.call(rbind, bloques))
  if (!nrow(dat)) stop("No GO enrichment results available to plot.")
  dat$Exp <- factor(dat$Exp, levels = nombres)

  ggballoonplot(dat, x = "Exp", y = "GODesc",
                size = "Enrichment", fill = "Log10Qvalue") +
    scale_fill_gradientn(colors = brewer.pal(8, "YlOrRd")) +
    guides(size = "none") +
    theme_minimal(base_size = 11) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 28)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title  = element_blank())
}


#' Run GO Enrichment Suite from a Differential Table
#'
#' Takes a binary differential-expression table, derives the background gene
#' universe automatically from the supplied OrgDb/keytype pair, runs full and
#' simplified enrichGO analyses, optionally prunes by GO level, and exports a
#' multi-page balloon-plot PDF.
#'
#' @param diff_table     Either a data.frame/tibble or a path to a TSV file. The
#'   first column must contain gene IDs; remaining columns must be binary 0/1.
#' @param output_dir     Directory where GO result tables and plots will be saved.
#' @param orgdb          OrgDb annotation object.
#' @param keytype        Key type matching the gene IDs in `diff_table`.
#' @param espacio        GO namespace: "BP", "MF", or "CC".
#' @param qvalue_cutoff  Q-value cutoff for enrichGO.
#' @param pvalue_cutoff  P-value cutoff for enrichGO.
#' @param simplify_cutoff Similarity cutoff for simplify().
#' @param go_level       GO level to retain in pruned outputs.
#' @param pdf_name       Output PDF filename for balloon plots.
#' @return Named list with total/simple and pruned GO result lists.
#' @export
run_go_enrichment_suite <- function(diff_table,
                                    output_dir,
                                    orgdb,
                                    keytype,
                                    espacio         = "BP",
                                    qvalue_cutoff   = 0.05,
                                    pvalue_cutoff   = 0.05,
                                    simplify_cutoff = 0.7,
                                    go_level        = 6,
                                    pdf_name        = "GO_enrichment.pdf") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.character(diff_table) && length(diff_table) == 1) {
    tabla_df <- read.table(diff_table, header = TRUE, sep = "\t", check.names = FALSE)
  } else {
    tabla_df <- as.data.frame(diff_table, check.names = FALSE)
  }

  if (ncol(tabla_df) < 2) {
    stop("Differential table must contain one gene-ID column plus at least one comparison column.")
  }

  gene_col <- colnames(tabla_df)[1]
  rownames(tabla_df) <- tabla_df[[gene_col]]
  tabla <- as.matrix(tabla_df[, -1, drop = FALSE])
  mode(tabla) <- "numeric"
  tabla <- tabla[, colSums(tabla != 0, na.rm = TRUE) > 0, drop = FALSE]

  if (!ncol(tabla)) {
    stop("Differential table contains no non-zero comparison columns.")
  }

  universo <- keys(orgdb, keytype = keytype)

  go_total <- correr_enriquecimiento_go(
    tabla            = tabla,
    universo         = universo,
    espacio          = espacio,
    orgdb            = orgdb,
    keytype          = keytype,
    qvalueCutoff     = qvalue_cutoff,
    pvalueCutoff     = pvalue_cutoff,
    simplificar      = FALSE,
    umbral_simply    = simplify_cutoff,
    output_dir       = output_dir
  )

  go_simple <- correr_enriquecimiento_go(
    tabla            = tabla,
    universo         = universo,
    espacio          = espacio,
    orgdb            = orgdb,
    keytype          = keytype,
    qvalueCutoff     = qvalue_cutoff,
    pvalueCutoff     = pvalue_cutoff,
    simplificar      = TRUE,
    umbral_simply    = simplify_cutoff,
    output_dir       = output_dir
  )

  go_total_podado <- podar_go(
    resuGO         = go_total,
    nivel          = go_level,
    espacio        = espacio,
    qvalueCutoff   = qvalue_cutoff,
    simplificar    = FALSE,
    output_dir     = output_dir
  )

  go_simple_podado <- podar_go(
    resuGO         = go_simple,
    nivel          = go_level,
    espacio        = espacio,
    qvalueCutoff   = qvalue_cutoff,
    simplificar    = TRUE,
    output_dir     = output_dir
  )

  pdf(file.path(output_dir, pdf_name), width = 18, height = 18)
  on.exit(dev.off(), add = TRUE)

  for (obj in list(go_total, go_simple, go_total_podado, go_simple_podado)) {
    tryCatch(print(graficar_go_balones(obj)),
             error = function(e) message("Skipping GO balloon plot: ", e$message))
  }

  invisible(list(
    tabla            = tabla,
    universo         = universo,
    total            = go_total,
    simple           = go_simple,
    total_podado     = go_total_podado,
    simple_podado    = go_simple_podado
  ))
}


#' Run Heatmap, Rank-Based Coexpression, TOM and GO per Cluster
#'
#' Builds a log2FC heatmap from a combined differential table, derives gene
#' clusters from dynamic tree cutting, computes a rank-based coexpression
#' network with TOM, and runs GO enrichment for each resulting cluster/module.
#'
#' @param diff_table       Path to a log2FC TSV file or a data.frame/tibble. The
#'   first column must contain gene IDs.
#' @param output_dir       Output directory for plots and tables.
#' @param selected_cols    Optional character vector of column names to retain.
#'   If NULL, all non-gene columns are used.
#' @param min_genes        Minimum cluster/module size for dynamic tree cut.
#' @param deepSplit_val    `cutreeDynamic` deepSplit parameter.
#' @param breaks           Two-element numeric vector controlling heatmap scale.
#' @param network_power    Soft-threshold power applied to correlation similarity.
#' @param network_type     Network type: "signed" or "unsigned".
#' @param cor_method       Correlation method; use "spearman" for rank-based
#'   coexpression.
#' @param go_orgdb         OrgDb object matching the organism.
#' @param go_keytype       Key type matching the gene IDs in `diff_table`.
#' @param go_space         GO namespace: "BP", "MF", or "CC".
#' @param go_qvalue_cutoff Q-value cutoff for enrichGO.
#' @param go_pvalue_cutoff P-value cutoff for enrichGO.
#' @param go_simplify_cutoff Similarity cutoff for simplify().
#' @param go_level         GO level used for pruning.
#' @param heatmap_pdf      Output PDF filename for the heatmap.
#' @param tom_pdf          Output PDF filename for the TOM heatmap.
#' @param go_pdf           Output PDF filename for GO balloon plots.
#' @return Named list with matrix, cluster assignments, TOM, and GO results.
#' @export
run_coexpression_cluster_suite <- function(diff_table,
                                           output_dir,
                                           selected_cols      = NULL,
                                           min_genes          = 1,
                                           deepSplit_val      = 0,
                                           breaks             = c(-5, 5),
                                           network_power      = 6,
                                           network_type       = c("signed", "unsigned"),
                                           cor_method         = "spearman",
                                           go_orgdb,
                                           go_keytype,
                                           go_space           = "BP",
                                           go_qvalue_cutoff   = 0.05,
                                           go_pvalue_cutoff   = 0.05,
                                           go_simplify_cutoff = 0.7,
                                           go_level           = 6,
                                           heatmap_pdf        = "coexpression_heatmap.pdf",
                                           tom_pdf            = "tom_heatmap.pdf",
                                           go_pdf             = "GO_clusters.pdf") {

  network_type <- match.arg(network_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.character(diff_table) && length(diff_table) == 1) {
    df <- read.table(diff_table, header = TRUE, sep = "\t", check.names = FALSE)
  } else {
    df <- as.data.frame(diff_table, check.names = FALSE)
  }

  if (ncol(df) < 2) {
    stop("diff_table must contain one gene-ID column plus at least one numeric column.")
  }

  gene_col <- colnames(df)[1]
  rownames(df) <- df[[gene_col]]

  if (is.null(selected_cols)) {
    selected_cols <- colnames(df)[-1]
  } else {
    selected_cols <- intersect(selected_cols, colnames(df))
  }

  if (!length(selected_cols)) {
    stop("No selected columns found in diff_table.")
  }

  Mz <- as.matrix(df[, selected_cols, drop = FALSE])
  mode(Mz) <- "numeric"
  Mz[is.na(Mz)] <- 0

  if (nrow(Mz) < 2 || ncol(Mz) < 2) {
    stop("Need at least two genes and two selected columns for coexpression analysis.")
  }

  dist_rows <- dist(Mz, method = "euclidean")
  hc_rows   <- hclust(dist_rows, method = "complete")

  clust <- cutreeDynamic(
    dendro            = hc_rows,
    distM             = as.matrix(dist_rows),
    deepSplit         = deepSplit_val,
    minClusterSize    = min_genes,
    pamRespectsDendro = FALSE
  )

  pca_res <- prcomp(t(Mz), scale. = FALSE)
  var_exp <- summary(pca_res)$importance[3, ]
  n_pcs   <- which(var_exp >= 0.90)[1]
  if (is.na(n_pcs)) n_pcs <- ncol(pca_res$x)
  hc_cols <- hclust(dist(pca_res$x[, seq_len(n_pcs), drop = FALSE]), method = "complete")

  clusters_unicos <- sort(unique(clust[clust > 0]))
  paleta <- colorRampPalette(brewer.pal(12, "Dark2"))(max(1, length(clusters_unicos)))
  names(paleta) <- if (length(clusters_unicos)) clusters_unicos else "0"

  annotation_rows <- data.frame(Cluster = as.factor(clust))
  rownames(annotation_rows) <- rownames(Mz)

  breaks_seq  <- seq(breaks[1], breaks[2], length.out = 80)
  color_scale <- colorRampPalette(c("blue", "black", "yellow"))(length(breaks_seq) - 1)

  pdf(file.path(output_dir, heatmap_pdf), width = 10, height = 20)
  heatmap_obj <- pheatmap(
    Mz,
    cluster_rows      = hc_rows,
    cluster_cols      = hc_cols,
    annotation_row    = annotation_rows,
    annotation_colors = list(Cluster = paleta),
    color             = color_scale,
    breaks            = breaks_seq,
    show_rownames     = TRUE,
    border_color      = NA,
    fontsize_row      = 1,
    fontsize_col      = 20,
    fontsize          = 22,
    use_raster        = FALSE,
    treeheight_row    = 50,
    treeheight_col    = 50,
    main              = sprintf("Heatmap (%d genes) - min %d genes por cluster",
                                nrow(Mz), min_genes)
  )
  dev.off()

  datExpr <- t(Mz)
  gene_cor <- suppressWarnings(cor(datExpr, method = cor_method, use = "pairwise.complete.obs"))
  gene_cor[is.na(gene_cor)] <- 0
  diag(gene_cor) <- 1

  adjacency_mat <- if (network_type == "signed") {
    ((1 + gene_cor) / 2) ^ network_power
  } else {
    abs(gene_cor) ^ network_power
  }
  diag(adjacency_mat) <- 1

  TOM <- TOMsimilarity(adjacency_mat, TOMType = network_type)
  rownames(TOM) <- rownames(Mz)
  colnames(TOM) <- rownames(Mz)

  dissTOM <- 1 - TOM
  gene_tree <- hclust(as.dist(dissTOM), method = "average")
  tom_clusters <- cutreeDynamic(
    dendro            = gene_tree,
    distM             = dissTOM,
    deepSplit         = deepSplit_val,
    minClusterSize    = min_genes,
    pamRespectsDendro = FALSE
  )

  tom_annotation <- data.frame(Module = as.factor(tom_clusters))
  rownames(tom_annotation) <- rownames(Mz)
  modules_unicos <- sort(unique(tom_clusters[tom_clusters > 0]))
  tom_paleta <- colorRampPalette(brewer.pal(12, "Set3"))(max(1, length(modules_unicos)))
  names(tom_paleta) <- if (length(modules_unicos)) modules_unicos else "0"

  tom_order <- gene_tree$order
  tom_plot_mat <- TOM[tom_order, tom_order, drop = FALSE]
  tom_plot_ann <- tom_annotation[tom_order, , drop = FALSE]
  pdf(file.path(output_dir, tom_pdf), width = 10, height = 10)
  pheatmap(
    tom_plot_mat,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    annotation_row    = tom_plot_ann,
    annotation_col    = tom_plot_ann,
    annotation_colors = list(Module = tom_paleta),
    color             = colorRampPalette(c("white", "steelblue", "navy"))(80),
    border_color      = NA,
    main              = sprintf("TOM Heatmap (%d genes)", nrow(Mz))
  )
  dev.off()

  cluster_assignments <- data.frame(
    gene_id         = rownames(Mz),
    heatmap_cluster = clust,
    tom_module      = tom_clusters,
    stringsAsFactors = FALSE
  )
  write.table(cluster_assignments,
              file.path(output_dir, "gene_cluster_assignments.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cluster_table <- data.frame(gene_id = rownames(Mz), stringsAsFactors = FALSE)
  for (cluster_id in sort(unique(tom_clusters[tom_clusters > 0]))) {
    cluster_name <- paste0("cluster_", cluster_id)
    cluster_table[[cluster_name]] <- as.integer(tom_clusters == cluster_id)
  }

  if (ncol(cluster_table) < 2) {
    stop("No positive TOM modules were detected for GO enrichment.")
  }

  go_results <- run_go_enrichment_suite(
    diff_table      = cluster_table,
    output_dir      = file.path(output_dir, "GO_clusters"),
    orgdb           = go_orgdb,
    keytype         = go_keytype,
    espacio         = go_space,
    qvalue_cutoff   = go_qvalue_cutoff,
    pvalue_cutoff   = go_pvalue_cutoff,
    simplify_cutoff = go_simplify_cutoff,
    go_level        = go_level,
    pdf_name        = go_pdf
  )

  invisible(list(
    matrix              = Mz,
    heatmap_clusters    = clust,
    tom_modules         = tom_clusters,
    cluster_assignments = cluster_assignments,
    tom                 = TOM,
    go_results          = go_results
  ))
}


#' Prepare Log2FC Matrix for Coexpression Analysis
#'
#' Loads a combined log2FC table, selects columns of interest, and returns a
#' numeric matrix with gene IDs as row names.
#'
#' @param diff_table    Path to a log2FC TSV file or a data.frame/tibble. The
#'   first column must contain gene IDs.
#' @param selected_cols Optional character vector of columns to retain. If NULL,
#'   all non-gene columns are used.
#' @return Numeric matrix with genes as rows and selected contrasts/cell types as columns.
#' @export
prepare_coexpression_matrix <- function(diff_table, selected_cols = NULL) {

  if (is.character(diff_table) && length(diff_table) == 1) {
    df <- read.table(diff_table, header = TRUE, sep = "\t", check.names = FALSE)
  } else {
    df <- as.data.frame(diff_table, check.names = FALSE)
  }

  if (ncol(df) < 2) {
    stop("diff_table must contain one gene-ID column plus at least one numeric column.")
  }

  gene_col <- colnames(df)[1]
  rownames(df) <- df[[gene_col]]

  if (is.null(selected_cols)) {
    selected_cols <- colnames(df)[-1]
  } else {
    selected_cols <- intersect(selected_cols, colnames(df))
  }

  if (!length(selected_cols)) {
    stop("No selected columns found in diff_table.")
  }

  Mz <- as.matrix(df[, selected_cols, drop = FALSE])
  mode(Mz) <- "numeric"
  Mz[is.na(Mz)] <- 0

  if (nrow(Mz) < 2 || ncol(Mz) < 2) {
    stop("Need at least two genes and two selected columns for coexpression analysis.")
  }

  Mz
}


#' Build Differential Gene Heatmap and Dynamic Clusters
#'
#' Generates a clustered heatmap from a log2FC matrix and returns the heatmap
#' gene cluster assignments.
#'
#' @param Mz            Numeric matrix with genes as rows.
#' @param output_dir    Output directory.
#' @param min_genes     Minimum cluster size for dynamic tree cut.
#' @param deepSplit_val `cutreeDynamic` deepSplit parameter.
#' @param breaks        Two-element numeric vector controlling heatmap scale.
#' @param heatmap_pdf   Output PDF filename.
#' @return Named list with matrix, dendrograms, cluster assignments, and pheatmap object.
#' @export
build_heatmap_clusters <- function(Mz,
                                   output_dir,
                                   min_genes     = 1,
                                   deepSplit_val = 0,
                                   breaks        = c(-5, 5),
                                   heatmap_pdf   = "coexpression_heatmap.pdf") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  dist_rows <- dist(Mz, method = "euclidean")
  hc_rows   <- hclust(dist_rows, method = "complete")

  clust <- cutreeDynamic(
    dendro            = hc_rows,
    distM             = as.matrix(dist_rows),
    deepSplit         = deepSplit_val,
    minClusterSize    = min_genes,
    pamRespectsDendro = FALSE
  )

  pca_res <- prcomp(t(Mz), scale. = FALSE)
  var_exp <- summary(pca_res)$importance[3, ]
  n_pcs   <- which(var_exp >= 0.90)[1]
  if (is.na(n_pcs)) n_pcs <- ncol(pca_res$x)
  hc_cols <- hclust(dist(pca_res$x[, seq_len(n_pcs), drop = FALSE]), method = "complete")

  clusters_unicos <- sort(unique(clust[clust > 0]))
  paleta <- colorRampPalette(brewer.pal(12, "Dark2"))(max(1, length(clusters_unicos)))
  names(paleta) <- if (length(clusters_unicos)) clusters_unicos else "0"

  annotation_rows <- data.frame(Cluster = as.factor(clust))
  rownames(annotation_rows) <- rownames(Mz)

  breaks_seq  <- seq(breaks[1], breaks[2], length.out = 80)
  color_scale <- colorRampPalette(c("blue", "black", "yellow"))(length(breaks_seq) - 1)

  pdf(file.path(output_dir, heatmap_pdf), width = 10, height = 20)
  pheatmap(
    Mz,
    cluster_rows      = hc_rows,
    cluster_cols      = hc_cols,
    annotation_row    = annotation_rows,
    annotation_colors = list(Cluster = paleta),
    color             = color_scale,
    breaks            = breaks_seq,
    show_rownames     = TRUE,
    border_color      = NA,
    fontsize_row      = 1,
    fontsize_col      = 20,
    fontsize          = 22,
    use_raster        = FALSE,
    treeheight_row    = 50,
    treeheight_col    = 50,
    main              = sprintf("Heatmap (%d genes) - min %d genes por cluster",
                                nrow(Mz), min_genes)
  )
  dev.off()

  cluster_assignments <- data.frame(
    gene_id         = rownames(Mz),
    heatmap_cluster = clust,
    stringsAsFactors = FALSE
  )
  write.table(cluster_assignments,
              file.path(output_dir, "heatmap_gene_clusters.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(
    matrix               = Mz,
    row_tree             = hc_rows,
    col_tree             = hc_cols,
    heatmap_clusters     = clust,
    cluster_assignments  = cluster_assignments
  ))
}


#' Build Rank-Based Coexpression Network and TOM Modules
#'
#' Computes a rank-based gene-gene correlation matrix, adjacency, TOM, and TOM
#' modules from a log2FC matrix.
#'
#' @param Mz             Numeric matrix with genes as rows.
#' @param output_dir     Output directory.
#' @param min_genes      Minimum module size for dynamic tree cut.
#' @param deepSplit_val  `cutreeDynamic` deepSplit parameter.
#' @param network_power  Soft-threshold power applied to similarity.
#' @param network_type   Network type: "signed" or "unsigned".
#' @param cor_method     Correlation method, e.g. "spearman".
#' @param tom_pdf        Output PDF filename.
#' @return Named list with correlation, TOM, gene tree, module assignments, and pheatmap object.
#' @export
build_coexpression_modules <- function(Mz,
                                       output_dir,
                                       min_genes     = 1,
                                       deepSplit_val = 0,
                                       network_power = 6,
                                       network_type  = c("signed", "unsigned"),
                                       cor_method    = "spearman",
                                       tom_pdf       = "tom_heatmap.pdf") {

  network_type <- match.arg(network_type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  datExpr <- t(Mz)
  gene_cor <- suppressWarnings(cor(datExpr, method = cor_method, use = "pairwise.complete.obs"))
  gene_cor[is.na(gene_cor)] <- 0
  diag(gene_cor) <- 1

  adjacency_mat <- if (network_type == "signed") {
    ((1 + gene_cor) / 2) ^ network_power
  } else {
    abs(gene_cor) ^ network_power
  }
  diag(adjacency_mat) <- 1

  TOM <- TOMsimilarity(adjacency_mat, TOMType = network_type)
  rownames(TOM) <- rownames(Mz)
  colnames(TOM) <- rownames(Mz)

  dissTOM <- 1 - TOM
  gene_tree <- hclust(as.dist(dissTOM), method = "average")
  tom_clusters <- cutreeDynamic(
    dendro            = gene_tree,
    distM             = dissTOM,
    deepSplit         = deepSplit_val,
    minClusterSize    = min_genes,
    pamRespectsDendro = FALSE
  )

  tom_annotation <- data.frame(Module = as.factor(tom_clusters))
  rownames(tom_annotation) <- rownames(Mz)
  modules_unicos <- sort(unique(tom_clusters[tom_clusters > 0]))
  tom_paleta <- colorRampPalette(brewer.pal(12, "Set3"))(max(1, length(modules_unicos)))
  names(tom_paleta) <- if (length(modules_unicos)) modules_unicos else "0"

  tom_order <- gene_tree$order
  tom_plot_mat <- TOM[tom_order, tom_order, drop = FALSE]
  tom_plot_ann <- tom_annotation[tom_order, , drop = FALSE]
  pdf(file.path(output_dir, tom_pdf), width = 10, height = 10)
  pheatmap(
    tom_plot_mat,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    annotation_row    = tom_plot_ann,
    annotation_col    = tom_plot_ann,
    annotation_colors = list(Module = tom_paleta),
    color             = colorRampPalette(c("white", "steelblue", "navy"))(80),
    border_color      = NA,
    main              = sprintf("TOM Heatmap (%d genes)", nrow(Mz))
  )
  dev.off()

  module_assignments <- data.frame(
    gene_id     = rownames(Mz),
    tom_module  = tom_clusters,
    stringsAsFactors = FALSE
  )
  write.table(module_assignments,
              file.path(output_dir, "tom_gene_modules.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(
    correlation        = gene_cor,
    adjacency          = adjacency_mat,
    tom                = TOM,
    gene_tree          = gene_tree,
    tom_modules        = tom_clusters,
    module_assignments = module_assignments,
    tom_heatmap        = tom_plot
  ))
}


#' Run GO Enrichment for Gene Clusters or Modules
#'
#' Converts a gene-to-cluster assignment table into a binary membership matrix
#' and runs GO enrichment for each positive cluster/module.
#'
#' @param assignments        Data frame with columns `gene_id` and one cluster/module column.
#' @param cluster_col        Column name holding the cluster/module IDs.
#' @param output_dir         Output directory.
#' @param orgdb              OrgDb object matching the organism.
#' @param keytype            Key type matching the gene IDs.
#' @param espacio            GO namespace.
#' @param qvalue_cutoff      Q-value cutoff for enrichGO.
#' @param pvalue_cutoff      P-value cutoff for enrichGO.
#' @param simplify_cutoff    Similarity cutoff for simplify().
#' @param go_level           GO level used for pruning.
#' @param pdf_name           Output PDF filename.
#' @return Output of `run_go_enrichment_suite()`.
#' @export
run_go_for_gene_clusters <- function(assignments,
                                     cluster_col,
                                     output_dir,
                                     orgdb,
                                     keytype,
                                     espacio         = "BP",
                                     qvalue_cutoff   = 0.05,
                                     pvalue_cutoff   = 0.05,
                                     simplify_cutoff = 0.7,
                                     go_level        = 6,
                                     pdf_name        = "GO_clusters.pdf") {

  if (!all(c("gene_id", cluster_col) %in% colnames(assignments))) {
    stop("assignments must contain 'gene_id' and the requested cluster column.")
  }

  cluster_ids <- sort(unique(assignments[[cluster_col]][assignments[[cluster_col]] > 0]))
  if (!length(cluster_ids)) {
    stop("No positive clusters/modules found for GO enrichment.")
  }

  cluster_table <- data.frame(gene_id = assignments$gene_id, stringsAsFactors = FALSE)
  for (cluster_id in cluster_ids) {
    cluster_name <- paste0("cluster_", cluster_id)
    cluster_table[[cluster_name]] <- as.integer(assignments[[cluster_col]] == cluster_id)
  }

  run_go_enrichment_suite(
    diff_table      = cluster_table,
    output_dir      = output_dir,
    orgdb           = orgdb,
    keytype         = keytype,
    espacio         = espacio,
    qvalue_cutoff   = qvalue_cutoff,
    pvalue_cutoff   = pvalue_cutoff,
    simplify_cutoff = simplify_cutoff,
    go_level        = go_level,
    pdf_name        = pdf_name
  )
}

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
# 8. PIPELINE SETUP HELPERS
# =============================================================================

#' Create All Pipeline Output Directories
#'
#' Creates the standard folder structure under base_dir and returns a named
#' list of paths. Call list2env() on the result to assign dir_01…dir_08 in
#' the global environment.
#'
#' @param base_dir Root results directory (e.g. file.path(DATA_DIR, "results")).
#' @return Named list with dir_01 through dir_08 plus dir_objects.
#' @export
create_pipeline_dirs <- function(base_dir) {
  mkd <- function(name) {
    d <- file.path(base_dir, name)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    d
  }
  list(
    dir_01      = mkd("01_qc"),
    dir_02      = mkd("02_clustering"),
    dir_03      = mkd("03_annotation"),
    dir_04      = mkd("04_expression"),
    dir_05 = mkd("05_curacion"),
    dir_06      = mkd("06_de_results"),
    dir_07      = mkd("07_go"),
    dir_08      = mkd("08_networks"),
    dir_objects = mkd("objects")
  )
}


# =============================================================================
# 9. PSEUDOBULK DE & GO ENRICHMENT PIPELINE
# =============================================================================

#' Run the Full Pseudobulk DE and GO Enrichment Pipeline (Chapter 2)
#'
#' Executes sections 14–20 in one call:
#'   1. Global pseudobulk replicate-correlation heatmap
#'   2. Per-cell-type subsets
#'   3. Pseudo-replicate assignment and pseudobulk count matrices
#'   4. DESeq2 pairwise contrasts
#'   5. Volcano plots
#'   6. DE gene tables and log2FC heatmaps
#'   7. GO enrichment (full + simplified, raw + level-pruned)
#'
#' Common OrgDb / keytype combinations:
#'   Arabidopsis thaliana : orgdb = org.At.tair.db, keytype = "TAIR"
#'   Homo sapiens         : orgdb = org.Hs.eg.db,   keytype = "ENSEMBL"
#'   Mus musculus         : orgdb = org.Mm.eg.db,   keytype = "ENSEMBL"
#'
#' @param obj           Integrated Seurat object (post-harmony).
#' @param comparaciones List of lists, each with conds = c("ref","treat") and tag.
#' @param orgdb         OrgDb annotation object for GO enrichment.
#' @param keytype       Key type matching gene IDs in the Seurat object.
#' @param annot_col     Metadata column with curated cell-type labels.
#' @param n_reps        Number of pseudo-replicates per condition (default 3).
#' @param padj_cut      Adjusted p-value cutoff for DE and volcano (default 0.05).
#' @param lfc_cut       Log2 fold-change cutoff (default 1).
#' @param espacio       GO namespace: "BP", "MF", or "CC" (default "BP").
#' @param qval          Q-value cutoff for GO enrichment (default 0.05).
#' @param nivel_poda    Maximum GO hierarchy depth for term pruning (default 6).
#' @param dir_pseudobulk Output directory for sections 14–16.
#' @param dir_deseq2    Output directory for section 17.
#' @param dir_volcano   Output directory for section 18.
#' @param dir_heatmaps  Output directory for section 19.
#' @param dir_go        Output directory for section 20.
#' @return Invisible NULL. All outputs are written to disk.
#' @export
run_pseudobulk_pipeline <- function(obj,
                                     comparaciones,
                                     orgdb,
                                     keytype,
                                     annot_col      = "celltype_grouped",
                                     n_reps         = 3,
                                     padj_cut       = 0.05,
                                     lfc_cut        = 1,
                                     espacio        = "BP",
                                     qval           = 0.05,
                                     nivel_poda     = 6,
                                     dir_pseudobulk,
                                     dir_deseq2,
                                     dir_volcano,
                                     dir_heatmaps,
                                     dir_go) {

  # ── [1/6] Global pseudobulk correlation ──────────────────────────────────────
  message("[1/6] Global pseudobulk replicate correlation...")
  pseudobulk <- generate_pseudobulk(obj, group_by = "orig.ident")
  ggsave(file.path(dir_pseudobulk, "pseudobulk_correlation.pdf"),
         plot_replicate_correlation(pseudobulk$by_sample),
         width = 8, height = 8, dpi = 300)

  # ── [2/6] Per-cell-type subsets + pseudo-replicates ──────────────────────────
  message("[2/6] Building cell-type subsets and pseudo-replicates...")
  cell_types <- unique(obj@meta.data[[annot_col]])
  cell_type_subsets <- setNames(
    lapply(cell_types, function(t)
      subset(obj, cells = colnames(obj)[obj@meta.data[[annot_col]] == t])),
    gsub("[^[:alnum:]_]", "_", cell_types)
  )

  cell_type_subsets_replicates <- Filter(Negate(is.null),
    lapply(cell_type_subsets, assign_pseudo_replicates, n_reps = n_reps))
  pseudobulk_list <- lapply(cell_type_subsets_replicates, run_pseudobulk)

  rep_dir <- file.path(dir_pseudobulk, "pseudobulk_replicas")
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  for (tipo in names(pseudobulk_list))
    write.csv(pseudobulk_list[[tipo]],
              file.path(rep_dir, paste0("Pseudobulk_", tipo, ".csv")),
              row.names = TRUE)

  # ── [3/6] DESeq2 ─────────────────────────────────────────────────────────────
  message("[3/6] Running DESeq2 differential expression...")
  for (tag in sapply(comparaciones, `[[`, "tag"))
    dir.create(file.path(dir_deseq2, tag), recursive = TRUE, showWarnings = FALSE)
  for (tipo in names(pseudobulk_list))
    correr_deseq2(as.matrix(pseudobulk_list[[tipo]]), comparaciones,
                  output_dir = dir_deseq2, tipo = tipo)

  # ── [4/6] Volcano plots ───────────────────────────────────────────────────────
  message("[4/6] Generating volcano plots...")
  for (comp in comparaciones) {
    tag       <- comp$tag
    csv_files <- list.files(file.path(dir_deseq2, tag),
                            pattern = "\\.csv$", full.names = TRUE)
    if (!length(csv_files)) { message("  No CSV for: ", tag); next }

    pdf(file.path(dir_volcano, paste0("VolcanoPlots_", tag, ".pdf")),
        width = 12, height = 6)
    plots <- list()
    for (f in csv_files) {
      plots <- c(plots, list(plot_volcano(f, padj_cut = padj_cut, lfc_cut = lfc_cut)))
      if (length(plots) == 2) { grid.arrange(grobs = plots, ncol = 2); plots <- list() }
    }
    if (length(plots)) grid.arrange(grobs = plots, ncol = 1)
    dev.off()
  }

  # ── [5/6] DE tables and heatmaps ─────────────────────────────────────────────
  message("[5/6] Building DE tables and heatmaps...")
  for (comp in comparaciones) {
    tag       <- comp$tag
    diff_dir  <- file.path(dir_heatmaps, tag)
    csv_files <- list.files(file.path(dir_deseq2, tag),
                            pattern = "\\.csv$", full.names = TRUE)

    dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)
    if (!length(csv_files)) { message("  No CSV for: ", tag); next }

    listas      <- lapply(csv_files, procesar_deseq2_resultado,
                          output_dir = diff_dir, padj_cut = padj_cut, lfc_cut = lfc_cut)
    tabla_class <- Reduce(function(x, y) full_join(x, y, by = "gene_id"),
                          lapply(listas, `[[`, "class"))
    tabla_logfc <- Reduce(function(x, y) full_join(x, y, by = "gene_id"),
                          lapply(listas, `[[`, "logfc"))

    tabla_class <- tabla_class %>% filter(apply(select(., -gene_id) != 0, 1, any))
    tabla_logfc <- tabla_logfc %>% filter(gene_id %in% tabla_class$gene_id)

    write_tsv(tabla_class, file.path(diff_dir, "tabla_diferenciales.tsv"))
    write_tsv(tabla_logfc, file.path(diff_dir, "tabla_log2FC.tsv"))

    matriz <- as.matrix(column_to_rownames(tabla_logfc, "gene_id"))
    matriz[is.na(matriz)] <- 0
    if (nrow(matriz) > 1) {
      pdf(file.path(diff_dir, paste0("heatmap_", tag, ".pdf")), width = 14, height = 18)
      tryCatch(plot_heatmap(matriz), error = function(e) message("Heatmap error: ", e$message))
      dev.off()
    }
  }

  # ── [6/6] GO enrichment ───────────────────────────────────────────────────────
  message("[6/6] GO enrichment analysis...")
  universo <- keys(orgdb, keytype = keytype)

  for (comp in comparaciones) {
    tag        <- comp$tag
    enr_dir    <- file.path(dir_go, tag)
    tabla_path <- file.path(dir_heatmaps, tag, "tabla_diferenciales.tsv")

    dir.create(enr_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(tabla_path)) { message("  No differential table for: ", tag); next }

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

  message("Part 2 complete. All outputs saved.")
  invisible(NULL)
}


# =============================================================================
# 10. PIPELINE WORKFLOW FIGURE
# =============================================================================

#' Plot Pipeline Workflow Overview
#'
#' Generates a visual diagram of the full scRNA-seq pipeline and saves it as a
#' PDF. Useful as a quick reference for what each chapter covers.
#'
#' @param outfile Full path for the output PDF file.
#' @return Invisible ggplot object. Side effect: writes PDF to disk.
#' @export
plot_pipeline_workflow <- function(outfile) {

  # \u2500\u2500 Node positions \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  # Layout: inverted-Y
  #   Stem  (Part 1) \u2014 center column x = 5, y = 9 \u2192 1
  #   Junction        \u2014 Export h5ad at (5, 1)
  #   Left  branch    \u2014 Pseudotime,  x = 2,  y = -1 \u2192 -3
  #   Right branch    \u2014 Part 2 DE/Networks, x = 8, y = -1 \u2192 -5

  trunk <- data.frame(
    x = 5,
    y = c(9, 8, 7, 6, 5, 4, 3, 2, 1),
    label = c(
      "FASTQ\n+ CellRanger", "CellBender\n(optional)",
      "Load &\nQC", "Filter +\nDoubletFinder",
      "Merge &\nNormalize", "Harmony\nIntegration",
      "Clustering", "Annotate &\nCurate",
      "Export\nh5ad"
    ),
    group    = c("Pre-processing", "Pre-processing", rep("Part 1", 6), "Junction"),
    optional = c(FALSE, TRUE, rep(FALSE, 7)),
    stringsAsFactors = FALSE
  )

  left_nodes <- data.frame(
    x = 2, y = c(-1, -2, -3),
    label = c("Load h5ad\n(Python)", "Trajectory\nInference", "Pseudotime\nPlots"),
    group = "Pseudotime", optional = FALSE,
    stringsAsFactors = FALSE
  )

  right_nodes <- data.frame(
    x = 8, y = c(-1, -2, -3, -4, -5),
    label = c("Cell-type\nSubsets", "Pseudobulk\n+ Replicates",
              "DESeq2\nDE", "Volcano +\nHeatmap + GO", "Network\nInference"),
    group = "Part 2", optional = FALSE,
    stringsAsFactors = FALSE
  )

  nodes <- rbind(trunk[, c("x","y","label","group","optional")],
                 left_nodes, right_nodes)

  # \u2500\u2500 Edges \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  trunk_edges <- data.frame(x1 = 5, y1 = 9:2, x2 = 5, y2 = 8:1)

  fork_edges  <- data.frame(      # Export \u2192 branch starts
    x1 = c(5, 5), y1 = c(1, 1),
    x2 = c(2, 8), y2 = c(-1, -1)
  )

  left_edges  <- data.frame(x1 = 2, y1 = c(-1,-2), x2 = 2, y2 = c(-2,-3))
  right_edges <- data.frame(x1 = 8, y1 = -1:-4,    x2 = 8, y2 = -2:-5)

  edges <- rbind(trunk_edges, fork_edges, left_edges, right_edges)

  # \u2500\u2500 Colors \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  group_colors <- c(
    "Pre-processing" = "#d9d9d9",
    "Part 1"         = "#b3cde3",
    "Junction"       = "#fdcdac",
    "Pseudotime"     = "#e0c8f0",
    "Part 2"         = "#ccebc5"
  )

  w <- 0.82; h <- 0.52

  req  <- nodes[!nodes$optional, ]
  opt  <- nodes[ nodes$optional, ]

  p <- ggplot() +
    # Section labels
    annotate("text", x = 5,   y = 9.65,
             label = "PART 1 \u2014 Single-Cell Analysis",
             size = 3.5, color = "#1565c0", fontface = "bold") +
    annotate("text", x = 2,   y = -0.3,
             label = "Pseudotime", size = 3.2, color = "#6a1b9a", fontface = "bold") +
    annotate("text", x = 8,   y = -0.3,
             label = "PART 2 \u2014 Pseudobulk DE & Networks",
             size = 3.2, color = "#2e7d32", fontface = "bold") +
    # Background bands
    annotate("rect", xmin = 3.9, xmax = 6.1, ymin = 0.35, ymax = 9.45,
             fill = "#f0f4ff", color = NA) +
    annotate("rect", xmin = 0.7, xmax = 3.3, ymin = -3.65, ymax = -0.45,
             fill = "#f8f0ff", color = NA) +
    annotate("rect", xmin = 6.7, xmax = 9.3, ymin = -5.65, ymax = -0.45,
             fill = "#f0fff4", color = NA) +
    # Edges
    geom_segment(data = edges,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 arrow     = arrow(length = unit(0.20, "cm"), type = "closed"),
                 linewidth = 0.55, color = "grey45", lineend = "round") +
    # Required nodes
    geom_rect(data = req,
              aes(xmin = x-w/2, xmax = x+w/2, ymin = y-h/2, ymax = y+h/2, fill = group),
              color = "grey35", linewidth = 0.4) +
    # Optional nodes (dashed border)
    geom_rect(data = opt,
              aes(xmin = x-w/2, xmax = x+w/2, ymin = y-h/2, ymax = y+h/2, fill = group),
              color = "grey55", linewidth = 0.4, linetype = "dashed") +
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              size = 2.6, lineheight = 0.9) +
    scale_fill_manual(values = group_colors, name = NULL,
                      guide  = guide_legend(nrow = 1)) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 9),
          plot.margin     = margin(15, 20, 15, 20)) +
    xlim(0, 10) + ylim(-6.3, 10.2)

  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  ggsave(outfile, p, width = 12, height = 20, dpi = 300)
  message("Workflow figure saved: ", outfile)
  invisible(p)
}


# =============================================================================
# END OF FUNCTIONS
# =============================================================================


# =============================================================================
# DATA LOADING
# =============================================================================

#' Load Seurat objects from samples (CellBender or CellRanger)
#'
#' @param samples List of sample configurations (file, label, condition)
#' @param DATA_DIR Root data directory path
#' @param USE_CELLBENDER Logical; if TRUE load CellBender HDF5, if FALSE load CellRanger
#' @param mt_pattern Regex pattern for mitochondrial genes (e.g., "^ATMG" for Arabidopsis)
#' @param cp_pattern Regex pattern for chloroplast genes (e.g., "^ATCG" for Arabidopsis)
#'
#' @return Named list of Seurat objects with QC metrics (percent.mt, percent.cp)
#'
#' @export
load_seurat_samples <- function(samples, DATA_DIR, USE_CELLBENDER, mt_pattern, cp_pattern) {
  
  if (USE_CELLBENDER) {
    # Load CellBender-filtered HDF5 files
    seurat_list <- lapply(samples, load_sample,
                          mt_pattern = mt_pattern,
                          cp_pattern = cp_pattern)
  } else {
    # Load CellRanger filtered_feature_bc_matrix directories
    seurat_list <- lapply(samples, function(s) {
      mat <- Read10X(data.dir = file.path(DATA_DIR, s$file))
      obj <- CreateSeuratObject(counts = mat, project = s$label,
                                min.cells = 3, min.features = 200)
      obj$condition <- s$condition
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
      if (!is.null(cp_pattern))
        obj[["percent.cp"]] <- PercentageFeatureSet(obj, pattern = cp_pattern)
      obj
    })
  }
  
  names(seurat_list) <- sapply(samples, `[[`, "label")
  return(seurat_list)
}



# =============================================================================
# CELL-TYPE SUBSETTING
# =============================================================================

#' Create Seurat subsets for each cell type
#'
#' @param seurat_obj Seurat object with cell-type annotations
#' @param annot_col Metadata column name containing cell-type labels
#'
#' @return Named list of Seurat subsets (one per cell type), with sanitized names
#'
#' @details
#' Cell-type names are sanitised (special characters replaced with "_") so they
#' can be used safely as list names and output filenames.
#'
#' @export
create_cell_type_subsets <- function(seurat_obj, annot_col = "celltype_grouped") {
  
  # Get unique cell types (excluding NA)
  cell_types <- sort(unique(na.omit(seurat_obj@meta.data[[annot_col]])))
  
  # Create subsets and sanitise names
  subsets <- setNames(
    lapply(cell_types, function(tipo) {
      subset(seurat_obj,
             cells = colnames(seurat_obj)[seurat_obj@meta.data[[annot_col]] == tipo])
    }),
    gsub("[^[:alnum:]_]", "_", cell_types)
  )
  
  # Print summary
  cell_counts <- setNames(
    vapply(subsets, function(x) as.integer(ncol(x)), integer(1)),
    names(subsets)
  )
  
  message("Cell-type subsets created:")
  print(cell_counts)
  
  return(subsets)
}



# =============================================================================
# PSEUDOBULK PREPARATION
# =============================================================================

#' Assign pseudo-replicates to cell-type subsets
#'
#' @param cell_type_subsets Named list of Seurat objects (one per cell type)
#' @param pseudobulk_conditions Optional character vector of conditions to retain (NULL = all)
#' @param n_pseudoreps Number of pseudo-replicates per condition
#'
#' @return Named list of Seurat subsets with pseudo-replicate assignments,
#'         filtered to include only subsets with ≥2 conditions
#'
#' @details
#' The random seed must be set globally before calling this function (e.g., set.seed(1807))
#' to ensure reproducibility.
#'
#' @export
assign_pseudoreplicates_batch <- function(cell_type_subsets,
                                          pseudobulk_conditions = NULL,
                                          n_pseudoreps = 3) {
  
  subsets_with_reps <- Filter(
    Negate(is.null),
    lapply(cell_type_subsets,
           assign_pseudo_replicates,
           conditions = pseudobulk_conditions,
           n_reps = n_pseudoreps)
  )
  
  message("Cell types retained for pseudobulk (≥2 conditions):")
  print(names(subsets_with_reps))
  
  return(subsets_with_reps)
}



# =============================================================================
# PSEUDOBULK AND DESEQ2 ANALYSIS
# =============================================================================

#' Run pseudobulk aggregation and DESeq2 analysis
#'
#' @param cell_type_subsets_replicates Named list of Seurat objects with pseudo-replicate assignments
#' @param comparisons List of contrast definitions, each with:
#'                    - conds: character vector c("reference", "treatment")
#'                    - tag: character label for output folder
#' @param output_dir Base output directory for results
#' @param cell_types Optional character vector of specific cell types to process
#'                   (NULL = process all)
#' @param pseudobulk_dir Directory to save pseudobulk count tables
#'
#' @return Named list of DESeq2 results per cell type
#'
#' @export
run_pseudobulk_deseq2_analysis <- function(cell_type_subsets_replicates,
                                           comparisons,
                                           output_dir,
                                           cell_types = NULL,
                                           pseudobulk_dir = NULL) {
  
  # Set default pseudobulk directory
  if (is.null(pseudobulk_dir))
    pseudobulk_dir <- file.path(dirname(output_dir), "pseudobulk_replicas")
  
  # Create output directories
  dir.create(pseudobulk_dir, recursive = TRUE, showWarnings = FALSE)
  for (tag in sapply(comparisons, `[[`, "tag")) {
    dir.create(file.path(output_dir, tag), recursive = TRUE, showWarnings = FALSE)
  }
  
  # If cell_types not specified, use all available
  if (is.null(cell_types))
    cell_types <- names(cell_type_subsets_replicates)
  
  # Generate pseudobulk count tables
  message("Generating pseudobulk count tables...")
  pseudobulk_list <- guardar_tablas_pseudobulk(
    cell_type_subsets_replicates,
    output_dir = pseudobulk_dir
  )
  
  # Filter to requested cell types
  pseudobulk_list <- pseudobulk_list[names(pseudobulk_list) %in% cell_types]
  
  message("Processing ", length(pseudobulk_list), " cell types:")
  print(names(pseudobulk_list))
  
  # Run DESeq2 for each cell type (shows all output)
  deseq2_results <- list()
  for (tipo in names(pseudobulk_list)) {
    message("
", strrep("─", 70))
    message("DESeq2 analysis for cell type: ", tipo)
    message(strrep("─", 70))
    
    deseq2_results[[tipo]] <- correr_deseq2(
      counts_mat = as.matrix(pseudobulk_list[[tipo]]),
      comparaciones = comparisons,
      output_dir = output_dir,
      tipo = tipo
    )
  }
  
  message("
", strrep("═", 70))
  message("✓ DESeq2 analysis complete for all cell types")
  message(strrep("═", 70))
  return(deseq2_results)
}



# =============================================================================
# GO ENRICHMENT (SIMPLE)
# =============================================================================

#' Simple GO enrichment analysis
#'
#' @param diff_table Path to differential expression table (gene IDs in first column)
#' @param output_dir Output directory for results
#' @param orgdb OrgDb object for the organism (e.g., org.At.tair.db)
#' @param keytype Key type matching gene IDs (e.g., "TAIR" for Arabidopsis)
#' @param go_space Ontology namespace: "BP" (biological process), "MF" (molecular function), "CC" (cellular component)
#' @param padj_cutoff Adjusted p-value threshold
#'
#' @return Enrichment results with GO terms, counts, and visualizations
#'
#' @export
run_simple_go_enrichment <- function(diff_table,
                                     output_dir,
                                     orgdb,
                                     keytype = "TAIR",
                                     go_space = "BP",
                                     padj_cutoff = 0.05,
                                     cell_type = NULL,
                                     contrast_tag = NULL) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Read table
  if (is.character(diff_table)) {
    tabla_df <- read.table(diff_table, header = TRUE, sep = "	", check.names = FALSE)
  } else {
    tabla_df <- as.data.frame(diff_table, check.names = FALSE)
  }

  gene_col <- colnames(tabla_df)[1]
  genes_all <- tabla_df[[gene_col]]

  # Get universe of genes
  universo <- keys(orgdb, keytype = keytype)

  # Run enrichment
  go_result <- enrichGO(
    gene = genes_all,
    universe = universo,
    OrgDb = orgdb,
    keyType = keytype,
    ont = go_space,
    pvalueCutoff = padj_cutoff,
    pAdjustMethod = "BH",
    readable = TRUE
  )

  # If no enrichment results, return early
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    return(NULL)
  }

  # Build filename prefix
  file_prefix <- "GO"
  if (!is.null(cell_type)) file_prefix <- paste0(file_prefix, "_", cell_type)
  if (!is.null(contrast_tag)) file_prefix <- paste0(file_prefix, "_", contrast_tag)
  file_prefix <- paste0(file_prefix, "_", go_space)

  # Save results
  write.table(go_result@result,
              file = file.path(output_dir, paste0(file_prefix, ".tsv")),
              sep = "	", quote = FALSE, row.names = FALSE)

  # Generate plots with title
  plot_title <- if (!is.null(cell_type)) paste0("GO Enrichment - ", cell_type) else "GO Enrichment"

  tryCatch({
    pdf(file.path(output_dir, paste0(file_prefix, "_bubble.pdf")), width = 12, height = 8)
    p <- dotplot(go_result, showCategory = 20) + ggtitle(plot_title)
    print(p)
    dev.off()
  }, error = function(e) message("Could not generate bubble plot: ", e$message))

  return(go_result)
}


# =============================================================================
# build_logfc_heatmap
# =============================================================================
# Builds a log2FC heatmap per cell type.
# Genes are ordered by the internal hierarchical clustering dendrogram
# (cluster_rows = TRUE) — no manual clustering or cluster assignment.
build_logfc_heatmap <- function(logfc_table,
                                contrast_tag,
                                output_dir,
                                limits = c(-5, 5)) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  mat <- as.matrix(logfc_table[, -1])
  rownames(mat) <- logfc_table[[1]]
  colnames(mat) <- gsub(paste0("_", contrast_tag, "$"), "", colnames(mat))
  mat[is.na(mat)] <- 0

  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "log2FC",
    col  = circlize::colorRamp2(c(limits[1], 0, limits[2]), c("blue", "black", "yellow")),

    cluster_rows    = TRUE,
    cluster_columns = FALSE,

    show_row_names    = FALSE,
    show_column_names = TRUE,
    column_names_gp   = grid::gpar(fontsize = 14, fontface = "bold"),
    column_names_rot  = 45,

    column_title    = sprintf("log2FC — %s  (%d genes)", contrast_tag, nrow(mat)),
    column_title_gp = gpar(fontsize = 15, fontface = "bold"),

    heatmap_legend_param = list(
      title         = "log2FC",
      title_gp      = gpar(fontsize = 12, fontface = "bold"),
      legend_height = unit(4, "cm")
    )
  )

  pdf(file.path(output_dir, paste0("heatmap_", contrast_tag, ".pdf")),
      width = 12, height = 16)
  ComplexHeatmap::draw(ht, merge_legend = TRUE, padding = grid::unit(c(2, 2, 2, 10), "mm"))
  dev.off()

  invisible(NULL)
}



# =============================================================================
# get_tfs_from_orgdb
# =============================================================================
# Generic helper: returns transcription factor IDs from any Bioconductor OrgDb
# using GO:0003700 (DNA-binding transcription factor activity).
# Works for any organism with GO annotation (Arabidopsis, human, mouse, etc.).
get_tfs_from_orgdb <- function(orgdb, keytype = "TAIR") {
  tryCatch({
    tfs <- AnnotationDbi::select(orgdb,
                                  keys    = "GO:0003700",
                                  columns = keytype,
                                  keytype = "GO")[[keytype]]
    unique(stats::na.omit(tfs))
  }, error = function(e) {
    message("Could not retrieve TFs via GO:0003700 (", e$message, ")")
    character(0)
  })
}


# =============================================================================
# load_pseudobulk_matrix
# =============================================================================
# Internal helper: loads and merges pseudobulk replicate counts from a directory
# of files named Pseudobulk_Reps_<celltype>.csv. Returns a CPM + log2 matrix.
load_pseudobulk_matrix <- function(pseudobulk_dir, normalize = TRUE) {
  pb_files <- list.files(pseudobulk_dir,
                         pattern = "^Pseudobulk_Reps_.*\\.csv$",
                         full.names = TRUE)
  if (!length(pb_files)) stop("No Pseudobulk_Reps_*.csv files in: ", pseudobulk_dir)
  pb_list <- lapply(pb_files, function(f) {
    df <- read.csv(f, row.names = 1, check.names = FALSE)
    ct <- gsub("Pseudobulk_Reps_|\\.csv$", "", basename(f))
    colnames(df) <- paste0(ct, "_", colnames(df))
    df
  })
  common <- Reduce(intersect, lapply(pb_list, rownames))
  mat    <- as.matrix(do.call(cbind, lapply(pb_list, function(d) d[common, ])))

  if (normalize) {
    lib_sizes <- colSums(mat)
    mat <- sweep(mat, 2, lib_sizes / 1e6, FUN = "/")
    mat <- log2(mat + 1)
  }
  mat
}


# =============================================================================
# run_genie3_per_cluster
# =============================================================================
# GENIE3-only analysis per cluster. Uses transcription factors as regulators
# (TF -> target directed edges) and filters edges by absolute Pearson correlation.
#
# Generic across organisms: provide an OrgDb (org.At.tair.db, org.Hs.eg.db, ...)
# and a matching keytype, or pass a custom_tfs vector to override.
run_genie3_per_cluster <- function(cluster_assignments,
                                   pseudobulk_dir,
                                   output_dir,
                                   orgdb,
                                   keytype        = "TAIR",
                                   custom_tfs     = NULL,
                                   n_top_clusters = 3,
                                   cor_min        = 0.90,
                                   genie3_ntrees  = 100,
                                   n_cores        = 4,
                                   min_var_filter = 0.01) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ── Pseudobulk: CPM + log2 ──────────────────────────────────────────────────
  exprMatr_full <- load_pseudobulk_matrix(pseudobulk_dir, normalize = TRUE)
  message("Pseudobulk matrix: ", nrow(exprMatr_full), " genes x ",
          ncol(exprMatr_full), " samples (CPM + log2)")

  # ── Transcription factor list ───────────────────────────────────────────────
  tf_list <- if (!is.null(custom_tfs)) {
    unique(custom_tfs)
  } else {
    get_tfs_from_orgdb(orgdb, keytype)
  }
  if (!length(tf_list)) stop("Empty TF list — provide custom_tfs or check OrgDb/keytype")
  message("TFs available (organism-wide): ", length(tf_list))

  # ── Top clusters ────────────────────────────────────────────────────────────
  ca <- cluster_assignments[cluster_assignments$cluster != "grey", ]
  cluster_sizes <- sort(table(ca$cluster), decreasing = TRUE)
  n_use <- if (is.null(n_top_clusters)) length(cluster_sizes) else min(n_top_clusters, length(cluster_sizes))
  top_clusters  <- names(cluster_sizes)[seq_len(n_use)]

  # ── Per-cluster GENIE3 ──────────────────────────────────────────────────────
  results <- list()
  for (clust_id in top_clusters) {
    message("\n── GENIE3 on cluster: ", clust_id, " ──")
    genes_ok <- intersect(ca$gene_id[ca$cluster == clust_id], rownames(exprMatr_full))
    expr     <- exprMatr_full[genes_ok, , drop = FALSE]
    expr     <- expr[apply(expr, 1, var) > min_var_filter, , drop = FALSE]

    TFreg   <- intersect(tf_list, rownames(expr))
    targets <- setdiff(rownames(expr), TFreg)
    message("  Genes: ", nrow(expr), " | TFs: ", length(TFreg), " | Targets: ", length(targets))

    if (length(TFreg) == 0 || length(targets) == 0) {
      message("  Skip (no TFs or no targets)"); next
    }

    set.seed(123)
    weightMat <- GENIE3(expr,
                        regulators = TFreg, targets = targets,
                        treeMethod = "RF", K = "sqrt",
                        nCores = n_cores, verbose = FALSE)

    linkList <- getLinkList(weightMat)
    colnames(linkList) <- c("source", "target", "weight")
    dt <- data.table::as.data.table(linkList)
    dt[, source := as.character(source)]
    dt[, target := as.character(target)]
    dt <- dt[weight > 0]

    cor_matrix <- cor(t(expr), method = "pearson")
    dt[, pearson_cor := abs(cor_matrix[cbind(source, target)])]

    final <- dt[pearson_cor >= cor_min]
    data.table::setorder(final, -pearson_cor, -weight)
    message("  All edges (weight > 0): ", nrow(dt),
            " | Filtered (|r| >= ", cor_min, "): ", nrow(final))

    # Save outputs
    write.table(dt[, .(source, target, weight, pearson_cor)],
                file.path(output_dir, paste0("GENIE3_", clust_id, "_all_edges.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(final[, .(source, target, weight, pearson_cor)],
                file.path(output_dir,
                          paste0("GENIE3_", clust_id,
                                 "_cor", round(cor_min * 100), ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    results[[clust_id]] <- list(
      cluster      = clust_id,
      n_genes      = nrow(expr),
      n_tfs        = length(TFreg),
      n_targets    = length(targets),
      all_edges    = dt,
      filtered     = final,
      n_all        = nrow(dt),
      n_filtered   = nrow(final)
    )
  }

  saveRDS(results, file.path(output_dir, "GENIE3_results.rds"))
  message("\n✓ GENIE3 outputs saved to: ", output_dir)

  invisible(results)
}


# =============================================================================
# run_wgcna_per_cluster
# =============================================================================
# WGCNA-only coexpression analysis per cluster. Builds a TOM-based undirected
# network and filters edges above a TOM threshold.
run_wgcna_per_cluster <- function(cluster_assignments,
                                  pseudobulk_dir,
                                  output_dir,
                                  n_top_clusters = 3,
                                  soft_power     = 6,
                                  network_type   = "signed",
                                  tom_threshold  = 0.10,
                                  min_var_filter = 0.01) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  exprMatr_full <- load_pseudobulk_matrix(pseudobulk_dir, normalize = TRUE)
  message("Pseudobulk matrix: ", nrow(exprMatr_full), " genes x ",
          ncol(exprMatr_full), " samples (CPM + log2)")

  ca <- cluster_assignments[cluster_assignments$cluster != "grey", ]
  cluster_sizes <- sort(table(ca$cluster), decreasing = TRUE)
  n_use <- if (is.null(n_top_clusters)) length(cluster_sizes) else min(n_top_clusters, length(cluster_sizes))
  top_clusters  <- names(cluster_sizes)[seq_len(n_use)]

  results <- list()
  for (clust_id in top_clusters) {
    message("\n── WGCNA on cluster: ", clust_id, " ──")
    genes_ok <- intersect(ca$gene_id[ca$cluster == clust_id], rownames(exprMatr_full))
    expr     <- exprMatr_full[genes_ok, , drop = FALSE]
    expr     <- expr[apply(expr, 1, var) > min_var_filter, , drop = FALSE]
    if (nrow(expr) < 10) { message("  Skip (<10 genes)"); next }

    adj <- adjacency(t(expr), power = soft_power, type = network_type)
    TOM <- TOMsimilarity(adj, TOMType = network_type, verbose = 0)
    rownames(TOM) <- colnames(TOM) <- rownames(expr)

    idx <- which(upper.tri(TOM), arr.ind = TRUE)
    edges <- data.frame(
      source = rownames(TOM)[idx[, 1]],
      target = colnames(TOM)[idx[, 2]],
      adjacency = adj[idx],
      TOM       = TOM[idx],
      stringsAsFactors = FALSE
    )
    final <- edges[edges$TOM >= tom_threshold, ]
    final <- final[order(-final$TOM), ]
    message("  Total pairs: ", nrow(edges),
            " | Filtered (TOM >= ", tom_threshold, "): ", nrow(final))

    write.table(edges,
                file.path(output_dir, paste0("WGCNA_", clust_id, "_all_edges.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(final,
                file.path(output_dir,
                          paste0("WGCNA_", clust_id,
                                 "_TOM", gsub("\\.", "", as.character(tom_threshold)),
                                 ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    results[[clust_id]] <- list(
      cluster    = clust_id,
      n_genes    = nrow(expr),
      all_edges  = edges,
      filtered   = final,
      n_all      = nrow(edges),
      n_filtered = nrow(final)
    )
  }

  saveRDS(results, file.path(output_dir, "WGCNA_results.rds"))
  message("\n✓ WGCNA outputs saved to: ", output_dir)

  invisible(results)
}


# =============================================================================
# run_synergistic_network
# =============================================================================
# SYNERGISTIC analysis (NOT a consensus / intersection).
#
# Each method contributes what the other cannot:
#   - GENIE3 provides DIRECTIONALITY (TF -> target) and predictive power
#   - WGCNA  provides COEXPRESSION ROBUSTNESS via TOM (shared neighbors)
#
# For every TF -> target edge proposed by GENIE3, the WGCNA TOM of that pair
# is looked up. The synergistic score combines both with a geometric mean of
# rank-normalized values, which penalizes any edge where one side is weak:
#
#   score_synergy = sqrt( rank_norm(weight_genie3) * rank_norm(TOM) )
#
# Output: directed edges (TF -> target) with both metrics + synergistic score.
run_synergistic_network <- function(cluster_assignments,
                                    pseudobulk_dir,
                                    output_dir,
                                    orgdb,
                                    keytype        = "TAIR",
                                    custom_tfs     = NULL,
                                    n_top_clusters = 3,
                                    soft_power     = 6,
                                    network_type   = "signed",
                                    genie3_ntrees  = 100,
                                    n_cores        = 4,
                                    min_var_filter = 0.01,
                                    cor_min        = 0.90,
                                    tom_min        = 0.15) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  exprMatr_full <- load_pseudobulk_matrix(pseudobulk_dir, normalize = TRUE)
  message("Pseudobulk matrix: ", nrow(exprMatr_full), " genes x ",
          ncol(exprMatr_full), " samples (CPM + log2)")

  tf_list <- if (!is.null(custom_tfs)) unique(custom_tfs) else get_tfs_from_orgdb(orgdb, keytype)
  if (!length(tf_list)) stop("Empty TF list — provide custom_tfs or check OrgDb/keytype")
  message("TFs available (organism-wide): ", length(tf_list))

  ca <- cluster_assignments[cluster_assignments$cluster != "grey", ]
  cluster_sizes <- sort(table(ca$cluster), decreasing = TRUE)
  n_use <- if (is.null(n_top_clusters)) length(cluster_sizes) else min(n_top_clusters, length(cluster_sizes))
  top_clusters  <- names(cluster_sizes)[seq_len(n_use)]

  rank_norm <- function(x) rank(x, ties.method = "average") / length(x)

  results <- list()
  for (clust_id in top_clusters) {
    message("\n── Synergistic analysis on cluster: ", clust_id, " ──")
    genes_ok <- intersect(ca$gene_id[ca$cluster == clust_id], rownames(exprMatr_full))
    expr     <- exprMatr_full[genes_ok, , drop = FALSE]
    expr     <- expr[apply(expr, 1, var) > min_var_filter, , drop = FALSE]
    if (nrow(expr) < 10) { message("  Skip (<10 genes)"); next }

    TFreg   <- intersect(tf_list, rownames(expr))
    targets <- setdiff(rownames(expr), TFreg)
    if (!length(TFreg) || !length(targets)) { message("  Skip (no TFs/targets)"); next }
    message("  Genes: ", nrow(expr), " | TFs: ", length(TFreg),
            " | Targets: ", length(targets))

    # WGCNA layer (coexpression backbone)
    adj <- adjacency(t(expr), power = soft_power, type = network_type)
    TOM <- TOMsimilarity(adj, TOMType = network_type, verbose = 0)
    rownames(TOM) <- colnames(TOM) <- rownames(expr)
    cor_matrix <- cor(t(expr), method = "pearson")

    # GENIE3 layer (directed TF -> target)
    set.seed(123)
    weightMat <- GENIE3(expr, regulators = TFreg, targets = targets,
                        treeMethod = "RF", K = "sqrt",
                        nCores = n_cores, verbose = FALSE)
    linkList <- getLinkList(weightMat)
    colnames(linkList) <- c("source", "target", "weight_genie3")
    dt <- data.table::as.data.table(linkList)
    dt[, source := as.character(source)]
    dt[, target := as.character(target)]
    dt <- dt[weight_genie3 > 0]

    # Lookup WGCNA metrics for each directed edge
    dt[, TOM         := TOM[cbind(source, target)]]
    dt[, adjacency   := adj[cbind(source, target)]]
    dt[, pearson_cor := abs(cor_matrix[cbind(source, target)])]

    # Synergistic score: geometric mean of rank-normalized values
    dt[, score_synergy := sqrt(rank_norm(weight_genie3) * rank_norm(TOM))]
    data.table::setorder(dt, -score_synergy)

    # Filter: BOTH layers must pass their own threshold
    final <- dt[pearson_cor >= cor_min & TOM >= tom_min]
    data.table::setorder(final, -score_synergy)
    message("  All TF->target edges: ", nrow(dt),
            " | Filtered (|r|>=", cor_min, " AND TOM>=", tom_min, "): ",
            nrow(final))

    write.table(dt[, .(source, target, weight_genie3, TOM, adjacency,
                       pearson_cor, score_synergy)],
                file.path(output_dir, paste0("SYNERGY_", clust_id, "_all_edges.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(final[, .(source, target, weight_genie3, TOM, adjacency,
                          pearson_cor, score_synergy)],
                file.path(output_dir, paste0("SYNERGY_", clust_id, "_filtered.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    results[[clust_id]] <- list(
      cluster    = clust_id,
      n_genes    = nrow(expr),
      n_tfs      = length(TFreg),
      n_targets  = length(targets),
      all_edges  = dt,
      filtered   = final,
      n_all      = nrow(dt),
      n_filtered = nrow(final)
    )
  }

  saveRDS(results, file.path(output_dir, "synergy_results.rds"))
  message("\n✓ Synergistic outputs saved to: ", output_dir)

  invisible(results)
}


# =============================================================================
# generate_network_pdf
# =============================================================================
# Generic PDF report for any of the three network methods (GENIE3, WGCNA, SYNERGY).
# Per cluster: network plot, top hubs bar chart, summary text.
#
# Expects results = list-of-clusters with each cluster providing:
#   $cluster, $n_genes, $filtered (data.frame/data.table with source, target,
#                                  + a "weight_col" used for edge width)
generate_network_pdf <- function(results,
                                 output_dir,
                                 method_name,
                                 weight_col,
                                 directed     = FALSE,
                                 n_top_hubs   = 10,
                                 max_nodes    = 80,
                                 edge_color   = "#1f78b4") {

  if (!length(results)) { message("No results to plot for ", method_name); return(invisible(NULL)) }

  pdf_path <- file.path(output_dir, paste0(method_name, "_network_report.pdf"))
  pdf(pdf_path, width = 12, height = 10)

  # ── Cover page with description ─────────────────────────────────────────────
  grid::grid.newpage()
  grid::grid.text(paste0(method_name, " - Network Analysis"),
                  y = 0.92, gp = grid::gpar(fontsize = 22, fontface = "bold"))
  grid::grid.text(paste0(length(results), " clusters analyzed"),
                  y = 0.86, gp = grid::gpar(fontsize = 12, fontface = "italic"))

  desc <- switch(method_name,
    "GENIE3" = paste(
      "GENIE3 — Gene Network Inference with Ensemble of Trees",
      "",
      "  • Random Forest based regulatory network inference",
      "  • Edges are DIRECTED: regulator (TF) -> target",
      "  • Edge score = variable importance from RF model",
      "  • Filtered by absolute Pearson correlation",
      "",
      "Each edge represents a putative TF -> target interaction.",
      "Hubs are TFs that regulate many targets in the cluster.",
      sep = "\n"
    ),
    "WGCNA" = paste(
      "WGCNA — Weighted Gene Co-expression Network Analysis",
      "",
      "  • Pearson correlation raised to soft-power (default 6)",
      "  • TOM (Topological Overlap Measure) refines using shared neighbors",
      "  • Edges are UNDIRECTED",
      "  • Filtered by TOM threshold",
      "",
      "Each edge represents two genes that co-vary AND share many neighbors.",
      "Hubs are genes highly connected within the co-expression module.",
      sep = "\n"
    ),
    "SYNERGY" = paste(
      "SYNERGY — Complementary GENIE3 + WGCNA",
      "",
      "  • Each method contributes what the other cannot:",
      "      GENIE3 -> directionality (TF -> target)",
      "      WGCNA  -> coexpression robustness via TOM",
      "  • score_synergy = sqrt(rank(weight_genie3) * rank(TOM))",
      "  • Edges retained only if BOTH thresholds are met",
      "",
      "Hubs are TFs whose targets are also strongly co-expressed.",
      "The strongest candidates for biological follow-up.",
      sep = "\n"
    ),
    paste("Method:", method_name)
  )
  grid::grid.text(desc, x = 0.05, y = 0.5, just = c("left", "center"),
                  gp = grid::gpar(fontsize = 11, fontfamily = "mono"))

  # ── Summary table ───────────────────────────────────────────────────────────
  summary_df <- do.call(rbind, lapply(results, function(r) data.frame(
    Cluster        = r$cluster,
    Genes          = r$n_genes,
    Edges_total    = r$n_all,
    Edges_filtered = r$n_filtered
  )))
  grid::grid.newpage()
  grid::grid.text("Summary across clusters",
                  y = 0.92, gp = grid::gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(y = 0.55, height = 0.7))
  grid::grid.draw(gridExtra::tableGrob(summary_df, rows = NULL,
                                       theme = gridExtra::ttheme_default(base_size = 12)))
  popViewport()

  # ── Per-cluster pages ───────────────────────────────────────────────────────
  for (r in results) {
    edges <- as.data.frame(r$filtered)
    if (nrow(edges) == 0) next

    g <- igraph::graph_from_data_frame(
      edges[, c("source", "target", weight_col)],
      directed = directed
    )

    if (length(igraph::V(g)) > max_nodes) {
      keep <- names(sort(igraph::degree(g), decreasing = TRUE))[seq_len(max_nodes)]
      g <- igraph::induced_subgraph(g, keep)
    }
    deg <- sort(igraph::degree(g), decreasing = TRUE)
    hubs <- names(deg)[seq_len(min(n_top_hubs, length(deg)))]

    # Page A: network
    igraph::V(g)$color <- ifelse(igraph::V(g)$name %in% hubs, "#ff7f0e", "#cccccc")
    igraph::V(g)$size  <- ifelse(igraph::V(g)$name %in% hubs, 8, 3)
    igraph::V(g)$label <- ifelse(igraph::V(g)$name %in% hubs, igraph::V(g)$name, "")
    par(mar = c(0, 0, 3, 0))
    plot(g, layout = igraph::layout_with_fr(g),
         vertex.label.cex   = 0.6,
         vertex.label.color = "black",
         vertex.frame.color = NA,
         edge.color = adjustcolor(edge_color, alpha.f = 0.4),
         edge.width = scales::rescale(igraph::edge_attr(g, weight_col), to = c(0.3, 1.8)),
         edge.arrow.size = if (directed) 0.3 else 0,
         main = sprintf("%s - Cluster '%s'  (%d nodes, %d edges)",
                        method_name, r$cluster,
                        length(igraph::V(g)), length(igraph::E(g))))

    # Page B: hub barchart
    hubs_df <- data.frame(gene = hubs, degree = deg[seq_len(length(hubs))])
    print(ggplot2::ggplot(hubs_df,
                          ggplot2::aes(x = stats::reorder(gene, degree), y = degree)) +
          ggplot2::geom_col(fill = edge_color) +
          ggplot2::coord_flip() +
          ggplot2::labs(title    = sprintf("Top %d hubs - cluster '%s'",
                                            n_top_hubs, r$cluster),
                        subtitle = paste0("Method: ", method_name),
                        x = NULL, y = "Degree") +
          ggplot2::theme_bw(base_size = 12))

    # Page C: top edges table + interpretation
    grid::grid.newpage()
    grid::grid.text(sprintf("Cluster '%s' - Details", r$cluster),
                    y = 0.97, gp = grid::gpar(fontsize = 16, fontface = "bold"))

    # Top 15 edges table (above)
    top15 <- head(edges, 15)
    pushViewport(viewport(y = 0.70, height = 0.45))
    grid::grid.text(paste0("Top 15 edges (sorted by ", weight_col, ")"),
                    y = 1.02, gp = grid::gpar(fontsize = 12, fontface = "bold"))
    grid::grid.draw(gridExtra::tableGrob(
      do.call(data.frame, lapply(top15, function(x)
        if (is.numeric(x)) round(x, 4) else x)),
      rows = NULL,
      theme = gridExtra::ttheme_default(base_size = 8)))
    popViewport()

    # Interpretation text (below)
    interp <- paste(
      sprintf("Genes in cluster (after variance filter): %d", r$n_genes),
      sprintf("Total edges (unfiltered): %d", r$n_all),
      sprintf("Filtered edges in this network: %d", r$n_filtered),
      "",
      sprintf("Top hubs by degree:"),
      paste0("  ", paste(hubs, collapse = ", ")),
      "",
      sprintf("Strongest edge: %s %s %s  (%s = %.4f)",
              top15$source[1],
              if (directed) "->" else "--",
              top15$target[1],
              weight_col,
              as.numeric(top15[[weight_col]][1])),
      sep = "\n"
    )
    grid::grid.text(interp, x = 0.05, y = 0.20, just = c("left", "center"),
                    gp = grid::gpar(fontsize = 10, fontfamily = "mono"))
  }

  dev.off()
  message("✓ PDF saved: ", pdf_path)
  invisible(pdf_path)
}


# =============================================================================
# visualize_network_per_cluster
# =============================================================================
# Create clean, force-directed network visualization for each cluster.
# Uses igraph with Fruchterman-Reingold layout.
visualize_network_per_cluster <- function(network_results,
                                          cluster_assignments,
                                          output_dir,
                                          method_name = "GENIE3") {

  weight_col <- if (method_name == "GENIE3") "weight" else "TOM"
  directed   <- method_name == "GENIE3"
  edge_color <- if (method_name == "GENIE3") "#2ca02c" else "#1f77b4"

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!length(network_results)) {
    message("No network results to visualize for ", method_name)
    return(invisible(NULL))
  }

  pdf_path <- file.path(output_dir, paste0(method_name, "_network_visualization.pdf"))
  pdf(pdf_path, width = 14, height = 11)

  grid::grid.newpage()
  grid::grid.text(paste0(method_name, " - Network Visualization"),
                  y = 0.92, gp = grid::gpar(fontsize = 22, fontface = "bold"))
  grid::grid.text("Force-directed layout (Fruchterman-Reingold)",
                  y = 0.86, gp = grid::gpar(fontsize = 12, fontface = "italic"))
  grid::grid.text(sprintf("Clusters analyzed: %d", length(network_results)),
                  y = 0.80, gp = grid::gpar(fontsize = 11))

  for (clust_id in names(network_results)) {
    res <- network_results[[clust_id]]

    if (is.null(res$filtered) || nrow(res$filtered) == 0) {
      grid::grid.newpage()
      grid::grid.text(paste0("Cluster ", clust_id, ": No edges after filtering"),
                      y = 0.5, gp = grid::gpar(fontsize = 14))
      next
    }

    edges_df <- as.data.frame(res$filtered[, c("source", "target", weight_col), with = FALSE])
    colnames(edges_df) <- c("source", "target", "weight")

    # Normalize to [0.1, 1] to ensure all weights are positive for FR layout
    w_min <- min(edges_df$weight, na.rm = TRUE)
    w_max <- max(edges_df$weight, na.rm = TRUE)
    if (w_min == w_max) {
      edges_df$weight_norm <- rep(0.5, nrow(edges_df))
    } else {
      edges_df$weight_norm <- 0.1 + 0.9 * (edges_df$weight - w_min) / (w_max - w_min)
    }
    edges_df$weight_norm <- pmax(edges_df$weight_norm, 0.01)  # Ensure > 0

    g <- igraph::graph_from_data_frame(edges_df[, c("source", "target")],
                                       directed = directed)
    # Assign normalized weights (FR layout needs positive weights)
    igraph::E(g)$weight <- edges_df$weight_norm
    igraph::E(g)$original_weight <- edges_df$weight

    node_degree <- igraph::degree(g)
    igraph::V(g)$size <- 4 + (node_degree / max(node_degree)) * 12
    igraph::V(g)$color <- colorRampPalette(c("#ffffcc", "#ff7f00"))(100)[
      ceiling(node_degree / max(node_degree) * 100)
    ]

    set.seed(123)
    layout <- igraph::layout_with_fr(g, niter = 500, dim = 2, weights = igraph::E(g)$weight)

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(x = 0.05, y = 0.05, width = 0.9, height = 0.85, just = c("left", "bottom")))

    plot(g,
         layout = layout,
         edge.width = edges_df$weight_norm * 3,
         edge.color = edge_color,
         edge.arrow.size = ifelse(directed, 0.3, 0),
         edge.curved = 0.2,
         vertex.label.cex = 0.8,
         vertex.label.dist = 1.5,
         asp = 0.8,
         margin = 0.1)

    grid::popViewport()

    grid::pushViewport(grid::viewport(x = 0.05, y = 0.88, width = 0.9, height = 0.1, just = c("left", "bottom")))
    grid::grid.text(paste0("Cluster ", clust_id),
                    x = 0.05, y = 0.8, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 14, fontface = "bold"))
    grid::grid.text(sprintf("Genes: %d | Edges: %d | Layout: Fruchterman-Reingold",
                            length(unique(c(res$filtered$source, res$filtered$target))),
                            nrow(res$filtered)),
                    x = 0.05, y = 0.3, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 10))
    grid::popViewport()
  }

  dev.off()
  message("✓ Network visualization PDF saved: ", pdf_path)
  invisible(pdf_path)
}


# =============================================================================
# generate_cluster_profile_report
# =============================================================================
# Generate cluster profile: expression heatmap, box plots, GO enrichment, stats.
generate_cluster_profile_report <- function(cluster_assignments,
                                            pseudobulk_matrix,
                                            output_dir,
                                            method_name = "MIXED") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  clusters_to_plot <- unique(cluster_assignments$cluster[cluster_assignments$cluster != "grey"])
  pdf_path <- file.path(output_dir, paste0(method_name, "_cluster_profiles.pdf"))

  pdf(pdf_path, width = 14, height = 11)

  grid::grid.newpage()
  grid::grid.text("Gene Cluster Profiles",
                  y = 0.92, gp = grid::gpar(fontsize = 22, fontface = "bold"))
  grid::grid.text(sprintf("Method: %s | Clusters: %d", method_name, length(clusters_to_plot)),
                  y = 0.86, gp = grid::gpar(fontsize = 12, fontface = "italic"))

  for (clust_id in clusters_to_plot) {
    genes_in_cluster <- cluster_assignments$gene_id[cluster_assignments$cluster == clust_id]
    genes_in_cluster <- intersect(genes_in_cluster, rownames(pseudobulk_matrix))

    if (length(genes_in_cluster) == 0) next

    expr_matrix <- pseudobulk_matrix[genes_in_cluster, , drop = FALSE]

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(x = 0.05, y = 0.35, width = 0.9, height = 0.6, just = c("left", "bottom")))

    expr_scaled <- t(scale(t(expr_matrix)))
    expr_scaled[is.na(expr_scaled)] <- 0
    expr_scaled <- pmin(pmax(expr_scaled, -3), 3)

    h <- ComplexHeatmap::Heatmap(
      expr_scaled,
      name = "log2(CPM)\n(scaled)",
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_names = length(genes_in_cluster) <= 50,
      show_column_names = TRUE,
      column_title = paste0("Cluster ", clust_id, " - Expression Heatmap"),
      col = circlize::colorRamp2(c(-3, 0, 3), c("#3182bd", "#fff7fb", "#e6550d")),
      width = grid::unit(10, "cm"),
      height = grid::unit(8, "cm")
    )

    ComplexHeatmap::draw(h, newpage = FALSE)
    grid::popViewport()

    grid::pushViewport(grid::viewport(x = 0.05, y = 0.05, width = 0.9, height = 0.25, just = c("left", "bottom")))
    stats_text <- sprintf(
      "Cluster: %s\nGenes: %d | Samples: %d\nExpression: mean(log2 CPM) = %.2f ± %.2f",
      clust_id, nrow(expr_matrix), ncol(expr_matrix),
      mean(expr_matrix), sd(as.vector(expr_matrix))
    )
    grid::grid.text(stats_text, x = 0.05, y = 0.9, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 10, fontfamily = "mono"))
    grid::popViewport()
  }

  dev.off()
  message("✓ Cluster profile PDF saved: ", pdf_path)

  summary_list <- list()
  for (clust_id in clusters_to_plot) {
    genes <- cluster_assignments$gene_id[cluster_assignments$cluster == clust_id]
    genes <- intersect(genes, rownames(pseudobulk_matrix))

    if (length(genes) == 0) next

    expr <- pseudobulk_matrix[genes, , drop = FALSE]
    summary_list[[clust_id]] <- data.frame(
      cluster = clust_id,
      n_genes = nrow(expr),
      mean_expr = round(mean(expr), 3),
      sd_expr = round(sd(as.vector(expr)), 3),
      max_expr = round(max(expr), 3),
      min_expr = round(min(expr), 3)
    )
  }

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  write.table(summary_df, file.path(output_dir, "cluster_summary_stats.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message("✓ Cluster profile report saved to: ", output_dir)
  invisible(list(pdf = pdf_path, stats = summary_df))
}


# =============================================================================
# run_network_inference_pipeline
# =============================================================================
# Encapsulates network inference (GENIE3, WGCNA, SYNERGY) with method selection.
# Returns all results in a named list for downstream use.
#
# Parameters:
#   heatmap_results      : cluster assignments from Section 20
#   pseudobulk_dir       : directory with pseudobulk replicas (Section 9)
#   output_base_dir      : base directory for results (dir_08/<contrast>)
#   methods              : vector of methods to run ("GENIE3", "WGCNA", "SYNERGY")
#   orgdb, keytype, custom_tfs : for GENIE3/SYNERGY
#   cor_min, genie3_ntrees, n_cores : GENIE3/SYNERGY parameters
#   soft_power, network_type, tom_threshold : WGCNA parameters
#   n_top_clusters, min_var_filter : common parameters
#
# Returns: list with $results and $pdfs (named by method)
#
run_network_inference_pipeline <- function(heatmap_results,
                                           pseudobulk_dir,
                                           output_base_dir,
                                           methods = c("GENIE3", "WGCNA", "SYNERGY"),
                                           orgdb = org.At.tair.db,
                                           keytype = "TAIR",
                                           custom_tfs = NULL,
                                           cor_min = 0.90,
                                           genie3_ntrees = 100,
                                           n_cores = 4,
                                           soft_power = 6,
                                           network_type = "signed",
                                           tom_threshold = 0.15,
                                           n_top_clusters = 3,
                                           min_var_filter = 0.01) {

  message("\n═══════════════════════════════════════════════════════════════════")
  message("NETWORK INFERENCE PIPELINE")
  message("Methods: ", paste(methods, collapse = ", "))
  message("═══════════════════════════════════════════════════════════════════\n")

  results <- list()
  pdf_paths <- list()

  # ── Run GENIE3 ──────────────────────────────────────────────────────────────
  if ("GENIE3" %in% methods) {
    message("\n▶ Running GENIE3...")
    results$GENIE3 <- run_genie3_per_cluster(
      cluster_assignments = heatmap_results,
      pseudobulk_dir      = pseudobulk_dir,
      output_dir          = file.path(output_base_dir, "GENIE3"),
      orgdb               = orgdb,
      keytype             = keytype,
      custom_tfs          = custom_tfs,
      n_top_clusters      = n_top_clusters,
      cor_min             = cor_min,
      genie3_ntrees       = genie3_ntrees,
      n_cores             = n_cores,
      min_var_filter      = min_var_filter
    )
    pdf_paths$GENIE3 <- generate_network_pdf(
      results[[1]],
      file.path(output_base_dir, "GENIE3"),
      method_name = "GENIE3",
      weight_col = "weight",
      directed = TRUE,
      edge_color = "#2ca02c"
    )
  }

  # ── Run WGCNA ───────────────────────────────────────────────────────────────
  if ("WGCNA" %in% methods) {
    message("\n▶ Running WGCNA...")
    results$WGCNA <- run_wgcna_per_cluster(
      cluster_assignments = heatmap_results,
      pseudobulk_dir      = pseudobulk_dir,
      output_dir          = file.path(output_base_dir, "WGCNA"),
      n_top_clusters      = n_top_clusters,
      soft_power          = soft_power,
      network_type        = network_type,
      tom_threshold       = tom_threshold,
      min_var_filter      = min_var_filter
    )
    pdf_paths$WGCNA <- generate_network_pdf(
      results$WGCNA,
      file.path(output_base_dir, "WGCNA"),
      method_name = "WGCNA",
      weight_col = "TOM",
      directed = FALSE,
      edge_color = "#1f77b4"
    )
  }

  # ── Run SYNERGY ─────────────────────────────────────────────────────────────
  if ("SYNERGY" %in% methods) {
    message("\n▶ Running SYNERGY...")
    results$SYNERGY <- run_synergistic_network(
      cluster_assignments = heatmap_results,
      pseudobulk_dir      = pseudobulk_dir,
      output_dir          = file.path(output_base_dir, "SYNERGY"),
      orgdb               = orgdb,
      keytype             = keytype,
      custom_tfs          = custom_tfs,
      n_top_clusters      = n_top_clusters,
      soft_power          = soft_power,
      network_type        = network_type,
      genie3_ntrees       = genie3_ntrees,
      n_cores             = n_cores,
      min_var_filter      = min_var_filter,
      cor_min             = cor_min,
      tom_min             = tom_threshold
    )
    pdf_paths$SYNERGY <- generate_network_pdf(
      results$SYNERGY,
      file.path(output_base_dir, "SYNERGY"),
      method_name = "SYNERGY",
      weight_col = "score_synergy",
      directed = TRUE,
      edge_color = "#d62728"
    )
  }

  message("\n═══════════════════════════════════════════════════════════════════")
  message("✓ NETWORK INFERENCE COMPLETE")
  message("  Methods run: ", paste(names(results), collapse = ", "))
  message("═══════════════════════════════════════════════════════════════════\n")

  invisible(list(results = results, pdfs = pdf_paths))
}


# =============================================================================
# test_network_thresholds
# =============================================================================
# Test multiple threshold combinations and generate comparison report.
# Helps identify optimal thresholds for network inference.
#
# Generates a PDF with:
#   - Table of thresholds and edge counts per cluster
#   - Recommendation based on cluster coverage
#
test_network_thresholds <- function(heatmap_results,
                                     pseudobulk_dir,
                                     output_dir,
                                     method = "SYNERGY",
                                     orgdb = org.At.tair.db,
                                     keytype = "TAIR",
                                     custom_tfs = NULL,
                                     genie3_ntrees = 100,
                                     n_cores = 4,
                                     soft_power = 6,
                                     network_type = "signed",
                                     n_top_clusters = 3,
                                     min_var_filter = 0.01) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Define 5 threshold combinations (increasing permissiveness)
  thresholds <- list(
    list(name = "Strict (top 5%)",     cor = 0.90, tom = 0.15),
    list(name = "Moderate (top 10%)",  cor = 0.80, tom = 0.08),
    list(name = "Exploratory (top 15%)", cor = 0.75, tom = 0.05),
    list(name = "Permissive (top 20%)", cor = 0.70, tom = 0.03),
    list(name = "Very Lax (top 25%)",  cor = 0.65, tom = 0.01)
  )

  results_summary <- list()

  # Test each threshold
  for (t in thresholds) {
    message("\nTesting: ", t$name)

    if (method == "SYNERGY") {
      res <- run_synergistic_network(
        cluster_assignments = heatmap_results,
        pseudobulk_dir      = pseudobulk_dir,
        output_dir          = file.path(output_dir, "temp"),
        orgdb               = orgdb,
        keytype             = keytype,
        custom_tfs          = custom_tfs,
        n_top_clusters      = n_top_clusters,
        soft_power          = soft_power,
        network_type        = network_type,
        genie3_ntrees       = genie3_ntrees,
        n_cores             = n_cores,
        min_var_filter      = min_var_filter,
        cor_min             = t$cor,
        tom_min             = t$tom
      )
    } else if (method == "GENIE3") {
      res <- run_genie3_per_cluster(
        cluster_assignments = heatmap_results,
        pseudobulk_dir      = pseudobulk_dir,
        output_dir          = file.path(output_dir, "temp"),
        orgdb               = orgdb,
        keytype             = keytype,
        custom_tfs          = custom_tfs,
        n_top_clusters      = n_top_clusters,
        cor_min             = t$cor,
        genie3_ntrees       = genie3_ntrees,
        n_cores             = n_cores,
        min_var_filter      = min_var_filter
      )
    } else if (method == "WGCNA") {
      res <- run_wgcna_per_cluster(
        cluster_assignments = heatmap_results,
        pseudobulk_dir      = pseudobulk_dir,
        output_dir          = file.path(output_dir, "temp"),
        n_top_clusters      = n_top_clusters,
        soft_power          = soft_power,
        network_type        = network_type,
        tom_threshold       = t$tom,
        min_var_filter      = min_var_filter
      )
    }

    # Summarize results
    if (is.null(res) || !length(res)) {
      message("  No results for this threshold")
      next
    }

    edge_counts <- sapply(res, function(x) {
      if (is.null(x$filtered)) return(0)
      nrow(x$filtered)
    })
    edge_counts <- as.numeric(edge_counts)  # Ensure numeric
    clusters_with_edges <- sum(edge_counts > 0, na.rm = TRUE)

    results_summary[[t$name]] <- list(
      threshold = t,
      n_clusters_with_edges = clusters_with_edges,
      total_edges = sum(edge_counts, na.rm = TRUE),
      edge_counts = edge_counts
    )

    message(sprintf("  Clusters with edges: %d/%d | Total edges: %d",
                    clusters_with_edges, length(res), sum(edge_counts, na.rm = TRUE)))
  }

  # Generate comparison PDF
  pdf_path <- file.path(output_dir, paste0(method, "_threshold_comparison.pdf"))
  pdf(pdf_path, width = 14, height = 11)

  # Title page
  grid::grid.newpage()
  grid::grid.text("Network Threshold Comparison",
                  y = 0.92, gp = grid::gpar(fontsize = 22, fontface = "bold"))
  grid::grid.text(paste0("Method: ", method, " | Testing 5 threshold combinations"),
                  y = 0.86, gp = grid::gpar(fontsize = 12, fontface = "italic"))

  # Summary table
  grid::grid.newpage()
  grid::grid.text("Threshold Summary Table",
                  y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))

  # Summary table as text instead of grid.table (simpler, more reliable)
  table_text <- sprintf(
    "Threshold Summary:\n\n%s",
    paste(sapply(names(results_summary), function(name) {
      r <- results_summary[[name]]
      sprintf(
        "%-30s | Cor=%.2f TOM=%.2f | Clusters=%d | Edges=%d",
        name, r$threshold$cor, r$threshold$tom,
        r$n_clusters_with_edges, r$total_edges
      )
    }), collapse = "\n")
  )

  grid::grid.text(table_text, x = 0.1, y = 0.8, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 10, fontfamily = "mono"))

  # Per-cluster detail pages
  for (name in names(results_summary)) {
    grid::grid.newpage()
    r <- results_summary[[name]]

    grid::grid.text(name, y = 0.95, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    grid::grid.text(sprintf("Pearson: %.2f | TOM: %.2f", r$threshold$cor, r$threshold$tom),
                    y = 0.90, gp = grid::gpar(fontsize = 10, fontface = "italic"))

    # Edge count per cluster
    cluster_names <- names(r$edge_counts)
    cluster_text <- sprintf("%s: %d edges\n", cluster_names, r$edge_counts)

    detail_text <- sprintf(
      "Total clusters tested: %d\nClusters with edges: %d\nTotal edges: %d\n\nPer-cluster breakdown:\n%s",
      length(r$edge_counts),
      r$n_clusters_with_edges,
      r$total_edges,
      paste(cluster_text, collapse = "")
    )

    grid::grid.text(detail_text, x = 0.1, y = 0.75, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 11, fontfamily = "mono"))
  }

  # Recommendation page
  grid::grid.newpage()
  grid::grid.text("Recommendation",
                  y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))

  # Guard against empty results
  if (!length(results_summary)) {
    rec_text <- "No results found for any threshold combination.\nPlease check the parameters and data."
  } else {
    best_idx <- which.max(sapply(results_summary, function(x) x$n_clusters_with_edges))
    best_name <- names(results_summary)[best_idx]
    best <- results_summary[[best_name]]

    rec_text <- sprintf(
      "RECOMMENDED THRESHOLD: %s\n\nReason:\n• Covers %d clusters (most coverage)\n• Total edges: %d\n• Balance: not too strict, not too lax\n\nYou can adjust based on your needs:\n- Need fewer but high-confidence edges → use stricter\n- Need more exploratory edges → use more permissive",
      best_name,
      best$n_clusters_with_edges,
      best$total_edges
    )
  }

  grid::grid.text(rec_text, x = 0.1, y = 0.8, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 12, fontfamily = "mono"))

  dev.off()

  message("\n✓ Threshold comparison PDF saved: ", pdf_path)
  if (length(results_summary)) {
    message("✓ Recommendation: ", best_name)
  } else {
    message("⚠ No results found for threshold testing")
    best_name <- NA
  }

  invisible(list(pdf = pdf_path, summary = results_summary, recommendation = best_name))
}
