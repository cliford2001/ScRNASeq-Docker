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
#     - resumen_nFeature_plot
#
#  2. PREPROCESSING AND DOUBLET DETECTION
#     - preprocesar_y_doubletfinder
#     - doubletfinder_pipeline
#     - load_sample          (load + annotate only, no filtering)
#     - filter_sample        (filter + DoubletFinder on annotated object)
#     - process_sample       (shortcut: load_sample + filter_sample)
#
#  3. BULK / PSEUDOBULK UTILITIES
#     - normalizar_bulk_pseudobulk
#     - clasificar_residuos
#     - generate_pseudobulk
#     - plot_replicate_correlation
#
#  4. SEURAT UTILITIES
#     - unificar_nombres
#     - mostrar_tabla
#     - exportar_para_scanpy
#     - safe_vln
#     - unir_layers_counts
#
#  5. ANNOTATION
#     - find_markers
#     - annotate_by_markers
#     - annotate_by_reference
#     - subclustar_tipo
#
#  6. PSEUDOBULK, DESEQ2, VOLCANO, HEATMAP
#     - asignar_pseudoreplicados
#     - hacer_pseudobulk
#     - correr_deseq2
#     - hacer_volcano
#     - procesar_deseq2_resultado
#     - hacer_heatmap
#     - hacer_dotplot_marcadores
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
#' @param etiquetas Labels for each object.
#' @param colores   Color vector (named or positional).
#' @export
resumen_nFeature_plot <- function(obj_list, etiquetas = NULL, colores = NULL) {

  if (is.null(etiquetas)) etiquetas <- paste0("Group", seq_along(obj_list))
  if (length(etiquetas) != length(obj_list)) stop("Labels must match objects.")

  if (is.null(colores)) {
    colores       <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")[1:length(obj_list)]
    names(colores) <- etiquetas
  }

  lista_df <- lapply(seq_along(obj_list), function(i) {
    obj <- obj_list[[i]]
    data.frame(nFeature_RNA = obj@meta.data$nFeature_RNA, grupo = etiquetas[i])
  })

  meta_comb       <- bind_rows(lista_df)
  meta_comb$grupo <- factor(meta_comb$grupo, levels = etiquetas)

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
    arrange(factor(grupo, levels = etiquetas))

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
    arrange(factor(grupo, levels = etiquetas))

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
normalizar_bulk_pseudobulk <- function(pseudobulk_counts, bulk_counts) {

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
  cat("Generando pseudobulk para:", paste(groups, collapse = ", "), "\n")

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
#'   hacer_pseudobulk().
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
exportar_para_scanpy <- function(seurat_obj,
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
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    message("Writing h5ad with zellkonverter...")
    zellkonverter::writeH5AD(sce, file = outfile, X_name = X_name)
    ok <- TRUE
  } else if (requireNamespace("SeuratDisk", quietly = TRUE)) {
    message("Using SeuratDisk (zellkonverter not available)...")
    tmp_h5seu <- file.path(tempdir(), paste0(basename(outfile), ".h5seurat"))
    if (file.exists(tmp_h5seu)) file.remove(tmp_h5seu)
    SeuratDisk::SaveH5Seurat(seurat_obj, filename = tmp_h5seu, overwrite = TRUE)
    SeuratDisk::Convert(source = tmp_h5seu, dest = "h5ad", overwrite = TRUE)
    gen_h5ad <- sub("\\.h5seurat$", ".h5ad", tmp_h5seu)
    if (!file.exists(gen_h5ad)) stop("Failed to generate h5ad.")
    file.rename(gen_h5ad, outfile)
    ok <- TRUE
  } else {
    stop("Install 'zellkonverter' or 'SeuratDisk' to export h5ad.")
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
    cat("Cargando marcadores existentes:", output_file, "\n")
    markers <- read.table(output_file, header = TRUE, sep = "\t", quote = "")
  } else {
    cat("Calculando marcadores...\n")
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
    reference_file <- file.choose(caption = "Selecciona archivo de referencia (gene | cell.types)")
  }

  cat("Usando referencia:", reference_file, "\n")

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
    ref_file      <- file.choose(caption = "Selecciona objeto Seurat de referencia (.rds)")
    cat("Cargando referencia:", ref_file, "\n")
    reference_obj <- readRDS(ref_file)
  }

  if (is.null(reference_col)) {
    cat("\nColumnas disponibles en referencia:\n")
    cols <- colnames(reference_obj@meta.data)
    for (i in seq_along(cols)) {
      cat(" ", i, "->", cols[i], "\n")
    }
    selection     <- as.integer(readline("Selecciona numero de columna: "))
    reference_col <- cols[selection]
  }

  cat("Usando columna:", reference_col, "\n")

  anchors <- FindTransferAnchors(
    reference = reference_obj,
    query     = seurat_obj,
    dims      = dims
  )

  predictions <- TransferData(
    anchorset = anchors,
    refdata   = reference_obj@meta.data[[reference_col]],
    dims      = dims
  )

  seurat_obj$celltype_reference <- predictions$predicted.id

  cat("\nAnotacion por referencia:\n")
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
subclustar_tipo <- function(obj, tipo, annot_col = "annotation_agrupada",
                            resolution = 0.3, dims = 1:20) {

  sub <- subset(obj, cells = colnames(obj)[obj@meta.data[[annot_col]] %in% tipo])
  sub <- sub %>%
    RunPCA() %>%
    RunUMAP(dims = dims) %>%
    FindNeighbors(dims = dims) %>%
    FindClusters(resolution = resolution)

  sub$cluster_subtipo <- as.character(sub$seurat_clusters)

  return(sub)
}


# =============================================================================
# 6. PSEUDOBULK, DESEQ2, VOLCANO, HEATMAP
# =============================================================================

#' Assign Pseudo-replicates
#'
#' Randomly assigns cells within each condition to pseudo-replicate groups.
#' Conditions are auto-detected from orig.ident_uni unless explicitly provided.
#'
#' @param obj         Seurat object with orig.ident_uni metadata.
#' @param condiciones Character vector of condition names to include. NULL
#'   (default) uses all conditions present in the data.
#' @param n_reps      Number of pseudo-replicates per condition.
#' @param seed        Random seed for reproducibility.
#' @return Seurat object with a replicate metadata column, or NULL if fewer
#'   than 2 conditions are present.
#' @export
asignar_pseudoreplicados <- function(obj,
                                     condiciones = NULL,
                                     n_reps      = 3,
                                     seed        = 1807) {

  set.seed(seed)

  # Auto-detect conditions from data if not provided
  all_conds <- unique(obj$orig.ident_uni)
  condiciones_presentes <- if (!is.null(condiciones)) intersect(all_conds, condiciones) else all_conds

  if (length(condiciones_presentes) < 2) return(NULL)

  obj$replicate <- NA
  for (cond in condiciones_presentes) {
    idx              <- obj$orig.ident_uni == cond
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
hacer_pseudobulk <- function(obj) {

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
      hacer_pseudobulk(obj),
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
hacer_volcano <- function(file, padj_cut = 0.05, lfc_cut = 1) {

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
    p <- hacer_volcano(file, padj_cut = padj_cut, lfc_cut = lfc_cut) +
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
hacer_heatmap <- function(matriz,
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
hacer_dotplot_marcadores <- function(seurat_obj,
                                     marks,
                                     annot_col       = "celltype_reference_curated",
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
                                nrow(Mz), min_genes),
    silent            = TRUE
  )
  ggsave(file.path(output_dir, heatmap_pdf), heatmap_obj$gtable, width = 10, height = 20, dpi = 300)

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
  tom_plot <- pheatmap(
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
    main              = sprintf("TOM Heatmap (%d genes)", nrow(Mz)),
    silent            = TRUE
  )
  ggsave(file.path(output_dir, tom_pdf), tom_plot$gtable, width = 10, height = 10, dpi = 300)

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
                                nrow(Mz), min_genes),
    silent            = TRUE
  )
  ggsave(file.path(output_dir, heatmap_pdf), heatmap_obj$gtable, width = 10, height = 20, dpi = 300)

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
    cluster_assignments  = cluster_assignments,
    heatmap              = heatmap_obj
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
  tom_plot <- pheatmap(
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
    main              = sprintf("TOM Heatmap (%d genes)", nrow(Mz)),
    silent            = TRUE
  )
  ggsave(file.path(output_dir, tom_pdf), tom_plot$gtable, width = 10, height = 10, dpi = 300)

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
#' list of paths. Call list2env() on the result to assign dir_00…dir_13 in
#' the global environment.
#'
#' @param base_dir Root results directory (e.g. file.path(DATA_DIR, "results")).
#' @return Named list with dir_00 through dir_13.
#' @export
create_pipeline_dirs <- function(base_dir) {
  mkd <- function(name) {
    d <- file.path(base_dir, name)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    d
  }
  list(
    dir_00 = mkd("00_workflow"),
    dir_01 = mkd("01_qc"),
    dir_02 = mkd("02_filtering"),
    dir_03 = mkd("03_integration"),
    dir_04 = mkd("04_clustering"),
    dir_05 = mkd("05_annotation"),
    dir_06 = mkd("06_expression"),
    dir_07 = mkd("07_curation"),
    dir_08 = mkd("08_export"),
    dir_09 = mkd("09_pseudobulk"),
    dir_10 = mkd("10_deseq2"),
    dir_11 = mkd("11_volcano"),
    dir_12 = mkd("12_heatmaps"),
    dir_13 = mkd("13_go")
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
                                     annot_col      = "celltype_reference_curated",
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
  celular_subsets <- setNames(
    lapply(cell_types, function(t)
      subset(obj, cells = colnames(obj)[obj@meta.data[[annot_col]] == t])),
    gsub("[^[:alnum:]_]", "_", cell_types)
  )

  celular_subsets_replicados <- Filter(Negate(is.null),
    lapply(celular_subsets, asignar_pseudoreplicados, n_reps = n_reps))
  pseudobulk_list <- lapply(celular_subsets_replicados, hacer_pseudobulk)

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
      plots <- c(plots, list(hacer_volcano(f, padj_cut = padj_cut, lfc_cut = lfc_cut)))
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
      tryCatch(hacer_heatmap(matriz), error = function(e) message("Heatmap error: ", e$message))
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

  top <- data.frame(
    x     = 1:10,
    y     = rep(2, 10),
    label = c("FASTQ\nFiles", "CellRanger\nCount", "CellBender\n(optional)",
              "Load &\nQC", "Filter +\nDoubletFinder",
              "Merge &\nNormalize", "Harmony\nIntegration",
              "Clustree +\nElbow plot", "Final\nClustering", "Annotate\n+ Export"),
    group = c(rep("Pre-processing", 3), rep("Chapter 1", 7)),
    stringsAsFactors = FALSE
  )

  bot <- data.frame(
    x     = 9:5,
    y     = rep(0, 5),
    label = c("Cell-type\nSubsets", "Pseudobulk\n+ Replicates",
              "DESeq2\nDE", "Volcano +\nHeatmap", "GO\nEnrichment"),
    group = rep("Chapter 2", 5),
    stringsAsFactors = FALSE
  )

  nodes <- rbind(top, bot)

  group_colors <- c(
    "Pre-processing" = "#d9d9d9",
    "Chapter 1"      = "#b3cde3",
    "Chapter 2"      = "#ccebc5"
  )

  edges_top <- data.frame(x1 = 1:9, y1 = rep(2, 9), x2 = 2:10, y2 = rep(2, 9))
  edge_down <- data.frame(x1 = 10,  y1 = 2,          x2 = 9,    y2 = 0)
  edges_bot <- data.frame(x1 = 9:6, y1 = rep(0, 4),  x2 = 8:5,  y2 = rep(0, 4))
  edges     <- rbind(edges_top, edge_down, edges_bot)

  w <- 0.85; h <- 0.50

  p <- ggplot() +
    annotate("rect", xmin = 0.4, xmax = 3.6, ymin = 1.4, ymax = 2.6,
             fill = "#f5f5f5", color = "grey75", linetype = "dashed", linewidth = 0.4) +
    annotate("text", x = 2,   y = 2.70,
             label = "Pre-processing  (bash)", size = 3, color = "grey55", fontface = "italic") +
    annotate("text", x = 6.5, y = 2.70,
             label = "CHAPTER 1 \u2014 Single-Cell Analysis",
             size = 3.5, color = "#1565c0", fontface = "bold") +
    annotate("text", x = 7,   y = -0.70,
             label = "CHAPTER 2 \u2014 Pseudobulk DE & GO Enrichment",
             size = 3.5, color = "#2e7d32", fontface = "bold") +
    geom_segment(data = edges,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 arrow     = arrow(length = unit(0.22, "cm"), type = "closed"),
                 linewidth = 0.55, color = "grey40", lineend = "round") +
    geom_rect(data = nodes,
              aes(xmin = x - w/2, xmax = x + w/2,
                  ymin = y - h/2, ymax = y + h/2, fill = group),
              color = "grey35", linewidth = 0.4) +
    scale_fill_manual(values = group_colors, name = NULL,
                      guide  = guide_legend(nrow = 1)) +
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              size = 2.7, lineheight = 0.9) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text      = element_text(size = 10),
          plot.margin      = margin(20, 10, 20, 10)) +
    coord_fixed(ratio = 1.5) +
    xlim(0.3, 10.7) + ylim(-1.1, 3.1)

  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  ggsave(outfile, p, width = 18, height = 7, dpi = 300)
  message("Workflow figure saved: ", outfile)
  invisible(p)
}


# =============================================================================
# END OF FUNCTIONS
# =============================================================================
