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

  obj    <- JoinLayers(obj)
  pseudo <- AggregateExpression(obj,
                                group.by     = "replicate",
                                assays       = "RNA",
                                return.seurat = FALSE,
                                layer        = "counts")

  counts          <- as.data.frame(pseudo$RNA)
  colnames(counts) <- sub("^g", "", colnames(counts))

  counts[, sort(colnames(counts))]
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
  condition <- gsub("-rep[0-9]+$", "", sub("^g", "", rep_names))

  if (length(unique(condition)) < 2) return(invisible(NULL))

  # Auto-detect condition levels from the data
  cond_levels <- unique(condition)

  colData <- data.frame(
    row.names = rep_names,
    condition = factor(condition, levels = cond_levels)
  )

  dds       <- DESeqDataSetFromMatrix(countData = counts_mat,
                                      colData   = colData,
                                      design    = ~ condition)
  dds       <- DESeq(dds)
  available <- levels(colData$condition)[levels(colData$condition) %in% unique(condition)]

  for (comp in comparaciones) {
    conds <- comp$conds
    tag   <- comp$tag
    if (!all(conds %in% available)) next
    res    <- results(dds, contrast = c("condition", conds[2], conds[1]))
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

  dat     <- na.omit(do.call(rbind, bloques))
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


# =============================================================================
# END OF FUNCTIONS
# =============================================================================
