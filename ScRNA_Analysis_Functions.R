# =============================================================================
# ScRNA-Seq Single-Cell Analysis - Custom Functions Library
# =============================================================================
# Author: Ellie Malcolm et al.
# Date: 2025-03
# Description: Reusable functions for QC, preprocessing, clustering, and DE analysis
# =============================================================================

# ============================================================================
# 1. QC AND VISUALIZATION FUNCTIONS
# ============================================================================

#' Load CellBender Filtered HDF5 Data
#' 
#' Reads filtered expression matrix from CellBender HDF5 output
#' @param h5_path Path to filtered HDF5 file
#' @param project Project name for metadata
#' @return Seurat object with raw counts
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
  
  message("Creating sparse gene × cell matrix...")
  mat <- new("dgCMatrix",
             x = as.numeric(data),
             i = indices,
             p = indptr,
             Dim = shape,
             Dimnames = list(gene_ids, barcodes))
  
  seu <- CreateSeuratObject(counts = mat, project = project)
  return(seu)
}

#' QC Violin Plot Grid
#' 
#' Visualizes nFeature, nCount, percent.mt, percent.cp by condition
#' @param obj1 Seurat object
#' @param label Condition label
#' @param color Color for plotting
#' @return ggplot object
#' @export
plot_qc_violin_grid <- function(obj1, label, color) {
  n1 <- ncol(obj1)
  obj1$cond <- label
  
  p1 <- VlnPlot(obj1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"),
                pt.size = 0.1, ncol = 4, group.by = "cond", cols = color) +
    ggtitle(paste0(label, " (", n1, " cells)")) +
    theme_minimal(base_size = 12)
  
  return(p1)
}

#' Summary of nFeature_RNA Distribution
#' 
#' Creates boxplot + quantile/quintile tables
#' @param obj_list List of Seurat objects
#' @param etiquetas Labels for each object
#' @param colores Color vector
#' @export
resumen_nFeature_plot <- function(obj_list, etiquetas = NULL, colores = NULL) {
  if (is.null(etiquetas)) etiquetas <- paste0("Group", seq_along(obj_list))
  if (length(etiquetas) != length(obj_list)) stop("Labels must match objects.")
  
  if (is.null(colores)) {
    colores <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")[1:length(obj_list)]
    names(colores) <- etiquetas
  }
  
  lista_df <- lapply(seq_along(obj_list), function(i) {
    obj <- obj_list[[i]]
    data.frame(nFeature_RNA = obj@meta.data$nFeature_RNA, grupo = etiquetas[i])
  })
  
  meta_comb <- bind_rows(lista_df)
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
      Min = quantile(nFeature_RNA, 0),
      Q1 = quantile(nFeature_RNA, 0.25),
      Median = quantile(nFeature_RNA, 0.5),
      Q3 = quantile(nFeature_RNA, 0.75),
      Max = quantile(nFeature_RNA, 1),
      .groups = "drop"
    ) %>%
    arrange(factor(grupo, levels = etiquetas))
  
  quintiles <- meta_comb %>%
    group_by(grupo) %>%
    summarise(
      `0%`  = quantile(nFeature_RNA, 0.0),
      `20%` = quantile(nFeature_RNA, 0.2),
      `40%` = quantile(nFeature_RNA, 0.4),
      `60%` = quantile(nFeature_RNA, 0.6),
      `80%` = quantile(nFeature_RNA, 0.8),
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
    ncol = 1, rel_heights = c(0.15, 1, 0.15, 1)
  )
  
  final_plot <- plot_grid(p_box, panel_tablas, ncol = 2, rel_widths = c(1.5, 1))
  print(final_plot)
}

#' Preprocessing + DoubletFinder Pipeline
#' 
#' Normalize, scale, PCA, and run DoubletFinder
#' @param seurat_obj Seurat object
#' @param pcs PCs to use
#' @param expected_doublet_rate Expected doublet rate
#' @param project_id Project ID
#' @return Seurat object with doublet classification
#' @export
preprocesar_y_doubletfinder <- function(seurat_obj, pcs = 1:20, 
                                       expected_doublet_rate = 0.075, 
                                       project_id = "sample") {
  message("Normalizing...")
  seurat_obj <- NormalizeData(seurat_obj)
  message("Finding variable features...")
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  message("Scaling...")
  seurat_obj <- ScaleData(seurat_obj)
  message("PCA...")
  seurat_obj <- RunPCA(seurat_obj, npcs = max(pcs))
  
  message("Running DoubletFinder...")
  sweep.res <- paramSweep(seurat_obj, PCs = pcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  nExp <- round(expected_doublet_rate * ncol(seurat_obj))
  
  seurat_obj <- doubletFinder(
    seurat_obj, PCs = pcs, pN = 0.25, pK = best.pK,
    nExp = nExp, reuse.pANN = FALSE, sct = FALSE
  )
  return(seurat_obj)
}

#' Full DoubletFinder Pipeline with Clustering
#' 
#' Comprehensive doublet detection with best parameter selection
#' @param obj Seurat object
#' @param etiqueta Sample label
#' @param PCs PCs for analysis
#' @param resolution Leiden resolution
#' @param return_singlets If TRUE, return only singlets
#' @param sct Whether to use SCT normalization
#' @return Seurat object, optionally filtered
#' @export
doubletfinder_pipeline <- function(obj, etiqueta = "Sample", PCs = 1:20, 
                                   resolution = 0.5, return_singlets = TRUE, 
                                   sct = FALSE) {
  message("Processing: ", etiqueta)
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = max(PCs))
  obj <- FindNeighbors(obj, dims = PCs)
  obj <- FindClusters(obj, resolution = resolution)
  
  sweep.res <- paramSweep(obj, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
  sweep.stats$pN <- as.numeric(as.character(sweep.stats$pN))
  
  best_row <- sweep.stats[which.max(sweep.stats$BCreal), ]
  best.pK <- best_row$pK
  best.pN <- best_row$pN
  nExp <- round(best.pN * ncol(obj))
  
  message("Best pK: ", best.pK, ", pN: ", best.pN, ", nExp: ", nExp)
  
  obj <- doubletFinder(
    obj, PCs = PCs, pN = best.pN, pK = best.pK,
    nExp = nExp, reuse.pANN = NULL, sct = sct
  )
  
  for (col in colnames(obj@meta.data)) {
    if (is.data.frame(obj@meta.data[[col]])) {
      message("Fixing column: ", col)
      obj@meta.data[[col]] <- obj@meta.data[[col]][, 1]
    }
  }
  
  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
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

#' Normalize Pseudobulk vs Bulk Counts
#' 
#' DESeq2-based normalization and log2 transformation
#' @param pseudobulk_counts Named vector of pseudobulk counts
#' @param bulk_counts Named vector of bulk counts
#' @return Data frame with normalized values
#' @export
normalizar_bulk_pseudobulk <- function(pseudobulk_counts, bulk_counts) {
  common_genes <- intersect(names(pseudobulk_counts), names(bulk_counts))
  
  if (length(common_genes) < 10) {
    stop("Too few common genes.")
  }
  
  counts_matrix <- data.frame(
    pseudobulk = round(pseudobulk_counts[common_genes]),
    bulk = round(bulk_counts[common_genes])
  )
  rownames(counts_matrix) <- common_genes
  
  condition <- factor(c("pseudobulk", "bulk"))
  col_data <- data.frame(condition = condition)
  
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = col_data,
                                design = ~ condition)
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  
  log_norm_counts <- log2(norm_counts + 1)
  
  df <- data.frame(
    gene = rownames(log_norm_counts),
    pseudobulk = log_norm_counts[, "pseudobulk"],
    bulk = log_norm_counts[, "bulk"]
  )
  
  return(df)
}

#' Classify Genes by Residuals
#' 
#' Linear model fitting and classification
#' @param df Data frame with pseudobulk and bulk columns
#' @param umbral Residual threshold
#' @return Data frame with status classification
#' @export
clasificar_residuos <- function(df, umbral = 5) {
  modelo <- lm(bulk ~ pseudobulk, data = df)
  df$residuals <- resid(modelo)
  
  df$status <- case_when(
    df$residuals >  umbral ~ "Upregulated",
    df$residuals < -umbral ~ "Downregulated",
    TRUE                   ~ "Consistent"
  )
  
  return(df)
}

#' Unify Ident Names
#' 
#' Remove numeric suffixes (e.g., ".1", "_2") from cluster names
#' @param obj Seurat object
#' @return Seurat object with updated Idents
#' @export
unificar_nombres <- function(obj) {
  old_levels <- levels(obj)
  new_levels <- gsub("[._][0-9]+$", "", old_levels)
  new_ids <- setNames(new_levels, old_levels)
  obj <- RenameIdents(obj, new_ids)
  return(obj)
}

#' Display Annotation Table as Grid
#' 
#' Creates and displays cell type counts table
#' @param filtered_vec Filtered vector
#' @param cellbender_vec CellBender vector
#' @param titulo Table title
#' @export
mostrar_tabla <- function(filtered_vec, cellbender_vec, titulo = "Annotations") {
  t1 <- table(filtered_vec)
  t2 <- table(cellbender_vec)
  all_types <- union(names(t1), names(t2))

  df <- data.frame(
    celltype = all_types,
    filtered = as.integer(t1[all_types]),
    cellbender = as.integer(t2[all_types]),
    stringsAsFactors = FALSE
  )
  df[is.na(df)] <- 0

  total_row <- data.frame(
    celltype = "Total",
    filtered = sum(df$filtered),
    cellbender = sum(df$cellbender),
    stringsAsFactors = FALSE
  )
  df <- rbind(df, total_row)

  grid.newpage()
  grid.draw(tableGrob(df, rows = NULL, theme = ttheme_minimal()))
}

#' Export Seurat to Scanpy h5ad Format
#' 
#' Converts Seurat to SingleCellExperiment and exports to h5ad
#' @param seurat_obj Seurat object
#' @param outfile Output file path
#' @param assay_name Assay to export
#' @param use_reduc Reductions to include
#' @param X_name Name for .X in Scanpy
#' @param overwrite Overwrite existing file
#' @return Invisible SingleCellExperiment
#' @export
exportar_para_scanpy <- function(
  seurat_obj,
  outfile,
  assay_name = "RNA",
  use_reduc  = c("pca", "umap", "harmony"),
  X_name     = "logcounts",
  overwrite  = TRUE
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)

  message("Converting Seurat → SCE…")
  sce <- as.SingleCellExperiment(seurat_obj, assay = assay_name)

  if (is.null(rownames(sce))) stop("Object missing rownames.")
  rownames(sce) <- make.unique(rownames(sce))
  if (is.null(colnames(sce)) || any(!nzchar(colnames(sce)))) {
    colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
  }
  stopifnot(!anyDuplicated(colnames(sce)))

  message("Extracting counts and logcounts…")
  assay(sce, "counts")    <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "counts")
  assay(sce, "logcounts") <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = "data")

  if (ncol(colData(sce)) == 0) {
    colData(sce) <- S4Vectors::DataFrame(row.names = colnames(sce))
  }
  colData(sce) <- cbind(
    colData(sce),
    seurat_obj@meta.data[colnames(sce), , drop = FALSE]
  )

  message("Exporting reductions…")
  reds <- Seurat::Reductions(seurat_obj)
  for (red in use_reduc) {
    if (red %in% reds) {
      message("   ✓ Including: ", red)
      reducedDims(sce)[[toupper(red)]] <- seurat_obj@reductions[[red]]@cell.embeddings
    }
  }

  if (file.exists(outfile) && overwrite) file.remove(outfile)

  ok <- FALSE
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    message("Writing h5ad with zellkonverter…")
    zellkonverter::writeH5AD(sce, file = outfile, X_name = X_name)
    ok <- TRUE
  } else if (requireNamespace("SeuratDisk", quietly = TRUE)) {
    message("Using SeuratDisk (zellkonverter not available)…")
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
#' Violin plot with custom colors
#' @param obj Seurat object
#' @param feature Gene to plot
#' @param colors Color palette
#' @return ggplot object
#' @export
safe_vln <- function(obj, feature, colors) {
  p <- VlnPlot(
    obj,
    features = feature,
    group.by = "orig.ident",
    pt.size = 0
  )
  p <- p + scale_fill_manual(values = colors)
  return(p)
}

#' Unify and Merge Seurat Layers
#' 
#' Combines multiple sparse count matrices
#' @param obj Seurat object
#' @param capas Layers to merge
#' @return Merged sparse matrix
#' @export
unir_layers_counts <- function(obj, capas) {
  if (length(capas) == 1) {
    return(GetAssayData(obj[["RNA"]], layer = capas))
  }
  
  message("Merging ", length(capas), " layers...")
  mats <- lapply(capas, function(x) GetAssayData(obj[["RNA"]], layer = x))
  merged <- Reduce(RowMergeSparseMatrices, mats)
  return(merged)
}

process_sample <- function(sample_info,
                           # Filtros nFeature_RNA
                           min_features = 0,
                           max_features = Inf,
                           # Filtros nCount_RNA
                           min_counts = 0,
                           max_counts = Inf,
                           # Filtros porcentajes
                           max_mt = 100,
                           max_cp = 100) {
  
  obj <- load_cellbender_filtered_h5(sample_info$file, sample_info$label)
  
  # Calcular porcentajes
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^ATMG")
  obj[["percent.cp"]] <- PercentageFeatureSet(obj, pattern = "^ATCG")
  
  # Filtrar
  obj <- subset(obj, subset = 
                  nFeature_RNA > min_features &
                  nFeature_RNA < max_features &
                  nCount_RNA > min_counts &
                  nCount_RNA < max_counts &
                  percent.mt < max_mt &
                  percent.cp < max_cp
  )
  
  # Doublets
  obj <- doubletfinder_pipeline(obj, etiqueta = sample_info$label)
  
  # Renombrar células y agregar condition
  obj <- RenameCells(obj, add.cell.id = sample_info$condition)
  obj$condition <- sample_info$condition
  
  return(obj)
}

# ============================================================
# PSEUDOBULK AUTOMÁTICO
# ============================================================

generate_pseudobulk <- function(seurat_obj, 
                                 group_by = "orig.ident",
                                 merge_replicates = TRUE) {
  
  # Obtener grupos únicos
  groups <- unique(seurat_obj@meta.data[[group_by]])
  cat("Generando pseudobulk para:", paste(groups, collapse = ", "), "\n")
  
  # Función para procesar un grupo

process_group <- function(group_name) {
    cells <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[group_by]] == group_name])
    
    # Detectar y unir layers de counts
    layers <- grep("^counts", Layers(cells[["RNA"]]), value = TRUE)
    
    if (length(layers) == 0) {
      counts <- GetAssayData(cells[["RNA"]], layer = "data")
    } else if (length(layers) == 1) {
      counts <- GetAssayData(cells[["RNA"]], layer = layers)
    } else {
      mats <- lapply(layers, function(x) GetAssayData(cells[["RNA"]], layer = x))
      counts <- Reduce(RowMergeSparseMatrices, mats)
    }
    
    # Sumar counts por gen
    gene_sums <- Matrix::rowSums(counts)
    cat(" ", group_name, "→", ncol(cells), "células,", length(gene_sums), "genes\n")
    return(gene_sums)
  }
  
  # Generar pseudobulk por grupo
  pseudobulk_list <- lapply(groups, process_group)
  names(pseudobulk_list) <- groups
  
  # Crear matriz unificada
  all_genes <- unique(unlist(lapply(pseudobulk_list, names)))
  
  pseudobulk_matrix <- sapply(pseudobulk_list, function(x) {
    v <- x[all_genes]
    v[is.na(v)] <- 0
    return(v)
  })
  rownames(pseudobulk_matrix) <- all_genes
  
  # Opción: merge réplicas (R1 + R2)
  if (merge_replicates) {
    conditions <- unique(seurat_obj$condition)
    
    merged_matrix <- sapply(conditions, function(cond) {
      # Encontrar columnas que pertenecen a esta condición
      cols <- grep(paste0("^", cond), colnames(pseudobulk_matrix), value = TRUE)
      if (length(cols) == 0) {
        cols <- colnames(pseudobulk_matrix)[grepl(cond, colnames(pseudobulk_matrix))]
      }
      if (length(cols) == 1) {
        return(pseudobulk_matrix[, cols])
      }
      return(rowSums(pseudobulk_matrix[, cols, drop = FALSE]))
    })
    colnames(merged_matrix) <- conditions
    
    cat("\nRéplicas fusionadas:\n")
    print(colnames(merged_matrix))
    
    return(list(
      by_sample = pseudobulk_matrix,
      by_condition = merged_matrix
    ))
  }
  
  return(pseudobulk_matrix)
}



# ============================================================
# PASO 1: Encontrar marcadores en tus datos
# ============================================================

find_markers <- function(seurat_obj, 
                         output_file = "results/FindAllMarkers.tsv",
                         only_pos = TRUE,
                         min_pct = 0.25,
                         logfc_threshold = 0.25,
                         force = FALSE) {
  
  seurat_obj <- JoinLayers(seurat_obj)
  Idents(seurat_obj) <- "seurat_clusters"
  
  if (file.exists(output_file) && !force) {
    cat("Cargando marcadores existentes:", output_file, "\n")
    markers <- read.table(output_file, header = TRUE, sep = "\t", quote = "")
  } else {
    cat("Calculando marcadores...\n")
    markers <- FindAllMarkers(
      seurat_obj, 
      only.pos = only_pos, 
      min.pct = min_pct, 
      logfc.threshold = logfc_threshold
    )
    write.table(markers, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
    cat("Guardado en:", output_file, "\n")
  }
  
  return(markers)
}

# ============================================================
# PASO 2: Anotar con tu lista de referencia
# ============================================================

annotate_by_markers <- function(seurat_obj,
                                markers,
                                reference_file = NULL) {
  
  # Interactivo si no se especifica
  if (is.null(reference_file)) {
    reference_file <- file.choose(caption = "Selecciona archivo de referencia (gene | cell.types)")
  }
  
  cat("Usando referencia:", reference_file, "\n")
  
  # Cargar tu lista de marcadores conocidos
  reference <- read.table(reference_file, header = TRUE, sep = "\t", quote = "")
  
  # Cruzar: tus marcadores vs referencia
  merged <- merge(markers, reference, by.x = "gene", by.y = "gene")
  merged <- merged[order(merged$cluster, merged$p_val_adj), ]
  merged <- merged[!duplicated(merged$cluster), ]
  
  cat("\nCoincidencias encontradas:\n")
  print(merged[, c("cluster", "gene", "cell.types")])
  
  # Asignar nombres a clusters
  Idents(seurat_obj) <- "seurat_clusters"
  new_ids <- merged$cell.types
  names(new_ids) <- merged$cluster
  seurat_obj <- RenameIdents(seurat_obj, new_ids)
  seurat_obj$celltype <- Idents(seurat_obj)
  
  cat("\nAnotación final:\n")
  print(table(seurat_obj$celltype))
  
  return(seurat_obj)
}


# --- Método 2: Por transferencia de referencia ---

annotate_by_reference <- function(seurat_obj,
                                  reference_obj = NULL,
                                  reference_col = NULL,
                                  dims = 1:30) {
  
  # Interactivo: cargar objeto de referencia
  if (is.null(reference_obj)) {
    ref_file <- file.choose(caption = "Selecciona objeto Seurat de referencia (.rds)")
    cat("Cargando referencia:", ref_file, "\n")
    reference_obj <- readRDS(ref_file)
  }
  
  # Interactivo: seleccionar columna de anotación
  if (is.null(reference_col)) {
    cat("\nColumnas disponibles en referencia:\n")
    cols <- colnames(reference_obj@meta.data)
    for (i in seq_along(cols)) {
      cat(" ", i, "→", cols[i], "\n")
    }
    selection <- as.integer(readline("Selecciona número de columna: "))
    reference_col <- cols[selection]
  }
  
  cat("Usando columna:", reference_col, "\n")
  
  # Transfer learning
  anchors <- FindTransferAnchors(
    reference = reference_obj, 
    query = seurat_obj, 
    dims = dims
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_obj@meta.data[[reference_col]],
    dims = dims
  )
  
  seurat_obj$celltype_reference <- predictions$predicted.id
  
  cat("\nAnotación por referencia:\n")
  print(table(seurat_obj$celltype_reference))
  
  return(seurat_obj)
}

# --- PASO 2: Función para subclustar un tipo celular ---
subclustar_tipo <- function(obj, tipo, resolution = 0.3, dims = 1:20) {
  sub <- subset(obj, subset = annotation_agrupada %in% tipo)
  sub <- sub %>%
    RunPCA() %>%
    RunUMAP(dims = dims) %>%
    FindNeighbors(dims = dims) %>%
    FindClusters(resolution = resolution)
  sub$cluster_subtipo <- as.character(sub$seurat_clusters)
  return(sub)
}
# ============================================================
# FUNCIONES: PSEUDOBULK, DESEQ2, VOLCANO, HEATMAP
# ============================================================

asignar_pseudoreplicados <- function(obj, condiciones = c("0N", "0.5N", "5N"), n_reps = 3, seed = 1807) {
  set.seed(seed)
  condiciones_presentes <- intersect(unique(obj$orig.ident_uni), condiciones)
  if (length(condiciones_presentes) < 2) return(NULL)
  obj$replicate <- NA
  for (cond in condiciones_presentes) {
    idx <- obj$orig.ident_uni == cond
    obj$replicate[idx] <- sample(paste0(cond, "_rep", 1:n_reps), sum(idx), replace = TRUE)
  }
  return(obj)
}

hacer_pseudobulk <- function(obj) {
  obj <- JoinLayers(obj)
  pseudo <- AggregateExpression(obj, group.by = "replicate", assays = "RNA",
                                return.seurat = FALSE, slot = "counts")
  counts <- as.data.frame(pseudo$RNA)
  colnames(counts) <- sub("^g", "", colnames(counts))
  counts[, sort(colnames(counts))]
}

correr_deseq2 <- function(counts_mat, comparaciones, output_dir) {
  rep_names   <- colnames(counts_mat)
  condition   <- gsub("-rep[0-9]+$", "", sub("^g", "", rep_names))
  if (length(unique(condition)) < 2) return(invisible(NULL))

  colData <- data.frame(row.names = rep_names,
                        condition = factor(condition, levels = c("0N", "0.5N", "5N")))
  if (any(is.na(colData$condition))) return(invisible(NULL))

  dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  available <- levels(colData$condition)[levels(colData$condition) %in% unique(condition)]

  for (comp in comparaciones) {
    conds <- comp$conds
    tag   <- comp$tag
    if (!all(conds %in% available)) next
    res <- results(dds, contrast = c("condition", conds[2], conds[1]))
    write.csv(as.data.frame(res), file = file.path(output_dir, tag, paste0("DESeq2_", tag, ".csv")))
  }
}

hacer_volcano <- function(file, output_dir, padj_cut = 0.05, lfc_cut = 1) {
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
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    labs(title = titulo, x = "Log2 Fold Change", y = "-Log10 adj p-value", color = "Significance") +
    theme_minimal()
}

procesar_deseq2_resultado <- function(file_path, output_dir, padj_cut = 0.05, lfc_cut = 1) {
  df          <- read_csv(file_path, show_col_types = FALSE)
  comparacion <- gsub("^DESeq2_(.*)\\.csv$", "\\1", basename(file_path))

  df_class <- df %>%
    mutate(AGI = `...1`,
           clasificacion = case_when(
             padj <= padj_cut & log2FoldChange >  lfc_cut ~  1,
             padj <= padj_cut & log2FoldChange < -lfc_cut ~ -1,
             TRUE ~ 0)) %>%
    dplyr::select(AGI, clasificacion) %>%
    setNames(c("AGI", comparacion))

  df_logfc <- df %>%
    mutate(AGI = `...1`,
           logfc = ifelse(padj <= padj_cut & abs(log2FoldChange) > lfc_cut, log2FoldChange, NA_real_)) %>%
    dplyr::select(AGI, logfc) %>%
    setNames(c("AGI", comparacion))

  df_filt <- df %>% filter(padj <= padj_cut, abs(log2FoldChange) > lfc_cut)
  write_csv(df_filt, file.path(output_dir, paste0(comparacion, "_filtrado.csv")))

  list(class = df_class, logfc = df_logfc)
}

hacer_heatmap <- function(matriz, min_genes = 1, deepSplit_val = 0, breaks = c(-5, 5)) {
  dist_rows <- dist(matriz, method = "euclidean")
  hc_rows   <- hclust(dist_rows, method = "complete")
  clust     <- cutreeDynamic(dendro = hc_rows, distM = as.matrix(dist_rows),
                             deepSplit = deepSplit_val, minClusterSize = min_genes,
                             pamRespectsDendro = FALSE)

  pca_res  <- prcomp(t(matriz), scale. = FALSE)
  var_exp  <- summary(pca_res)$importance[3, ]
  n_pcs    <- which(var_exp >= 0.90)[1]
  hc_cols  <- hclust(dist(pca_res$x[, 1:n_pcs]), method = "complete")

  paleta         <- colorRampPalette(brewer.pal(12, "Dark2"))(length(unique(clust[clust > 0])))
  annotation_row <- data.frame(Cluster = as.factor(clust))
  rownames(annotation_row) <- rownames(matriz)

  breaks_seq  <- seq(breaks[1], breaks[2], length.out = 80)
  color_scale <- colorRampPalette(c("blue", "black", "yellow"))(length(breaks_seq) - 1)

  pheatmap(matriz,
           cluster_rows = hc_rows, cluster_cols = hc_cols,
           annotation_row = annotation_row,
           annotation_colors = list(Cluster = setNames(paleta, sort(unique(clust[clust > 0])))),
           color = color_scale, breaks = breaks_seq,
           show_rownames = TRUE, border_color = NA,
           fontsize_row = 1, fontsize_col = 20, fontsize = 22,
           main = sprintf("Heatmap (%d genes)", nrow(matriz)))
}

# ============================================================
# FUNCIONES: GO ENRICHMENT
# ============================================================

correr_enriquecimiento_go <- function(tabla, universo, espacio,
                                       qvalueCutoff = 0.05, pvalueCutoff = 0.05,
                                       simplificar = FALSE, umbral_simply = 0.7,
                                       output_dir = "results/Enrichment") {
  salida <- vector("list", ncol(tabla))
  names(salida) <- colnames(tabla)

  for (n in seq_len(ncol(tabla))) {
    gene <- unique(trimws(gsub("\\..*", "", rownames(tabla)[tabla[, n] == 1])))
    if (length(gene) == 0) { message("Sin genes: ", colnames(tabla)[n]); next }

    enri <- tryCatch(
      enrichGO(gene = gene, universe = universo, OrgDb = org.At.tair.db,
               keyType = "TAIR", ont = espacio, pAdjustMethod = "BH",
               pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, readable = FALSE),
      error = function(e) NULL
    )

    if (is.null(enri) || nrow(enri@result) == 0) { message("Sin GO: ", colnames(tabla)[n]); next }

    # Guardar crudo y legible
    sufijo <- paste(colnames(tabla)[n], espacio, qvalueCutoff, sep = ".")
    write.table(as.data.frame(enri),
                file.path(output_dir, paste0(sufijo, ".txt")),
                sep = "\t", col.names = NA, quote = FALSE)
    write.table(as.data.frame(setReadable(enri, OrgDb = org.At.tair.db)),
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

podar_go <- function(resuGO, nivel, espacio, qvalueCutoff,
                     simplificar = FALSE, output_dir = "results/Enrichment") {
  salida <- vector("list", length(resuGO))
  names(salida) <- names(resuGO)

  for (k in seq_along(resuGO)) {
    if (is.null(resuGO[[k]])) next
    res <- tryCatch(gofilter(resuGO[[k]], nivel), error = function(e) NULL)
    if (is.null(res) || nrow(res@result) == 0) next
    salida[[k]] <- res

    sufijo <- paste(names(resuGO)[k], espacio, qvalueCutoff,
                    if (simplificar) "simply" else "total",
                    paste0("nivel_", nivel), "txt", sep = ".")
    write.table(as.data.frame(res), file.path(output_dir, sufijo),
                sep = "\t", col.names = NA, quote = FALSE)
  }
  return(salida)
}

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
    data.frame(Exp = nombres[k], GOid = df$ID, GODesc = df$Description,
               Log10Qvalue = -log10(df$qvalue), Enrichment = gr / br)
  })

  dat <- na.omit(do.call(rbind, bloques))
  dat$Exp <- factor(dat$Exp, levels = nombres)

  ggballoonplot(dat, x = "Exp", y = "GODesc", size = "Enrichment", fill = "Log10Qvalue") +
    scale_fill_gradientn(colors = brewer.pal(8, "YlOrRd")) +
    guides(size = "none") +
    theme_minimal(base_size = 11) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 28)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank())
}

# =============================================================================
# END OF FUNCTIONS
# =============================================================================
