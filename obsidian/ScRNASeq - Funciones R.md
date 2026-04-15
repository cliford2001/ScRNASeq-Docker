---
tags: [scrnaseq, R, funciones, referencia]
type: reference
relacionado: [[ScRNASeq-Docker MOC]]
---

# ScRNASeq — Funciones R

Referencia de las ~50 funciones en `ScRNA_Analysis_Functions.R` (2.764 líneas).  
Todas se cargan con: `source("ScRNA_Analysis_Functions.R")`

---

## A. QC y Visualización

### `load_cellbender_filtered_h5(h5_path, project="Sample")`
- **Input**: ruta al HDF5 de CellBender filtrado
- **Output**: objeto Seurat con conteos crudos
- **Nota**: Parsea manualmente la estructura CSR sparse del HDF5 (data, indices, indptr, shape)

### `plot_qc_violin_grid(obj, label, color)`
- **Input**: Seurat object, etiqueta de muestra, color hex
- **Output**: ggplot2 de 4 paneles (nFeature_RNA, nCount_RNA, percent.mt, percent.cp)
- **Uso**: QC pre- y post-filtro

### `resumen_nFeature_plot(obj_list, etiquetas=NULL, colores=NULL)`
- **Input**: lista de Seurat objects
- **Output**: boxplot + tablas de cuartiles/quintiles comparando nFeature entre muestras

---

## B. Preprocesamiento y Detección de Dobletes

### `preprocesar_y_doubletfinder(seurat_obj, pcs=1:20, expected_doublet_rate=0.075, project_id="sample")`
- Normalización → VariableFeatures → ScaleData → PCA → DoubletFinder sweep
- **Output**: Seurat con clasificación de dobletes en metadata

### `doubletfinder_pipeline(obj, etiqueta="Sample", PCs=1:20, resolution=0.5, return_singlets=TRUE, sct=FALSE)`
- Como `preprocesar_y_doubletfinder` pero incluye clustering antes del sweep
- `return_singlets=TRUE` devuelve solo células clasificadas como singlets
- Metadata añadida: `doublet_class`

### `load_sample(sample_info, mt_pattern="^ATMG", cp_pattern="^ATCG")`
- **Input**: lista con `$file`, `$label`, `$condition`
- **Output**: Seurat anotado con percent.mt, percent.cp y prefijo de condición en barcode
- **Nota**: No aplica filtros — solo carga y agrega métricas

### `filter_sample(obj, min_features=200, max_features=Inf, min_counts=0, max_counts=Inf, max_mt=5, max_cp=100, run_doubletfinder=TRUE)`
- Aplica umbrales QC + opcionalmente DoubletFinder
- **Parámetro clave para plantas**: `max_cp` (en animales, pasar `cp_pattern=NULL` a `load_sample`)

### `process_sample(sample_info, mt_pattern, cp_pattern, min_features, ..., run_doubletfinder=TRUE)`
- Atajo que llama `load_sample()` + `filter_sample()` en una sola línea
- Usar cuando no se necesita inspección intermedia

---

## C. Utilidades Pseudobulk vs. Bulk

### `normalizar_bulk_pseudobulk(pseudobulk_counts, bulk_counts)`
- **Input**: dos vectores numéricos con nombres de genes
- **Output**: data frame con columnas `gene`, `pseudobulk`, `bulk` (log2-normalizados con DESeq2)
- **Uso**: validar agregados single-cell contra RNA-seq bulk

### `clasificar_residuos(df, umbral=5)`
- **Input**: data frame de `normalizar_bulk_pseudobulk()`
- **Output**: clasificación por residuo (Upregulated / Downregulated / Consistent)
- Usa regresión lineal `bulk ~ pseudobulk`

### `generate_pseudobulk(seurat_obj, group_by="orig.ident", merge_replicates=TRUE)`
- **Output**: lista con `$by_sample` y `$by_condition` (matrices sparse)
- Maneja la estructura multi-capa de Seurat 5

### `plot_replicate_correlation(pseudobulk_mat, main="Replicate Correlation")`
- **Output**: pheatmap con correlación de Pearson entre réplicas

---

## D. Utilidades Seurat

### `unificar_nombres(obj)`
- Limpia sufijos numéricos en niveles de identidad (e.g., `.1`, `_2`)

### `mostrar_tabla(filtered_vec, cellbender_vec, titulo="Annotations")`
- Tabla comparativa side-by-side de anotaciones (filtrado vs. CellBender)

### `exportar_para_scanpy(seurat_obj, outfile, assay_name="RNA", use_reduc=c("pca","umap","harmony"), X_name="logcounts", overwrite=TRUE)`
- Convierte a SingleCellExperiment → escribe `.h5ad` via zellkonverter
- Fallback a SeuratDisk si falla zellkonverter

### `safe_vln(obj, feature, colors)`
- VlnPlot seguro para RMarkdown (agrupado por orig.ident, sin overlay de puntos)

### `unir_layers_counts(obj, capas)`
- Concatena matrices sparse de múltiples capas de Seurat 5
- Necesario antes de ciertos análisis downstream en datasets integrados

---

## E. Anotación de Tipos Celulares

### `find_markers(seurat_obj, output_file="results/FindAllMarkers.tsv", only_pos=TRUE, min_pct=0.25, logfc_threshold=0.25, force=FALSE)`
- `FindAllMarkers()` con test Wilcoxon
- **Cacheo automático**: si el TSV existe, lo lee en lugar de recalcular
- `force=TRUE` fuerza recalcular
- ⚠️ Computacionalmente costoso en datasets grandes

### `annotate_by_markers(seurat_obj, markers, reference_file)`
- **Input**: markers de `find_markers()` + TSV bibliográfico (`gene | cell.types`)
- **Output**: Seurat con columna `celltype` en metadata

### `annotate_by_reference(seurat_obj, reference_obj, reference_col, dims=1:30)`
- `FindTransferAnchors()` + `TransferData()` para transferencia de etiquetas desde atlas
- **Output**: Seurat con columna `celltype_reference` en metadata
- Atlas usado: GSE273033 (Arabidopsis root)

### `subclustar_tipo(obj, tipo, annot_col="annotation_agrupada", resolution=0.3, dims=1:20)`
- Sub-clustering de un tipo celular específico
- Re-ejecuta PCA, UMAP, neighbors y clustering a menor resolución
- **Output**: Seurat del tipo celular con sub-clusters

---

## F. Pseudobulk, DESeq2, Volcanes, Heatmaps

### `asignar_pseudoreplicados(obj, condiciones=NULL, n_reps=3, seed=1807)`
- Asigna células a `n_reps` pseudo-réplicas por condición al azar
- **Columna**: `replicate` en metadata
- **Uso crítico**: permite DESeq2 en experimentos con una sola réplica biológica

### `hacer_pseudobulk(obj)`
- `AggregateExpression()` por pseudo-réplica
- **Output**: data frame genes × réplicas (formato listo para DESeq2)

### `correr_deseq2(counts_mat, comparaciones, output_dir, tipo=NULL)`
- **Input comparaciones**: `list(list(conds=c("ref", "treat"), tag="label"))`
- Construye DESeqDataSet, auto-detecta niveles de condición, escribe CSVs
- **Output**: `output_dir/tag/DESeq2_tag.csv` por comparación

### `hacer_volcano(file, padj_cut=0.05, lfc_cut=1)`
- Lee CSV de DESeq2, clasifica genes, genera volcano plot con colores
- **Output**: ggplot2

### `procesar_deseq2_resultado(file_path, output_dir, padj_cut=0.05, lfc_cut=1)`
- Clasifica genes (up=1, down=-1, unchanged=0)
- Escribe `*_significant_genes.csv`
- **Output**: lista con `$class` y `$logfc`
- **Uso**: preparar matriz binaria para enriquecimiento GO

### `hacer_heatmap(matriz, min_genes=1, deepSplit_val=0, breaks=c(-5, 5))`
- **Input**: matriz genes × condiciones (típicamente log2FC)
- Clustering jerárquico filas (Euclidiana, dynamic tree cut) y columnas (distancia PCA)
- **Escala de color**: azul-negro-amarillo
- `deepSplit_val`: 0 = conservador, mayor = más splits

### `hacer_dotplot_marcadores(seurat_obj, marks, annot_col, cell_order, clusters_remove=NULL, rename_map=NULL, outfile=NULL, width=20, height=10, dot_scale=12, base_size=18)`
- DotPlot coordinado con orden personalizado de tipos celulares y genes
- **Output**: ggplot2; opcionalmente guarda en PDF
- Los genes deben formar patrón diagonal para validar anotación

---

## G. Enriquecimiento GO

### `correr_enriquecimiento_go(tabla, universo, espacio, orgdb=org.At.tair.db, keytype="TAIR", qvalueCutoff=0.05, pvalueCutoff=0.05, simplificar=FALSE, umbral_simply=0.7, output_dir="results/Enrichment")`
- **Input tabla**: matriz binaria genes × comparaciones (1=upregulated)
- **Input espacio**: "BP" / "MF" / "CC"
- Escribe TSVs con IDs y con símbolos de genes
- `simplificar=TRUE` elimina términos GO redundantes (umbral Jaccard 0.7)

### `podar_go(resuGO, nivel, espacio, qvalueCutoff, simplificar=FALSE, output_dir="results/Enrichment")`
- Filtra términos GO por nivel de jerarquía (`gofilter()`)
- **Uso**: nivel 4 para términos más amplios

### `graficar_go_balones(resuGO)`
- Gráfico de burbujas: tamaño = fold enrichment, color = -log10(qvalue)
- **Output**: ggplot2

---

## H. Guardar Figuras

### `save_pdf(plot, file, w=10, h=8)`
- Guarda en `output_dir/file` a 300 dpi
- ⚠️ Requiere variable global `output_dir` definida

### `save_vln(plot, file, n=1)`
- Guarda VlnPlot con altura auto-escalada (14 × 6n pulgadas)

### `save_qc(plot_list, file)`
- Guarda lista de plots apilados verticalmente (patchwork), ancho 14, alto 6 por panel

---

## I. Infraestructura del Pipeline

### `create_pipeline_dirs(base_dir)`
- Crea subdirectorios `dir_00` a `dir_12` bajo `base_dir`
- **Output**: lista nombrada de rutas por sección

### `plot_pipeline_workflow(outfile)`
- Genera PDF con diagrama visual del flujo completo del pipeline

---

## Scripts auxiliares

### `load_libraries.R`
Carga ~60 paquetes en el orden correcto de dependencias. Llamar al inicio de cada análisis:
```r
source("load_libraries.R")
# → imprime: "All libraries loaded successfully."
```

### `custom_seurat.R`
Función única: `plot_integrated_clusters(srat)`
- **Output**: patchwork de 2 paneles:
  - Izquierda: stacked bar chart (fracción de células por dataset dentro de cada cluster)
  - Derecha: log10 de conteo de células por cluster con etiquetas

---

## Backlinks

- [[ScRNASeq-Docker MOC]]
- [[ScRNASeq - Flujo de Datos]]
- [[ScRNASeq - Configuración]]
