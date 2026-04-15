---
tags: [scrnaseq, pipeline, workflow]
type: note
relacionado: [[ScRNASeq-Docker MOC]]
---

# ScRNASeq — Flujo de Datos

## Diagrama general

```
FASTQ (lecturas crudas)
        ↓
[CellRanger 9.0.1]  ← alineamiento + conteo de UMIs
        ↓
raw_feature_bc_matrix.h5
        ↓
[CellBender 0.3.0]  ← remoción RNA ambiente (recomendado)
        ↓
cellbender_filtered.h5
        ↓
────────────────────────────── Pipeline R ──────────────────────────────
SECCIÓN 1  →  Carga + QC pre-filtro
SECCIÓN 2  →  Filtrado + DoubletFinder
SECCIÓN 3  →  Merge + Normalización + PCA + UMAP
SECCIÓN 4  →  Harmony (corrección de batch)
SECCIÓN 5  →  Selección de resolución (elbow + clustree)
SECCIÓN 6  →  Clustering final
SECCIÓN 7  →  DotPlot de marcadores (pre-anotación)
SECCIÓN 8  →  Anotación de tipos celulares
SECCIÓN 9  →  Clustree anotado
SECCIÓN 10 →  Visualización de expresión génica
SECCIÓN 11 →  Agrupamiento de anotaciones
SECCIÓN 12 →  Curación manual interactiva
        ↓
Objeto final Seurat (pbmc_harmony)
        ↓
    ┌───────────────────────────────────────┐
    ↓                                       ↓
Pseudobulk + DESeq2                  exportar_para_scanpy()
    ↓                                       ↓
Volcanos + Heatmaps              pbmc_harmony.h5ad
    ↓                            (para trajectories en Python)
Enriquecimiento GO
```

---

## Secciones del pipeline en detalle

### SECTION 0 — Diagrama del flujo
- Genera PDF con el diagrama visual del pipeline completo.
- Función: `plot_pipeline_workflow(outfile)`

---

### SECTION 1 — Carga y QC pre-filtro

**Entradas**: archivos HDF5 por muestra (CellBender o CellRanger)

**Funciones**:
- `load_sample(sample_info, mt_pattern, cp_pattern)` → objeto Seurat con métricas QC

**Métricas calculadas**:
- `percent.mt` → contenido mitocondrial (patrón `^ATMG` en Arabidopsis)
- `percent.cp` → contenido cloroplástico (patrón `^ATCG` en Arabidopsis)
- `nFeature_RNA` → genes detectados por célula
- `nCount_RNA` → UMIs totales por célula

**Salida**: lista `seurat_list_raw` + gráficos QC en violín (pre-filtro)

---

### SECTION 2 — Filtrado y detección de dobletes

**Función principal**: `filter_sample()` (llama a `doubletfinder_pipeline()` si `run_doubletfinder=TRUE`)

**Parámetros de filtrado**:
| Parámetro | Default | Descripción |
|---|---|---|
| `min_features` | 200 | Mínimo de genes detectados |
| `max_features` | Inf | Máximo de genes (por defecto sin límite) |
| `max_mt` | 5% | Máximo contenido mitocondrial |
| `max_cp` | 100% | Máximo contenido cloroplástico |
| `expected_doublet_rate` | 0.075 | Tasa esperada de dobletes (7.5%) |

> ⚠️ La tasa de dobletes de 7.5% es superior al estándar 10x (~0.8% per 1.000 células). Puede sobre-filtrar en algunos datasets.

**Salida**: lista `seurat_list` (singlets filtrados) + plots QC post-filtro

**Checkpoint**: `dir_02/seurat_list_postfilter.rds`

---

### SECTION 3 — Merge e integración inicial

**Pasos**:
1. `merge()` de todos los Seurat filtrados
2. `NormalizeData()` (log-normalización, factor 10.000)
3. `FindVariableFeatures()` (2.000 genes, método VST)
4. `ScaleData()`
5. `RunPCA()` (30 PCs)
6. `RunUMAP()` (en reducción PCA)

**Metadata creada**: `orig.ident_uni` → condición de cada célula

**Objeto**: `pbmc_harmony_preBatch` (con efectos de batch visibles en UMAP)

---

### SECTION 4 — Corrección de batch con Harmony

**Función**: `RunHarmony("orig.ident")`

> **Importante**: Todos los pasos downstream usan `reduction = "harmony"` (no PCA).

**Salida**: `pbmc_harmony_postBatch`

**Checkpoint**: `dir_04/pbmc_harmony_postharmony.rds`

Para restaurar desde checkpoint:
```r
pbmc_harmony <- readRDS("dir_04/pbmc_harmony_postharmony.rds")
```

---

### SECTION 5 — Selección de resolución

**Herramientas**:
- **Elbow plot**: k-means WSS para k=2-40 → orientativo para número de clusters
- **Clustree**: sweep de múltiples resoluciones → visualizar estabilidad de clusters

**Salida**: resolución óptima seleccionada manualmente (e.g., `0.3`)

---

### SECTION 6 — Clustering final

- `FindNeighbors()` + `FindClusters()` con resolución seleccionada
- `RunUMAP()` sobre reducción Harmony
- Guarda plots UMAP por cluster
- Imprime conteo de células por cluster

---

### SECTION 7 — DotPlot pre-anotación

- Lee archivo de referencia de marcadores bibliográficos (TSV: `gene | cell.types`)
- Genera dotplot de genes marcadores conocidos por cluster
- Guía la asignación manual de tipos celulares

---

### SECTION 8 — Anotación de tipos celulares

**8a — Por marcadores bibliográficos**:
```r
annotate_by_markers(seurat_obj, markers, reference_file)
# → metadata: celltype
```

**8b — Por transferencia de etiquetas (atlas de referencia)**:
```r
annotate_by_reference(seurat_obj, reference_obj, reference_col, dims=1:30)
# Usa atlas GSE273033 (Arabidopsis root)
# → metadata: celltype_reference
```

Ambos métodos se comparan mediante dotplots.

---

### SECTION 9 — Clustree anotado

- Re-ejecuta sweep de resoluciones con etiquetas de tipo celular superpuestas
- Confirma que la resolución elegida captura los tipos celulares conocidos

---

### SECTION 10 — Visualización de expresión

- Feature plots y violin plots para genes de interés
- Dentro de tipos celulares específicos
- ⚠️ En Seurat 5, requiere `JoinLayers()` antes de subsetear:
  ```r
  obj <- JoinLayers(obj)
  ```

---

### SECTION 11 — Agrupamiento de anotaciones

- Colapsa etiquetas finas en categorías amplias
- Crea columna `annotation_agrupada` en metadata
- Usa mapa de recodificación `grouping` definido por el usuario

---

### SECTION 12 — Curación interactiva ⚠️

> **Solo ejecutar en modo interactivo** — no hacer `source()` completo del script.

- Sub-clustering de poblaciones heterogéneas con `subclustar_tipo()`
- Inspección manual y reasignación de células
- Correcciones aplicadas al objeto global

---

## Análisis downstream

### Pseudobulk + Expresión Diferencial

```r
pseudobulk <- generate_pseudobulk(obj, group_by="orig.ident")
obj <- asignar_pseudoreplicados(obj, n_reps=3, seed=1807)
counts_mat <- hacer_pseudobulk(obj)
correr_deseq2(counts_mat, comparaciones, output_dir)
hacer_volcano(deseq2_file, padj_cut=0.05, lfc_cut=1)
```

### Enriquecimiento GO

```r
correr_enriquecimiento_go(tabla_binaria, universo, espacio="BP",
                          orgdb=org.At.tair.db, keytype="TAIR")
graficar_go_balones(resuGO)
```

### Exportar a Python/Scanpy

```r
exportar_para_scanpy(obj, "output.h5ad",
                     use_reduc=c("pca","umap","harmony"))
```

---

## Archivos de salida importantes

| Archivo | Descripción |
|---|---|
| `dir_02/seurat_list_postfilter.rds` | Checkpoint post-filtrado |
| `dir_04/pbmc_harmony_postharmony.rds` | Checkpoint post-Harmony |
| `results/FindAllMarkers.tsv` | Marcadores por cluster (cacheado) |
| `results/DESeq2_*/DESeq2_*.csv` | Resultados DE por comparación |
| `results/Enrichment/*.tsv` | Tablas de enriquecimiento GO |
| `pbmc_harmony.h5ad` | AnnData para scanpy/Python |
| `*.pdf` | Figuras del análisis |

---

## Backlinks

- [[ScRNASeq-Docker MOC]]
- [[ScRNASeq - Funciones R]]
- [[ScRNASeq - Configuración]]
