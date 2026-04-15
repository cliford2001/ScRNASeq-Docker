---
tags: [scrnaseq, configuracion, parametros]
type: note
relacionado: [[ScRNASeq-Docker MOC]]
---

# ScRNASeq — Configuración

Todos los parámetros se definen en el bloque **CONFIGURATION** de `scrnaseq_pipeline.R` (líneas 29-70). No hay archivo de configuración externo.

---

## Bloque de configuración (scrnaseq_pipeline.R)

```r
# ── PATHS ──────────────────────────────────────────────
PIPELINE_DIR <- "/ruta/a/scripts"      # Dónde están los .R helpers
DATA_DIR     <- "/ruta/al/proyecto"    # Directorio raíz del proyecto
base_dir     <- file.path(DATA_DIR, "results")  # Subdirectorio de resultados

# ── TOGGLES ────────────────────────────────────────────
USE_CELLBENDER <- TRUE   # TRUE = input desde CellBender HDF5
                          # FALSE = input desde CellRanger matrix

# ── MUESTRAS ───────────────────────────────────────────
samples <- list(
  list(file      = "sample1/cellbender_filtered.h5",
       label     = "Sample1",
       condition = "control"),
  list(file      = "sample2/cellbender_filtered.h5",
       label     = "Sample2",
       condition = "nitrogen_low"),
  # ... más muestras
)

# ── COLORES (nombrados por condición) ───────────────────
colors <- c(
  "control"      = "#1f77b4",
  "nitrogen_low" = "#ff7f0e",
  "nitrogen_high"= "#2ca02c"
)

# ── SEMILLA ─────────────────────────────────────────────
set.seed(1807)   # Reproducibilidad global
```

---

## Parámetros de QC (Sección 2)

| Parámetro | Default | Descripción |
|---|---|---|
| `min_features` | 200 | Mínimo de genes detectados por célula |
| `max_features` | Inf | Máximo de genes (omitir límite superior) |
| `min_counts` | 0 | Mínimo de UMIs |
| `max_counts` | Inf | Máximo de UMIs |
| `max_mt` | 5 | % máximo de contenido mitocondrial |
| `max_cp` | 100 | % máximo de contenido cloroplástico |
| `run_doubletfinder` | TRUE | Activar detección de dobletes |
| `expected_doublet_rate` | 0.075 | Tasa esperada de dobletes (7.5%) |

---

## Parámetros de integración (Secciones 3-6)

| Parámetro | Default | Descripción |
|---|---|---|
| `nVariableFeatures` | 2000 | Genes variables para PCA |
| `nPCs` | 30 | Número de PCs calculados |
| `dims_use` | 1:20 | PCs usados en vecinos/UMAP |
| `k.param` | 20 | k-vecinos más cercanos |
| `cluster_resolution` | 0.3 | Resolución de clustering Leiden |

---

## Adaptación a otros organismos

| Parámetro | Arabidopsis | Humano | Ratón | Arroz |
|---|---|---|---|---|
| `mt_pattern` | `^ATMG` | `^MT-` | `^mt-` | `^OsMG` |
| `cp_pattern` | `^ATCG` | `NULL` | `NULL` | `^OsCG` |
| `orgdb` | `org.At.tair.db` | `org.Hs.eg.db` | `org.Mm.eg.db` | `org.Os.eg.db` |
| `keytype` | `"TAIR"` | `"ENSEMBL"` | `"ENSEMBL"` | `"ENTREZID"` |

Para desactivar QC cloroplástico (animales):
```r
load_sample(sample_info, mt_pattern="^MT-", cp_pattern=NULL)
```

---

## Parámetros de expresión diferencial

| Parámetro | Default | Descripción |
|---|---|---|
| `n_reps` | 3 | Número de pseudo-réplicas por condición |
| `padj_cut` | 0.05 | Umbral FDR para genes significativos |
| `lfc_cut` | 1 | Umbral log2FC para clasificación |

---

## Parámetros de enriquecimiento GO

| Parámetro | Default | Descripción |
|---|---|---|
| `espacio` | "BP" | Namespace GO (BP / MF / CC) |
| `qvalueCutoff` | 0.05 | Umbral q-value |
| `pvalueCutoff` | 0.05 | Umbral p-value |
| `simplificar` | FALSE | Eliminar términos redundantes |
| `umbral_simply` | 0.7 | Umbral Jaccard para simplificación |

---

## Archivo de referencia de marcadores

Requerido por `annotate_by_markers()`. Formato TSV:

```
gene	cell.types
AT1G01010	Epidermis
AT2G03000	Cortex
AT3G05020	Endodermis
...
```

- La pipeline **no funciona** sin este archivo
- Puede construirse desde literatura o bases de datos públicas
- El atlas GSE273033 provee validación independiente

---

## Variable global `output_dir`

Las funciones `save_pdf`, `save_vln`, `save_qc` usan una variable global:
```r
output_dir <- dir_XX   # Reasignar al inicio de cada sección
save_pdf(plot, "figura.pdf")  # Guarda en dir_XX/figura.pdf
```
> ⚠️ Olvidar reasignar `output_dir` guarda figuras en el directorio equivocado.

---

## Backlinks

- [[ScRNASeq-Docker MOC]]
- [[ScRNASeq - Flujo de Datos]]
- [[ScRNASeq - Decisiones de Diseño]]
