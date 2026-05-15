# Capitulo 3 — Pseudotime Trajectory Analysis

Pseudotime analysis of single-cell RNA-seq data using scFates and Palantir.
Requires the curated AnnData object produced by `capitulo1_single_cell.R`.

---

## Overview

This module orders cells along a developmental trajectory and identifies which
genes change progressively as cells differentiate. Instead of asking "what is
different between two conditions?", pseudotime asks "in what order do cells
change, and what genes drive that change?".

---

## How to run

Open the notebook in Jupyter and run each section in order:

```
workflow/capitulo3_pseudotime.ipynb
```

The Python environment must be `scrna_seba`. This is set automatically by
`.Rprofile` when opening RStudio, and by the shebang line when running as a
script. In Jupyter, select the kernel **Python (scrna_seba)**.

---

## Sections

### SECTION 24 — Setup
Loads Python libraries and defines output directories. No parameters to change
unless the folder structure differs from capitulo1.

**Output:** creates `09_pseudotime/` inside the results folder.

---

### SECTION 25 — Load Data
Reads `pbmc_harmony_curated.h5ad` exported by capitulo1 and converts coordinate
keys from R/Seurat format to scanpy format. Prints available cell types and
plots a full UMAP as a visual check.

**Parameter to set:**
- `ANNOTATION_COL` — metadata column with the cell type labels to use

**Output:**
```
09_pseudotime/umap_overview.png   — UMAP of all cells coloured by cell type
```

---

### SECTION 26 — Cell Type Selection
Displays all available cell types with their cell counts and lets the user
define which ones to include in the trajectory. Plots a UMAP of the subset
to verify the selection before proceeding.

**Parameter to set:**
- `TRAJECTORY_CLUSTERS` — list of cell types to include, e.g.
  `["Pavement Cell", "Guard Cell", "Stomatal Line"]`

**Output:**
```
09_pseudotime/umap_selection.png  — UMAP of the selected cell type subset
```

---

### SECTION 27 — Trajectory Inference
Builds a tree-shaped trajectory through the selected cell types:

1. **Palantir diffusion maps** — reduces high-dimensional gene expression into
   a space that captures developmental relationships better than PCA alone.
2. **scFates PPT** (Principal Polynomial Tree) — fits a tree onto the diffusion
   space that traces the most likely differentiation path.
3. **Root selection** — the most-connected cell in `ROOT_CLUSTER` is set as
   the start of pseudotime (t = 0).
4. **Pseudotime** — each cell receives a value from 0 (root) to 1 (most
   differentiated), representing its position in the differentiation process.

**Parameters to set:**
- `ROOT_CLUSTER` — least-differentiated cell type; pseudotime starts here
- `NODES` — tree resolution (50–200); more nodes = finer branch detail
- `SIGMA` — smoothing (0.1–0.5); lower = tree follows cells more tightly
- `PPT_LAMBDA` — branch complexity penalty; higher = simpler, fewer branches
- `N_EIGS` — diffusion map dimensions used (15–25, must be < 50)

**Output:**
```
09_pseudotime/trajectory/trajectory_annotation.pdf  — tree coloured by cell type
09_pseudotime/trajectory/trajectory_leiden.pdf      — tree coloured by Leiden cluster
```

---

### SECTION 28 — Dendrogram and Milestones
Organises the trajectory branches into a dendrogram and labels each branch
endpoint as a "milestone". Each milestone represents a terminal differentiation
state (e.g. Guard Cell, Pavement Cell).

Run this section to discover the milestone names, which you need for Section 29.

**Output:**
```
09_pseudotime/umap_milestones.pdf  — dendrogram coloured by milestone
```
Console prints the available milestone names — copy them into Section 29.

---

### SECTION 29 — Milestone Analysis
For each selected milestone (branch endpoint), runs two analyses:

1. **Association test** — identifies which genes change significantly along
   that differentiation branch. Uses a non-linear association score; genes
   below the `A_CUT` threshold are discarded.
2. **Trend fitting** — fits smooth expression curves over pseudotime for the
   significant genes, producing the fitted values used in Section 30.

Saves intermediate checkpoints so results can be reloaded without re-running.

**Parameters to set:**
- `MILESTONES_TO_ANALYZE` — branch endpoints to analyse (from Section 28 output)
- `A_CUT` — association threshold (0–1); lower = more genes, higher = more selective
- `P_VAL_CUT` — p-value threshold for gene–pseudotime significance

**Output:**
```
09_pseudotime/milestones/adata_pseudotime_<milestone>_association.h5ad
09_pseudotime/milestones/adata_pseudotime_<milestone>_fitted.h5ad
```

---

### SECTION 30 — Gene Expression Trends
For each milestone, generates:

1. **Trends heatmap** — genes ordered by the time point of peak expression,
   showing which genes are active early, mid, or late in the trajectory.
   Known marker genes can be highlighted.
2. **Peak expression table** — CSV ranking all significant genes by their
   pseudotime peak, including Pearson correlation with pseudotime and the
   dominant cell cluster at peak.

**Parameters to set:**
- `HIGHLIGHT_GENES` — gene IDs to highlight in the heatmap (leave `[]` if none)
  For Arabidopsis, use TAIR IDs (e.g. `"AT3G24140"`).

**Output:**
```
09_pseudotime/trends/<milestone>_gene_trends.pdf    — expression heatmap
09_pseudotime/tables/<milestone>_genes_by_peak.csv  — ranked gene table
```

---

### SECTION 31 — Module Scores *(optional)*
Projects user-defined gene lists onto the trajectory as a normalised per-cell
score. Useful for overlaying published gene signatures or custom gene modules
onto the differentiation map.

Skip this section entirely if no custom gene lists are available.

**Parameters to set:**
- `MODULE_GENE_FILES` — dict of `{label: path_to_file}`; leave `{}` to skip
- `MODULE_ID_COL` — column with gene IDs in each file (e.g. `"ID"` for Arabidopsis)

**Output:**
```
09_pseudotime/draw_graph_module_scores.png  — trajectory coloured by module score
```

---

## Output directory structure

```
09_pseudotime/
├── umap_overview.png          Section 25 — all cells
├── umap_selection.png         Section 26 — selected subset
├── trajectory/
│   ├── trajectory_annotation.pdf
│   └── trajectory_leiden.pdf
├── umap_milestones.pdf        Section 28
├── milestones/
│   ├── adata_pseudotime_<milestone>_association.h5ad
│   └── adata_pseudotime_<milestone>_fitted.h5ad
├── trends/
│   └── <milestone>_gene_trends.pdf
└── tables/
    └── <milestone>_genes_by_peak.csv
```

---

## Dependencies

All managed by the `scrna_seba` conda environment:

| Package | Role |
|---------|------|
| `scanpy` | Data handling, UMAP visualisation |
| `scFates` | Trajectory inference (PPT), pseudotime, milestone analysis |
| `palantir` | Diffusion maps for trajectory embedding |
| `matplotlib` / `seaborn` | Plotting |
| `pandas` / `numpy` | Data manipulation |

---

## Key concepts

**Pseudotime** — a continuous value assigned to each cell representing its
position along a differentiation path. It is inferred from gene expression
similarity, not from real time.

**Milestone** — a branch endpoint in the trajectory tree. Each milestone
corresponds to a terminal cell state (a fully differentiated cell type).

**Diffusion map** — a dimensionality reduction that captures the global
geometry of cell differentiation better than PCA, used internally by Palantir
to compute the developmental distances between cells.

**PPT (Principal Polynomial Tree)** — the algorithm used by scFates to fit
a tree-shaped path through the cells in diffusion map space.
