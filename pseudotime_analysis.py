"""
pseudotime_analysis.py
======================
Análisis de pseudotiempo sobre el objeto AnnData exportado desde el pipeline R.

Requiere:
    scanpy, scFates, palantir, matplotlib, seaborn, pandas, numpy

Input esperado:
    pbmc_harmony.h5ad  (generado con exportar_para_scanpy() en scrnaseq_pipeline.R)

Uso:
    python3 pseudotime_analysis.py
    # o dentro del contenedor Docker:
    docker exec -it r45 bash -c "source /opt/venv/bin/activate && python3 /workspace/pseudotime_analysis.py"
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

# Manipulación de datos y cálculos numéricos
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import scale
from sklearn.metrics import pairwise_distances_argmin

# Análisis de Single-Cell
import scanpy as sc
import scFates as scf
import palantir

# Visualización
import matplotlib
matplotlib.use("Agg")  # Backend no interactivo (para servidor/Docker)
import matplotlib.pyplot as plt
import seaborn as sns

# Configuración
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout

# ─────────────────────────────────────────────
# FUNCIONES
# ─────────────────────────────────────────────

def scfates_trajectories_alignment2(nombre, adata, clusters, root_cluster, nodes, sigma, eigen, seeed, plambda):
    """
    Construye una trayectoria usando scFates + Palantir + ElPiGraph
    en un subconjunto de células definidas por sus anotaciones agrupadas.
    """

    # --------------------------
    # 1. Subset por clusters
    # --------------------------
    if clusters is None or len(clusters) == 0:
        adata = adata.copy()
    elif len(clusters) == 1:
        adata = adata[adata.obs['annotation_curada_esp'] == clusters[0]].copy()
    else:
        adata = adata[adata.obs['annotation_curada_esp'].isin(clusters)].copy()

    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=0.5)

    # --------------------------
    # 2. Palantir → Diffusion Maps
    # --------------------------
    pca_projections = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)
    dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=50)
    ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs=eigen)

    adata.obsm["X_palantir"] = ms_data.values

    # --------------------------
    # 3. Vecindarios + ForceAtlas2
    # --------------------------
    sc.pp.neighbors(adata, n_neighbors=50, use_rep="X_palantir", method='umap')
    adata.obsm["X_pca2d"] = adata.obsm["X_pca"][:, :2]
    sc.tl.draw_graph(adata, init_pos='X_pca2d')

    title = f"{nombre}_{nodes}_nodes_{sigma}_sigma_{eigen}_eigen_{plambda}_lambda"
    sc.set_figure_params(figsize=(10, 10), dpi=100, dpi_save=100)

    # Visualización 1
    sc.pl.draw_graph(
        adata,
        color=['annotation_curada_esp'],
        add_outline=True,
        legend_fontsize=10, legend_fontoutline=2, size=50,
        title="ESP_FA_clusters",
        palette=sns.color_palette('colorblind'),
        save=f"ESP_FA_clusters.png"
    )

    # Visualización 2 (comentada — descomentár si se necesita)
    # sc.pl.draw_graph(
    #     adata,
    #     color=['annotation_chinos_uni'],
    #     add_outline=True,
    #     legend_fontsize=10, legend_fontoutline=2, size=50,
    #     title="CN_FA_clusters",
    #     palette=sns.color_palette('colorblind'),
    #     save=f"CN_FA_clusters2.png"
    # )

    # --------------------------
    # 4. ElPiGraph con scFates
    # --------------------------
    scf.tl.tree(
        adata,
        method="ppt",
        Nodes=nodes,
        use_rep="palantir",
        plot=False,
        device="cpu",
        seed=seeed,
        ppt_lambda=plambda,
        ppt_nsteps=200,
        ppt_sigma=sigma
    )

    scf.pl.graph(
        adata,
        title=title,
        color_cells='annotation_curada_esp',
        palette=sns.color_palette('colorblind'),
        save=f"_{title}_NODES_FA_diffusion_trajectory_clusters.png"
    )

    # --------------------------
    # 5. Elegir raíz (root)
    # --------------------------
    mask = (adata.obs['annotation_curada_esp'] == root_cluster).to_numpy()
    cluster_cells = adata.obs_names[mask]

    sub_conn = adata.obsp['connectivities'][mask][:, mask]
    degrees = np.array(sub_conn.sum(axis=1)).flatten()

    root_cell = cluster_cells[degrees.argmax()]
    print("Célula representativa:", root_cell)

    cell_id = root_cell

    # Resaltar célula raíz en el grafo
    adata.obs["highlight"] = adata.obs_names.isin([cell_id])
    adata_sorted = adata[adata.obs.sort_values("highlight").index, :]

    sc.pl.draw_graph(
        adata_sorted,
        color="highlight",
        palette=["lightgray", "red"],
        size=20,
        title=f"Célula resaltada: {cell_id}",
        frameon=False
    )

    # Marcar raíz y calcular pseudotiempo
    adata.obs["is_root"] = adata.obs_names == root_cell
    scf.tl.root(adata, root="is_root")

    # --------------------------
    # 6. Pseudotiempo
    # --------------------------
    scf.tl.pseudotime(adata)

    # Exportar pseudotiempo por célula
    cell_time = adata.obs[["t"]].copy()
    cell_time.index.name = "cell_id"

    if "annotation_curada_esp" in adata.obs.columns:
        cell_time["annotation_curada_esp"] = adata.obs["annotation_curada_esp"]

    out_file = f"{nombre}_pseudotime_por_celula.tsv"
    cell_time.to_csv(out_file, sep="\t")
    print(f"Pseudotiempo exportado a {out_file}")

    # --------------------------
    # 7. Trayectoria coloreada por pseudotiempo
    # --------------------------
    scf.pl.trajectory(
        adata,
        color_seg="t",
        basis="draw_graph_fa",
        frameon=False,
        s=50,
        scale_path=0.6,
        save=f"_{title}_pseudotime.png"
    )

    sc.pl.umap(adata, color="leiden")

    return adata


def scfates_trajectories_dendogram(adata):
    """
    Genera pseudotiempo, dendrograma y trayectoria segmentada.
    Requiere que scf.tl.tree() haya sido ejecutado previamente.
    """

    if "seg" not in adata.obs:
        raise ValueError("ERROR: adata.obs['seg'] no existe. Primero corre scf.tl.tree()")

    # -------------------------------
    # 1. Grafo coloreado por milestones
    # -------------------------------
    sc.pl.draw_graph(
        adata,
        color=["milestones", "leiden"],
        palette=sns.color_palette("colorblind"),
        add_outline=True,
        legend_fontsize=10,
        legend_fontoutline=2
    )

    # -------------------------------
    # 2. Dendrograma scFates
    # -------------------------------
    scf.tl.dendrogram(adata)

    scf.pl.dendrogram(
        adata,
        color="seg",
        palette=sns.color_palette("colorblind"),
        legend_fontoutline=True,
        legend_loc="on data",
        save="_seg.pdf"
    )

    scf.pl.dendrogram(
        adata,
        color="milestones",
        palette=sns.color_palette("colorblind"),
        legend_fontoutline=True,
        legend_loc="on data",
        save="_milestones.pdf"
    )

    print("Dendrograma generado")

    return adata


# ─────────────────────────────────────────────
# CONFIGURACIÓN
# ─────────────────────────────────────────────

INPUT_H5AD      = "/home/mvergara/projects2/eleo/ScRNA/metodologia/resultados/08_export/pbmc_harmony_curated.h5ad"
BASE_DIR        = "/home/mvergara/projects2/eleo/ScRNA/results/pseudotime"
ROOT_CELL_TYPE  = "Stem"                           # Tipo celular raíz del árbol (ajustar)
ANNOTATION_COL  = "celltype_reference_curated"     # Columna de anotación exportada desde R
N_NEIGHBORS     = 30                               # k-vecinos para el grafo
N_PCS           = 30                               # PCs a usar
RANDOM_SEED     = 1807

# ── Crear estructura de carpetas ──────────────────────────
os.makedirs(BASE_DIR, exist_ok=True)

OUTPUT_DIR     = os.path.join(BASE_DIR, "figures")
OUTPUT_TABLAS  = os.path.join(BASE_DIR, "tablas")

os.makedirs(OUTPUT_DIR,    exist_ok=True)
os.makedirs(OUTPUT_TABLAS, exist_ok=True)

sc.settings.figdir = OUTPUT_DIR

# ─────────────────────────────────────────────
# 1. CARGA DEL OBJETO (generado por pipeline R)
# ─────────────────────────────────────────────

print("── 1. Cargando AnnData ──")
Pc_esp = sc.read_h5ad(INPUT_H5AD)

# Recuperar embeddings exportados desde Seurat
# (exportar_para_scanpy() los guarda como 'UMAP' y 'PCA')
Pc_esp.obsm['X_umap'] = Pc_esp.obsm['UMAP'].values
Pc_esp.obsm['X_pca']  = Pc_esp.obsm['PCA'].values

sc.set_figure_params(figsize=(10, 10), dpi=300, dpi_save=300)

print(Pc_esp)

# UMAP de referencia coloreado por tipo celular
sc.pl.umap(
    Pc_esp,
    color=ANNOTATION_COL,
    show=False,
    save=f"_{ANNOTATION_COL}.png"
)

adata = Pc_esp  # alias para el resto del pipeline

# Verificar que existe la columna de anotación
if ANNOTATION_COL not in adata.obs.columns:
    available = list(adata.obs.columns)
    raise ValueError(
        f"Columna '{ANNOTATION_COL}' no encontrada en obs.\n"
        f"Columnas disponibles: {available}"
    )

# ─────────────────────────────────────────────
# 2. PREPROCESAMIENTO BÁSICO
# ─────────────────────────────────────────────

print("── 2. Preprocesamiento ──")

# Usar reducción Harmony si está disponible, si no usar PCA
if "X_harmony" in adata.obsm:
    reduction_key = "X_harmony"
    print("   Usando reducción: Harmony")
elif "X_pca" in adata.obsm:
    reduction_key = "X_pca"
    print("   Usando reducción: PCA")
else:
    print("   Calculando PCA desde cero...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, n_comps=N_PCS)
    reduction_key = "X_pca"

# Reconstruir grafo de vecinos sobre la reducción disponible
sc.pp.neighbors(adata, use_rep=reduction_key, n_neighbors=N_NEIGHBORS, random_state=RANDOM_SEED)

# Recalcular UMAP si no existe
if "X_umap" not in adata.obsm:
    print("   Calculando UMAP...")
    sc.tl.umap(adata, random_state=RANDOM_SEED)

# ─────────────────────────────────────────────
# 3. VISUALIZACIÓN INICIAL
# ─────────────────────────────────────────────

print("── 3. UMAP de referencia ──")

sc.pl.umap(
    adata,
    color=ANNOTATION_COL,
    title="Tipos celulares",
    save="_cell_types.pdf",
    show=False
)

# ─────────────────────────────────────────────
# 4. PSEUDOTIEMPO CON PALANTIR
# ─────────────────────────────────────────────

print("── 4. Pseudotiempo con Palantir ──")

# Diffusion map (requerido por Palantir)
sc.tl.diffmap(adata, n_comps=10)

# Seleccionar célula raíz: la primera célula del tipo celular raíz
root_cells = adata.obs.index[adata.obs[ANNOTATION_COL] == ROOT_CELL_TYPE]
if len(root_cells) == 0:
    raise ValueError(
        f"No se encontraron células de tipo '{ROOT_CELL_TYPE}' en '{ANNOTATION_COL}'.\n"
        f"Tipos disponibles: {adata.obs[ANNOTATION_COL].unique().tolist()}"
    )

# Usar el centroide en diffusion map como célula raíz
dm_coords = adata.obsm["X_diffmap"][adata.obs[ANNOTATION_COL] == ROOT_CELL_TYPE]
centroid   = dm_coords.mean(axis=0)
dists      = np.linalg.norm(adata.obsm["X_diffmap"] - centroid, axis=1)
root_cell  = adata.obs.index[np.argmin(dists)]
print(f"   Célula raíz seleccionada: {root_cell}")

# Correr Palantir
palantir_result = palantir.core.run_palantir(
    adata,
    early_cell=root_cell,
    use_early_cell_as_start=True,
    num_waypoints=500,
    n_jobs=4,
)

# Añadir pseudotiempo y entropía al objeto
adata.obs["palantir_pseudotime"] = palantir_result.pseudotime
adata.obs["palantir_entropy"]    = palantir_result.entropy

# Graficar pseudotiempo sobre UMAP
sc.pl.umap(
    adata,
    color=["palantir_pseudotime", "palantir_entropy"],
    cmap="viridis",
    save="_palantir_pseudotime.pdf",
    show=False
)

# ─────────────────────────────────────────────
# 5. ÁRBOL DE TRAYECTORIAS CON scFates
# ─────────────────────────────────────────────

print("── 5. Árbol de trayectorias con scFates ──")

# Aprender el árbol principal sobre el diffusion map
scf.tl.tree(
    adata,
    method="ppt",           # Principal Progression Tree
    Nodes=20,               # Número de nodos del árbol (ajustar según complejidad)
    use_rep="X_diffmap",
    seed=RANDOM_SEED,
)

# Calcular pseudotiempo de scFates desde la raíz
scf.tl.pseudotime(
    adata,
    n_jobs=4,
    seed=RANDOM_SEED,
)

# Detectar bifurcaciones
scf.tl.test_fork(adata, n_jobs=4)

# Graficar árbol sobre UMAP
scf.pl.graph(
    adata,
    basis="umap",
    color_cells=ANNOTATION_COL,
    save=f"{OUTPUT_DIR}/scfates_tree.pdf",
    show=False
)

# Graficar pseudotiempo de scFates
scf.pl.trajectory(
    adata,
    basis="umap",
    color_cells="t",   # 't' = pseudotiempo calculado por scFates
    save=f"{OUTPUT_DIR}/scfates_pseudotime.pdf",
    show=False
)

# ─────────────────────────────────────────────
# 6. GENES ASOCIADOS AL PSEUDOTIEMPO
# ─────────────────────────────────────────────

print("── 6. Genes asociados al pseudotiempo (scFates) ──")

# Test de asociación gen ~ pseudotiempo (GAM)
scf.tl.test_association(
    adata,
    n_jobs=4,
    A_cut=0.5,    # Amplitud mínima de cambio
)

# Filtrar genes significativos
sig_genes = adata.var[adata.var["signi"]].index.tolist()
print(f"   Genes significativos: {len(sig_genes)}")

# Guardar tabla de genes asociados
assoc_df = adata.var[adata.var["signi"]].sort_values("A", ascending=False)
assoc_df.to_csv(f"{OUTPUT_DIR}/pseudotime_genes_associated.csv")

# Heatmap de genes significativos a lo largo del pseudotiempo
if len(sig_genes) > 0:
    top_genes = assoc_df.head(50).index.tolist()
    scf.pl.trends(
        adata,
        features=top_genes,
        save=f"{OUTPUT_DIR}/pseudotime_gene_trends.pdf",
        show=False
    )

# ─────────────────────────────────────────────
# 7. GENES EN BIFURCACIONES
# ─────────────────────────────────────────────

print("── 7. Genes diferenciales en bifurcaciones ──")

# Solo si se detectaron bifurcaciones
if "fork" in adata.uns:
    scf.tl.test_fork(adata, n_jobs=4)

    fork_genes = adata.var[adata.var["signi_fk"]].index.tolist()
    print(f"   Genes en bifurcación: {len(fork_genes)}")

    fork_df = adata.var[adata.var["signi_fk"]].sort_values("A_fk", ascending=False)
    fork_df.to_csv(f"{OUTPUT_DIR}/pseudotime_fork_genes.csv")

    if len(fork_genes) > 0:
        scf.pl.trends(
            adata,
            features=fork_df.head(30).index.tolist(),
            save=f"{OUTPUT_DIR}/fork_gene_trends.pdf",
            show=False
        )
else:
    print("   No se detectaron bifurcaciones en el árbol.")

# ─────────────────────────────────────────────
# 8. GUARDAR OBJETO FINAL
# ─────────────────────────────────────────────

print("── 8. Guardando objeto final ──")

output_h5ad = f"{OUTPUT_DIR}/pbmc_pseudotime.h5ad"
adata.write_h5ad(output_h5ad)
print(f"   Guardado en: {output_h5ad}")

# ─────────────────────────────────────────────
# RESUMEN
# ─────────────────────────────────────────────

print("\n══════════════════════════════════════")
print("  PSEUDOTIEMPO COMPLETADO")
print("══════════════════════════════════════")
print(f"  Input:              {INPUT_H5AD}")
print(f"  Output dir:         {OUTPUT_DIR}")
print(f"  Células analizadas: {adata.n_obs:,}")
print(f"  Genes totales:      {adata.n_vars:,}")
print(f"  Genes en pseudot.:  {len(sig_genes)}")
print(f"  Objeto final:       {output_h5ad}")
print("══════════════════════════════════════\n")
