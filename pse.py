

```{r setup-python, include=FALSE}
library(reticulate)

# Seleccionar Python del entorno scrna_seba
use_python("/home/mvergara/projects3/app/miniconda/envs/scrna_seba/bin/python", required = TRUE)

# Confirmar entorno cargado
py_config()
```

```{python carga_python}
#os.environ["R_HOME"] = "/home/mvergara/projects3/app/miniconda/envs/scrna_seba/lib/R"
import os, sys
import scFates as scf
import scanpy as sc
import warnings
import numpy as np
import pandas as pd
import palantir
import seaborn as sns
from sklearn.metrics import pairwise_distances_argmin
warnings.filterwarnings("ignore")
import sys
# Manipulación de datos y cálculos numéricos
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.preprocessing import scale

# Análisis de Single-Cell
import scanpy as sc
import scFates as scf

# Visualización
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

# Sistema y archivos
import os

# Configuración específica (opcional, según tu entorno)
matplotlib.use("Agg")
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
```

```{python}
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


    sc.pp.neighbors(adata)    # si no existe
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

    # Visualización 2
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
    #cell_id = "5N_GGTTACGCAGGGACTC_1_5"
    #root_cell = cell_id
    cell_id = root_cell
    # Crear columna booleana para resaltar la célula
    adata.obs["highlight"] = adata.obs_names.isin([cell_id])

    # Reordenar para que la célula roja quede sobre todas
    adata_sorted = adata[adata.obs.sort_values("highlight").index, :]

    # Graficar UMAP
    sc.pl.draw_graph(
        adata_sorted,
        color="highlight",
        palette=["lightgray", "red"],  # colores fijos
        size=20,
        title=f"Célula resaltada: {cell_id}",
        frameon=False
    )
  

    # Marcar la célula raíz
    adata.obs["is_root"] = adata.obs_names == root_cell
    scf.tl.root(adata, root="is_root")

    # --------------------------
    # 6. Pseudotiempo
    # --------------------------
    scf.tl.pseudotime(adata)

    # Exportar pseudotiempo
    cell_time = adata.obs[["t"]].copy()
    cell_time.index.name = "cell_id"

    if "annotation_esp" in adata.obs.columns:
        cell_time["annotation_curada_esp"] = adata.obs["annotation_curada_esp"]

    out_file = f"{nombre}_pseudotime_por_celula.tsv"
    cell_time.to_csv(out_file, sep="\t")
    print(f"Pseudotiempo exportado a {out_file}")

    # --------------------------
    # 7. Graficar trayectoria coloreada por pseudotiempo
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

```

```{python}
def scfates_trajectories_dendogram(adata):
    """
    Genera pseudotiempo, dendrograma y trayectoria segmentada.
    """

    if "seg" not in adata.obs:
        raise ValueError("ERROR: adata.obs['seg'] no existe. Primero corre scf.tl.tree()")
      
    # -------------------------------
    # 3. Grafo coloreado por milestones
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
    # 4. Dendrograma scFates
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
```

```{python}

os.makedirs("/home/mvergara/projects2/eleo/ScRNA/results/pseudotime/", exist_ok=True)
os.makedirs("/home/mvergara/projects2/eleo/ScRNA/results/pseudotime/figures2", exist_ok=True)
os.makedirs("/home/mvergara/projects2/eleo/ScRNA/results/pseudotime/tablas", exist_ok=True)
outdir = "/home/mvergara/projects2/eleo/ScRNA/results/pseudotime/figures2"
os.makedirs(outdir, exist_ok=True)
sc.settings.figdir = outdir
```

```{python, fig.width=10, fig.height=10}

Pc_esp =  sc.read_h5ad('/home/mvergara/projects2/eleo/ScRNA/metodologia/resultados/08_export/pbmc_harmony_curated.h5ad')

# Convert the DataFrame to a NumPy array and assign to X_umap
Pc_esp.obsm['X_umap'] = Pc_esp.obsm['UMAP'].values
Pc_esp.obsm['X_pca'] = Pc_esp.obsm['PCA'].values

sc.set_figure_params(figsize=(10, 10), dpi=300, dpi_save=300)


sc.pl.umap(Pc_esp, color='celltype_reference_curated', show=False, save='.png')

```

```{python}

### Epidermal trajectory
sel = ['Pavement Cell', 'Guard Cell', 'Meristemoid']

for eigen in [51, 52, 53, 54]:
   for nodes in [150]:
     for sigma in [0.2]:
       for plambda in [60, 100]: #100 y 60
             Pc_esp_2 = scfates_trajectories_alignment2(
                  nombre="Pc_esp",
                  adata=Pc_esp,
                  clusters=sel,                 # solo estos clusters
                  root_cluster="Meristemoid",    # raíz en procambium
                  nodes=nodes,
                  sigma=sigma,
                  eigen=eigen,
                  seeed=3,
                  plambda=plambda)

```

```{python, fig.width=10, fig.height=10}

fig, ax = plt.subplots(figsize=(10, 10))

scf.pl.graph(
    Pc_esp_2,
    title="Pavement - Meristemoid - Guard Cell",
    color_cells='annotation_curada_esp',
    palette=sns.color_palette('colorblind'),
    ax=ax,             # <-- clave
    show=True         # <-- no mostrar en pantalla
)


plt.tight_layout()
plt.subplots_adjust(right=0.85)   # deja espacio para la leyenda
#fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/Tree_Pv_Ms_Gc_60.pdf")


scf.pl.trajectory(
    Pc_esp_2,
    color_seg="t",
    basis="draw_graph_fa",
    frameon=False,
    s=50,
    scale_path=0.6,
    ax=ax,             # <-- clave
    show=True         # <-- no mostrar en pantalla
)

plt.tight_layout()
plt.subplots_adjust(right=0.85)   # deja espacio para la leyenda
#fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/Trayectory_Pv_Ms_Gc_60.pdf")
```

```{python, fig.width=15, fig.height=5}

# -------------------------------------------------------------
# 1. Cargar listas de genes
# -------------------------------------------------------------
cluster2 = pd.read_csv("/home/projects2/mvergara/eleo/ScRNA/recursos/karin/Cluster2.txt", sep="\t")
cluster4 = pd.read_csv("/home/projects2/mvergara/eleo/ScRNA/recursos/karin/Cluster4.txt", sep="\t")
cluster4v2 = pd.read_csv("/home/projects2/mvergara/eleo/ScRNA/recursos/karin/cluster4v2.txt", sep="\t")

print(cluster4)
vector_cluster2 = cluster2["Cluster2"].dropna().unique().tolist()
vector_cluster4 = cluster4["Cluster4"].dropna().unique().tolist()
vector_cluster4v2 = cluster4v2["ID"].dropna().unique().tolist()

# -------------------------------------------------------------
# 2. Definir objeto principal
# -------------------------------------------------------------
# Aquí tu objeto AnnData equivalente a pbmc_harmony
adata = Pc_esp_2  # si ya está cargado

# -------------------------------------------------------------
# 3. Filtrar genes que estén presentes en la matriz
# -------------------------------------------------------------
genes2 = [g for g in vector_cluster2 if g in adata.var_names]
genes4 = [g for g in vector_cluster4 if g in adata.var_names]
genes4v2 = [g for g in vector_cluster4v2 if g in adata.var_names]

# -------------------------------------------------------------
# 4. Función para calcular el módulo
# -------------------------------------------------------------
def compute_module_score(adata, gene_list, prefix):
    X = adata[:, gene_list].X  # matriz de expresión para esos genes
    if not isinstance(X, np.ndarray):
        X = X.toarray()  # convertir de sparse matrix a array
    
    # a. número de genes detectados (>0)
    detected = (X > 0).sum(axis=1)
    
    # b. z-score del número detectado
    detected_z = scale(detected)
    
    # c. promedio de expresión
    mean_expr = X.mean(axis=1)
    
    # d. score ponderado (z × promedio)
    module_score = detected_z * mean_expr

    # e. guardar en .obs
    adata.obs[f"detected_genes_{prefix}"] = detected
    adata.obs[f"detected_genes_{prefix}_z"] = detected_z
    adata.obs[f"module_expr_{prefix}"] = mean_expr
    adata.obs[f"module_score_{prefix}"] = module_score

    print(f"Calculado módulo: {prefix} ({len(gene_list)} genes)")

# -------------------------------------------------------------
# 5. Calcular para cada cluster
# -------------------------------------------------------------
compute_module_score(adata, genes2, "cluster2")
compute_module_score(adata, genes4, "cluster4")
compute_module_score(adata, genes4v2, "cluster4v2")

# -------------------------------------------------------------
# 6. Visualización opcional
# -------------------------------------------------------------
sc.set_figure_params(figsize=(15, 1), dpi=300, dpi_save=300)
sc.pl.draw_graph(adata, color=["module_score_cluster2", "module_score_cluster4v2", "module_score_cluster4"], cmap="viridis", save="karin.png")
sc.pl.draw_graph(adata, color="AT3G10525", cmap="viridis", save="lgo.png")
```

```{python}


def module_sum_score(adata, gene_list, prefix):
    # Filtrar genes válidos
    genes = [g for g in gene_list if g in adata.var_names]
    if len(genes) == 0:
        print(f"No se encontraron genes del módulo {prefix}")
        return
    
    # Matriz de expresión células × genes
    X = adata[:, genes].X
    if not isinstance(X, np.ndarray):
        X = X.toarray()

    # Suma de expresión del módulo por célula
    raw_score = np.asarray(X.sum(axis=1)).flatten()

    # Normalización por library size (total UMIs por célula) + log-transform
    library_size = np.asarray(adata.X.sum(axis=1)).flatten()
    norm_score = np.log1p(raw_score / library_size)

    # Z-score para comparabilidad entre células
    z_score = scale(norm_score)

    # Almacenar score final
    adata.obs[f"{prefix}_module_score"] = z_score

    print(f"✓ Score de módulo generado ({prefix}) → {len(genes)} genes")
    
module_sum_score(adata, genes2, "cluster2")
module_sum_score(adata, genes4, "cluster4")
module_sum_score(adata, genes4v2, "cluster4v2")
sc.pl.draw_graph(
    adata,
    color=["cluster2_module_score", "cluster4_module_score", "cluster4v2_module_score"],
    cmap="viridis"
)
```





```{python}
Pc_esp_3 = scfates_trajectories_dendogram(Pc_esp_2)
```

```{python, fig.width=15, fig.height=15}
sc.set_figure_params(figsize=(10, 10))  # ancho=20, alto=6
sc.pl.draw_graph(Pc_esp_3,color='t',color_map = 'viridis', add_outline=True, size =200, legend_fontsize=10, legend_fontoutline=2)


# sc.pl.draw_graph(
#     Pc_esp_3,
#     color=["AT3G17820", "AT1G76650", "AT5G53460", 
#                         "AT5G16570", "AT4G35270", "AT5G40850", 
#                         "AT4G24020", "AT5G24400", "AT1G66200", 
#                         "AT5G60410", "AT2G15620", "AT5G50950"],
#     add_outline=True,
#     legend_fontsize=10,
#     legend_fontoutline=2,
#     size=200,
#     cmap='viridis',
#     title='AT3G10525')
    
sc.pl.draw_graph(
    Pc_esp_3,
    color=["AT4G21750", "AT2G26250"],
    add_outline=True,
    legend_fontsize=10,
    legend_fontoutline=2,
    size=200,
    cmap='viridis',
    title='Epidermis')
    
sc.pl.draw_graph(
    Pc_esp_3,
    color=["AT5G38420", "AT1G29910"],
    add_outline=True,
    legend_fontsize=10,
    legend_fontoutline=2,
    size=200,
    cmap='viridis',
    title='Mesophyll')
```

```{python}
### chage milestones names
milestone = pd.DataFrame(Pc_esp_3.obs.groupby('milestones'))[0].tolist()
milestone
new_names_milestone = ['0', 'Guard_Cell', '16', '22', 'Meristemoid', '3', '34', '38', '43', '46', '6']

#
scf.tl.rename_milestones(Pc_esp_3,new_names_milestone)

# fig, ax = plt.subplots(figsize=(10, 10))
# 
# scf.pl.graph(
#     Pc_esp_3,
#     title="Milestone",
#     color_cells="milestones",
#     palette=sns.color_palette('colorblind'),
#     ax=ax,             # <-- clave
#     show=True)
# 
# 
# plt.tight_layout()
# plt.subplots_adjust(right=0.65)   # deja espacio para la leyenda
# fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_eleo/Pseudotime/Milestones_Pv_Ms.pdf")
```

```{python}
### REPLACE MILESTONE FOR SOMETHING MORE MEANINGFUL
seg = Pc_esp_3.obs['milestones'].tolist()
seg = ['Guard_Cell' if cell == 'guard_cell' else cell for cell in seg]
seg = ['Meristemoid' if cell == 'meristemoid' else cell for cell in seg]
seg = ['Pavement_1' if cell == 'pavement_1' else cell for cell in seg]
seg = ['Pavement_2' if cell == 'pavement_2' else cell for cell in seg]
seg = ['Pavement_3' if cell == 'pavement_3' else cell for cell in seg]
seg = ['Pavement_4' if cell == 'pavement_4' else cell for cell in seg]
seg = ['Brake' if cell == 'brake' else cell for cell in seg]
Pc_esp_3.obs['milestones'] = seg
```

```{python}
leiden = Pc_esp_3.obs['leiden'].tolist()

leiden = ['Guard_Cell' if cell == 'guard_cell' else cell for cell in leiden]
leiden = ['Meristemoid' if cell == 'meristemoid' else cell for cell in leiden]
leiden = ['Pavement_1' if cell == 'pavement_1' else cell for cell in leiden]
leiden = ['Pavement_2' if cell == 'pavement_2' else cell for cell in leiden]
leiden = ['Pavement_3' if cell == 'pavement_3' else cell for cell in leiden]
leiden = ['Pavement_4' if cell == 'pavement_4' else cell for cell in leiden]
leiden = ['Brake' if cell == 'brake' else cell for cell in leiden]

Pc_esp_3.obs['leiden'] = leiden
```

```{python}

# ============================
# CONFIGURACIÓN GENERAL
# ============================

path = "/home/mvergara/projects2/eleo/ScRNA/results/pseudotime/tablas"
name_file = "all_internal"

os.makedirs(path, exist_ok=True)

# Asegurar que milestones sea categoría — SCFates lo exige
if Pc_esp_3.obs["milestones"].dtype.name != "category":
    Pc_esp_3.obs["milestones"] = Pc_esp_3.obs["milestones"].astype("category")

# Milestones que quieres recorrer
milestones_to_analyze = ['Guard_Cell']


# ============================
# LOOP PRINCIPAL
# ============================

for milestone in milestones_to_analyze:

    print(f"\n===== Procesando milestone: {milestone} =====")

    sc.set_figure_params(figsize=(3, 4), dpi_save=600, frameon=False)

    # ----------------------------
    # 1) SUBSET DEL ÁRBOL
    # ----------------------------
    try:
        adata_Ic = scf.tl.subset_tree(
            Pc_esp_3,
            root_milestone="Meristemoid",
            milestones=[milestone],
            copy=True
        )
    except Exception as e:
        print(f"ERROR subsetting tree en milestone {milestone}: {e}")
        continue

    # ----------------------------
    # 2) PLOTS OPCIONALES
    # ----------------------------
    sc.pl.draw_graph(
        adata_Ic,
        color="leiden",
        frameon=True,
        palette=sns.color_palette("Dark2"),
        add_outline=True,
        legend_fontsize=10,
        legend_fontoutline=2,
        show=False
    )

    sc.pl.draw_graph(
        adata_Ic,
        color="milestones",
        palette=sns.color_palette("colorblind"),
        add_outline=True,
        legend_fontsize=10,
        legend_fontoutline=2,
        show=False
    )

    # ----------------------------
    # 3) TEST DE ASOCIACIÓN
    # ----------------------------
    scf.tl.test_association(
        adata_Ic,
        n_jobs=80,
        A_cut=0.3
    )

    # Guardar archivo intermedio
    assoc_file = f"{path}/adata_scfates_{name_file}_{milestone}_association.h5ad"
    adata_Ic.write_h5ad(assoc_file)

    # ----------------------------
    # 4) RECARGAR ARCHIVO + FILTRO DE SIGNIFICANCIA
    # ----------------------------
    adata_Ic = sc.read_h5ad(assoc_file)
    adata_Ic.var["signi"] = adata_Ic.var["p_val"] < 0.001

    # ----------------------------
    # 5) AJUSTE DE TENDENCIAS
    # ----------------------------
    scf.tl.fit(adata_Ic, n_jobs=80)

    fitted_file = f"{path}/adata_scfates_{name_file}_{milestone}_fitted.h5ad"
    adata_Ic.write_h5ad(fitted_file)

    print(f"✔ Milestone '{milestone}' procesado y guardado correctamente.")
```

```{python}
guard_cell = sc.read_h5ad(f"{path}/adata_scfates_{name_file}_Guard_Cell_fitted.h5ad")

```

```{python}

# --- Configuración general ---
milestones = ['Guard_cell']   
output_dir = "./results/pseudotime/tablas" 
os.makedirs(output_dir, exist_ok=True)

# --- Bucle principal ---
for milestone in milestones:
    print(f"\n🔹 Analizando milestone: {milestone}")

    # --- Seleccionar el objeto correspondiente dinámicamente ---
    # Debes tener objetos llamados 'pavement_1' y 'guardcell' ya cargados en memoria
    # Ejemplo: pavement_1 = adata[adata.obs['milestone'] == 'Pavement_1'].copy()
    milestone_obj = globals()[milestone.lower()]   # obtiene el objeto según su nombre en minúsculas

    # --- Calcular correlación de cada gen con pseudotiempo ---
    milestone_obj.var["corr"] = list(map(
        lambda g: pearsonr(
            milestone_obj.obs.t,
            milestone_obj[:, g].layers["fitted"].flatten()
        )[0],
        milestone_obj.var_names
    ))

    # --- Marcar los genes "up" (positivamente correlacionados) ---
    milestone_obj.var["up"] = milestone_obj.var["corr"] > 0

    # --- Guardar tabla ordenada por correlación ---
    corr_out = f"{output_dir}/{milestone}_DEG_ordered_peak_expression.csv"
    milestone_obj.var.sort_values(by="corr", ascending=True).to_csv(corr_out, index=False)
    print(f"Guardado: {corr_out}")

    # --- Calcular fitted expression normalizado ---
    fitted = pd.DataFrame(
        milestone_obj[:, milestone_obj.var_names].layers["fitted"],
        index=milestone_obj.obs_names,
        columns=milestone_obj.var_names
    ).T.copy(deep=True)

    # --- Ordenar genes por punto medio de máxima expresión (pseudotiempo) ---
    feature_order = (
        fitted.apply(
            lambda x: milestone_obj.obs.t[
                ((x - x.min()) / (x.max() - x.min())) > 0.7
            ].mean(),
            axis=1
        ).sort_values().index
    )

    # --- Exportar tabla ordenada por tiempo de pico ---
    df_peak = milestone_obj.var.loc[feature_order]
    peak_out = f"{output_dir}/{milestone}_DEG_ordered_peak_expression2.csv"
    df_peak.to_csv(peak_out)
    print(f"💾 Guardado: {peak_out}")

print("\nAnálisis completado para todos los milestones.")
```

```{python}

base = "/home/mvergara/projects2/eleo/ScRNA"
os.makedirs(f"{base}/figures/trends", exist_ok=True)

# ============================================================
# Guard Cell
# ============================================================
sc.set_figure_params(figsize=(6,20), dpi_save=600, frameon=False)

guard_cell.var["gene"] = guard_cell.var_names
guard_cell.var.index = guard_cell.var['gene']
guard_cell.var_names = guard_cell.var_names.astype(str)
guard_cell.var_names_make_unique()

# scFates retorna LISTA de Axes → capturamos
axes_list = scf.pl.trends(
    guard_cell,
    highlight_features=['AT5G53210', 'AT3G06120','AT3G24140'],
    style="italic",
    add_outline=True,
    basis="dendro",
    show_segs=True,
    fontsize=10,
    figsize=(3,5),
    ordering="max",
    show=False,
    title = "Guard Cell"
)

# Obtener la figura real desde el primer eje
fig = axes_list[0].get_figure()

# Ajustes opcionales
plt.tight_layout()

# Guardar en PDF
fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/guad_heatmap_alone.pdf")
print("Guardado GC OK")

fig, ax = plt.subplots(figsize=(10, 10))

scf.pl.graph(
    guard_cell,
    title="Pavement - Meristemoid - Guard Cell",
    color_cells='annotation_curada_esp',
    palette=sns.color_palette('colorblind'),
    ax=ax,             # <-- clave
    show=True         # <-- no mostrar en pantalla
)


plt.tight_layout()
plt.subplots_adjust(right=0.85)   # deja espacio para la leyenda
fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/Tree_Ms_Gc_60_alone.pdf")

scf.pl.trajectory(
    guard_cell,
    color_seg="t",
    basis="draw_graph_fa",
    frameon=False,
    s=50,
    scale_path=0.6,
    ax=ax,             # <-- clave
    show=True         # <-- no mostrar en pantalla
)

plt.tight_layout()
plt.subplots_adjust(right=0.85)   # deja espacio para la leyenda
fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/Trayectory_Ms_Gc_60_alone.pdf")
```

#### epidermis

```{python, fig.width=10, fig.height=10}

Pc_esp =  sc.read_h5ad('/home/mvergara/projects2/eleo/ScRNA/results/objs/subtipo_epi_ms.h5ad')

# Convert the DataFrame to a NumPy array and assign to X_umap
Pc_esp.obsm['X_umap'] = Pc_esp.obsm['UMAP'].values
Pc_esp.obsm['X_pca'] = Pc_esp.obsm['PCA'].values

sc.set_figure_params(figsize=(10, 10), dpi=300, dpi_save=300)


sc.pl.umap(Pc_esp,  color = 'ident',  add_outline=False, alpha=0.6, size=18,
           frameon=False, legend_fontsize=7, legend_fontoutline=2,
           palette='tab20')

```

```{python}

### Epidermal trajectory
sel = ['Pavement Cell', 'Meristemoid']

for eigen in [8]:
   for nodes in [15]:
     for sigma in [0.2]:
       for plambda in [100]:
             Pc_esp_2 = scfates_trajectories_alignment2(
                  nombre="Pc_esp",
                  adata=Pc_esp,
                  clusters=sel,                 # solo estos clusters
                  root_cluster="Meristemoid",    # raíz en procambium
                  nodes=nodes,
                  sigma=sigma,
                  eigen=eigen,
                  seeed=3,
                  plambda=plambda)

```

```{python, fig.width=15, fig.height=5}

# -------------------------------------------------------------
# 1. Cargar listas de genes
# -------------------------------------------------------------
cluster2 = pd.read_csv("/home/projects2/mvergara/eleo/ScRNA/recursos/karin/Cluster2.txt", sep="\t")
cluster4 = pd.read_csv("/home/projects2/mvergara/eleo/ScRNA/recursos/karin/Cluster4.txt", sep="\t")
cluster4v2 = pd.read_csv("/home/projects2/mvergara/eleo/ScRNA/recursos/karin/cluster4v2.txt", sep="\t")

print(cluster4)
vector_cluster2 = cluster2["Cluster2"].dropna().unique().tolist()
vector_cluster4 = cluster4["Cluster4"].dropna().unique().tolist()
vector_cluster4v2 = cluster4v2["ID"].dropna().unique().tolist()

# -------------------------------------------------------------
# 2. Definir objeto principal
# -------------------------------------------------------------
# Aquí tu objeto AnnData equivalente a pbmc_harmony
adata = obj1  # si ya está cargado

# -------------------------------------------------------------
# 3. Filtrar genes que estén presentes en la matriz
# -------------------------------------------------------------
genes2 = [g for g in vector_cluster2 if g in adata.var_names]
genes4 = [g for g in vector_cluster4 if g in adata.var_names]
genes4v2 = [g for g in vector_cluster4v2 if g in adata.var_names]

# -------------------------------------------------------------
# 4. Función para calcular el módulo
# -------------------------------------------------------------
def compute_module_score(adata, gene_list, prefix):
    X = adata[:, gene_list].X  # matriz de expresión para esos genes
    if not isinstance(X, np.ndarray):
        X = X.toarray()  # convertir de sparse matrix a array
    
    # a. número de genes detectados (>0)
    detected = (X > 0).sum(axis=1)
    
    # b. z-score del número detectado
    detected_z = scale(detected)
    
    # c. promedio de expresión
    mean_expr = X.mean(axis=1)
    
    # d. score ponderado (z × promedio)
    module_score = detected_z * mean_expr

    # e. guardar en .obs
    adata.obs[f"detected_genes_{prefix}"] = detected
    adata.obs[f"detected_genes_{prefix}_z"] = detected_z
    adata.obs[f"module_expr_{prefix}"] = mean_expr
    adata.obs[f"module_score_{prefix}"] = module_score

    print(f"Calculado módulo: {prefix} ({len(gene_list)} genes)")

# -------------------------------------------------------------
# 5. Calcular para cada cluster
# -------------------------------------------------------------
compute_module_score(adata, genes2, "cluster2")
compute_module_score(adata, genes4, "cluster4")
compute_module_score(adata, genes4v2, "cluster4v2")
```

```{python}

def module_sum_counts(adata, gene_list, prefix, use_log=True, lib_norm=True):
    # 1) Filtrar genes presentes
    genes = [g for g in gene_list if g in adata.var_names]
    if not genes:
        print(f"No genes for {prefix}")
        return
    
    # 2) Usar SIEMPRE conteos crudos
    if "counts" in adata.layers:
        M = adata[:, genes].layers["counts"]
        lib = adata.layers["counts"].sum(axis=1)
    else:
        M = adata[:, genes].X
        lib = adata.X.sum(axis=1)

    # Convertir a ndarray denso siempre (cubre sparse y np.matrix)
    if not isinstance(M, np.ndarray):
        M = np.asarray(M.todense())
    
    # FIX PRINCIPAL: asarray + flatten garantiza vector 1D sin importar el tipo
    lib = np.asarray(lib).flatten()

    # 3) Suma del módulo por célula
    raw = np.asarray(M.sum(axis=1)).flatten()

    # 4) Normalización opcional por library size
    if lib_norm:
        # Evitar división por cero
        lib = np.where(lib == 0, 1, lib)
        score = raw / lib
    else:
        score = raw

    # 5) log opcional
    if use_log:
        score = np.log1p(score)

    adata.obs[f"{prefix}_module_score"] = score
    print(f"{prefix}: {len(genes)} genes, score ≥ 0 listo.")


module_sum_counts(adata, genes2, "cluster2")
module_sum_counts(adata, genes4, "cluster4")
module_sum_counts(adata, genes4v2, "cluster4v2")

import matplotlib
matplotlib.use("Agg")  # backend sin pantalla, ponlo al inicio de tu sesión
import matplotlib.pyplot as plt

sc.pl.draw_graph(
    adata,
    color=["cluster2_module_score", "cluster4_module_score", "cluster4v2_module_score"],
    cmap="viridis",
    show=False,   # <-- no intenta plt.show()
    save="_module_scores.png"  # guarda en figures/draw_graph_module_scores.png
)
plt.close()
```

```{python, fig.width=10, fig.height=10}

fig, ax = plt.subplots(figsize=(10, 10))

# scf.pl.graph(
#     Pc_esp_2,
#     title="Pavement - Meristemoid",
#     color_cells='annotation_curada_esp',
#     palette=sns.color_palette('colorblind'),
#     ax=ax,             # <-- clave
#     show=True         # <-- no mostrar en pantalla
# )
# 
# 
# plt.tight_layout()
# plt.subplots_adjust(right=0.85)   # deja espacio para la leyenda
# fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/Tree_Pv_Ms.pdf")
# 
# 
# scf.pl.trajectory(
#     Pc_esp_2,
#     color_seg="t",
#     basis="draw_graph_fa",
#     frameon=False,
#     s=50,
#     scale_path=0.6,
#     ax=ax,             # <-- clave
#     show=True         # <-- no mostrar en pantalla
# )
# 
# plt.tight_layout()
# plt.subplots_adjust(right=0.85)   # deja espacio para la leyenda
# fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/Trayectory_Pv_Ms.pdf")
# 
# 
# sc.set_figure_params(figsize=(5, 5), dpi=300, dpi_save=300)
# fig = sc.pl.draw_graph(
#     adata,
#     color="module_score_cluster4",
#     cmap="viridis",
#     return_fig=True,   # <- clave
#     show=False         # <- evita mostrarla en pantalla
# )
# 
# fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/Cluster4_Pv_Ms.pdf", bbox_inches="tight")

fig = sc.pl.draw_graph(
    adata,
    color=["cluster4_module_score"],
    cmap="viridis",
    return_fig=True,   # <- clave
    show=False  
)
fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/Cluster44_Pv_Ms.pdf", bbox_inches="tight")

# 
# fig = sc.pl.draw_graph(
#     adata,
#     color="AT3G10525",
#     cmap="viridis",
#     return_fig=True,
#     show=True
# )
# 
# fig.axes[0].set_title("LGO", fontsize=16)
# fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/LGO_Pv_Ms.pdf",bbox_inches="tight")
# 
# fig = sc.pl.draw_graph(
#     adata,
#     color="AT4G21750",
#     cmap="viridis",
#     return_fig=True,
#     show=True
# )
# 
# fig.axes[0].set_title("ATML1", fontsize=16)
# fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_22nov/Pseudotime/ATML1_Pv_Ms.pdf",bbox_inches="tight")


```


```{python}

# -------------------------------------------------------------
# 1. Cargar listas de genes
# -------------------------------------------------------------
base_path = "/home/mvergara/projects2/eleo/ScRNA/recursos/karin/Karin_publicados"

files = {
    "NyP_cluster1":    f"{base_path}/GO_anova2vias_NyP_cluster1.txt",
    "NyP_cluster2":    f"{base_path}/GO_anova2vias_NyP_cluster2.txt",
    "NyP_cluster3":    f"{base_path}/GO_anova2vias_NyP_cluster3.txt",
    "NyP_cluster4":    f"{base_path}/GO_anova2vias_NyP_cluster4.txt",
    "P_cluster1":      f"{base_path}/GO_anova2vias_P_cluster1.txt",
    "P_cluster2":      f"{base_path}/GO_anova2vias_P_cluster2.txt",
    "NyP_cluster2_2C": f"{base_path}/NyP_cluster2_2C.txt",
    "NyP_cluster4_4C": f"{base_path}/NyP_cluster4_4C.txt",
    "P_cluster1_4C":   f"{base_path}/P_cluster1_4C.txt",
    "P_cluster2_2C":   f"{base_path}/P_cluster2_2C.txt",
}

gene_lists = {}
for name, path in files.items():
    df = pd.read_csv(path, sep="\t")
    gene_lists[name] = df["ID"].dropna().unique().tolist()
    print(f"{name}: {len(gene_lists[name])} genes cargados")

# -------------------------------------------------------------
# 2. Objeto AnnData
# -------------------------------------------------------------
adata = obj1

# -------------------------------------------------------------
# 3. Filtrar genes presentes en la matriz
# -------------------------------------------------------------
genes_filtered = {}
for name, genes in gene_lists.items():
    filtered = [g for g in genes if g in adata.var_names]
    genes_filtered[name] = filtered
    print(f"{name}: {len(filtered)}/{len(genes)} genes encontrados en adata")

# -------------------------------------------------------------
# 4. Función módulo
# -------------------------------------------------------------
def compute_module_score(adata, gene_list, prefix):
    if not gene_list:
        print(f"Sin genes para {prefix}, saltando.")
        return

    X = adata[:, gene_list].X
    if not isinstance(X, np.ndarray):
        X = np.asarray(X.todense())

    detected     = np.asarray((X > 0).sum(axis=1)).flatten()
    detected_z   = scale(detected.reshape(-1, 1)).flatten()
    mean_expr    = np.asarray(X.mean(axis=1)).flatten()
    module_score = detected_z * mean_expr

    adata.obs[f"detected_genes_{prefix}"]   = detected
    adata.obs[f"detected_genes_{prefix}_z"] = detected_z
    adata.obs[f"module_expr_{prefix}"]      = mean_expr
    adata.obs[f"module_score_{prefix}"]     = module_score

    print(f"Calculado módulo: {prefix} ({len(gene_list)} genes)")

# -------------------------------------------------------------
# 5. Calcular para todos
# -------------------------------------------------------------
for name, genes in genes_filtered.items():
    compute_module_score(adata, genes, name)

# -------------------------------------------------------------
# 6. Visualizar — una figura por grupo (NyP y P)
# -------------------------------------------------------------
groups = {
    "NyP": [k for k in genes_filtered if k.startswith("NyP")],
    "P":   [k for k in genes_filtered if k.startswith("P")],
}

for group_name, keys in groups.items():
    cols = [f"module_score_{k}" for k in keys]
    n = len(cols)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
    if n == 1:
        axes = [axes]
    for ax, col, key in zip(axes, cols, keys):
        sc.pl.draw_graph(adata, color=col, cmap="viridis", ax=ax, show=False, title=key)
    fig.savefig(f"module_scores_{group_name}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Figura guardada: module_scores_{group_name}.png")
    
    
    # -------------------------------------------------------------
# 6. Visualizar — grilla de 3 columnas con todos los módulos
# -------------------------------------------------------------
all_keys = list(genes_filtered.keys())
all_cols = [f"module_score_{k}" for k in all_keys]

n_total = len(all_cols)
n_cols  = 3
n_rows  = int(np.ceil(n_total / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
axes = axes.flatten()

for i, (col, key) in enumerate(zip(all_cols, all_keys)):
    sc.pl.draw_graph(adata, color=col, cmap="viridis", ax=axes[i], show=False, title=key)

# Apagar ejes sobrantes si n_total no es múltiplo de 3
for j in range(n_total, len(axes)):
    axes[j].set_visible(False)

fig.tight_layout()
fig.savefig("module_scores_all.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Figura guardada: module_scores_all.png ({n_rows} filas x {n_cols} columnas)")
```








```{python}
Pc_esp_3 = scfates_trajectories_dendogram(Pc_esp_2)
```

```{python}
### chage milestones names
milestone = pd.DataFrame(Pc_esp_3.obs.groupby('milestones'))[0].tolist()
milestone
new_names_milestone = ['Pavement_3', 'Pavement_2', 'Pavement_4', '12', 'Meristemoid', 'Pavement_1', '8']

#
scf.tl.rename_milestones(Pc_esp_3,new_names_milestone)

fig, ax = plt.subplots(figsize=(10, 10))

scf.pl.graph(
    Pc_esp_3,
    title="Milestone",
    color_cells="milestones",
    palette=sns.color_palette('colorblind'),
    ax=ax,             # <-- clave
    show=True)


plt.tight_layout()
plt.subplots_adjust(right=0.65)   # deja espacio para la leyenda
fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_eleo/Pseudotime/Milestones_Pv_Ms.pdf")
```

```{python}
### REPLACE MILESTONE FOR SOMETHING MORE MEANINGFUL
seg = Pc_esp_3.obs['milestones'].tolist()
seg = ['Guard_Cell' if cell == 'guard_cell' else cell for cell in seg]
seg = ['Meristemoid' if cell == 'meristemoid' else cell for cell in seg]
seg = ['Pavement_1' if cell == 'pavement_1' else cell for cell in seg]
seg = ['Pavement_2' if cell == 'pavement_2' else cell for cell in seg]
seg = ['Pavement_3' if cell == 'pavement_3' else cell for cell in seg]
seg = ['Pavement_4' if cell == 'pavement_4' else cell for cell in seg]
seg = ['Brake' if cell == 'brake' else cell for cell in seg]
Pc_esp_3.obs['milestones'] = seg
```

```{python}
leiden = Pc_esp_3.obs['leiden'].tolist()

leiden = ['Guard_Cell' if cell == 'guard_cell' else cell for cell in leiden]
leiden = ['Meristemoid' if cell == 'meristemoid' else cell for cell in leiden]
leiden = ['Pavement_1' if cell == 'pavement_1' else cell for cell in leiden]
leiden = ['Pavement_2' if cell == 'pavement_2' else cell for cell in leiden]
leiden = ['Pavement_3' if cell == 'pavement_3' else cell for cell in leiden]
leiden = ['Pavement_4' if cell == 'pavement_4' else cell for cell in leiden]
leiden = ['Brake' if cell == 'brake' else cell for cell in leiden]

Pc_esp_3.obs['leiden'] = leiden
```

```{python}

# ============================
# CONFIGURACIÓN GENERAL
# ============================

path = "/home/mvergara/projects2/eleo/ScRNA/results/pseudotime/tablas"
name_file = "all_internal"

os.makedirs(path, exist_ok=True)

# Asegurar que milestones sea categoría — SCFates lo exige
if Pc_esp_3.obs["milestones"].dtype.name != "category":
    Pc_esp_3.obs["milestones"] = Pc_esp_3.obs["milestones"].astype("category")

# Milestones que quieres recorrer
milestones_to_analyze = ['Pavement_1', 'Pavement_2', 'Pavement_3', 'Pavement_4']


# ============================
# LOOP PRINCIPAL
# ============================

for milestone in milestones_to_analyze:

    print(f"\n===== Procesando milestone: {milestone} =====")

    sc.set_figure_params(figsize=(3, 4), dpi_save=600, frameon=False)

    # ----------------------------
    # 1) SUBSET DEL ÁRBOL
    # ----------------------------
    try:
        adata_Ic = scf.tl.subset_tree(
            Pc_esp_3,
            root_milestone="Meristemoid",
            milestones=[milestone],
            copy=True
        )
    except Exception as e:
        print(f"ERROR subsetting tree en milestone {milestone}: {e}")
        continue

    # ----------------------------
    # 2) PLOTS OPCIONALES
    # ----------------------------
    sc.pl.draw_graph(
        adata_Ic,
        color="leiden",
        frameon=True,
        palette=sns.color_palette("Dark2"),
        add_outline=True,
        legend_fontsize=10,
        legend_fontoutline=2,
        show=False
    )

    sc.pl.draw_graph(
        adata_Ic,
        color="milestones",
        palette=sns.color_palette("colorblind"),
        add_outline=True,
        legend_fontsize=10,
        legend_fontoutline=2,
        show=False
    )

    # ----------------------------
    # 3) TEST DE ASOCIACIÓN
    # ----------------------------
    scf.tl.test_association(
        adata_Ic,
        n_jobs=80,
        A_cut=0.3
    )

    # Guardar archivo intermedio
    assoc_file = f"{path}/adata_scfates_{name_file}_{milestone}_association.h5ad"
    adata_Ic.write_h5ad(assoc_file)

    # ----------------------------
    # 4) RECARGAR ARCHIVO + FILTRO DE SIGNIFICANCIA
    # ----------------------------
    adata_Ic = sc.read_h5ad(assoc_file)
    adata_Ic.var["signi"] = adata_Ic.var["p_val"] < 0.001

    # ----------------------------
    # 5) AJUSTE DE TENDENCIAS
    # ----------------------------
    scf.tl.fit(adata_Ic, n_jobs=80)

    fitted_file = f"{path}/adata_scfates_{name_file}_{milestone}_fitted.h5ad"
    adata_Ic.write_h5ad(fitted_file)

    print(f"✔ Milestone '{milestone}' procesado y guardado correctamente.")
```

```{python}
pavement_1 = sc.read_h5ad("/home/mvergara/projects2/eleo/ScRNA/results/archivos_def/adata_scfates_all_internal_Pavement_1_fitted.h5ad")
pavement_2 = sc.read_h5ad(f"{path}/adata_scfates_{name_file}_Pavement_2_fitted.h5ad")
pavement_3 = sc.read_h5ad(f"{path}/adata_scfates_{name_file}_Pavement_3_fitted.h5ad")
pavement_4 = sc.read_h5ad(f"{path}/adata_scfates_{name_file}_Pavement_4_fitted.h5ad")
```

```{python}

# --- Configuración general ---
milestones = ['Pavement_1', 'Pavement_2', 'Pavement_3', 'Pavement_4']   # lista de milestones a analizar
output_dir = "./results/pseudotime/tablas"  # carpeta de salida
os.makedirs(output_dir, exist_ok=True)

# --- Bucle principal ---
for milestone in milestones:
    print(f"\nAnalizando milestone: {milestone}")

    # --- Seleccionar el objeto correspondiente dinámicamente ---
    # Debes tener objetos llamados 'pavement_1' y 'guardcell' ya cargados en memoria
    # Ejemplo: pavement_1 = adata[adata.obs['milestone'] == 'Pavement_1'].copy()
    milestone_obj = globals()[milestone.lower()]   # obtiene el objeto según su nombre en minúsculas

    # --- Calcular correlación de cada gen con pseudotiempo ---
    milestone_obj.var["corr"] = list(map(
        lambda g: pearsonr(
            milestone_obj.obs.t,
            milestone_obj[:, g].layers["fitted"].flatten()
        )[0],
        milestone_obj.var_names
    ))

    # --- Marcar los genes "up" (positivamente correlacionados) ---
    milestone_obj.var["up"] = milestone_obj.var["corr"] > 0

    # --- Guardar tabla ordenada por correlación ---
    corr_out = f"{output_dir}/{milestone}_DEG_ordered_peak_expression.csv"
    milestone_obj.var.sort_values(by="corr", ascending=True).to_csv(corr_out, index=False)
    print(f"Guardado: {corr_out}")

    # --- Calcular fitted expression normalizado ---
    fitted = pd.DataFrame(
        milestone_obj[:, milestone_obj.var_names].layers["fitted"],
        index=milestone_obj.obs_names,
        columns=milestone_obj.var_names
    ).T.copy(deep=True)

    # --- Ordenar genes por punto medio de máxima expresión (pseudotiempo) ---
    feature_order = (
        fitted.apply(
            lambda x: milestone_obj.obs.t[
                ((x - x.min()) / (x.max() - x.min())) > 0.7
            ].mean(),
            axis=1
        ).sort_values().index
    )

    # --- Exportar tabla ordenada por tiempo de pico ---
    df_peak = milestone_obj.var.loc[feature_order]
    peak_out = f"{output_dir}/{milestone}_DEG_ordered_peak_expression2.csv"
    df_peak.to_csv(peak_out)
    print(f"Guardado: {peak_out}")

print("\nAnálisis completado para todos los milestones.")
```

```{python}

base = "/home/mvergara/projects2/eleo/ScRNA"
os.makedirs(f"{base}/figures/trends", exist_ok=True)

# ============================================================
# PAVEMENT 1
# ============================================================
sc.set_figure_params(figsize=(6,20), dpi_save=600, frameon=False)

pavement_1.var["gene"] = pavement_1.var_names
pavement_1.var.index = pavement_1.var['gene']
pavement_1.var_names = pavement_1.var_names.astype(str)
pavement_1.var_names_make_unique()

# scFates retorna LISTA de Axes → capturamos
axes_list = scf.pl.trends(
    pavement_1,
    style="italic",
    add_outline=True,
    basis="dendro",
    show_segs=True,
    fontsize=10,
    figsize=(3,5),
    ordering="max",
    show=False,
    title = "Pavement 1"
)

# Obtener la figura real desde el primer eje
fig = axes_list[0].get_figure()

# Ajustes opcionales
plt.tight_layout()

# Guardar en PDF
fig.savefig(f"{base}/figures/trends/pavement_1_heatmap.pdf")
print("Guardado P1 OK")
```

```{python}

base = "/home/mvergara/projects2/eleo/ScRNA"
os.makedirs(f"{base}/figures/trends", exist_ok=True)

# ============================================================
# PAVEMENT 2
# ============================================================
sc.set_figure_params(figsize=(6,20), dpi_save=600, frameon=False)

pavement_2.var["gene"] = pavement_2.var_names
pavement_2.var.index = pavement_2.var['gene']
pavement_2.var_names = pavement_2.var_names.astype(str)
pavement_2.var_names_make_unique()

# scFates retorna LISTA de Axes → capturamos
axes_list = scf.pl.trends(
    pavement_2,
    style="italic",
    add_outline=True,
    basis="dendro",
    show_segs=True,
    fontsize=10,
    figsize=(3,5),
    ordering="max",
    show=False,
    title = "Pavement 2"
)

# Obtener la figura real desde el primer eje
fig = axes_list[0].get_figure()

# Ajustes opcionales
plt.tight_layout()

# Guardar en PDF
fig.savefig(f"{base}/figures/trends/pavement_2_heatmap.pdf")
print("Guardado P2 OK")
```

```{python}

base = "/home/mvergara/projects2/eleo/ScRNA"
os.makedirs(f"{base}/figures/trends", exist_ok=True)

# ============================================================
# PAVEMENT 3
# ============================================================
sc.set_figure_params(figsize=(6,20), dpi_save=600, frameon=False)

pavement_3.var["gene"] = pavement_3.var_names
pavement_3.var.index = pavement_3.var['gene']
pavement_3.var_names = pavement_3.var_names.astype(str)
pavement_3.var_names_make_unique()

# scFates retorna LISTA de Axes → capturamos
axes_list = scf.pl.trends(
    pavement_3,
    style="italic",
    add_outline=True,
    basis="dendro",
    show_segs=True,
    fontsize=10,
    figsize=(3,5),
    ordering="max",
    show=False,
    title = "Pavement 3"
)

# Obtener la figura real desde el primer eje
fig = axes_list[0].get_figure()

# Ajustes opcionales
plt.tight_layout()

# Guardar en PDF
fig.savefig(f"{base}/figures/trends/pavement_3_heatmap.pdf")
print("Guardado P3 OK")
```

```{python}

base = "/home/mvergara/projects2/eleo/ScRNA"
os.makedirs(f"{base}/figures/trends", exist_ok=True)

# ============================================================
# PAVEMENT 4
# ============================================================
sc.set_figure_params(figsize=(6,20), dpi_save=600, frameon=False)

pavement_4.var["gene"] = pavement_4.var_names
pavement_4.var.index = pavement_4.var['gene']
pavement_4.var_names = pavement_4.var_names.astype(str)
pavement_4.var_names_make_unique()

# scFates retorna LISTA de Axes → capturamos
axes_list = scf.pl.trends(
    pavement_4,
    style="italic",
    add_outline=True,
    basis="dendro",
    show_segs=True,
    fontsize=10,
    figsize=(3,5),
    ordering="max",
    show=False,
    title = "Pavement 4"
)

# Obtener la figura real desde el primer eje
fig = axes_list[0].get_figure()

# Ajustes opcionales
plt.tight_layout()

# Guardar en PDF
fig.savefig(f"{base}/figures/trends/pavement_4_heatmap.pdf")
print("Guardado P4 OK")
```

########### Guard Cells.

```{python, fig.width=10, fig.height=10}

Pc_esp =  sc.read_h5ad('/home/mvergara/projects2/eleo/ScRNA/results/objs/subtipo_gc_ms.h5ad')

# Convert the DataFrame to a NumPy array and assign to X_umap
Pc_esp.obsm['X_umap'] = Pc_esp.obsm['UMAP'].values
Pc_esp.obsm['X_pca'] = Pc_esp.obsm['PCA'].values

sc.set_figure_params(figsize=(10, 10), dpi=300, dpi_save=300)


sc.pl.umap(Pc_esp,  color = 'ident',  add_outline=False, alpha=0.6, size=18,
           frameon=False, legend_fontsize=7, legend_fontoutline=2,
           palette='tab20')

```

```{python}

### Epidermal trajectory
sel = ['Meristemoid', 'Guard Cell']
for eigen in [8]:
   for nodes in [15]:
     for sigma in [0.2]:
       for plambda in [100]:
             Pc_esp_2 = scfates_trajectories_alignment2(
                  nombre="Pc_esp",
                  adata=Pc_esp,
                  clusters=sel,                 # solo estos clusters
                  root_cluster="Meristemoid",    # raíz en procambium
                  nodes=nodes,
                  sigma=sigma,
                  eigen=eigen,
                  seeed=3,
                  plambda=plambda)

```

```{python, fig.width=10, fig.height=10}

fig, ax = plt.subplots(figsize=(10, 10))

scf.pl.graph(
    Pc_esp_2,
    title="Guard Cell - Meristemoid",
    color_cells='annotation_curada_esp',
    palette=sns.color_palette('colorblind'),
    ax=ax,             # <-- clave
    show=True         # <-- no mostrar en pantalla
)


plt.tight_layout()
plt.subplots_adjust(right=0.85)   # deja espacio para la leyenda
fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_eleo/Pseudotime/Tree_Gc_Ms.pdf")


scf.pl.trajectory(
    Pc_esp_2,
    color_seg="t",
    basis="draw_graph_fa",
    frameon=False,
    s=50,
    scale_path=0.6,
    ax=ax,             # <-- clave
    show=True         # <-- no mostrar en pantalla
)

plt.tight_layout()
plt.subplots_adjust(right=0.85)   # deja espacio para la leyenda
fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_eleo/Pseudotime/Trayectory_Gc_Ms.pdf")





### FAMA
fig = sc.pl.draw_graph(
    Pc_esp_2,
    color="AT3G24140",
    title="FMA",
    cmap="viridis",
    return_fig=True,
    show=True
)

fig.axes[0].set_title("FMA", fontsize=16)
fig.savefig(
    "/home/mvergara/projects2/eleo/ScRNA/figuras_eleo/Pseudotime/FMA_Gc_Ms.pdf",
    bbox_inches="tight"
)

### MUTE

fig = sc.pl.draw_graph(
    Pc_esp_2,
    color="AT3G06120",
    title="MUTE",
    cmap="viridis",
    return_fig=True,
    show=True
)

fig.axes[0].set_title("MUTE", fontsize=16)
fig.savefig(
    "/home/mvergara/projects2/eleo/ScRNA/figuras_eleo/Pseudotime/MUTE_Gc_Ms.pdf",
    bbox_inches="tight"
)

### SPCH

fig = sc.pl.draw_graph(
    Pc_esp_2,
    color="AT5G53210",
    title="SPEECHLESS",
    cmap="viridis",
    return_fig=True,
    show=True
)

fig.axes[0].set_title("SPEECHLESS", fontsize=16)
fig.savefig(
    "/home/mvergara/projects2/eleo/ScRNA/figuras_eleo/Pseudotime/SPEECHLESS_Gc_Ms.pdf",
    bbox_inches="tight"
)

```

```{python}
Pc_esp_3 = scfates_trajectories_dendogram(Pc_esp_2)
```

```{python}
### chage milestones names
milestone = pd.DataFrame(Pc_esp_3.obs.groupby('milestones'))[0].tolist()
milestone
new_names_milestone = ['Guard_Cell', 'Meristemoid']

#
scf.tl.rename_milestones(Pc_esp_3,new_names_milestone)

fig, ax = plt.subplots(figsize=(10, 10))

scf.pl.graph(
    Pc_esp_3,
    title="Milestone",
    color_cells="milestones",
    palette=sns.color_palette('colorblind'),
    ax=ax,             # <-- clave
    show=True)


#plt.tight_layout()
#plt.subplots_adjust(right=0.65)   # deja espacio para la leyenda
#fig.savefig("/home/mvergara/projects2/eleo/ScRNA/figuras_eleo/Pseudotime/Milestones_Gc_Ms.pdf")


```

```{python}
### REPLACE MILESTONE FOR SOMETHING MORE MEANINGFUL
seg = Pc_esp_3.obs['milestones'].tolist()
seg = ['Guard_Cell' if cell == 'guard_cell' else cell for cell in seg]
seg = ['Meristemoid' if cell == 'meristemoid' else cell for cell in seg]
seg = ['Pavement_1' if cell == 'pavement_1' else cell for cell in seg]
seg = ['Pavement_2' if cell == 'pavement_2' else cell for cell in seg]
seg = ['Pavement_3' if cell == 'pavement_3' else cell for cell in seg]
seg = ['Pavement_4' if cell == 'pavement_4' else cell for cell in seg]
seg = ['Brake' if cell == 'brake' else cell for cell in seg]
Pc_esp_3.obs['milestones'] = seg
```

```{python}
leiden = Pc_esp_3.obs['leiden'].tolist()

leiden = ['Guard_Cell' if cell == 'guard_cell' else cell for cell in leiden]
leiden = ['Meristemoid' if cell == 'meristemoid' else cell for cell in leiden]
leiden = ['Pavement_1' if cell == 'pavement_1' else cell for cell in leiden]
leiden = ['Pavement_2' if cell == 'pavement_2' else cell for cell in leiden]
leiden = ['Pavement_3' if cell == 'pavement_3' else cell for cell in leiden]
leiden = ['Pavement_4' if cell == 'pavement_4' else cell for cell in leiden]
leiden = ['Brake' if cell == 'brake' else cell for cell in leiden]

Pc_esp_3.obs['leiden'] = leiden
```

```{python}

path = "/home/mvergara/projects2/eleo/ScRNA/results/pseudotime/tablas/"
name_file = "Gc_ms_FULL"

scf.tl.test_association(Pc_esp_3, n_jobs=80, A_cut=0.3)
Pc_esp_3.var["signi"] = Pc_esp_3.var["p_val"] < 0.001
scf.tl.fit(Pc_esp_3, n_jobs=80)
```

```{python}

adata = Pc_esp_3
output_dir = "./results/pseudotime/tablas"
os.makedirs(output_dir, exist_ok=True)

adata.var["corr"] = [
    pearsonr(adata.obs.t, adata[:, g].layers["fitted"].flatten())[0]
    for g in adata.var_names
]

adata.var["up"] = adata.var["corr"] > 0
adata.var.sort_values("corr").to_csv(f"{output_dir}/Guard_Cell_DEG_ordered_correlation.csv")

fitted = pd.DataFrame(
    adata.layers["fitted"].T,
    index=adata.var_names,
    columns=adata.obs_names
)

feature_order = (
    fitted.apply(
        lambda x: adata.obs.loc[x.index, "t"][((x - x.min()) / (x.max() - x.min())) > 0.7].mean()
        if x.max() != x.min() else np.nan,
        axis=1
    ).sort_values().index
)

adata.var.loc[feature_order].to_csv(f"{output_dir}/Guard_Cell_DEG_ordered_peak_expression.csv")
```

```{python}
sc.set_figure_params(figsize=(6,20),dpi_save=600,frameon=False)
adata.var["gene"] = adata.var_names
adata.var.index = adata.var['gene']
adata.var_names = adata.var_names.astype(str)
adata.var_names_make_unique()

# scFates retorna LISTA de Axes → capturamos
axes_list = scf.pl.trends(
    adata,
    highlight_features=['AT5G53210', 'AT3G06120','AT3G24140'],
    style="italic",
    add_outline=True,
    basis="dendro",
    show_segs=True,
    fontsize=10,
    figsize=(3,5),
    ordering="max",
    show=False,
    title = "Pavement 4"
)

# Obtener la figura real desde el primer eje
fig = axes_list[0].get_figure()

# Ajustes opcionales
plt.tight_layout()

# Guardar en PDF
fig.savefig(f"{base}/figures/trends/Guard_Cell_heatmap.pdf")
print("Guardado GC OK")




```

















```{python}

# ============================
# CONFIGURACIÓN GENERAL
# ============================

path = "/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime"
name_file = "all_internal"

os.makedirs(path, exist_ok=True)

# Asegurar que milestones sea categoría — SCFates lo exige
if obj1.obs["milestones"].dtype.name != "category":
    obj1.obs["milestones"] = obj1.obs["milestones"].astype("category")

# Milestones que quieres recorrer
milestones_to_analyze = ['Guard_Cell']


# ============================
# LOOP PRINCIPAL
# ============================

for milestone in milestones_to_analyze:

    print(f"\n===== Procesando milestone: {milestone} =====")

    sc.set_figure_params(figsize=(10, 10), dpi_save=600, frameon=False)

    # ----------------------------
    # 1) SUBSET DEL ÁRBOL
    # ----------------------------
    try:
        adata_Ic = scf.tl.subset_tree(
            obj1,
            root_milestone="Meristemoid",
            milestones=[milestone],
            copy=True
        )
    except Exception as e:
        print(f"ERROR subsetting tree en milestone {milestone}: {e}")
        continue

    sc.pl.draw_graph(
        adata_Ic,
        color="milestones",
        palette=sns.color_palette("colorblind"),
        add_outline=True,
        legend_fontsize=10,
        legend_fontoutline=2,
        show=True
    )

    # ----------------------------
    # 3) TEST DE ASOCIACIÓN
    # ----------------------------
    scf.tl.test_association(
        adata_Ic,
        n_jobs=80,
        A_cut=0.3
    )

    # Guardar archivo intermedio
    assoc_file = f"{path}/adata_scfates_{name_file}_{milestone}_association.h5ad"
    adata_Ic.write_h5ad(assoc_file)

    # ----------------------------
    # 4) RECARGAR ARCHIVO + FILTRO DE SIGNIFICANCIA
    # ----------------------------
    adata_Ic = sc.read_h5ad(assoc_file)
    adata_Ic.var["signi"] = adata_Ic.var["p_val"] < 0.001

    # ----------------------------
    # 5) AJUSTE DE TENDENCIAS
    # ----------------------------
    scf.tl.fit(adata_Ic, n_jobs=80)

    fitted_file = f"{path}/adata_scfates_{name_file}_{milestone}_fitted.h5ad"
    adata_Ic.write_h5ad(fitted_file)

    print(f"✔ Milestone '{milestone}' procesado y guardado correctamente.")
```



```{python, fig.width=10, fig.height=10}


pv1 = sc.read_h5ad("/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime/adata_scfates_all_internal_Pavement_1_fitted.h5ad")
pv2 = sc.read_h5ad("/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime/adata_scfates_all_internal_Pavement_2_fitted.h5ad")
pv3 = sc.read_h5ad("/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime/adata_scfates_all_internal_Pavement_3_fitted.h5ad")
pv4 = sc.read_h5ad("/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime/adata_scfates_all_internal_Pavement_4_fitted.h5ad")
gc = sc.read_h5ad("/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime/adata_scfates_all_internal_Guard_Cell_fitted.h5ad")
b1 = sc.read_h5ad("/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime/adata_scfates_all_internal_Brake_1_fitted.h5ad")
b2 = sc.read_h5ad("/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime/adata_scfates_all_internal_Brake_2_fitted.h5ad")


def graficos(adata, nombre, outdir="/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime"):
    """
    Genera y guarda dos gráficos independientes del árbol:
    - coloreado por milestones
    - coloreado por leiden
    """

    # ============================
    # 1) MILESTONES
    # ============================
    fig, ax = plt.subplots(figsize=(10, 10))

    scf.pl.graph(
        adata,
        color_cells="milestones",
        ax=ax,
        show=True
    )

    plt.tight_layout()
    #plt.subplots_adjust(right=0.65)

    fig.savefig(f"{outdir}/Milestones_{nombre}.pdf")
    plt.close(fig)

    # ============================
    # 2) LEIDEN
    # ============================
    fig, ax = plt.subplots(figsize=(10, 10))

    scf.pl.graph(
        adata,
        color_cells="leiden",
        ax=ax,
        show=True
    )

    plt.tight_layout()
    #plt.subplots_adjust(right=0.65)

    fig.savefig(f"{outdir}/Leiden_{nombre}.pdf")
    plt.close(fig)


  
graficos(pv1, "Pavement_1")
graficos(pv2, "Pavement_2")
graficos(pv3, "Pavement_3")
graficos(pv4, "Pavement_4")
graficos(gc, "Guard_Cell")
graficos(b1, "Brake_1")
graficos(b2, "Brake_2")




```

```{python}


def genes_by_pseudotime_peak(
    adata,
    milestone_name,
    output_dir,
    t_key="t",
    leiden_key="leiden",
    layer_key="fitted",
    peak_threshold=0.7
):
    """
    Calcula correlación con pseudotime, tiempo de peak de expresión
    y Leiden dominante en el peak. Exporta tabla ordenada por peak_t.

    Parameters
    ----------
    adata : AnnData
        Objeto con pseudotime y expresión fitted
    milestone_name : str
        Nombre del milestone (ej. pv1, gc, b1, etc.)
    output_dir : str
        Directorio de salida
    t_key : str
        Columna en adata.obs con pseudotime
    leiden_key : str
        Columna en adata.obs con clusters Leiden
    layer_key : str
        Capa con expresión fitted (cells × genes)
    peak_threshold : float
        Umbral (0–1) para definir peak de expresión

    Returns
    -------
    df_final : pandas.DataFrame
        Tabla ordenada por peak_t
    """

    os.makedirs(output_dir, exist_ok=True)

    # ============================
    # 1) MATRICES BASE
    # ============================
    t = adata.obs[t_key].values
    leiden = adata.obs[leiden_key].values
    X = adata.layers[layer_key]

    # ============================
    # 2) CORRELACIÓN CON PSEUDOTIME
    # ============================
    t_c = t - t.mean()
    X_c = X - X.mean(axis=0)

    corr = (t_c @ X_c) / (
        np.sqrt((t_c ** 2).sum()) *
        np.sqrt((X_c ** 2).sum(axis=0))
    )

    adata.var["corr"] = corr
    adata.var["up"] = corr > 0

    # ============================
    # 3) NORMALIZAR EXPRESIÓN (0–1)
    # ============================
    X_min = X.min(axis=0)
    X_max = X.max(axis=0)
    X_norm = (X - X_min) / (X_max - X_min + 1e-9)

    # ============================
    # 4) PEAK TIME
    # ============================
    mask = X_norm > peak_threshold

    peak_t = np.full(X.shape[1], np.nan)
    for g in range(X.shape[1]):
        if mask[:, g].any():
            peak_t[g] = t[mask[:, g]].mean()

    adata.var["peak_t"] = peak_t

    # ============================
    # 5) LEIDEN DOMINANTE DEL PEAK
    # ============================
    peak_leiden = np.full(X.shape[1], np.nan, dtype=object)
    for g in range(X.shape[1]):
        if mask[:, g].any():
            peak_leiden[g] = (
                pd.Series(leiden[mask[:, g]])
                .value_counts()
                .idxmax()
            )

    adata.var["peak_leiden"] = peak_leiden

    # ============================
    # 6) ORDEN FINAL
    # ============================
    df_final = (
        adata.var
        .assign(_order=adata.var["peak_t"].fillna(2))
        .sort_values("_order")
        .drop(columns="_order")
    )

    # ============================
    # 7) EXPORTAR
    # ============================
    out_file = os.path.join(
        output_dir,
        f"{milestone_name}_genes_ordered_by_peak_with_leiden.csv"
    )
    df_final.to_csv(out_file)

    print(f"✔ Archivo generado: {out_file}")

    return df_final
  
  
import scanpy as sc

output_dir = "/home/mvergara/projects2/eleo/ScRNA/results/files_eleo_22dic_2025/Pseudotime"

genes_by_pseudotime_peak(pv1, "pv1", output_dir)
genes_by_pseudotime_peak(pv2, "pv2", output_dir)
genes_by_pseudotime_peak(pv3, "pv3", output_dir)
genes_by_pseudotime_peak(pv4, "pv4", output_dir)
genes_by_pseudotime_peak(gc,  "gc",  output_dir)
genes_by_pseudotime_peak(b1,  "b1",  output_dir)
genes_by_pseudotime_peak(b2,  "b2",  output_dir)
```





