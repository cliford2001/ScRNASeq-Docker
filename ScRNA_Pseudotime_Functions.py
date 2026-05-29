# =============================================================================
# ScRNA_Pseudotime_Functions.py — Pseudotime analysis helper functions
# =============================================================================
# Companion to capitulo3_pseudotime.py.
# Loaded automatically at the start of that script.


def trajectory_run(name, nodes=50, sigma=0.2, lambda_value=60, eigs=20):
    return {
        "name": name,
        "nodes": nodes,
        "sigma": sigma,
        "lambda": lambda_value,
        "eigs": eigs,
    }


# =============================================================================
# load_curated_object
# =============================================================================
# Notebook wrapper for Step 25. Loads the full curated object, fixes Seurat-style
# coordinate names, plots the overview UMAP, and prints available cell types.
def load_curated_object(input_h5ad, dir_pseudotime, annotation_col, n_jobs=4):

    print(f"Loading object: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    adata.obsm["X_umap"] = adata.obsm["UMAP"].values
    adata.obsm["X_pca"] = adata.obsm["PCA"].values

    sc.settings.figdir = dir_pseudotime
    sc.set_figure_params(figsize=(10, 8), dpi=80, dpi_save=300)

    fig = sc.pl.umap(
        adata,
        color              = annotation_col,
        legend_loc         = "on data",
        legend_fontsize    = 9,
        legend_fontoutline = 3,
        frameon            = False,
        show               = False,
        return_fig         = True,
    )
    fig.savefig(
        os.path.join(dir_pseudotime, "umap_overview.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(fig)

    counts = adata.obs[annotation_col].value_counts()
    print(f"\nCell types in '{annotation_col}':\n")
    for i, ct in enumerate(sorted(adata.obs[annotation_col].unique()), 1):
        print(f"  {i:2d}.  {ct:<30s}  {counts[ct]:>5d} cells")

    print("\nSTEP 25 COMPLETE: full object loaded")
    return adata, n_jobs


# =============================================================================
# preview_trajectory_selection
# =============================================================================
# Notebook wrapper for Step 26. Subsets the selected cell types and displays a
# UMAP preview so the user can decide whether the selection makes biological sense.
def preview_trajectory_selection(adata, clusters, annotation_col, dir_pseudotime):
    missing = [x for x in clusters if x not in set(adata.obs[annotation_col])]
    if missing:
        raise ValueError(
            "These cell types were not found in "
            f"'{annotation_col}': {missing}. Copy names exactly from Step 25."
        )

    adata_sub = adata[adata.obs[annotation_col].isin(clusters)].copy()
    print(f"Selected {len(adata_sub)} cells: {clusters}")

    fig = sc.pl.umap(
        adata_sub,
        color              = annotation_col,
        legend_loc         = "on data",
        legend_fontsize    = 9,
        legend_fontoutline = 3,
        frameon            = False,
        show               = False,
        return_fig         = True,
    )
    fig.savefig(
        os.path.join(dir_pseudotime, "umap_selection.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(fig)

    print("\nSTEP 26 COMPLETE: continue only if this UMAP looks biologically sensible")
    return adata_sub


# =============================================================================
# build_pseudotime_trajectory
# =============================================================================
# Full pipeline: subset cells → Palantir diffusion maps → scFates tree →
# root selection → pseudotime assignment.
# Returns adata with trajectory and pseudotime stored.
def build_pseudotime_trajectory(
    adata,
    clusters,           # list of cell types to include (values from annotation_col)
    root_cluster,       # cell type where pseudotime = 0 (the biological progenitor)
    annotation_col,     # adata.obs column containing cell type labels
    nodes       = 150,  # tree nodes: more = finer branches, slower to compute
    sigma       = 0.2,  # smoothing: lower = tree follows cells more tightly
    ppt_lambda  = 60,   # complexity: higher = simpler tree with fewer branches
    n_components= 50,   # total diffusion map dimensions computed by Palantir
    n_eigs      = 20,   # dimensions retained (must be strictly < n_components)
    n_neighbors = 50,   # neighbors for final graph (higher = smoother layout)
    seed        = 3,    # random seed for reproducibility across runs
):
    # Subset to selected cell types
    if clusters:
        adata = adata[adata.obs[annotation_col].isin(clusters)].copy()

    # Initial neighbors + Leiden (used internally by scFates for graph structure)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=0.5)

    # Palantir diffusion maps: captures global developmental geometry
    pca_proj = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)
    dm_res   = palantir.utils.run_diffusion_maps(pca_proj, n_components=n_components)
    ms_data  = palantir.utils.determine_multiscale_space(
        dm_res, n_eigs=min(n_eigs, n_components - 1)
    )
    adata.obsm["X_palantir"] = ms_data.values

    # Graph in Palantir space for tree construction
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_palantir", method="umap")
    adata.obsm["X_pca2d"] = adata.obsm["X_pca"][:, :2]
    sc.tl.draw_graph(adata, init_pos="X_pca2d")

    # Build principal tree (PPT = Principal Polynomial Tree)
    scf.tl.tree(
        adata,
        method     = "ppt",
        Nodes      = nodes,
        use_rep    = "palantir",
        plot       = False,
        device     = "cpu",
        seed       = seed,
        ppt_lambda = ppt_lambda,
        ppt_nsteps = 200,
        ppt_sigma  = sigma,
    )

    # Root selection: most-connected cell within root_cluster
    mask          = (adata.obs[annotation_col] == root_cluster).to_numpy()
    cluster_cells = adata.obs_names[mask]
    sub_conn      = adata.obsp["connectivities"][mask][:, mask]
    degrees       = np.array(sub_conn.sum(axis=1)).flatten()
    root_cell     = cluster_cells[degrees.argmax()]

    adata.obs["is_root"] = adata.obs_names == root_cell
    scf.tl.root(adata, root="is_root")
    scf.tl.pseudotime(adata)

    print(f"✓ Trajectory built — root cell: {root_cell}")
    return adata


# =============================================================================
# build_dendrogram
# =============================================================================
# Adds a dendrogram to an adata object that already has a scFates tree.
# Call this after build_pseudotime_trajectory.
def build_dendrogram(adata):
    scf.tl.dendrogram(adata)
    return adata


# =============================================================================
# plot_trajectory_graphs
# =============================================================================
# Saves force-directed tree graphs colored by annotation and leiden clusters.
# When show_inline=True, displays only the annotation graph in Jupyter so the
# notebook stays readable while all graph files are still saved.
def plot_trajectory_graphs(
    adata,
    name,
    output_dir,
    annotation_col,
    show_inline=False,
    inline_color=None,
):
    os.makedirs(output_dir, exist_ok=True)
    inline_color = inline_color or annotation_col

    for color_by, suffix in [(annotation_col, "annotation"), ("leiden", "leiden")]:
        fig, ax = plt.subplots(figsize=(10, 10))
        scf.pl.graph(adata, color_cells=color_by, ax=ax, show=False)
        plt.tight_layout()

        for ext in ["pdf", "png"]:
            fig.savefig(
                os.path.join(output_dir, f"{name}_{suffix}.{ext}"),
                dpi=300,
                bbox_inches="tight",
            )

        if show_inline and color_by == inline_color:
            plt.show()

        plt.close(fig)

    print(f"✓ Trajectory graphs saved to {output_dir}")


# =============================================================================
# export_pseudotime_table
# =============================================================================
def export_pseudotime_table(adata, name, output_dir, annotation_col=None):
    os.makedirs(output_dir, exist_ok=True)
    cell_time = adata.obs[["t"]].copy()
    cell_time.index.name = "cell_id"
    if annotation_col is not None and annotation_col in adata.obs.columns:
        cell_time[annotation_col] = adata.obs[annotation_col]
    out_file = os.path.join(output_dir, f"{name}_pseudotime_by_cell.tsv")
    cell_time.to_csv(out_file, sep="	")
    print(f"Pseudotime table saved: {out_file}")
    return out_file


# =============================================================================
# plot_pseudotime_trajectory
# =============================================================================
# Overlay the principal graph and pseudotime trajectory on the FA layout.
def plot_pseudotime_trajectory(
    adata,
    name,
    output_dir,
    annotation_col,
    show_inline=False,
):
    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 10))
    scf.pl.graph(
        adata,
        basis       = "draw_graph_fa",
        color_cells = annotation_col,
        ax          = ax,
        show        = False,
    )
    scf.pl.trajectory(
        adata,
        color_seg  = "t",
        basis      = "draw_graph_fa",
        frameon    = False,
        s          = 50,
        scale_path = 0.6,
        ax         = ax,
        show       = False,
    )
    plt.tight_layout()
    fig.savefig(
        os.path.join(output_dir, f"{name}_pseudotime_trajectory.png"),
        dpi=600,
        bbox_inches="tight",
    )

    if show_inline:
        plt.show()

    plt.close(fig)
    print(f"Pseudotime trajectory saved to {output_dir}")


# =============================================================================
# plot_root_cell
# =============================================================================
# QC plot: show the automatically selected root cell as a red point on FA1/FA2.
def plot_root_cell(adata, name, output_dir, show_inline=False):
    os.makedirs(output_dir, exist_ok=True)
    if "is_root" not in adata.obs.columns:
        print("Root-cell plot skipped: is_root not found.")
        return None

    root_cells = adata.obs_names[adata.obs["is_root"]]
    if len(root_cells) == 0:
        print("Root-cell plot skipped: no root cell marked.")
        return None

    root_cell = root_cells[0]
    plot_data = adata.copy()
    plot_data.obs["root_cell"] = np.where(plot_data.obs["is_root"], "Root cell", "Other cells")
    plot_data = plot_data[plot_data.obs.sort_values("root_cell").index, :]
    plot_data.uns["root_cell_colors"] = ["lightgray", "red"]

    fig = sc.pl.embedding(
        plot_data,
        basis              = "draw_graph_fa",
        color              = "root_cell",
        legend_loc         = "right margin",
        title              = f"Root cell: {root_cell}",
        frameon            = False,
        size               = 26,
        alpha              = 0.9,
        show               = False,
        return_fig         = True,
    )
    out_file = os.path.join(output_dir, f"{name}_root_cell.png")
    fig.savefig(out_file, dpi=600, bbox_inches="tight")

    if show_inline:
        plt.show()

    plt.close(fig)
    print(f"Root-cell plot saved: {out_file}")
    return root_cell


# =============================================================================
# run_trajectory_runs
# =============================================================================
# Convenience wrapper for the notebook. Runs one or more trajectory parameter
# sets, saves each result in a separate folder, shows the graph, and returns the
# selected trajectory for downstream steps.
def run_trajectory_runs(
    adata,
    clusters,
    root_cluster,
    annotation_col,
    output_base_dir,
    runs,
    selected_run=None,
    show_inline=True,
):
    trajectory_runs = {}

    for run in runs:
        run_name = run.get("name", run.get("id"))
        if not run_name:
            raise ValueError("Each trajectory run needs a 'name'.")

        run_dir = os.path.join(output_base_dir, "trajectory", run_name)
        os.makedirs(run_dir, exist_ok=True)

        nodes       = run.get("nodes", 50)
        sigma       = run.get("sigma", 0.2)
        ppt_lambda  = run.get("lambda", run.get("ppt_lambda", 60))
        n_eigs      = run.get("eigs", run.get("n_eigs", 20))
        n_neighbors = run.get("n_neighbors", 50)
        seed        = run.get("seed", 3)

        print(f"\nRunning: {run_name}")
        print(f"Results folder: {run_dir}")

        adata_run = build_pseudotime_trajectory(
            adata          = adata,
            clusters       = clusters,
            root_cluster   = root_cluster,
            annotation_col = annotation_col,
            nodes          = nodes,
            sigma          = sigma,
            ppt_lambda     = ppt_lambda,
            n_eigs         = n_eigs,
            n_components   = 50,
            n_neighbors    = n_neighbors,
            seed           = seed,
        )

        params = dict(run)
        params.update({"root_cluster": root_cluster})
        adata_run.write_h5ad(os.path.join(run_dir, f"{run_name}_trajectory.h5ad"))
        pd.Series(params).to_csv(os.path.join(run_dir, "parameters.tsv"), sep="\t", header=False)

        plot_trajectory_graphs(
            adata          = adata_run,
            name           = run_name,
            output_dir     = run_dir,
            annotation_col = annotation_col,
            show_inline    = show_inline,
        )

        export_pseudotime_table(
            adata          = adata_run,
            name           = run_name,
            output_dir     = run_dir,
            annotation_col = annotation_col,
        )

        plot_pseudotime_trajectory(
            adata          = adata_run,
            name           = run_name,
            output_dir     = run_dir,
            annotation_col = annotation_col,
            show_inline    = show_inline,
        )

        plot_root_cell(
            adata       = adata_run,
            name        = run_name,
            output_dir  = run_dir,
            show_inline = show_inline,
        )

        trajectory_runs[run_name] = {
            "adata": adata_run,
            "output_dir": run_dir,
            "params": params,
        }

    if selected_run is None:
        selected_run = runs[-1].get("name", runs[-1].get("id"))

    if selected_run not in trajectory_runs:
        available = list(trajectory_runs)
        raise ValueError(f"selected_run '{selected_run}' not found. Available: {available}")

    selected = trajectory_runs[selected_run]
    print(f"\nSTEP 27 COMPLETE")
    print(f"Selected run: {selected_run}")
    print(f"Selected folder: {selected['output_dir']}")

    return selected["adata"], selected["output_dir"], trajectory_runs


# =============================================================================
# plot_genes_on_trajectory
# =============================================================================
# Step 28 helper: plot pseudotime and selected marker genes on the FA graph.
def plot_genes_on_trajectory(
    adata,
    output_dir,
    gene_sets=None,
    cmap="viridis",
    size=200,
):
    os.makedirs(output_dir, exist_ok=True)
    sc.settings.figdir = output_dir
    sc.set_figure_params(figsize=(10, 10), dpi_save=300)

    # Always plot pseudotime first.
    if "t" in adata.obs.columns:
        sc.pl.draw_graph(
            adata,
            color="t",
            color_map=cmap,
            add_outline=True,
            size=size,
            legend_fontsize=10,
            legend_fontoutline=2,
            show=True,
            save="_pseudotime_t.png",
        )
    else:
        print("Pseudotime column 't' was not found; skipping pseudotime plot.")

    gene_sets = gene_sets or {}
    for label, genes in gene_sets.items():
        genes_found = [g for g in genes if g in adata.var_names]
        genes_missing = [g for g in genes if g not in adata.var_names]

        if genes_missing:
            print(f"{label}: genes not found and skipped: {genes_missing}")
        if not genes_found:
            print(f"{label}: no valid genes to plot.")
            continue

        sc.pl.draw_graph(
            adata,
            color=genes_found,
            add_outline=True,
            legend_fontsize=10,
            legend_fontoutline=2,
            size=size,
            cmap=cmap,
            title=label,
            show=True,
            save=f"_{label}_genes.png",
        )

    print("\nSTEP 28 COMPLETE: gene plots saved")


# =============================================================================
# _ensure_milestones_category
# =============================================================================
def _ensure_milestones_category(adata):
    if "milestones" in adata.obs.columns:
        if not pd.api.types.is_categorical_dtype(adata.obs["milestones"]):
            adata.obs["milestones"] = adata.obs["milestones"].astype("category")


# =============================================================================
# review_available_branches
# =============================================================================
# Builds the scFates dendrogram, saves diagnostic plots, and writes a compact
# branch table for choosing milestones in the next step. It does not rename or
# reinterpret milestones.
def review_available_branches(
    adata,
    annotation_col,
    output_dir,
    root_cluster=None,
    show_plots=True,
):
    os.makedirs(output_dir, exist_ok=True)
    sc.settings.figdir = output_dir

    if annotation_col not in adata.obs.columns:
        available = list(adata.obs.columns)
        raise ValueError(
            f"Annotation column '{annotation_col}' was not found. Available columns: {available}"
        )

    if "seg" not in adata.obs or "milestones" not in adata.obs:
        raise ValueError(
            "Step 29 needs the trajectory from Step 27. Run Step 27 first, then Step 29 again."
        )
    _ensure_milestones_category(adata)

    # Diagnostic plots matching the old workflow.
    sc.pl.draw_graph(
        adata,
        color              = ["milestones", "leiden"],
        palette            = sns.color_palette("colorblind"),
        add_outline        = True,
        legend_fontsize    = 10,
        legend_fontoutline = 2,
        show               = show_plots,
        save               = "_milestones_leiden.png",
    )

    scf.tl.dendrogram(adata)

    scf.pl.dendrogram(
        adata,
        color              = "seg",
        palette            = sns.color_palette("colorblind"),
        legend_fontoutline = True,
        legend_loc         = "on data",
        show               = show_plots,
        save               = "_seg.pdf",
    )

    scf.pl.dendrogram(
        adata,
        color              = "milestones",
        palette            = sns.color_palette("colorblind"),
        legend_fontoutline = True,
        legend_loc         = "on data",
        show               = show_plots,
        save               = "_milestones.pdf",
    )

    # Compact branch table: one row per milestone, oriented to user decisions.
    rows = []
    root_milestone = None
    if root_cluster is not None and "is_root" in adata.obs.columns:
        root_cells = adata.obs_names[adata.obs["is_root"]]
        if len(root_cells) > 0:
            root_milestone = str(adata.obs.loc[root_cells[0], "milestones"])

    for milestone in list(adata.obs["milestones"].cat.categories):
        mask = adata.obs["milestones"] == milestone
        cell_counts = adata.obs.loc[mask, annotation_col].value_counts()
        leiden_counts = adata.obs.loc[mask, "leiden"].value_counts() if "leiden" in adata.obs.columns else pd.Series(dtype=int)
        n_cells = int(mask.sum())
        main_cell_type = str(cell_counts.index[0]) if len(cell_counts) else "NA"
        main_percent = round(float(cell_counts.iloc[0] / n_cells * 100), 1) if n_cells else 0.0
        note = "candidate endpoint"
        if root_milestone is not None and str(milestone) == root_milestone:
            note = "root-enriched; usually not endpoint"
        elif n_cells < 30:
            note = "small branch; review carefully"
        elif main_percent < 50:
            note = "mixed branch; review carefully"

        rows.append({
            "branch_id": str(milestone),
            "n_cells": n_cells,
            "main_cell_type": main_cell_type,
            "main_cell_type_percent": main_percent,
            "celltype_counts": "; ".join(f"{k}:{v}" for k, v in cell_counts.items()),
            "leiden_counts": "; ".join(f"{k}:{v}" for k, v in leiden_counts.items()),
            "note": note,
        })

    branch_table = pd.DataFrame(rows).sort_values(["note", "n_cells"], ascending=[True, False])
    out_table = os.path.join(output_dir, "milestone_branches_for_step30.tsv")
    branch_table.to_csv(out_table, sep="\t", index=False)

    summarize_milestones_by_celltype(
        adata          = adata,
        annotation_col = annotation_col,
        output_dir     = output_dir,
        root_cluster   = root_cluster,
    )

    print("\nAvailable branches for Step 30:")
    print(branch_table[["branch_id", "n_cells", "main_cell_type", "main_cell_type_percent", "note"]].to_string(index=False))
    print(f"\nBranch table saved: {out_table}")
    print("\nCopy branch_id values into Step 30, for example:")
    print('MILESTONES_TO_ANALYZE = ["7", "34"]')

    return adata, branch_table


# =============================================================================
# run_step28_dendrogram
# =============================================================================
# Backward-compatible wrapper for older notebooks that still call the previous
# helper name. New notebooks should call review_available_branches in Step 29.
def run_step28_dendrogram(
    adata,
    annotation_col,
    output_dir,
    root_cluster=None,
):
    """Backward-compatible wrapper for the old notebook name."""
    adata, _branch_table = review_available_branches(
        adata          = adata,
        annotation_col = annotation_col,
        output_dir     = output_dir,
        root_cluster   = root_cluster,
        show_plots     = True,
    )
    return adata


# =============================================================================
# resolve_milestone_value
# =============================================================================
def resolve_milestone_value(adata, value):
    _ensure_milestones_category(adata)
    categories = list(adata.obs["milestones"].cat.categories)
    for category in categories:
        if str(category) == str(value):
            return category
    available = [str(x) for x in categories]
    raise ValueError(f"milestone {value!r} not found. Available: {available}")


# =============================================================================
# run_milestone_analysis
# =============================================================================
# For a single branch endpoint (milestone): subsets the tree, tests which genes
# change significantly along that branch, and fits smooth expression curves.
# Saves two h5ad checkpoints: *_association.h5ad and *_fitted.h5ad.
def run_milestone_analysis(
    adata,
    milestone,          # name of the branch endpoint (from adata.obs["milestones"])
    root_milestone,     # name of the root (starting point of the tree)
    output_dir,         # folder for intermediate h5ad checkpoints
    n_jobs    = 4,      # parallel CPUs (reduce to 4-8 on a laptop)
    a_cut     = 0.3,    # association threshold (0-1): lower = more genes retained
    p_val_cut = 0.001,  # p-value cutoff for gene-pseudotime significance
    name_file = "pseudotime",
):
    os.makedirs(output_dir, exist_ok=True)

    try:
        root_milestone = resolve_milestone_value(adata, root_milestone)
    except ValueError:
        if "is_root" in adata.obs.columns:
            root_cells = adata.obs_names[adata.obs["is_root"]]
            if len(root_cells) > 0:
                root_milestone = adata.obs.loc[root_cells[0], "milestones"]
                print(f"Root milestone auto-detected: '{root_milestone}'")
            else:
                raise
        else:
            raise

    milestone = resolve_milestone_value(adata, milestone)

    print(f"\n{'='*50}\nProcessing milestone: {milestone}\n{'='*50}")

    adata_branch = scf.tl.subset_tree(
        adata, root_milestone=root_milestone, milestones=[milestone], copy=True
    )

    scf.tl.test_association(adata_branch, n_jobs=n_jobs, A_cut=a_cut)
    assoc_path = os.path.join(output_dir, f"adata_{name_file}_{milestone}_association.h5ad")
    adata_branch.write_h5ad(assoc_path)

    adata_branch = sc.read_h5ad(assoc_path)
    adata_branch.var["signi"] = adata_branch.var["p_val"] < p_val_cut

    scf.tl.fit(adata_branch, n_jobs=n_jobs)
    fitted_path = os.path.join(output_dir, f"adata_{name_file}_{milestone}_fitted.h5ad")
    adata_branch.write_h5ad(fitted_path)

    print(f"✓ Milestone '{milestone}' complete → {fitted_path}")
    return adata_branch


# =============================================================================
# plot_gene_trends
# =============================================================================
# Heatmap of the most dynamically expressed genes along a branch,
# with optional markers of interest highlighted.
def plot_gene_trends(adata_fitted, milestone_name, output_dir, highlight_genes=None):
    os.makedirs(output_dir, exist_ok=True)

    adata_fitted.var["_gene"] = adata_fitted.var_names
    adata_fitted.var.index    = adata_fitted.var["_gene"]
    adata_fitted.var_names    = adata_fitted.var_names.astype(str)
    adata_fitted.var_names_make_unique()

    sc.set_figure_params(figsize=(6, 20), dpi_save=600, frameon=False)
    axes_list = scf.pl.trends(
        adata_fitted,
        highlight_features = highlight_genes or [],
        style              = "italic",
        add_outline        = True,
        basis              = "dendro",
        show_segs          = True,
        fontsize           = 10,
        figsize            = (3, 5),
        ordering           = "max",
        show               = False,
        title              = milestone_name,
    )
    fig      = axes_list[0].get_figure()
    out_path = os.path.join(output_dir, f"{milestone_name}_gene_trends.pdf")
    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"✓ Gene trends saved: {out_path}")


# =============================================================================
# genes_by_pseudotime_peak
# =============================================================================
# For each gene: computes Pearson correlation with pseudotime, the time point
# of peak expression, and the dominant cell cluster at that peak.
# Exports a CSV table ranked by peak time (useful for finding wave-like genes).
def genes_by_pseudotime_peak(
    adata,
    milestone_name,
    output_dir,
    t_key          = "t",
    leiden_key     = "leiden",
    layer_key      = "fitted",
    peak_threshold = 0.7,  # expression must reach this fraction of max to count as "peak"
):
    os.makedirs(output_dir, exist_ok=True)

    t      = adata.obs[t_key].values
    leiden = adata.obs[leiden_key].values
    X      = adata.layers[layer_key]

    # Pearson correlation with pseudotime (vectorized)
    t_c  = t - t.mean()
    X_c  = X - X.mean(axis=0)
    corr = (t_c @ X_c) / (
        np.sqrt((t_c ** 2).sum()) * np.sqrt((X_c ** 2).sum(axis=0))
    )
    adata.var["corr"] = corr
    adata.var["up"]   = corr > 0

    # Peak expression time and dominant cluster at peak
    X_norm    = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0) + 1e-9)
    mask      = X_norm > peak_threshold
    peak_t    = np.full(X.shape[1], np.nan)
    peak_leid = np.full(X.shape[1], np.nan, dtype=object)

    for g in range(X.shape[1]):
        if mask[:, g].any():
            peak_t[g]    = t[mask[:, g]].mean()
            peak_leid[g] = pd.Series(leiden[mask[:, g]]).value_counts().idxmax()

    adata.var["peak_t"]      = peak_t
    adata.var["peak_leiden"] = peak_leid

    df_out   = (adata.var
                .assign(_order=adata.var["peak_t"].fillna(2))
                .sort_values("_order")
                .drop(columns="_order"))
    out_file = os.path.join(output_dir, f"{milestone_name}_genes_by_peak.csv")
    df_out.to_csv(out_file)
    print(f"✓ Gene peak table saved: {out_file}")
    return df_out


# =============================================================================
# compute_module_score
# =============================================================================
# Projects a user-defined gene list onto the trajectory as a per-cell score
# (normalized by library size, log-scaled, z-scored across cells).
# Useful for overlaying published signatures or custom gene modules.
def compute_module_score(adata, gene_list, prefix):
    genes = [g for g in gene_list if g in adata.var_names]
    if not genes:
        print(f"  No genes found for module '{prefix}' — skipping.")
        return

    M   = adata[:, genes].X
    lib = adata.X.sum(axis=1)

    if not isinstance(M, np.ndarray):
        M = np.asarray(M.todense())
    lib = np.asarray(lib).flatten()
    lib = np.where(lib == 0, 1, lib)  # avoid division by zero

    raw  = np.asarray(M.sum(axis=1)).flatten()
    norm = scale(np.log1p(raw / lib))
    adata.obs[f"{prefix}_module_score"] = norm

    print(f"✓ Module '{prefix}': {len(genes)}/{len(gene_list)} genes scored.")


# =============================================================================
# rename_milestones_by_celltype
# =============================================================================
# Renames scFates milestone IDs (numbers) to the dominant cell type at each
# branch endpoint. Kept for optional cosmetic renaming; it is not called
# automatically by the current notebook.
def rename_milestones_by_celltype(adata, annotation_col, priority_celltypes=None):
    categories = list(adata.obs["milestones"].cat.categories)
    priority_celltypes = priority_celltypes or []

    # Force biologically important cell types, such as the root, to appear in
    # the dendrogram at the milestone where most of their cells are located.
    forced_names = {}
    counts = pd.crosstab(adata.obs[annotation_col], adata.obs["milestones"])
    for cell_type in priority_celltypes:
        if cell_type in counts.index:
            forced_names[str(counts.loc[cell_type].idxmax())] = cell_type

    seen = {}
    new_names = []
    for ms in categories:
        ms_key = str(ms)
        if ms_key in forced_names:
            base_name = forced_names[ms_key]
        else:
            mask = adata.obs["milestones"] == ms
            values = adata.obs.loc[mask, annotation_col]
            if priority_celltypes:
                values_no_priority = values[~values.isin(priority_celltypes)]
                if len(values_no_priority) > 0:
                    values = values_no_priority
            base_name = values.value_counts().index[0]

        if base_name in seen:
            seen[base_name] += 1
            name = f"{base_name}_{seen[base_name]}"
        else:
            seen[base_name] = 0
            name = base_name
        new_names.append(name)

    scf.tl.rename_milestones(adata, new_names)
    print("Milestones renamed:")
    for old, new in zip(categories, new_names):
        print(f"  {str(old):>4s}  ->  {new}")
    return adata

# =============================================================================
# summarize_milestones_by_celltype
# =============================================================================
# Shows where each selected cell type falls across the inferred milestones.
# This is useful because the biological root can be present in the trajectory
# without becoming a terminal milestone in the dendrogram.
def summarize_milestones_by_celltype(
    adata,
    annotation_col,
    output_dir,
    root_cluster=None,
):
    os.makedirs(output_dir, exist_ok=True)

    counts = pd.crosstab(adata.obs[annotation_col], adata.obs["milestones"])
    percents = counts.div(counts.sum(axis=1), axis=0).fillna(0) * 100

    summary = pd.DataFrame({
        "cell_type": counts.index,
        "n_cells": counts.sum(axis=1).values,
        "main_milestone": counts.idxmax(axis=1).astype(str).values,
        "main_milestone_percent": percents.max(axis=1).round(1).values,
    })

    counts.to_csv(os.path.join(output_dir, "milestone_counts_by_celltype.tsv"), sep="\t")
    percents.round(1).to_csv(os.path.join(output_dir, "milestone_percent_by_celltype.tsv"), sep="\t")
    summary.to_csv(os.path.join(output_dir, "milestone_summary_by_celltype.tsv"), sep="\t", index=False)

    print("\nCell types across milestones:")
    print(summary.to_string(index=False))

    if root_cluster is not None and root_cluster in list(summary["cell_type"]):
        root_row = summary.loc[summary["cell_type"] == root_cluster].iloc[0]
        print(
            f"\nRoot '{root_cluster}' is present mainly in milestone "
            f"'{root_row['main_milestone']}' "
            f"({root_row['main_milestone_percent']}% of its cells)."
        )
    elif root_cluster is not None:
        print(f"\nRoot '{root_cluster}' was not found in column '{annotation_col}'.")

    print(f"\nTables saved to: {output_dir}")
    return summary, counts, percents
