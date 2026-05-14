# =============================================================================
# ScRNA_Pseudotime_Functions.py — Pseudotime analysis helper functions
# =============================================================================
# Companion to capitulo3_pseudotime.py.
# Loaded automatically at the start of that script.


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
def plot_trajectory_graphs(adata, name, output_dir, annotation_col):
    os.makedirs(output_dir, exist_ok=True)
    for color_by, suffix in [(annotation_col, "annotation"), ("leiden", "leiden")]:
        fig, ax = plt.subplots(figsize=(10, 10))
        scf.pl.graph(adata, color_cells=color_by, ax=ax, show=False)
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, f"{name}_{suffix}.pdf"))
        plt.close(fig)
    print(f"✓ Trajectory graphs saved to {output_dir}")


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
