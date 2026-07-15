import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import valency_anndata as val
    import numpy as np
    import pandas as pd
    import jscatter


    return jscatter, mo, np, pd, val


@app.cell
def adata_rev_state(mo):
    get_adata_rev, set_adata_rev = mo.state(0)
    return get_adata_rev, set_adata_rev


@app.cell
def _(val):
    adata = val.datasets.japanchoice(topic="2026_economy_taxation_employment", translate_to="en")
    return (adata,)


@app.cell
def _(adata, get_adata_rev, val):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.schematic_diagram(adata)
    return


@app.cell
def _(adata, get_adata_rev):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    adata.var
    return


@app.cell
def _(adata, set_adata_rev):
    if "is_meta" not in adata.var.columns:
        adata.var["is_meta"] = False
    elif adata.var["is_meta"].isna().any():
        adata.var["is_meta"] = adata.var["is_meta"].fillna(False).astype(bool)
    set_adata_rev(lambda v: v + 1)
    return


@app.cell
def _(adata, set_adata_rev, val):
    val.tools.recipe_polis(adata)
    set_adata_rev(lambda v: v + 1)
    return


@app.cell
def _(adata, set_adata_rev, val):
    val.preprocessing.impute(adata, strategy="knn", source_layer="X_masked", target_layer="X_masked_imputed_knn5")
    set_adata_rev(lambda v: v + 1)
    val.viz.schematic_diagram(adata)
    return


@app.cell
def _(adata, set_adata_rev, val):
    # PaCMAP-based recipe, modelled after val.tools.recipe_polis:
    # same masked/imputed layer and participant vote-count mask, but
    # embedding via PaCMAP instead of PCA + sparsity-aware scaling.
    val.tools.pacmap(
        adata,
        layer="X_masked_imputed_knn5",
        key_added="X_pacmap_polis",
    )
    val.tools.kmeans(
        adata,
        use_rep="X_pacmap_polis",
        k_bounds=(2, 7),
        init="polis",
        mask_obs="cluster_mask",
        key_added="kmeans_pacmap",
    )
    set_adata_rev(lambda v: v + 1)
    return


@app.cell
def _(adata, set_adata_rev, val):
    # LocalMAP-based recipe, modelled after val.tools.recipe_polis:
    # same participant vote-count mask, but embedding via LocalMAP on the
    # KNN-imputed layer instead of PCA + sparsity-aware scaling.
    val.tools.localmap(
        adata,
        layer="X_masked_imputed_knn5",
        key_added="X_localmap_polis",
    )
    val.tools.kmeans(
        adata,
        use_rep="X_localmap_polis",
        k_bounds=(2, 7),
        init="polis",
        mask_obs="cluster_mask",
        key_added="kmeans_localmap",
    )
    set_adata_rev(lambda v: v + 1)
    return


@app.cell
def _(adata, np, pd, set_adata_rev, val):
    # Leiden clustering on the PCA embedding, restricted to cluster_mask
    # (mirrors the kmeans_pacmap masking) since leiden itself has no mask_obs arg.
    # Uses the resolution=1 default; see the leiden_n_neighbors_slider/
    # leiden_resolution_slider sweep below for how cluster count varies with
    # these params, e.g. to match kmeans's k_bounds=(2, 7) range.
    _mask = adata.obs["cluster_mask"].to_numpy()
    _sub = adata[_mask].copy()
    val.preprocessing.neighbors(_sub, use_rep="X_pca_polis")
    val.tools.leiden(_sub, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)

    _labels = np.full(adata.n_obs, np.nan, dtype=object)
    _labels[_mask] = _sub.obs["leiden"].to_numpy()
    adata.obs["leiden"] = pd.Categorical(_labels)
    set_adata_rev(lambda v: v + 1)
    return


@app.cell
def leiden_sweep_sliders(mo):
    leiden_n_neighbors_slider = mo.ui.slider(5, 50, step=5, value=15, label="n_neighbors")
    leiden_resolution_slider = mo.ui.slider(0.1, 2.0, step=0.05, value=1.0, label="resolution")
    mo.hstack([leiden_n_neighbors_slider, leiden_resolution_slider])
    return leiden_n_neighbors_slider, leiden_resolution_slider


@app.cell(hide_code=True)
def leiden_sweep_compute(
    adata,
    leiden_n_neighbors_slider,
    leiden_resolution_slider,
    mo,
    np,
    pd,
    set_adata_rev,
    val,
):
    # Interactive leiden sweep: re-run neighbors + leiden on the PCA embedding
    # with the slider params, to see how cluster count/sizes vary.
    _mask = adata.obs["cluster_mask"].to_numpy()
    _sub = adata[_mask].copy()
    val.preprocessing.neighbors(
        _sub, n_neighbors=leiden_n_neighbors_slider.value, use_rep="X_pca_polis"
    )
    val.tools.leiden(
        _sub,
        resolution=leiden_resolution_slider.value,
        key_added="leiden_sweep",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    _labels = np.full(adata.n_obs, np.nan, dtype=object)
    _labels[_mask] = _sub.obs["leiden_sweep"].to_numpy()
    adata.obs["leiden_sweep"] = pd.Categorical(_labels)
    set_adata_rev(lambda v: v + 1)

    leiden_sweep_n_clusters = _sub.obs["leiden_sweep"].nunique()
    leiden_sweep_sizes = _sub.obs["leiden_sweep"].value_counts()
    mo.md(
        f"**{leiden_sweep_n_clusters} clusters** "
        f"(sizes {leiden_sweep_sizes.min()}-{leiden_sweep_sizes.max()}) "
        f"at n_neighbors={leiden_n_neighbors_slider.value}, "
        f"resolution={leiden_resolution_slider.value:.2f}"
    )
    return


@app.cell
def leiden_sweep_viz(adata, get_adata_rev, val):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.embedding(adata, basis="pca_polis", color="leiden_sweep")
    return


@app.cell
def _(adata, get_adata_rev, val):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.embedding(adata, basis="pca_polis", color="kmeans_polis")
    return


@app.cell
def _(adata, get_adata_rev, val):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.embedding(adata, basis="pca_polis", color="leiden")
    return


@app.cell
def _(adata, get_adata_rev, strict_mod_in_mask, val):
    # Using original data
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="kmeans_polis",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


@app.cell
def _(adata, get_adata_rev, val):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.embedding(adata, basis="pacmap_polis", color="kmeans_pacmap")
    return


@app.cell
def _(adata, get_adata_rev, val):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.embedding(adata, basis="localmap_polis", color="kmeans_localmap")
    return


@app.cell
def _(adata, get_adata_rev, val):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.embedding(adata, basis="pacmap_polis", color="leiden")
    return


@app.cell
def _(adata, get_adata_rev, val):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.embedding(adata, basis="localmap_polis", color="leiden")
    return


@app.cell
def jscatter_setup(adata, get_adata_rev, jscatter, mo, pd):
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere

    _mask = adata.obs["cluster_mask"].to_numpy()
    _sub = adata[_mask]

    projection_df = pd.DataFrame({
        "pca_x": _sub.obsm["X_pca_polis"][:, 0],
        "pca_y": _sub.obsm["X_pca_polis"][:, 1],
        "pacmap_x": _sub.obsm["X_pacmap_polis"][:, 0],
        "pacmap_y": _sub.obsm["X_pacmap_polis"][:, 1],
        "localmap_x": _sub.obsm["X_localmap_polis"][:, 0],
        "localmap_y": _sub.obsm["X_localmap_polis"][:, 1],
        "leiden": _sub.obs["leiden"].astype(str).to_numpy(),
    })

    _PROJECTIONS = {
        "PCA": ("pca_x", "pca_y"),
        "PaCMAP": ("pacmap_x", "pacmap_y"),
        "LocalMAP": ("localmap_x", "localmap_y"),
    }

    projection_scatter = jscatter.Scatter(
        data=projection_df, x="pca_x", y="pca_y", color_by="leiden"
    )
    projection_scatter.options(transition_points_duration=1500)
    projection_widget = mo.ui.anywidget(projection_scatter.widget)


    def _on_projection_change(name):
        x_col, y_col = _PROJECTIONS[name]
        projection_scatter.xy(x=x_col, y=y_col, animate=True)


    projection_radio = mo.ui.radio(
        options=list(_PROJECTIONS.keys()),
        value="PCA",
        label="Projection",
        on_change=_on_projection_change,
    )

    return projection_radio, projection_widget


@app.cell
def jscatter_display(mo, projection_radio, projection_widget):
    mo.vstack([projection_radio, projection_widget])
    return


@app.cell
def _(adata, get_adata_rev):
    # This conversation used strict moderation, so unless explicitly voted in, it wasn't shown.
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    strict_mod_in_mask = (adata.var["moderation_state"] > 0).to_numpy()
    return (strict_mod_in_mask,)


@app.cell
def _(adata, get_adata_rev, strict_mod_in_mask, val):
    # Using original data
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="kmeans_localmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


@app.cell
def _(adata, get_adata_rev, strict_mod_in_mask, val):
    # Using original data
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="leiden",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


@app.cell
def _(adata, get_adata_rev, strict_mod_in_mask, val):
    # Using original data
    _ = get_adata_rev()  # force rerun whenever adata is mutated elsewhere
    val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="leiden",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


if __name__ == "__main__":
    app.run()
