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
def jscatter_setup(adata, get_adata_rev, jscatter, mo, np, pd, val):
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
        "leiden_sweep": _sub.obs["leiden_sweep"].astype(str).to_numpy(),
        "kmeans_polis": _sub.obs["kmeans_polis"].astype(str).to_numpy(),
        "kmeans_pacmap": _sub.obs["kmeans_pacmap"].astype(str).to_numpy(),
        "kmeans_localmap": _sub.obs["kmeans_localmap"].astype(str).to_numpy(),
    })

    _PROJECTIONS = {
        "PCA": ("pca_x", "pca_y"),
        "PaCMAP": ("pacmap_x", "pacmap_y"),
        "LocalMAP": ("localmap_x", "localmap_y"),
    }

    _COLOR_OPTIONS = {
        "leiden (canonical)": "leiden",
        "leiden_sweep (slider)": "leiden_sweep",
        "kmeans (PCA)": "kmeans_polis",
        "kmeans (PaCMAP)": "kmeans_pacmap",
        "kmeans (LocalMAP)": "kmeans_localmap",
    }

    # Precompute a grid of leiden clusterings at different n_neighbors/resolution
    # combos, so switching colors via the dropdown is instant instead of
    # recomputing on the fly. Neighbors are computed once per n_neighbors value
    # and reused across resolutions, since only resolution affects the (cheap)
    # partitioning step.
    _GRID_N_NEIGHBORS = [10, 15, 30]
    _GRID_RESOLUTIONS = [0.25, 0.5, 1.0]
    for _nn in _GRID_N_NEIGHBORS:
        _sub_nn = _sub.copy()
        val.preprocessing.neighbors(_sub_nn, n_neighbors=_nn, use_rep="X_pca_polis")
        for _res in _GRID_RESOLUTIONS:
            _key = f"leiden_nn{_nn}_res{_res:.2f}"
            val.tools.leiden(
                _sub_nn,
                resolution=_res,
                key_added=_key,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            projection_df[_key] = _sub_nn.obs[_key].astype(str).to_numpy()
            _COLOR_OPTIONS[f"leiden (nn={_nn}, res={_res:.2f})"] = _key

    # Per-statement vote coloring: agree/disagree/pass/no-vote for each comment,
    # selectable via a separate "Statement" dropdown that only appears once
    # "Votes" is chosen as the color mode.
    _VOTE_LABELS = {1: "agree", -1: "disagree", 0: "pass"}
    _VOTE_COLOR_MAP = {"agree": "#2ca02c", "disagree": "#d62728", "pass": "#f1c40f", "no vote": "#d3d3d3"}
    _VOTE_OPTIONS = {}
    projection_vote_full_text = {}
    _raw_votes = _sub.layers["raw_sparse"]
    for _i, (_comment_id, _content) in enumerate(zip(adata.var.index, adata.var["content"])):
        _vote_col = _raw_votes[:, _i]
        _vote_str = np.full(_vote_col.shape, "no vote", dtype=object)
        for _v, _label in _VOTE_LABELS.items():
            _vote_str[_vote_col == _v] = _label
        _key = f"vote_{_comment_id}"
        projection_df[_key] = _vote_str
        _short_content = _content if len(_content) <= 60 else _content[:57] + "..."
        _VOTE_OPTIONS[f"#{_comment_id}: {_short_content}"] = _key
        projection_vote_full_text[_key] = f"#{_comment_id}: {_content}"

    _COLOR_OPTIONS["Votes..."] = "__votes__"

    projection_scatter = jscatter.Scatter(
        data=projection_df, x="pca_x", y="pca_y", color_by="leiden"
    )
    projection_scatter.options(transition_points_duration=1500)
    projection_widget = mo.ui.anywidget(projection_scatter.widget)


    def _on_projection_change(name):
        x_col, y_col = _PROJECTIONS[name]
        projection_scatter.xy(x=x_col, y=y_col, animate=True)


    def _on_color_change(column):
        if column == "__votes__":
            projection_scatter.color(by=projection_vote_dropdown.value, map=_VOTE_COLOR_MAP)
        else:
            projection_scatter.color(by=column, map="auto")


    def _on_vote_change(column):
        projection_scatter.color(by=column, map=_VOTE_COLOR_MAP)


    projection_dropdown = mo.ui.dropdown(
        options=list(_PROJECTIONS.keys()),
        value="PCA",
        label="Projection",
        on_change=_on_projection_change,
    )

    projection_color_dropdown = mo.ui.dropdown(
        options=_COLOR_OPTIONS,
        value="leiden (canonical)",
        label="Color by",
        on_change=_on_color_change,
    )

    projection_vote_dropdown = mo.ui.dropdown(
        options=_VOTE_OPTIONS,
        value=list(_VOTE_OPTIONS.keys())[0],
        label="Statement",
        on_change=_on_vote_change,
    )

    return (
        projection_color_dropdown,
        projection_dropdown,
        projection_vote_dropdown,
        projection_vote_full_text,
        projection_widget,
    )


@app.cell
def jscatter_display(
    mo,
    projection_color_dropdown,
    projection_dropdown,
    projection_vote_dropdown,
    projection_vote_full_text,
    projection_widget,
):
    _controls = [projection_dropdown, projection_color_dropdown]
    _extra = []
    if projection_color_dropdown.value == "__votes__":
        _controls.append(projection_vote_dropdown)
        _extra.append(mo.md(f"**{projection_vote_full_text[projection_vote_dropdown.value]}**"))

    mo.vstack(
        [
            mo.hstack(_controls, justify="start"),
            *_extra,
            projection_widget,
        ]
    )

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
