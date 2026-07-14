import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import valency_anndata as val
    import numpy as np
    import pandas as pd
    import scanpy as sc

    return np, pd, sc, val


@app.cell
def _(val):
    # adata = val.datasets.chile_protest(translate_to="en")
    adata = val.datasets.japanchoice(topic="2026_economy_taxation_employment", translate_to="en")

    # work-around for valency-anndata/scanpy doing in-memory edits that aren't reative.
    _ = True
    return (adata,)


@app.cell
def _(adata, val):
    _ = val.viz.schematic_diagram(adata)
    return


@app.cell
def _(adata):
    adata.var
    return


@app.cell
def _(adata):
    if "is_meta" not in adata.var.columns:
        adata.var["is_meta"] = False
    elif adata.var["is_meta"].isna().any():
        adata.var["is_meta"] = adata.var["is_meta"].fillna(False).astype(bool)
    return


@app.cell
def _(adata, val):
    _ = val.tools.recipe_polis(adata)
    return


@app.cell
def _(adata, val):
    _ = val.preprocessing.impute(adata, strategy="knn", source_layer="X_masked", target_layer="X_masked_imputed_knn5")
    _ = val.viz.schematic_diagram(adata)
    return


@app.cell
def _(adata, val):
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
    return


@app.cell
def _(adata, val):
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
    return


@app.cell
def _(adata, np, pd, sc, val):
    # Leiden clustering on the PaCMAP embedding, restricted to cluster_mask
    # (mirrors the kmeans_pacmap masking) since leiden itself has no mask_obs arg.
    _mask = adata.obs["cluster_mask"].to_numpy()
    _sub = adata[_mask].copy()
    sc.pp.neighbors(_sub, use_rep="X_pacmap_polis")
    val.tools.leiden(_sub, key_added="leiden_pacmap")

    _labels = np.full(adata.n_obs, np.nan, dtype=object)
    _labels[_mask] = _sub.obs["leiden_pacmap"].to_numpy()
    adata.obs["leiden_pacmap"] = pd.Categorical(_labels)
    return


@app.cell
def _(adata, np, pd, sc, val):
    # Leiden clustering on the LocalMAP embedding, restricted to cluster_mask.
    _mask = adata.obs["cluster_mask"].to_numpy()
    _sub = adata[_mask].copy()
    sc.pp.neighbors(_sub, use_rep="X_localmap_polis")
    val.tools.leiden(_sub, key_added="leiden_localmap")

    _labels = np.full(adata.n_obs, np.nan, dtype=object)
    _labels[_mask] = _sub.obs["leiden_localmap"].to_numpy()
    adata.obs["leiden_localmap"] = pd.Categorical(_labels)
    return


@app.cell
def _(adata, val):
    val.viz.embedding(adata, basis="pca_polis", color="kmeans_polis")
    return


@app.cell
def _(adata, val):
    val.viz.embedding(adata, basis="pacmap_polis", color="kmeans_pacmap")
    return


@app.cell
def _(adata, val):
    val.viz.embedding(adata, basis="localmap_polis", color="kmeans_localmap")
    return


@app.cell
def _(adata, val):
    val.viz.embedding(adata, basis="pacmap_polis", color="leiden_pacmap")
    return


@app.cell
def _(adata, val):
    val.viz.embedding(adata, basis="localmap_polis", color="leiden_localmap")
    return


@app.cell
def _(adata):
    # This conversation used strict moderation, so unless explicitly voted in, it wasn't shown.
    strict_mod_in_mask = (adata.var["moderation_state"] > 0).to_numpy()
    return (strict_mod_in_mask,)


@app.cell
def _(adata, strict_mod_in_mask, val):
    # Using original data
    _ = val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="kmeans_polis",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using mean imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_mean",
        groupby="kmeans_polis",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using k-nearest neighbors (N=5) imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_knn5",
        groupby="kmeans_polis",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


@app.cell
def _(adata, strict_mod_in_mask, val):
    # Using original data
    _ = val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="kmeans_pacmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using mean imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_mean",
        groupby="kmeans_pacmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using k-nearest neighbors (N=5) imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_knn5",
        groupby="kmeans_pacmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


@app.cell
def _(adata, strict_mod_in_mask, val):
    # Using original data
    _ = val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="kmeans_localmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using mean imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_mean",
        groupby="kmeans_localmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using k-nearest neighbors (N=5) imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_knn5",
        groupby="kmeans_localmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


@app.cell
def _(adata, strict_mod_in_mask, val):
    # Using original data
    _ = val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="leiden_pacmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using mean imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_mean",
        groupby="leiden_pacmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using k-nearest neighbors (N=5) imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_knn5",
        groupby="leiden_pacmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


@app.cell
def _(adata, strict_mod_in_mask, val):
    # Using original data
    _ = val.viz.heatmap(
        adata,
        discrete=True,
        layer="raw_sparse",
        groupby="leiden_localmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using mean imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_mean",
        groupby="leiden_localmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )

    # Using k-nearest neighbors (N=5) imputed data
    _ = val.viz.heatmap(
        adata,
        layer="X_masked_imputed_knn5",
        groupby="leiden_localmap",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


if __name__ == "__main__":
    app.run()
