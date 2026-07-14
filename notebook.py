import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import valency_anndata as val

    return (val,)


@app.cell
def _(val):
    adata = val.datasets.chile_protest(translate_to="en")

    # work-around for valency-anndata/scanpy doing in-memory edits that aren't reative.
    _ = True
    return (adata,)


@app.cell
def _(adata, val):
    _ = val.viz.schematic_diagram(adata)
    return


@app.cell
def _(adata, val):
    _ = val.preprocessing.impute(adata, strategy="knn", source_layer="X_masked", target_layer="X_masked_imputed_knn5")
    _ = val.viz.schematic_diagram(adata)
    return


@app.cell
def _(adata):
    adata.var
    return


@app.cell
def _(adata, val):
    _ = val.tools.recipe_polis(adata)
    return


@app.cell
def _(adata, val):
    _ = val.viz.schematic_diagram(adata)
    return


@app.cell
def _(adata):
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
    return


@app.cell
def _(adata, strict_mod_in_mask, val):
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


if __name__ == "__main__":
    app.run()
