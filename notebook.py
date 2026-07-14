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
def _(adata, val):
    _ = val.tools.recipe_polis(adata)
    return


@app.cell
def _(adata, val):
    _ = val.viz.schematic_diagram(adata)
    return


@app.cell
def _(adata, val):
    strict_mod_in_mask = (adata.var["moderation_state"] > 0).to_numpy()
    _ = val.viz.heatmap(
        adata,
        discrete=True,
        groupby="kmeans_polis",
        mask_obs="cluster_mask",
        mask_var=strict_mod_in_mask,
    )
    return


if __name__ == "__main__":
    app.run()
