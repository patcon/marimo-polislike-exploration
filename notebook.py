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
    val.viz.schematic_diagram(adata)
    return


if __name__ == "__main__":
    app.run()
