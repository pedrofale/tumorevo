"""
Create a cartoon of a tumor given the frequencies of different genotypes.
"""
from .util import *

import pandas as pd
import matplotlib.pyplot as plt

import click
import os
from pathlib import Path

from pymuller import muller


@click.command(help="Plot the evolution of a tumor.")
@click.argument(
    "genotype-counts",
    type=click.Path(exists=True, dir_okay=False),
)
@click.argument(
    "genotype-parents",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option("-c", "--cells", default=100, help="Number of cells in slice plot.")
@click.option(
    "-r",
    "--average-radius",
    default=10,
    help="Average radius of circles in slice plot.",
)
@click.option("--grid-file", default="", help="Path to grid file.")
@click.option("--colormap", default="gnuplot", help="Colormap for genotypes.")
@click.option("--dpi", default=100, help="DPI for figures.")
@click.option("--plot", is_flag=True, help="Plot all the figures.")
@click.option("--do-muller", is_flag=True, help="Make a Muller plot.")
@click.option("--do-slice", is_flag=True, help="Make a slice plot.")
@click.option("--do-tree", is_flag=True, help="Make a clone tree plot.")
@click.option(
    "--normalize", is_flag=True, help="Normalize the abundances in the Muller plot."
)
@click.option("--labels", is_flag=True, help="Annotate the clone tree plot.")
@click.option(
    "--remove", is_flag=True, help="Remove empty clones in the clone tree plot."
)
@click.option(
    "-o", "--output-path", default="./", help="Directory to write figures into."
)
def main(
    genotype_counts,
    genotype_parents,
    cells,
    average_radius,
    grid_file,
    colormap,
    dpi,
    plot,
    do_muller,
    do_slice,
    do_tree,
    normalize,
    labels,
    remove,
    output_path,
):
    genotype_counts = pd.read_csv(genotype_counts, index_col=0)
    genotype_parents = pd.read_csv(genotype_parents, index_col=0, dtype=str)
    if grid_file != "":
        grid = pd.read_csv(grid_file, index_col=0, dtype=str)

    pop_df, anc_df, color_by = prepare_plots(genotype_counts, genotype_parents)
    cmap, genotypes = get_colormap(pop_df, anc_df, color_by, colormap)

    if plot:
        fig, ax_list = plt.subplots(ncols=3, sharex=False, dpi=dpi, figsize=(8, 2))
        muller(
            pop_df,
            anc_df,
            color_by,
            ax=ax_list[0],
            colorbar=False,
            colormap=colormap,
            normalize=normalize,
            background_strain=False,
        )
        plt.axis("off")

        if grid_file == "":
            plot_deme(
                cells,
                genotype_counts.iloc[-1],
                pop_df,
                anc_df,
                color_by,
                average_radius=average_radius,
                colormap=colormap,
                ax=ax_list[1],
            )
        else:
            plot_grid(grid, cmap, genotypes, ax=ax_list[1])

        plot_tree(
            genotype_parents,
            pop_df,
            anc_df,
            color_by,
            genotype_counts=genotype_counts.iloc[-1],
            filter_clones=remove,
            labels=labels,
            colormap=colormap,
            ax=ax_list[2],
        )
        plt.show()
    else:
        Path(output_path).mkdir(parents=True, exist_ok=True)
        if do_muller:
            ax = muller(
                pop_df,
                anc_df,
                color_by,
                colorbar=False,
                colormap=colormap,
                normalize=normalize,
            )
            plt.axis("off")
            plt.savefig(
                os.path.join(output_path, "muller.pdf"), dpi=dpi, bbox_inches="tight"
            )

        if do_slice:
            if grid_file == "":
                ax = plot_deme(
                    cells,
                    genotype_counts.iloc[-1],
                    pop_df,
                    anc_df,
                    color_by,
                    average_radius=average_radius,
                    colormap=colormap,
                )
            else:
                ax = plot_grid(grid, cmap)
            plt.savefig(
                os.path.join(output_path, "slice.pdf"), dpi=dpi, bbox_inches="tight"
            )

        if do_tree:
            ax = plot_tree(
                genotype_parents,
                pop_df,
                anc_df,
                color_by,
                genotype_counts=genotype_counts.iloc[-1],
                filter_clones=remove,
                labels=labels,
                colormap=colormap,
            )
            plt.savefig(
                os.path.join(output_path, "tree.pdf"), dpi=dpi, bbox_inches="tight"
            )


if __name__ == "__main__":
    main()
