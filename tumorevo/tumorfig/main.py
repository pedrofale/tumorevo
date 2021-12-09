"""
Create a cartoon of a tumor given the frequencies of different genotypes.
"""
from .util import *

import pandas as pd
import matplotlib.pyplot as plt

import click
import os

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
@click.option(
	"-c",
	"--cells",
	default=100,
	help="Number of cells in slice plot."	
)
@click.option(
	"-r",
	"--average-radius",
	default=10,
	help="Average radius of circles in slice plot."
)
@click.option(
	"-m",
	"--colormap",
	default="terrain",
	help="Colormap for genotypes."
)
@click.option(
	"-p",
	"--plot",
	is_flag=True,
	help="Plot all the figures."
)
@click.option(
	"-o",
	"--output-path",
	default="./",
	help="Directory to write figures into."
)
def main(genotype_counts, genotype_parents, cells, average_radius, colormap, plot, output_path):
	genotype_counts = pd.read_csv(genotype_counts, index_col=0)
	genotype_parents = pd.read_csv(genotype_parents, index_col=0, dtype=str)
	
	pop_df, anc_df, color_by = prepare_plots(genotype_counts, genotype_parents)

	if plot:
		fig, ax_list = plt.subplots(ncols=2, sharex=False)
		muller(pop_df, anc_df, color_by, ax=ax_list[0], colorbar=False)
		plot_deme(cells, genotype_counts.iloc[-1], pop_df, anc_df, color_by, 
				average_radius=average_radius, colormap=colormap, ax=ax_list[1])
		plt.show()
	else:
		Path(output_path).mkdir(parents=True, exist_ok=True)
		ax = muller(pop_df, anc_df, color_by, colorbar=False)
		plt.savefig(os.path.join(output_path, 'muller.png'), bbox_inches='tight')

		ax = plot_deme(cells, genotype_counts.iloc[-1], pop_df, anc_df, color_by,
				average_radius=average_radius, colormap=colormap)
		plt.savefig(os.path.join(output_path, 'slice.png'), bbox_inches='tight')	

if __name__ == "__main__":
	main()
