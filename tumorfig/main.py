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
@click.option(
	"-f1",
	"--genotype-counts",
	help="CSV file containing counts of genotypes at each time step."
)
@click.option(
	"-f2",
	"--genotype-parents",
	help="CSV file containing the parents of each genotype across all time steps."
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
def main(genotype_counts, genotype_parents, cells, average_radius, colormap):
	genotype_counts = pd.read_csv(genotype_counts, index_col=0)
	genotype_parents = pd.read_csv(genotype_parents, index_col=0, dtype=str)
	
	pop_df, anc_df, color_by = prepare_plots(genotype_counts, genotype_parents)
	
	fig, ax_list = plt.subplots(ncols=2, sharex=False)
	muller(pop_df, anc_df, color_by, ax=ax_list[0], colorbar=False)
	plot_slice(cells, genotype_counts.iloc[-1], pop_df, anc_df, color_by, 
			average_radius=average_radius, colormap=colormap, ax=ax_list[1])
	plt.show()
		

if __name__ == "__main__":
	main()
