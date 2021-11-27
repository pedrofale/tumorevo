"""
Simulate tumor growth under different spatial models. Inspired by Noble et al, 2019.
"""
from .cell import TumorCell
from .modes import *
from .util import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import click
import os

from pymuller import muller


MODE_LIST = [
	simulate_nonspatial,
	simulate_invasion,
	simulate_fission,
	simulate_boundary, 
]

@click.command(help="Simulate tumor evolution under different spatial constraints.")
@click.option("-m", "--mode", default=0, help="Spatial structure.")
@click.option(
    "-k",
    "--carrying-capacity",
    default=100,
    help="Deme carrying capacity.",
)
@click.option(
    "-g",
    "--genes",
    default=100,
    help="Number of genes.",
)
@click.option(
    "-s",
    "--steps",
    default=1000,
    help="Number of steps in simulation.",
)
@click.option(
   "-r",
   "--grid-side",
   default=10,
   help="Number of units in grid side.",
)
@click.option(
   "-o",
   "--output-path",
   default="./",
   help="Output path."
)
@click.option(
   "--plot",
   is_flag=True,
   default=False,
   help="Show the results in a figure."
)
def main(mode, carrying_capacity, genes, steps, grid_side, output_path, plot):
	tumor_cell = TumorCell(n_genes=genes)
	env, traces = MODE_LIST[mode](steps, tumor_cell)
	print(env.get_genotype_frequencies())
	print(len(env.cells))
	#print(env.get_diversity())

	populations_df, adjacency_df, color_by = prepare_muller([t['genotypes_counts'] for t in traces], [t['genotypes_parents'] for t in traces])
	
	adjacency_df.to_csv(os.path.join(output_path, 'adjacency.csv'))
	
	#env.save_grid(os.path.join(output_path, f'mode_{mode}_grid.txt'))
	#env.save_event_tree(os.path.join(output_path, f'mode_{mode}_tree.txt'))
	if plot:
		ax = muller(populations_df, adjacency_df, color_by, colorbar=False)
		plt.show()
		#env.plot_grid(os.path.join(output_path, f'mode_{mode}_grid.png'))
		#env.plot_tree(os.path.join(output_path, f'mode_{mode}_tree.png'))	

if __name__ == '__main__':
	main()
