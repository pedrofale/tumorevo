"""
Simulate tumor growth under different spatial models. Inspired by Noble et al, 2019.
"""
from .cell import TumorCell
from .modes import *
from .util import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import argparse
import os

from pymuller import muller


parser = argparse.ArgumentParser()
parser.add_argument('--mode', default=0)
parser.add_argument('--carrying_capacity', default=10)
parser.add_argument('--n_genes', default=100)
parser.add_argument('--n_steps', default=10)
parser.add_argument('--grid_side', default=10)
parser.add_argument('--output_path', default='./')
parser.add_argument('--plot', default=False)

MODE_LIST = [
	simulate_nonspatial,
	simulate_invasion,
	simulate_fission,
	simulate_boundary, 
]

def main():
	args = parser.parse_args()

	mode = args.mode
	output_path = args.output_path
	plot = args.plot
	carrying_capacity = int(args.carrying_capacity)
	n_genes=int(args.n_genes)
	n_steps=int(args.n_steps)

	tumor_cell = TumorCell(n_genes=n_genes)
	env, traces = MODE_LIST[mode](n_steps, tumor_cell)
	print(env.get_genotype_frequencies())
	print(len(env.cells))
	#print(env.get_diversity())
	
	#env.save_grid(os.path.join(output_path, f'mode_{mode}_grid.txt'))
	#env.save_event_tree(os.path.join(output_path, f'mode_{mode}_tree.txt'))
	if plot:
		populations_df, adjacency_df, color_by = prepare_muller([t['genotypes_counts'] for t in traces], [t['genotypes_parents'] for t in traces])
		ax = muller(populations_df, adjacency_df, color_by, colorbar=False)
		plt.show()
		#env.plot_grid(os.path.join(output_path, f'mode_{mode}_grid.png'))
		#env.plot_tree(os.path.join(output_path, f'mode_{mode}_tree.png'))	

if __name__ == '__main__':
	main()
