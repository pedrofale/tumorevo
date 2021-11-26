"""
Simulate tumor growth under different spatial models. Inspired by Noble et al, 2019.
"""
from cell import TumorCell
import modes as modes
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
	modes.simulate_nonspatial,
	modes.simulate_invasion,
	modes.simulate_fission,
	modes.simulate_boundary, 
]

def main(mode, output_path, plot, n_steps=1000, n_genes=1000, **kwargs):
	tumor_cell = TumorCell(n_genes=n_genes)
	env, traces = MODE_LIST[mode](n_steps, tumor_cell, **kwargs)
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

def prepare_muller(genotype_counts, genotype_parents):
	pop_df = pd.DataFrame(genotype_counts).fillna(0)	
	pop_df['Generation'] = np.arange(pop_df.shape[0])
	pop_df = pop_df.melt(id_vars=['Generation'], value_name='Population', var_name='Identity').sort_values('Generation').reset_index(drop=True)

	anc_df = pd.DataFrame([genotype_parents[-1]]).melt(var_name='Identity', value_name='Parent')

	color_by = pd.Series(np.arange(anc_df.shape[0]+1), index=pop_df['Identity'].unique())

	return pop_df, anc_df, color_by

if __name__ == '__main__':
	args = parser.parse_args()
	main(args.mode, args.output_path, args.plot, carrying_capacity=int(args.carrying_capacity), n_genes=int(args.n_genes), n_steps=int(args.n_steps))
