import numpy as np
import pandas as pd
import packcircles as pc
import pymuller
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

def prepare_plots(genotype_counts, genotype_parents):
        pop_df = genotype_counts
        pop_df['Generation'] = np.arange(pop_df.shape[0])
        pop_df = pop_df.melt(id_vars=['Generation'], value_name='Population', var_name='Identity').sort_values('Generation').reset_index(drop=True)

        anc_df = genotype_parents.melt(var_name='Identity', value_name='Parent').astype(str)

        color_by = pd.Series(np.arange(anc_df.shape[0]+1), index=pop_df['Identity'].unique())

        return pop_df, anc_df, color_by


def plot_deme(n_cells, genotype_counts, populations_df, adjacency_df, color_by, average_radius=10, colormap="terrain", ax=None):
	"""Create a circle for each cell and color it by genotype.
	The colors are consistent with the Muller plots.
	"""
	# Re-use PyMuller to get the same genotype colors
	x = populations_df['Generation'].unique()
	y_table = pymuller.logic._get_y_values(populations_df, adjacency_df, 10)
	final_order = y_table.columns.values
	Y = y_table.to_numpy().T
	cmap = plt.get_cmap(colormap)
	color_by = color_by.copy()
	ordered_colors = color_by.loc[final_order]
	norm = matplotlib.colors.Normalize(vmin=np.min(ordered_colors), vmax=np.max(ordered_colors))
	color_map = cmap(norm(ordered_colors.values))

	radii = []
	colors = []
	n_total_cells = sum(genotype_counts.values)
	min_r = max(1, average_radius-10)
	max_r = average_radius+10
	for i, genotype in enumerate(final_order):
		n = max(1, int(n_cells * genotype_counts[genotype]/n_total_cells))
		g_radii = np.random.randint(min_r, max_r, size=n).tolist()
		radii = radii + g_radii
		colors = colors + [color_map[i]]*n

	# Randomly permute cells to model well mixed population in deme
	perm = np.random.choice(len(radii), size=len(radii), replace=False)
	radii = np.array(radii)[perm]
	colors = np.array(colors)[perm]

	# Pack circles
	circles = pc.pack(radii)
	
	if ax is None:
		fig, ax = plt.subplots(figsize=(5,5))
	else:
		plt.sca(ax) 
	for i, (x,y,radius) in enumerate(circles):
	    patch = plt.Circle(
		(x,y),
		radius,
		color=colors[i],
		alpha=0.65
	    )
	    ax.add_patch(patch)
	ax.set(xlim=(-150, 140), ylim=(-180, 170))
	plt.axis('off')

