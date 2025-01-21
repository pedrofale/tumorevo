import numpy as np

def sample(cell_ids, fraction):
    n_cells = len(cell_ids)

    # Take random subset
    sampled_cells = np.random.choice(n_cells, size=int(n_cells*fraction), replace=False)

    return cell_ids[sampled_cells]

def sample_region(grid):
    pass

def sample_regions(grid):
    pass
