from .tumor import Tumor
from .deme import Deme
from .cell import Cell, TumorCell

from tqdm import tqdm
from copy import deepcopy


def simulate_nonspatial(n_steps, cell, **kwargs):
    # Prevent cell from migrating
    cell.dispersal_rate = 0

    # Initialise tumor with single cell
    tumor = Tumor(cell, **kwargs)

    traces = []
    traces.append(dict(genotypes_counts=deepcopy(tumor.genotypes_counts)))
    # Simulate within-deme dynamics
    for step in tqdm(range(n_steps - 1)):
        tumor.update()
        traces.append(dict(genotypes_counts=deepcopy(tumor.genotypes_counts)))

    # Return tumor
    return tumor, traces


def simulate_invasion(n_steps, cell, **kwargs):
    # Initialise tumor with single cell
    tumor = Tumor(cell, **kwargs)

    traces = []
    traces.append(dict(genotypes_counts=deepcopy(tumor.genotypes_counts)))
    # Simulate tumor growth
    for step in tqdm(range(n_steps - 1)):
        tumor.update()
        traces.append(dict(genotypes_counts=deepcopy(tumor.genotypes_counts)))

    # Return tumor
    return tumor, traces


def simulate_fission():
    raise NotImplementedError


def simulate_boundary():
    raise NotImplementedError
