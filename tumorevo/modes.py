from .tumor import Tumor
from .deme import Deme
from .cell import Cell, TumorCell

from tqdm import tqdm
from copy import deepcopy

def simulate_nonspatial(n_steps, cell, **kwargs):
	# Create Deme
	deme = Deme(cell, **kwargs)
	
	traces = []
	traces.append(dict(genotypes_counts=deepcopy(deme.genotypes_counts)))
	# Simulate within-deme dynamics
	for step in tqdm(range(n_steps-1)):
		deme.update()
		traces.append(dict(genotypes_counts=deepcopy(deme.genotypes_counts)))
	
	# Return Deme
	return deme, traces

def simulate_invasion():
	raise NotImplementedError
	
def simulate_fission():
	raise NotImplementedError

def simulate_boundary():
	raise NotImplementedError
