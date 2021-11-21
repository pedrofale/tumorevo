"""
Simulate tumor growth under different spatial models according to Noble et al, 2019.
"""
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--mode', default=0)
parser.add_argument('--max_n_cells', default=100)
parser.add_argument('--n_genes', default=10)
parser.add_argument('--n_steps', default=10)
parser.add_argument('--grid_side', default=10)

class Cell(object):
	def __init__(self, cell_id=0, n_genes=1000, parent=None, division_rate=.1, death_rate=.1):
		self.cell_id = cell_id
		self.parent = parent
		self.n_genes = n_genes
		self.snvs = np.zeros((n_genes,)) # 0 is no mutation
		self.exp = np.zeros((n_genes,)) # Over/under expressed level
		
		self.division_rate = division_rate
		self.death_rate = death_rate

	def divide(self, new_cell_id=0):
		new_cell = deepcopy(self)
		new_cell.parent = self
		new_cell.cell_id = new_cell_id
		return new_cell

class TumorCell(Cell):
	def __init__(self, snv_prob=.1, dispersal_rate=.1):
		super(TumorCell, self).__init__()
		
		self.snv_probs = np.ones((self.n_genes,)) * snv_prob		

		# Infinite sites assumption
		self.snv_probs[np.where(self.snvs != 0)] = 0.

		self.dispersal_rate = dispersal_rate

	def mutate(self, n_events=1):
		is_mutated = np.random.binomial(1, self.n_probs)
		mutated_genes = np.random.choice(np.where(is_mutated), replace=False, n=min(len(is_mutated), n_events)
		snvs[mutated_genes] = 1


class Deme(object):
	def __init__(self, carrying_capacity=1):
		self.carrying_capacity = carrying_capacity
Driver mutation rate per cell division 10−5
Passenger mutation rate per cell division 0.1
Normal cell relative division rate 0.9
Mean value of driver effect on cell division rate 0.1
Maximum relative cell division rate 10
Maximum relative dispersal rate 10
Dispersal rate Conditional


class Tumor(object):
	def __init__(self, n_healthy_cells=1, grid_size=10, carrying_capacity=1, seed=42):
		n_cells = n_cancer_cells + n_healthy_cells
		assert n_cells <= grid_size*grid_size
		self.n_cancer_cells = n_cancer_cells
		self.n_healthy_cells = n_healthy_cells
		self.grid_size = grid_size
		self.carrying_capacity = carrying_capacity
		self.seed = seed

		self.grid = np.zeros((self.grid_size, self.grid_size)) # Contains the cell IDs

		self.cells = dict()

		init_cell_positions = np.random.choice(n_cells, replace=False, n=n_cells)

		for i in range(self.n_cancer_cells):
			self.cells[i] = TumorCell(cell_id=i)
			self.grid[

		for j in range(self.n_healthy_cells):
			i = self.n_cancer_cells + j
			self.cells[j] == Cell(cell_id=i)			
		
Tumour cells undergo stochastic division, death, dispersal, and mutation events, whereas normal cells undergo
only division and death
	def update(self):
		for cell in self.cells:
			# Apply rules to each cell, which depends on the neighborhood
			move_probs = cell.move_probs()
			move = np.random.choice(
			
				 		

	def get_event_tree(self):
		 

def main():
	tumor = Tumor()
	for step in range(n_steps):
		tumor.step()

if __name__ == 'main':
	main()
