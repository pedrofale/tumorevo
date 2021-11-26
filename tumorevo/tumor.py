class Tumor(object):
        def __init__(self, n_healthy_cells=1, grid_size=10, carrying_capacity=1, seed=42):
                n_cells = n_cancer_cells + n_healthy_cells
                assert n_cells <= grid_size*grid_size
                self.n_cancer_cells = n_cancer_cells
                self.n_healthy_cells = n_healthy_cells
                self.grid_size = grid_size
                self.carrying_capacity = carrying_capacity
                self.seed = seed

                self.grid = np.zeros((self.grid_size, self.grid_size)) # Contains the deme IDs

                self.cells = dict()

                init_cell_positions = np.random.choice(n_cells, replace=False, n=n_cells)

                for i in range(self.n_cancer_cells):
                        self.cells[i] = TumorCell(cell_id=i)

                for j in range(self.n_healthy_cells):
                        i = self.n_cancer_cells + j
                        self.cells[j] == Cell(cell_id=i)


        def update(self):
                deme = np.random.choice(demes)
                celltype = np.random.choice(['normal', 'tumor'])
                genotype = np.random.choice(deme.genotypes)
                cell = np.random.choice(deme.cells.genotypes)
                event = np.random.choice(cell)
                for cell in self.cells:
                        # Apply rules to each cell, which depends on the neighborhood
                        move_probs = cell.move_probs()

