import numpy as np


class Deme(object):
    def __init__(self, cell, carrying_capacity=1, mutation_rate=0.05, death_rate=0.98):
        self.carrying_capacity = carrying_capacity
        self.mutation_rate = mutation_rate
        self.initial_death_rate = cell.division_rate * death_rate
        self.maximum_death_rate = min(cell.max_birth_rate * 10, 0.2)
        self.death_rate = self.initial_death_rate
        self.init_cell = cell
        self.init_cell.genotype_id = str(0) + str(cell.genotype_id)
        self.cells = set()
        self.cells.add(self.init_cell)
        self.genotypes_counts = dict()
        self.genotypes_counts[cell.genotype_id] = 1
        self.genotypes_parents = dict()

    def update(self):
        # Choose subset of cells randomly
        cells = np.random.choice(
            list(self.cells), size=min(10, len(self.cells)), replace=False
        )

        # Try to apply events
        events = ["death", "division"]
        for i, cell in enumerate(cells):
            rates = [self.death_rate * (len(self.cells) > 1), cell.division_rate]
            rates_dict = dict(zip(events, rates))
            event = np.random.choice(events, p=np.array(rates) / np.sum(rates))
            success = np.random.binomial(1, rates_dict[event])
            if success:
                if event == "death":
                    self.cells.remove(cell)
                    self.genotypes_counts[cell.genotype_id] -= 1
                    del cell
                elif event == "division":
                    new_cell = cell.divide()
                    new_cell.set_params()
                    if cell.type == "tumor":
                        mutate = np.random.binomial(1, self.mutation_rate)

                        if mutate:
                            new_cell.mutate()
                            new_cell.genotype_id = str(i) + str(new_cell.genotype_id)
                            self.genotypes_counts[new_cell.genotype_id] = 1
                            self.genotypes_parents[
                                new_cell.genotype_id
                            ] = new_cell.parent.genotype_id
                            if (
                                self.genotypes_parents[new_cell.genotype_id]
                                == new_cell.genotype_id
                            ):
                                raise Exception(
                                    f"Oh no! genotype is its own parent?!: {new_cell.genotype_id}"
                                )
                        else:
                            self.genotypes_counts[new_cell.genotype_id] += 1
                    self.cells.add(new_cell)

        # Update death rate
        self.update_death_rate()

    def update_death_rate(self):
        if len(self.cells) <= self.carrying_capacity:
            self.death_rate = self.initial_death_rate
        else:
            self.death_rate = self.maximum_death_rate

    def get_genotype_frequencies(self):
        # Get unique genotypes and their frequencies
        genotypes = np.vstack([cell.snvs for cell in self.cells])
        unique, counts = np.unique(genotypes, return_counts=True, axis=0)
        freqs = counts / np.sum(counts)
        return genotypes, freqs

    def get_diversity(self):
        # Get Simpson index
        raise NotImplementedError

    def plot_grid(self):
        raise NotImplementedError

    def plot_tree(self):
        raise NotImplementedError
