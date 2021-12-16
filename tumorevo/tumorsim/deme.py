import numpy as np

from collections import Counter
import logging


class Deme(object):
    def __init__(
        self,
        cell=None,
        carrying_capacity=1,
        division_rate=0.1,
        max_birth_rate=0.5,
        mutation_rate=0.05,
        dispersal_rate=0.01,
        death_rate=0.98,
        tumor=None,
        x=None,
        y=None,
    ):
        if cell is None and tumor is None:
            raise ValueError(
                "Must initialise Deme with either a Cell or a Tumor object."
            )
        self.carrying_capacity = carrying_capacity
        self.division_rate = division_rate
        self.mutation_rate = mutation_rate
        self.dispersal_rate = dispersal_rate
        self.initial_death_rate = division_rate * death_rate
        self.death_rate = self.initial_death_rate
        self.maximum_death_rate = min(max_birth_rate * 10, 0.2)
        self.tumor = tumor
        self.x = x
        self.y = y
        self.genotypes_counts = Counter()
        self.genotypes_parents = dict()
        self.cells = set()

        if cell is not None:
            self.add_cell(cell)

    def add_cell(self, cell):
        cell.genotype_id = str(0) + str(cell.genotype_id)
        self.cells.add(cell)
        self.genotypes_counts[cell.genotype_id] = 1
        self.initial_death_rate = cell.division_rate * self.death_rate
        self.maximum_death_rate = min(cell.max_birth_rate * 10, 0.2)
        self.death_rate = self.initial_death_rate

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
                            new_cell.genotype_id = str(i) + str(
                                new_cell.genotype_id
                            )  # TODO: Maybe not a great idea to depend on number of cells
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
                            self.cells.add(new_cell)
                        else:
                            disperse = np.random.binomial(1, self.dispersal_rate)
                            if disperse:
                                possible_demes = self.tumor.get_neighboring_demes(self)
                                target_deme = np.random.choice(possible_demes)
                                target_deme.cells.add(new_cell)
                                if new_cell.genotype_id in target_deme.genotypes_counts:
                                    target_deme.genotypes_counts[
                                        new_cell.genotype_id
                                    ] += 1
                                else:
                                    target_deme.genotypes_counts[
                                        new_cell.genotype_id
                                    ] = 1
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

    def get_most_frequent_genotype(self):
        return self.genotypes_counts.most_common(1)[0][0]

    def get_diversity(self):
        # Get Simpson index
        raise NotImplementedError

    def plot_grid(self):
        raise NotImplementedError

    def plot_tree(self):
        raise NotImplementedError
