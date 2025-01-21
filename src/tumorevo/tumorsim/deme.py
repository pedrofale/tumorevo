import numpy as np

from collections import Counter
import logging


class Deme(object):
    def __init__(
        self,
        cell=None,
        carrying_capacity=1,
        initial_death_rate=0.1,
        maximum_death_rate=0.5,
        tumor=None,
        row=None,
        col=None,
    ):
        if cell is None and tumor is None:
            raise ValueError(
                "Must initialise Deme with either a Cell or a Tumor object."
            )
        self.carrying_capacity = carrying_capacity
        self.initial_death_rate = initial_death_rate
        self.maximum_death_rate = maximum_death_rate
        self.tumor = tumor
        self.row = row
        self.col = col
        self.types_counts = Counter()
        self.genotypes_counts = Counter()
        self.genotypes_parents = dict()
        self.cells = set()

        if cell is not None:
            self.add_cell(cell)

    def add_cell(self, cell, genotype_id=None):
        if genotype_id is None:
            if cell.type == 'cancer':
                cell.genotype_id = str(cell.genotype_id)
            else:
                cell.genotype_id = cell.type    
        else:
            cell.genotype_id = genotype_id
        self.cells.add(cell)
        if cell.genotype_id in self.genotypes_counts:
            self.genotypes_counts[cell.genotype_id] += 1
        else:
            self.genotypes_counts[cell.genotype_id] = 1
        
        if cell.type in self.types_counts:
            self.types_counts[cell.type] += 1
        else:
            self.types_counts[cell.type] = 1            

    def update(self, treat=False, treatment_target=None, rng=None):
        cells_killed = 0
        # Choose subset of cells randomly
        l = list(self.cells)
        # l.sort(key=lambda x: x.genotype_id)
        cells = rng.choice(
            l, size=min(5, len(self.cells)), replace=False
        )

        # Try to apply events
        for i, cell in enumerate(cells):
            if cell.type == 'cancer':
                events = ["death", "division"]
                if treat and cell.snv[treatment_target] == 1:
                    rates = [cell.treatment_effectiveness, cell.division_rate]
                else:
                    rates = [self.get_death_rate(cell.death_rate), cell.division_rate]
                event = rng.choice(events, p=np.array(rates) / np.sum(rates))                    
            else: # assume non-cancer cells don't divide
                events = ['death']
                rates = [self.get_death_rate(cell.death_rate)]
                event = 'death'
                
            rates_dict = dict(zip(events, rates))
            success = rng.binomial(1, rates_dict[event])
            if cell.viability == 0:
                event = "death"
                success = 1
            if success:
                if event == "death":
                    if treat and cell.snv[treatment_target] == 1:
                        cells_killed += 1
                    pre_death_count = self.genotypes_counts[cell.genotype_id]
                    self.cells.remove(cell)
                    self.genotypes_counts[cell.genotype_id] -= 1
                    self.types_counts[cell.type] -= 1
                    
                    if self.genotypes_counts[cell.genotype_id] < 0:
                        print(cell.genotype_id)
                        print(cell.type)

                    del cell    
                elif event == "division":
                    new_cell = cell.divide()
                    new_cell.set_params()
                    if cell.type == "cancer":
                        mutate = rng.binomial(1, cell.mutation_rate)

                        if mutate:
                            new_cell.mutate(rng, self.tumor.selection.update_dict)
                            # new_cell.genotype_id = f"{i:03}" + str(
                            #     new_cell.genotype_id
                            # )  # TODO: Maybe not a great idea to depend on number of cells
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
                            self.add_cell(new_cell, genotype_id=new_cell.genotype_id)
                        else:
                            disperse = rng.binomial(1, new_cell.dispersal_rate)
                            if disperse:
                                possible_demes = self.tumor.get_neighboring_demes(self)
                                target_deme = rng.choice(possible_demes)                                    
                                target_deme.add_cell(new_cell, genotype_id=new_cell.genotype_id)
                            else:
                                self.add_cell(new_cell, genotype_id=new_cell.genotype_id)

        # Update death rate
        self.update_death_rate()

        return cells_killed

    def get_death_rate(self, cell_death_rate):
        """Update prob of each cell dieing based on its own 
        rate and the deme's carrying capacity"""
        if len(self.cells) <= self.carrying_capacity:
            return cell_death_rate
        else:
            return min(cell_death_rate * self.carrying_capacity, self.maximum_death_rate)

    def update_death_rate(self):
        if len(self.cells) <= self.carrying_capacity:
            self.death_rate = self.initial_death_rate
        else:
            self.death_rate = self.maximum_death_rate

    def get_genotype_frequencies(self, normalize=True):
        # Get unique genotypes and their frequencies
        genotypes = np.vstack([cell.snv for cell in self.cells])
        unique, counts = np.unique(genotypes, return_counts=True, axis=0)
        freqs = np.array(counts)
        if normalize:
            freqs = freqs / np.sum(freqs)
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
