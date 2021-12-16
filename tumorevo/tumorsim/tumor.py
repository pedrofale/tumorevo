from .deme import Deme

import numpy as np

from collections import Counter
import logging


class Tumor(object):
    def __init__(self, cell, grid_size=10, carrying_capacity=1, seed=42, **kwargs):
        self.grid_size = grid_size
        self.carrying_capacity = carrying_capacity
        self.seed = seed

        # Initialize grid of empty demes
        self.positions = []
        self.grid = []
        self.deme_list = []
        for x in range(self.grid_size):
            row = []
            for y in range(self.grid_size):
                deme = Deme(
                    tumor=self,
                    x=x,
                    y=y,
                    division_rate=cell.division_rate,
                    max_birth_rate=cell.max_birth_rate,
                    carrying_capacity=self.carrying_capacity,
                    **kwargs
                )
                self.deme_list.append(deme)
                row.append(deme)
            self.grid.append(row)

        # Put tumor cell in center deme
        center = int(self.grid_size / 2)
        self.grid[center][center].add_cell(cell)
        self.genotypes_parents = dict()
        self.genotypes_counts = Counter()

    def get_genotype_frequencies(self):
        # Get unique genotypes and their frequencies
        genotypes = list(self.genotypes_counts.keys())
        snvs = []
        counts = []
        for genotype in genotypes:
            snvs.append(np.array(list(genotype), dtype=int)[1:])
            counts.append(self.genotypes_counts[genotype])
        freqs = np.array(counts) / np.sum(counts)
        return snvs, freqs

    def get_genotype_matrix(self):
        # Create grid containing most frequent genotype at each deme
        grid = []
        for x in range(self.grid_size):
            row = []
            for y in range(self.grid_size):
                if len(self.grid[x][y].cells) > 0:
                    row.append(self.grid[x][y].get_most_frequent_genotype())
                else:
                    row.append("")
            grid.append(row)
        return grid

    def update_genotype_parents(self):
        self.genotypes_parents = dict()
        for deme in self.deme_list:
            self.genotypes_parents.update(deme.genotypes_parents)

    def update_genotype_counts(self):
        self.genotypes_counts = Counter()
        for deme in self.deme_list:
            self.genotypes_counts = self.genotypes_counts + deme.genotypes_counts

    def get_neighboring_demes(self, deme):
        x = deme.x
        y = deme.y

        possible_demes = []
        pos = []
        for tup in [(x - 1, y), (x, y + 1), (x + 1, y), (x, y - 1)]:
            if tup[0] > 0 and tup[0] < self.grid_size:
                if tup[1] > 0 and tup[1] < self.grid_size:
                    pos.append(tup)
                    possible_demes.append(self.grid[tup[0]][tup[1]])

        return possible_demes

    def update(self):
        deme_list = [deme for deme in self.deme_list if len(deme.cells) > 0]
        demes = np.random.choice(deme_list, size=min(10, len(deme_list)), replace=False)
        for deme in demes:
            deme.update()

        self.update_genotype_parents()
        self.update_genotype_counts()
