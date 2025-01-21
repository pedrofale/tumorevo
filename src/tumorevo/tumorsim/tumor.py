from .deme import Deme
from .cell import EpithelialCell, StromalCell
from ..constants import *

import numpy as np
import pandas as pd

from collections import Counter
import logging

def bresenham_circumference(x0, y0, radius):
    x = radius
    y = 0
    err = 0

    points = []

    while x >= y:
        points.append((x0 + x, y0 + y))
        points.append((x0 + y, y0 + x))
        points.append((x0 - y, y0 + x))
        points.append((x0 - x, y0 + y))
        points.append((x0 - x, y0 - y))
        points.append((x0 - y, y0 - x))
        points.append((x0 + y, y0 - x))
        points.append((x0 + x, y0 - y))

        y += 1
        err += 1 + 2*y
        if 2*(err-x) + 1 > 0:
            x -= 1
            err += 1 - 2*x

    points = pd.DataFrame(points).drop_duplicates().values
    points = [(p[0], p[1]) for p in points]
    return points

def get_inside(border):
    points = []
    for (x,y) in border:
        # Add points inside -- general for any shape
        # find another x with this y
        for b in border:
            if b[1] == y:
                if b[0]-x > 1: # if to the right
                    points.extend([(x_,y) for x_ in range(x+1,b[0])])
                if b[0]-x > 1: # if to the left
                    points.extend([(x_,y) for x_ in range(b[0],x-1)])
    points = pd.DataFrame(points).drop_duplicates().values
    points = [(p[0], p[1]) for p in points]
    return points

class Tumor(object):
    def __init__(self, cancer_cell, selection, n_structures=1, structure_radius=0, grid_size=10, epithelial_cell_params=dict(), stromal_cell_params=dict(), immune_cell_params=dict(), deme_params=dict()):
        self.grid_size = grid_size
        self.selection = selection
        self.celltype_exps = dict()
        self.n_genes = self.selection.n_segments * self.selection.segment_size
        self.celltypes = ['cancer', 'epithelial', 'stromal', 'immune']
        self.make_celltype_exps()
            
        # Initialize grid of empty demes
        self.positions = []
        self.grid = []
        self.deme_list = []
        for grid_row in range(self.grid_size):
            row = []
            for grid_col in range(self.grid_size):
                deme = Deme(
                    tumor=self,
                    row=grid_row,
                    col=grid_col,
                    **deme_params
                )
                self.deme_list.append(deme)
                row.append(deme)
            self.grid.append(row)

        # Initialize cell positions
        center = int(self.grid_size / 2)
        if structure_radius <= 0:
            # Put cancer cell in center deme    
            self.grid[center][center].add_cell(cancer_cell)
        else:
            structure_borders = []
            structure_in_borders = []
            structure_circles = []
            for s_idx in range(n_structures):
                # Get border
                border = bresenham_circumference(center, center, structure_radius)
                structure_borders.append(border)

                # add healthy epithelial cells in border of shape
                for (row,col) in border:
                    # Fill them up to carrying capacity
                    i = 0
                    while i < self.grid[row][col].carrying_capacity:
                        epithelial_cell = EpithelialCell(**epithelial_cell_params)
                        self.grid[row][col].add_cell(epithelial_cell)
                        i += 1

                # Get inside
                circle = get_inside(border)
                structure_circles.append(circle)

                in_border = bresenham_circumference(center, center, structure_radius-1)                        
                structure_in_borders.append(in_border)
                if s_idx == 0:
                    # add a cancer cell inside the border 
                    pos = np.random.choice(len(in_border))
                    self.grid[in_border[pos][0]][in_border[pos][1]].add_cell(cancer_cell)

            # add healthy stromal cells outside of border
            for row in range(grid_size):
                for col in range(grid_size):
                    if (row,col) not in structure_circles and (row,col) not in structure_borders:
                        # Fill them up to carrying capacity
                        i = 0
                        while i < self.grid[row][col].carrying_capacity:
                            stromal_cell = StromalCell(**stromal_cell_params)
                            self.grid[row][col].add_cell(stromal_cell)
                            i += 1                    

        self.genotypes_parents = dict()
        self.genotypes_counts = Counter()

        self.update_genotype_counts()

    def make_celltype_exps(self):
        self.celltype_exps = dict()
        for celltype in self.celltypes:
            exp = np.random.beta(.1, 1., size=self.n_genes)
            exp[np.where(self.selection.get_tsgs())] = 0.8
            exp[np.where(self.selection.get_oncogenes())] = 0.01 
            self.celltype_exps[celltype] = exp

    def get_genotype_frequencies(self, normalize=True):
        # Get unique genotypes and their frequencies
        genotypes = list(self.genotypes_counts.keys())
        snvs = []
        counts = []
        for genotype in genotypes:
            if genotype not in normal_names:
                snvs.append(np.array(list(genotype), dtype=int)[1:])
                counts.append(self.genotypes_counts[genotype])
        freqs = np.array(counts)
        if normalize:
            freqs = freqs / np.sum(freqs)
        return snvs, freqs

    def get_deme_genotype_frequencies(self, normalize=True):
        # Get unique genotypes and their frequencies per deme
        # Create grid containing genotype freqs at each deme
        grid = []
        for grid_row in range(self.grid_size):
            row = []
            for grid_col in range(self.grid_size):
                if len(self.grid[grid_row][grid_col].cells) > 0:
                    row.append(self.grid[grid_row][grid_col].get_genotype_frequencies(normalize=normalize))
                else:
                    row.append("")
            grid.append(row)
        return grid

    def get_genotype_matrix(self):
        # Create grid containing most frequent genotype at each deme
        grid = []
        for grid_row in range(self.grid_size):
            row = []
            for grid_col in range(self.grid_size):
                if len(self.grid[grid_row][grid_col].cells) > 0:
                    row.append(self.grid[grid_row][grid_col].get_most_frequent_genotype())
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
        grid_row = deme.row
        grid_col = deme.col

        possible_demes = []
        pos = []
        # Von Neumann neighborhood
        for tup in [(grid_row - 1, grid_col), (grid_row, grid_col + 1), (grid_row + 1, grid_col), (grid_row, grid_col - 1)]:
            if tup[0] > 0 and tup[0] < self.grid_size:
                if tup[1] > 0 and tup[1] < self.grid_size:
                    pos.append(tup)
                    possible_demes.append(self.grid[tup[0]][tup[1]])

        # Other structure
        # for in_border in self.structure_in_borders:
        #     for tup in in_border:
        #         possible_demes.append(self.grid[tup[0]][tup[1]])

        return possible_demes

    def update(self, treat=False, treatment_target=None, rng=None):
        cells_killed = 0

        deme_list = [deme for deme in self.deme_list if deme.types_counts['cancer'] > 0]
        demes = rng.choice(deme_list, size=min(10, len(deme_list)), replace=False)
        for deme in demes:
            cells_killed += deme.update(treat=treat, treatment_target=treatment_target, rng=rng)

        self.update_genotype_parents()
        self.update_genotype_counts()

        return cells_killed

    def set_cell_exps(self):
        for deme in self.deme_list:
            for cell in deme.cells:        
                if cell.type == 'cancer':
                    cell.baseline_exp = self.celltype_exps['cancer']
                    cell.exp = cell.get_exp(self.selection.update_exp)
                else:
                    cell.exp = self.celltype_exps[cell.type]

    def get_cell_data(self,):
        n_cells = sum(self.genotypes_counts.values())
        cell_gen = []
        cell_exp = np.zeros((n_cells, self.selection.segment_size*self.selection.n_segments))
        cell_crd = np.zeros((n_cells, 2))
        cell_names = []
        cell_ids = []
        i = 0
        for deme in self.deme_list:
            for cell in deme.cells:
                cell_name = f'C{i}'
                cell_gen[i] = cell.get_genome_df()
                cell_gen[i]['cell'] = cell_name
                cell_exp[i] = cell.exp
                cell_crd[i][0], cell_crd[i][1] = deme.row, deme.col # should be after expanding demes
                cell_ids.append(cell.genotype_id)
                cell_names.append(cell_name)
                i += 1

        gene_names = []
        for segment in range(self.selection.n_segments):
            gn = [f'G{segment}_{i}' for i in range(self.selection.segment_size)]
            gene_names.extend(gn)

        cell_gn_df = pd.concat(cell_gen)

        return dict(cell_gen=cell_gn_df, 
                    cell_exp=pd.DataFrame(cell_exp.astype(float), index=cell_names, columns=gene_names), 
                    cell_crd=pd.DataFrame(cell_crd.astype(int), index=cell_names, columns=['row', 'col']),
                    cell_ids=pd.DataFrame(cell_ids, index=cell_names, columns=['cell_id']))

    def get_gene_data(self):
        gene_names = [f'G{i}' for i in range(self.n_genes)]
        return dict(driver_types=pd.DataFrame(self.selection.driver_types, index=gene_names), 
                    confers_resistance=pd.DataFrame(self.selection.confers_resistance, index=gene_names))