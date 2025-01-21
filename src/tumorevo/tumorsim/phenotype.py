import numpy as np

class Phenotype(object):
    def __init__(self, n_genes,  driver_effects=1.1, resistant_effects=0.1):
        self.n_genes = n_genes

        self.make_drivers()
        self.make_resistant()
        self.make_expmap()

    def make_expmap(self):
        # Effect of snv on exp: up or down depending on wether tsg or og
        self.mul_effect = np.ones((self.n_genes,))
        self.mul_effect[np.where(self.driver_types == 1)] = 2.
        self.mul_effect[np.where(self.driver_types == -1)] = 0.5
