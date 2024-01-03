import numpy as np
from copy import copy


class Cell(object):
    def __init__(
        self,
        n_genes=1000,
        parent=None,
        division_rate=0.1,
        death_rate=0.1,
        max_birth_rate=0.3,
        dispersal_rate=0.1,
    ):
        self.type = "healthy"
        self.parent = parent
        self.n_genes = n_genes
        if parent is None:
            self.snvs = np.zeros((n_genes,))  # 0 is no mutation
            self.set_genotype_id()
        else:
            self.snvs = np.array(parent.snvs)
            self.genotype_id = self.parent.genotype_id
        self.exp = np.zeros((n_genes,))  # Over/under expressed level

        self.division_rate = division_rate
        self.max_birth_rate = max_birth_rate
        self.death_rate = death_rate
        self.dispersal_rate = dispersal_rate

    def set_genotype_id(self):
        self.genotype_id = "".join(self.snvs.astype(int).astype(str))

    def divide(self, new_cell_id=0):
        new_cell = copy(self)
        new_cell.parent = self
        return new_cell


class TumorCell(Cell):
    def __init__(self, snv_prob=0.1, prop_driver=0.1, prop_resistance=0.1, genotype_id=None, seed=42, **cell_kwargs):
        super(TumorCell, self).__init__(**cell_kwargs)
        self.type = "tumor"
        self.prop_driver = prop_driver
        self.prop_resistance = prop_resistance
        self.snv_prob = snv_prob
        self.driver_effects = 0.2
        self.treatment_effectiveness = 0.8
        self.seed  = seed
        self.set_params()

    def set_params(self):
        if self.parent is None:
            rng = np.random.default_rng(self.seed)
            self.is_driver = rng.binomial(1, self.prop_driver, size=self.n_genes)
            self.snv_probs = np.ones((self.n_genes,)) * self.snv_prob
            # Infinite sites assumption
            self.snv_probs[np.where(self.snvs != 0)] = 0.0
            # Treatment resistance -- can't overlap with drivers
            rng = np.random.default_rng(self.seed+1)
            self.confers_resistance = rng.binomial(1, self.prop_resistance, size=self.n_genes)
            self.confers_resistance[np.where(self.is_driver)] = 0.
        else:
            self.is_driver = np.array(self.parent.is_driver)
            self.snv_probs = np.array(self.parent.snv_probs)
            self.genotype_id = self.parent.genotype_id
            self.confers_resistance = np.array(self.parent.confers_resistance)

    def mutate(self, rng, n_events=1):
        mutated_genes = rng.choice(
            self.n_genes,
            p=self.snv_probs / np.sum(self.snv_probs),
            replace=False,
            size=n_events,
        )
        self.snvs[mutated_genes] = 1
        self.snv_probs[np.where(self.snvs != 0)] = 0.0

        # Apply driver and resistance effects
        for gene in mutated_genes:
            if self.is_driver[gene]:
                print("Cell has driver mutation!")
                self.division_rate += self.division_rate * self.driver_effects
            if self.confers_resistance[gene]:
                print("Cell has resistance-inducing mutation!")
                self.treatment_effectiveness *= 0.1

        self.division_rate = min(self.max_birth_rate, self.division_rate)
        self.set_genotype_id()
