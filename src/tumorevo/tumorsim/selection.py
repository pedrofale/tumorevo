import numpy as np

class Selection(object):
    def __init__(self, n_segments=10, segment_size=1000, prop_driver=0.1, prop_resistance=0.1,
                 driver_effects=1.1, resistant_effects=1.1,
                 max_ploidy=6, max_cn=12, max_nullisomy=2, max_mut_drivers=1000, ):
        # Fixed about the genome
        self.n_segments = n_segments
        self.segment_size = segment_size
        self.prop_driver = prop_driver
        self.prop_resistance = prop_resistance
        
        # Fixed about fitness
        self.driver_effects = driver_effects
        self.resistant_effects = resistant_effects

        # Fixed about viability
        self.max_ploidy = max_ploidy
        self.max_cn = max_cn
        self.max_nullisomy = max_nullisomy
        self.max_mut_drivers = max_mut_drivers

        self.drivers = []
        self.passengers = []

        # Put drivers and passengers in position
        self.make_drivers()
        self.make_expmap()
        self.update_dict = {'viability': self.update_viability,
                            'division_rate': self.update_division_rate,
                            'dispersal_rate': self.update_dispersal_rate,
                            'treatment_effectiveness': self.update_treatment_effectiveness}

    def make_drivers(self): 
        # supressor, notdriver, oncogene
        self.drivers = []
        self.passengers = []
        self.driver_types = []
        for _ in range(self.n_segments):
            driver_types = np.random.choice([-1,0,1], p=[self.prop_driver/2, 1.-self.prop_driver, self.prop_driver/2], 
                                                size=self.segment_size)
            self.drivers.append(np.where(driver_types!=0)[0])
            self.passengers.append(np.where(driver_types==0)[0])
            self.driver_types.append(driver_types)
    
    def make_expmap(self):
        # Effect of snv on exp: up or down depending on wether tsg or og
        self.mul_effect = []
        for _ in range(self.n_segments):
            mul_effect = np.ones((self.segment_size,))
            mul_effect[np.where(self.driver_types == 1)] = 2.
            mul_effect[np.where(self.driver_types == -1)] = 0.5
            self.mul_effect.append(mul_effect)

    def update_exp(self, seg_baseline_exp, seg_idx, allele_muts):
        return seg_baseline_exp * self.mul_effect[seg_idx] * allele_muts

    def get_tsgs(self):
        tsgs = []
        for seg in range(self.n_segments):
            idx = np.where(self.driver_types[seg] == -1)[0]
            tsgs.append(idx + seg*self.segment_size)
        tsgs = np.concatenate(tsgs)
        return tsgs

    def get_oncogenes(self):
        tsgs = []
        for seg in range(self.n_segments):
            idx = np.where(self.driver_types[seg] == 1)[0]
            tsgs.append(idx + seg*self.segment_size)
        tsgs = np.concatenate(tsgs)
        return tsgs

    def make_resistant(self):
        self.confers_resistance = np.random.binomial(1, self.prop_resistance, size=self.n_genes)
        self.confers_resistance[self.drivers] = 0.

    def update_viability(self, genome):
        avg_ploidy = 1
        if avg_ploidy > self.max_ploidy:
            return 0
        highest_cn = 1
        if highest_cn > self.max_cn:
            return 0 
        nullisomy_count = 1
        if nullisomy_count > self.max_nullisomy:
            return 0
        n_mutated_drivers = 1 
        if n_mutated_drivers > self.max_mut_drivers:
            return 0
        return 1

    def update_division_rate(self, baseline_fitness, genome):
        # Make it sensible to copy number
        fitness = baseline_fitness
        for seg in range(len(genome)):
            for hap in genome[seg]:
                for all in genome[seg][hap]:
                    fitness += self.driver_effects * len(all.intersection(self.drivers[seg])) # more mutated copies of a gene shouldn't make a big difference though
        fitness = np.min([1., fitness])                    
        return fitness

    def update_dispersal_rate(self, baseline_fitness, genome):
        # Make it sensible to copy number
        fitness = baseline_fitness
        for seg in range(len(genome)):
            for hap in genome[seg]:
                for all in genome[seg][hap]:
                    fitness += self.driver_effects * len(all.intersection(self.drivers[seg])) # more mutated copies of a gene shouldn't make a big difference though
        fitness = np.min([1., fitness])
        return fitness
    
    def update_treatment_effectiveness(self, baseline_fitness, genome):
        # Make it sensible to copy number
        fitness = baseline_fitness
        for seg in range(len(genome)):
            for hap in genome[seg]:
                for all in genome[seg][hap]:
                    fitness += self.driver_effects * len(all.intersection(self.drivers[seg])) # more mutated copies of a gene shouldn't make a big difference though
        fitness = np.min([1., fitness])                    
        return fitness
    
