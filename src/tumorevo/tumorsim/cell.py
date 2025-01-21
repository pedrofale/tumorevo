import numpy as np
import pandas as pd
from copy import copy


class Cell(object):
    def __init__(
        self,
        n_segments=10,
        segment_size=1000,
        parent=None,
        division_rate=0.1,
        death_rate=0.1,
        max_birth_rate=0.3,
        dispersal_rate=0.1,
    ):
        self.n_segments = n_segments
        self.segment_size = segment_size
        self.type = "healthy"
        self.parent = parent
        if parent is None:
            self.genome = [{'p':[set()], 'm':[set()]}] * n_segments # copy number = 2
            self.set_genotype_id()
        else:
            self.genome = list(parent.genome)
            self.genotype_id = self.parent.genotype_id
        # self.exp = np.random.beta(.1, 1, size=self.n_genes)  # Gene activity probability (each gene ranges from 0 to 1 indicating its prob of expression. transcripts will be sampled binomially)

        self.division_rate = division_rate
        self.max_birth_rate = max_birth_rate
        self.death_rate = death_rate
        self.dispersal_rate = dispersal_rate
        self.baseline_treatment_effectiveness = self.death_rate

        self.viability = 1.
        self.treatment_effectiveness = self.baseline_treatment_effectiveness

    def update_evolutionary_parameters(self, update_dict):
        self.viability = update_dict['viability'](self.genome)
        self.division_rate = update_dict['division_rate'](self.baseline_division_rate, self.genome)
        self.dispersal_rate = update_dict['dispersal_rate'](self.baseline_dispersal_rate, self.genome)
        self.treatment_effectiveness = update_dict['treatment_effectiveness'](self.baseline_treatment_effectiveness, self.genome)
        self.division_rate = min(self.max_birth_rate, self.division_rate)            

    def set_genotype_id(self):
        self.genotype_id = id(self)#"".join(self.snv.astype(int).astype(str))

    def divide(self, new_cell_id=0):
        new_cell = copy(self)
        new_cell.parent = self
        return new_cell
    
    def get_genome_df(self):
        rows = []
        for seg in range(self.n_segments):
            for hap in self.genome[seg]:
                for all in self.genome[seg][hap]:
                    for pos in range(self.segment_size):
                        mut = 0
                        if pos in self.genome[seg][hap][all]:
                            mut = 1
                        row = [seg, hap, pos, mut]
                        rows.append(row)
        df = pd.DataFrame.from_records(rows, columns=['seg', 'hap', 'pos', 'mut'])
        return df

    def set_baseline_exp(self):
        self.baseline_exp = np.random.beta(.1, 1, size=self.n_segments * self.segment_size) 

    def get_exp(self, make_exp_fun):
        exp = np.array(self.baseline_exp)
        for seg in range(self.n_segments):
            seg_baseline = self.baseline_exp[seg*self.segment_size:(seg+1)*self.segment_size]
            seg_exp = np.array(seg_baseline)
            if len(self.genome[seg]['p']) + len(self.genome[seg]['m']) == 0:
                seg_exp *= 0.
            else:
                for hap in self.genome[seg]:
                    for all in range(len(self.genome[seg][hap])):
                        allele_contrib = make_exp_fun(seg_baseline, seg, self.genome[seg][hap][all]) # depending on muts
                        seg_exp += allele_contrib
            exp[seg*self.segment_size:(seg+1)*self.segment_size] = seg_exp
        return exp

class EpithelialCell(Cell):
    def __init__(self, **cell_kwargs):
        super(EpithelialCell, self).__init__(**cell_kwargs)
        self.type = "epithelial"
        self.genotype_id = self.type
        # self.exp = np.random.beta(.1, 1, size=self.n_genes)  # Gene activity probability (each gene ranges from 0 to 1 indicating its prob of expression. transcripts will be sampled binomially)

class StromalCell(Cell):
    def __init__(self, **cell_kwargs):
        super(StromalCell, self).__init__(**cell_kwargs)
        self.type = "stromal"        
        self.genotype_id = self.type
        # self.exp = np.random.beta(.1, 1, size=self.n_genes)  # Gene activity probability (each gene ranges from 0 to 1 indicating its prob of expression. transcripts will be sampled binomially)

class ImmuneCell(Cell):
    def __init__(self, prob_kill=.01, **cell_kwargs):
        super(ImmuneCell, self).__init__(**cell_kwargs)
        self.prob_kill = prob_kill # probability of killing a neighboring cancer cell
        self.type = "immune"
        self.genotype_id = self.type
        # self.exp = np.random.beta(.1, 1, size=self.n_genes)  # Gene activity probability (each gene ranges from 0 to 1 indicating its prob of expression. transcripts will be sampled binomially)

class CancerCell(Cell):
    def __init__(self, mutation_rate=0.1, snv_prob=0.1, genotype_id=None, seed=42, **cell_kwargs):
        super(CancerCell, self).__init__(**cell_kwargs)
        self.type = "cancer"
        self.mutation_rate = mutation_rate
        self.seed  = seed
        self.baseline_division_rate = self.division_rate
        self.baseline_dispersal_rate = self.dispersal_rate
        self.set_params()

    def set_params(self):
        if self.parent is not None:
            self.genotype_id = self.parent.genotype_id
            
    def mutate(self, rng, update_dict, n_events=5, mut_prob=.1, cnv_prob=.1):
        event = np.random.choice(['cnv', 'mut'], p=np.array([mut_prob, cnv_prob])/sum([mut_prob, cnv_prob])) # add WGDs too...
        if event == 'mut':
            # Sample segment
            segment_probs = np.array([len(self.genome[seg]['p'] + self.genome[seg]['m']) for seg in range(self.n_segments)]) # can't select empty segment
            segment_probs = segment_probs / np.sum(segment_probs)
            seg = np.random.choice(range(self.n_segments), p=segment_probs)
            # Sample haplotype
            n_p, n_m = len(self.genome[seg]['p']), len(self.genome[seg]['m'])
            hap = np.random.choice(['p', 'm'], p=np.array([n_p, n_m])/(n_p + n_m))
            # Sample allele
            all = np.random.choice(range(len(self.genome[seg][hap])))
            # Sample number of mutations
            n_mutations = np.min([np.random.poisson(n_events)+1, self.segment_size])
            # choose mutations and reject the ones which are already in the allele (force ISA)
            muts = set(np.random.choice(self.segment_size, size=n_mutations, replace=False))
            while len(muts.intersection(self.genome[seg][hap][all])) == len(muts):
                muts = set(np.random.choice(self.segment_size, size=n_mutations, replace=False))
            self.genome[seg][hap][all].update(muts)
        elif event == 'cnv':
            # Sample segment
            segment_probs = np.array([len(self.genome[seg]['p'] + self.genome[seg]['m']) for seg in range(self.n_segments)]) # can't select empty segment
            segment_probs = segment_probs / np.sum(segment_probs)
            seg = np.random.choice(range(self.n_segments), p=segment_probs)
            # Sample haplotype
            n_p, n_m = len(self.genome[seg]['p']), len(self.genome[seg]['m'])
            hap = np.random.choice(['p', 'm'], p=np.array([n_p, n_m])/(n_p + n_m))
            # Sample allele
            all = np.random.choice(range(len(self.genome[seg][hap])))
            # Decide wether to delete or copy
            evt = np.random.choice(['del', 'amp'])
            if evt == 'amp':
                self.genome[seg][hap].append(self.genome[seg][hap][all]) # add a copy
            elif evt == 'del':
                if len(self.genome[seg][hap]) == 1:
                    self.genome[seg][hap][all] = set()
                else:
                    del self.genome[seg][hap][all] # remove

        self.update_evolutionary_parameters(update_dict)
        self.set_genotype_id()

    # def mutate(self, rng, update_dict, n_events=1):
    #     mutated_genes = rng.choice(
    #         self.n_genes,
    #         p=self.snv_probs / np.sum(self.snv_probs),
    #         replace=False,
    #         size=n_events,
    #     )
    #     # Randomly pick an allele
    #     allele = np.random.choice(len(self.genome))
    #     self.genome[allele][mutated_genes] = 1
    #     self.snv_probs[np.where(self.snv != 0)] = 0.0
    #     self.division_rate = update_dict['division_rate'](self.baseline_division_rate, self.genome)
    #     self.treatment_effectiveness = update_dict['treatment_effectiveness'](self.baseline_treatment_effectiveness, self.genome)
    #     self.exp = update_dict['exp'](self.baseline_exp, self.genome)

    #     self.division_rate = min(self.max_birth_rate, self.division_rate)
    #     self.set_genotype_id()
