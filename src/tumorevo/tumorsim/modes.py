from .tumor import Tumor

import numpy as np
from tqdm import tqdm
from copy import deepcopy


def simulate_nonspatial(n_steps, tumor, traces=None, seed=42, **kwargs):
    if traces is None:
        traces = []
        traces.append(dict(genotypes_counts=deepcopy(tumor.genotypes_counts)))

    # Simulate within-deme dynamics
    for step in tqdm(range(n_steps - 1)):
        rng = np.random.default_rng(seed + step)
        tumor.update(rng=rng)
        traces.append(dict(genotypes_counts=deepcopy(tumor.genotypes_counts)))

    # Return tumor
    return tumor, traces, -1


def simulate_invasion(n_steps, tumor, traces=None, treatment_duration=10, treatment_iteration=-1, treatment_target=-1, cells_killed=0, seed=42, **kwargs):
    if traces is None:
        traces = []
        traces.append(dict(genotypes_counts=deepcopy(tumor.genotypes_counts)))

    # Simulate tumor growth
    for step in tqdm(range(n_steps - 1)):

        treat = False
        if len(traces) == treatment_iteration:
            treat = True
            # For all genotypes which are alive
            active_genotypes = np.concatenate([np.tile(np.array(list(genotype)).astype(int), (tumor.genotypes_counts[genotype], 1)) for genotype in tumor.genotypes_counts.keys() if tumor.genotypes_counts[genotype] > 0])
            print(active_genotypes)
            print(active_genotypes.shape)
            # TODO: treatment types - all cells, most common mutation, immunotherapy
            # Most common mutation
            treatment_target = np.argmax(np.sum(active_genotypes[:,3:] == 1, axis=0)) # Ignore first indices which is the generation, used to order the genotypes for the plots
            prevalence = np.sum(active_genotypes[:,treatment_target+3]==1)
            print(f"Starting treatment with target {treatment_target}, which is present in {prevalence}/{active_genotypes.shape[0]} of cells")
        if len(traces) > treatment_iteration and len(traces) <= treatment_iteration + treatment_duration:
            treat = True
        if len(traces) > treatment_iteration + treatment_duration and treatment_target != -1:
            print(f"Stopped treating with target {treatment_target}, killed {cells_killed} cells")
            treatment_target = -1

        rng = np.random.default_rng(seed + step)
        cells_killed += tumor.update(treat=treat, treatment_target=treatment_target, rng=rng)
        traces.append(dict(genotypes_counts=deepcopy(tumor.genotypes_counts)))

    # Return tumor
    return tumor, traces, treatment_target, cells_killed

def simulate_fission():
    raise NotImplementedError


def simulate_boundary():
    raise NotImplementedError
