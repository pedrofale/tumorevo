"""
Simulate tumor growth under different spatial models. Inspired by Noble et al, 2019.
"""
from .cell import TumorCell
from .modes import *

import numpy as np
import pandas as pd

import click
import os
from pathlib import Path
import logging

MODE_LIST = [
    simulate_nonspatial,
    simulate_invasion,
    simulate_fission,
    simulate_boundary,
]


@click.command(help="Simulate tumor evolution under different spatial constraints.")
@click.option("-m", "--mode", default=0, help="Spatial structure.")
@click.option(
    "-k",
    "--carrying-capacity",
    default=100,
    help="Deme carrying capacity.",
)
@click.option(
    "-g",
    "--genes",
    default=100,
    help="Number of genes.",
)
@click.option(
    "-s",
    "--steps",
    default=1000,
    help="Number of steps in simulation.",
)
@click.option("--grid-size", default=50, help="Grid size.")
@click.option("--division-rate", default=0.1, help="Divison rate.")
@click.option("--mutation-rate", default=0.01, help="Mutation rate.")
@click.option("--dispersal-rate", default=0.05, help="Dispersal rate.")
@click.option(
    "-r",
    "--random_seed",
    default=42,
    help="Random seed for the pseudo random number generator.",
)
@click.option(
    "--log", default=0, help="Logging level. 0 for no logging, 1 for info, 2 for debug."
)
@click.option("-o", "--output-path", default="./out", help="Output directory")
def main(
    mode,
    carrying_capacity,
    genes,
    steps,
    grid_size,
    division_rate,
    mutation_rate,
    dispersal_rate,
    random_seed,
    log,
    output_path,
):
    if log == 0:
        log = logging.CRITICAL
    elif log == 1:
        log = logging.INFO
    elif log == 2:
        log == logging.DEBUG
    logging.basicConfig(level=log)

    np.random.seed(random_seed)
    tumor_cell = TumorCell(
        n_genes=genes,
        division_rate=division_rate,
    )
    env, traces = MODE_LIST[mode](
        steps,
        tumor_cell,
        grid_size=grid_size,
        mutation_rate=mutation_rate,
        dispersal_rate=dispersal_rate,
    )
    genotypes, _ = env.get_genotype_frequencies()
    parents = env.genotypes_parents

    print(f"Simulation in mode {mode} finished.")

    Path(output_path).mkdir(parents=True, exist_ok=True)
    pd.DataFrame([t["genotypes_counts"] for t in traces]).fillna(0).to_csv(
        os.path.join(output_path, "trace_counts.csv")
    )
    pd.DataFrame([parents]).to_csv(os.path.join(output_path, "parents.csv"))
    pd.DataFrame(genotypes).to_csv(os.path.join(output_path, "genotypes.csv"))

    if mode > 0:
        genotype_matrix = env.get_genotype_matrix()
        pd.DataFrame(genotype_matrix).to_csv(os.path.join(output_path, "grid.csv"))

    print(f"Saved results to {output_path}.")


if __name__ == "__main__":
    main()
