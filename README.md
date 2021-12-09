# tumorevo

[![PyPI](https://img.shields.io/pypi/v/tumorevo.svg?style=flat)](https://pypi.python.org/pypi/tumorevo)
[![Tests](https://github.com/pedrofale/tumorevo/actions/workflows/main.yaml/badge.svg)](https://github.com/pedrofale/tumorevo/actions/workflows/main.yaml)

Simulate tumor evolution under different spatial constraints according to Noble et al (2019).
`tumorevo` produces a cartoon of the 2D spatial organization of the tumor cells, a clone tree and a Muller plot.

## Installation

```bash
$ pip install tumorevo
```

## Usage

`tumorevo` contains two command line utilities: `tumorsim` and `tumorfig`.

`tumorsim` can be used to simulate the evolution of a tumor according to a specified spatial structure.
```bash
$ tumorsim --n_cells 2000 --n_genes 1000 --mode 0
```

This will create a folder containing:
* `parents.csv`: file indicating each clones's parent;
* `trace_counts.csv`: file indicating the number of cells of each clone at each time step;
* `genotypes.csv`: file containing the genotypes of each clone.

Full overview:
```bash
$ tumorsim --help
Usage: tumorsim [OPTIONS]

  Simulate tumor evolution under different spatial constraints.

Options:
  -m, --mode INTEGER              Spatial structure.
  -k, --carrying-capacity INTEGER
                                  Deme carrying capacity.
  -g, --genes INTEGER             Number of genes.
  -s, --steps INTEGER             Number of steps in simulation.
  -d, --division-rate FLOAT       Divison rate.
  -r, --random_seed INTEGER       Random seed for the pseudo random number
                                  generator.
  -o, --output-path TEXT          Output directory
  --help                          Show this message and exit.
```

`tumorfig` can be used to create a Muller plot of the tumor's evolution, the 2D spatial organization of the tumor cells, and a clone tree.
```bash
$ tumorfig trace_counts.csv parents.csv
```

Full overview:
```bash
Usage: tumorfig [OPTIONS] GENOTYPE_COUNTS GENOTYPE_PARENTS

  Plot the evolution of a tumor.

Options:
  -c, --cells INTEGER           Number of cells in slice plot.
  -r, --average-radius INTEGER  Average radius of circles in slice plot.
  -m, --colormap TEXT           Colormap for genotypes.
  --help                        Show this message and exit.
```
