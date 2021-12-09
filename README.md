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

### `tumorsim`
This can be used to simulate the evolution of a tumor according to a specified spatial structure.
```bash
$ tumorsim --n_cells 2000 --n_genes 1000 --mode 0
```

This will create a folder containing:
* `parents.csv`: file indicating each clones's parent;
* `trace_counts.csv`: file indicating the number of cells of each clone at each time step;
* `genotypes.csv`: file containing the genotypes of each clone.

### `tumorfig`
This can be used to create a Muller plot of the tumor's evolution, the 2D spatial organization of the tumor cells, and a clone tree.
```bash
$ tumorfig -f1 trace_counts.csv -f2 parents.csv
```

