# tumorevo

[![PyPI](https://img.shields.io/pypi/v/tumorevo.svg?style=flat)](https://pypi.python.org/pypi/tumorevo)
[![Tests](https://github.com/pedrofale/tumorevo/actions/workflows/main.yaml/badge.svg)](https://github.com/pedrofale/tumorevo/actions/workflows/main.yaml)

Simulate tumor evolution under different spatial constraints according to Noble et al (2019).
`tumorevo` produces a cartoon of the 2D spatial organization of the tumor cells, a phylogenetic tree and a Muller plot, using [pymuller](https://github.com/boaz85/pymuller).

## Installation

```bash
$ pip install tumorevo
```

## Usage

`tumorevo` can be used from the command line:

```bash
$ tumorevo --n_cells 2000 --n_genes 1000 --mode 0
```

This will generate two figures: a picture and a phylogeny.
