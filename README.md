# tumorevo

[![PyPI](https://img.shields.io/pypi/v/tumorevo.svg?style=flat)](https://pypi.python.org/pypi/tumorevo)
[![Tests](https://github.com/pedrofale/tumorevo/actions/workflows/main.yaml/badge.svg)](https://github.com/pedrofale/tumorevo/actions/workflows/main.yaml)

Simulate tumor evolution under different spatial constraints. This package aims to be as awesome as [demon](https://github.com/robjohnnoble/demon_model).
`tumorevo` simulates tumor growth and and produces a Muller plot, a cartoon of the 2D spatial organization of the tumor cells, and a clone tree.

## Installation

```bash
$ pip install tumorevo
```

## Usage

`tumorevo` contains two command line utilities: `tumorsim` and `tumorfig`.

### Simulating tumor evolution
`tumorsim` can be used to simulate the evolution of a tumor according to a specified spatial structure.
```bash
$ tumorsim --mode 1 --steps 2000 --genes 20 --carrying-capacity 5 --grid-size 20 --division-rate 0.2 --dispersal-rate 0.1
100%|████████████████████| 1999/1999 [00:07<00:00, 251.69it/s]
```

This will create a folder containing:
* `parents.csv`: file indicating each clones's parent;
* `trace_counts.csv`: file indicating the number of cells of each clone at each time step;
* `genotypes.csv`: file containing the genotypes of each clone;
* `grid.csv`: file containing the regular grid of genotypes if `mode` > 0.

Full overview:
```
$ tumorsim --help
Usage: tumorsim [OPTIONS]

  Simulate tumor evolution under different spatial constraints.

Options:
  -m, --mode INTEGER              Spatial structure.
  -k, --carrying-capacity INTEGER
                                  Deme carrying capacity.
  -g, --genes INTEGER             Number of genes.
  -s, --steps INTEGER             Number of steps in simulation.
  --grid-size INTEGER             Grid size.
  --division-rate FLOAT           Divison rate.
  --mutation-rate FLOAT           Mutation rate.
  --dispersal-rate FLOAT          Dispersal rate.
  -r, --random_seed INTEGER       Random seed for the pseudo random number
                                  generator.
  --log INTEGER                   Logging level. 0 for no logging, 1 for info,
                                  2 for debug.
  -o, --output-path TEXT          Output directory
  --help                          Show this message and exit.
```

### Plotting tumor evolution
`tumorfig` can be used to create a Muller plot of the tumor's evolution, the 2D spatial organization of the tumor cells, and a clone tree.
```bash
$ tumorfig out/trace_counts.csv out/parents.csv --plot --grid-file out/grid.csv --normalize --remove
```

This will open a figure like this:
<div align="center">
  <img src="https://github.com/pedrofale/tumorevo/raw/main/figures/example.png", width="700px">
</div>

Full overview:
```
$ tumorfig --help
Usage: tumorfig [OPTIONS] GENOTYPE_COUNTS GENOTYPE_PARENTS

  Plot the evolution of a tumor.

Options:
  -c, --cells INTEGER           Number of cells in slice plot.
  -r, --average-radius INTEGER  Average radius of circles in slice plot.
  --grid-file TEXT              Path to grid file.
  --colormap TEXT               Colormap for genotypes.
  --dpi INTEGER                 DPI for figures.
  --plot                        Plot all the figures.
  --do-muller                   Make a Muller plot.
  --do-slice                    Make a slice plot.
  --do-tree                     Make a clone tree plot.
  --normalize                   Normalize the abundances in the Muller plot.
  --labels                      Annotate the clone tree plot.
  --remove                      Remove empty clones in the clone tree plot.
  -o, --output-path TEXT        Directory to write figures into.
  --help                        Show this message and exit.
```
