import numpy as np

import pytest
from click.testing import CliRunner
import os

from tumorevo.tumorsim import main
from tumorevo.tumorsim.modes import *
from tumorevo.tumorsim.cell import TumorCell

MODE_LIST = [
    simulate_nonspatial,
    simulate_invasion,
    # simulate_fission,
    # simulate_boundary,
]


@pytest.mark.parametrize("mode", list(range(len(MODE_LIST))))
def test_simulation(mode):
    tumor_cell = TumorCell(n_genes=100)
    env, traces = MODE_LIST[mode](1000, tumor_cell)

    assert len(traces) == 1000
    assert isinstance(env.get_genotype_frequencies(), tuple)


@pytest.mark.parametrize("mode", list(range(len(MODE_LIST))))
def test_cli(mode):
    runner = CliRunner()

    # run program
    with runner.isolated_filesystem():
        result = runner.invoke(main, ["--mode", mode, "-o", "./outs"])

        # test output
        assert result.exit_code == 0
        assert os.path.isfile("./outs/parents.csv")
        assert os.path.isfile("./outs/trace_counts.csv")
        assert os.path.isfile("./outs/genotypes.csv")
