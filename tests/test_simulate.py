import pytest
from sstar import simulate

@pytest.fixture
def simulator(test_paths):
    return simulate.MSPrimeSimulator(
        demes_file = test_paths.model_file,
    )

def test_define_sample(simulator):

    assert len(simulator.samples) == 0

    simulator.define_sample(
        role = "tgt",
        population = "Bonobo",
        num_samples = 2,
        ploidy = 2,
    )
    assert len(simulator.samples) == 1


def test_define_missing_sample(simulator):

    assert len(simulator.samples) == 0

    with pytest.raises(ValueError, match="not found"):
        simulator.define_sample(
            role = "tgt",
            population = "NonexistentPopulation",
            num_samples = 2,
            ploidy = 2,
        )

def test_define_params(simulator):

    simulator.define_params(
        mut_rate = 1.4e-8,
        recomb_rate = 1e-8,
        seq_length = 1_000_000,
    )

    assert simulator.mut_rate == 1.4e-8
    assert simulator.recomb_rate == 1e-8
    assert simulator.seq_length == 1_000_000


def test_simulate(simulator):

    simulator.define_sample(
        role = "tgt",
        population = "Bonobo",
        num_samples = 2,
        ploidy = 2,
    )

    simulator.define_params(
        mut_rate = 1e-2,
        recomb_rate = 1e-8,
        seq_length = 1000,
    )

    simulator.simulate(random_seed = 42)

    assert simulator.ts is not None
    assert simulator.ts.num_samples == 4  # 2 diploid individuals

def test_save(tmp_path, simulator):

    simulator.define_sample(
        role = "tgt",
        population = "Bonobo",
        num_samples = 2,
        ploidy = 2,
    )

    simulator.define_params(
        mut_rate = 1e-2,
        recomb_rate = 1e-8,
        seq_length = 1000,
    )

    simulator.simulate(random_seed = 42)

    output_dir = tmp_path / "simulation_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    simulator.save(output_dir)

    # TODO change to define output file names as argument
    # TODO this is very fragile
    assert (output_dir / "sim_input.trees").exists()
    assert (output_dir / "sim_input.vcf").exists()
    assert (output_dir / "sim_input.tgt.ind.list").exists()