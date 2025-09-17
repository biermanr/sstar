import demes
import msprime
import pathlib

class MSPrimeSimulator:
    """
    A class to simulate genetic data using msprime based on a given demes file.
    """

    def __init__(self, demes_file: pathlib.Path):
        self.demes_file = demes_file
        self.graph = demes.load(demes_file)
        self.demography = msprime.Demography.from_demes(self.graph)
        self.roles = {}
        self.samples = []
        self.ts = None

    def define_sample(self, role: str, population: str, num_samples: int, time: int = None, ploidy: int = 2) -> None:
        """
        Define the sample sets for the simulation. Follows closely the msprime.SampleSet class.

        Parameters:
            num_samples (int): Number of samples to draw.
            population (str): Population name as defined in the demography.
            time (int | None): Time in generations before present to sample. None means present.
            ploidy (int): Ploidy of the samples (default 2).
        """
        demography_names = [pop.name for pop in self.demography.populations]
        if population not in demography_names:
            raise ValueError(f"Population '{population}' not found in the demography.")

        self.roles[role] = msprime.SampleSet(num_samples, ploidy=ploidy, population=population, time=time)
        self.samples.append(self.roles[role])

    def define_params(self, mut_rate: float = 1.4e-8, recomb_rate: float = 1e-8, seq_length: int = 1_000_000) -> None:
        """
        Define the simulation parameters.

        Parameters:
            mut_rate (float): Mutation rate per base per generation (default 1.4e-8).
            recomb_rate (float): Recombination rate per base per generation (default 1e-8).
            seq_length (int): Length of the sequence to simulate (default 1,000,000).
        """
        self.mut_rate = mut_rate
        self.recomb_rate = recomb_rate
        self.seq_length = seq_length

    def simulate(self, random_seed: int = None) -> None:
        """
        Run the simulation with the defined samples and parameters.

        Parameters:
            random_seed (int | None): Random seed for reproducibility (default None).
        """
        ts = msprime.sim_ancestry(
            recombination_rate=self.recomb_rate,
            sequence_length=self.seq_length,
            samples=self.samples,
            demography=self.demography,
            record_migrations=True,
            random_seed=random_seed,
        )

        self.ts = msprime.sim_mutations(ts, rate=self.mut_rate, random_seed=random_seed)

    def save(self, output_dir: pathlib.Path) -> None:
        """
        Save the simulated tree sequence and VCF to the specified output directory.

        Parameters:
            output_dir (pathlib.Path): Directory to save the output files.
        """
        if self.ts is None:
            raise RuntimeError("Simulation not yet run. Please call simulate() before saving.")

        output_dir.mkdir(parents=True, exist_ok=True)
        self.ts.dump(output_dir / "sim_input.trees")

        with open(output_dir / "sim_input.vcf", 'w') as f_out:
            self.ts.write_vcf(f_out, allow_position_zero=True)

        # Write the role sample lists to {role}.ind.list
        i = 0
        for role, sample_set in self.roles.items():
            with open(output_dir / f"sim_input.{role}.ind.list", 'w') as f_out:
                for _ in range(sample_set.num_samples):
                    f_out.write(f"tsk_{i}\n")
                    i += 1