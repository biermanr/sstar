import pathlib
import logging
from multiprocessing import Pool
import os

from .simulate import MSPrimeSimulator
from .cal_s_star import _cal_score_worker
from .utils import read_data

def archaic_matchrate_pvalue(*,
                             demes_file: pathlib.Path,
                             obs_matchrate_path: pathlib.Path, 
                             output_dir: pathlib.Path,
                             num_sims: int,
                             ref_pop: str,
                             tgt_pop: str,
                             src_pop: str,
                             ref_size: int,
                             tgt_size: int,
                             src_size: int,
                             src_sample_gen: int,
                             mut_rate: float,
                             rec_rate: float,
                             seq_length: int,
                             threads: int,
                             ):
    """archaic_matchrate_pvalue: calculate p-values for archaic match rates

    Performs the following steps once per simulation:
    1) Simulate genetic data without introgression using msprime based on a given demes file
    2) Calculate scores
    3) Calculate match rates

    Then combine the results from all simulations to calculate p-values for observed match rates.
    """
    # NOTE SETUP LOGGING SOMEWHERE ELSE
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO,
    )

    logging.info("Calculating null match rates from simulated data without introgression...")
    logging.info(f"Using observed match rate file: {obs_matchrate_path}")
    logging.info(f"Using demographic file: {demes_file}")
    logging.info(f"Using number of simulations: {num_sims}")
    logging.info(f"Reference population: {ref_pop}, size: {ref_size}")
    logging.info(f"Target population: {tgt_pop}, size: {tgt_size}")
    logging.info(f"Source population: {src_pop}, size: {src_size}, sampled {src_sample_gen} generations ago")
    logging.info(f"Mutation rate: {mut_rate}, recombination rate: {rec_rate}")
    logging.info(f"Output directory: {output_dir}, threads: {threads}")

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create an MSPrimeSimulator instance and define samples and parameters
    msprime_simulator = MSPrimeSimulator(demes_file=demes_file)
    msprime_simulator.define_sample(role="ref", population=ref_pop, num_samples=ref_size)
    msprime_simulator.define_sample(role="tgt", population=tgt_pop, num_samples=tgt_size)
    msprime_simulator.define_sample(role="nean_src", population=src_pop, num_samples=src_size, time=src_sample_gen)
    msprime_simulator.define_params(mut_rate=mut_rate, recomb_rate=rec_rate, seq_length=seq_length)

    processes = max(1, min(os.cpu_count()-1, threads))

    # Run simulations with a multiprocessing pool
    with Pool(processes=processes) as p:
        sim_vcf_paths = p.starmap(
            _run_sstar_sim_vcf,
            [(msprime_simulator, output_dir, sim_num) for sim_num in range(num_sims)]
        )

    # TODO Run score calculations with a multiprocessing pool
    with Pool(processes=processes) as p:
        score_paths = p.map(_run_sstar_score, sim_vcf_paths)

    # TODO Run match rate calculations with a multiprocessing pool

    # TODO Read all matchrate results and calculate p-values for observed match rates


def _run_sstar_sim_vcf(simulator: MSPrimeSimulator, out_simulation_dir: pathlib.Path, sim_num: int) -> pathlib.Path:
    """Run msprime simulation and save VCF and sample lists to the specified directory.

    Parameters:
        simulator (MSPrimeSimulator): An instance of the MSPrimeSimulator class with defined samples and parameters.
        out_simulation_dir (pathlib.Path): Directory to save the output files.
    """
    out_simulation_dir.mkdir(parents=True, exist_ok=True)
    out_simulation_dir = out_simulation_dir / f"sim{sim_num}"
    out_simulation_dir.mkdir(parents=True, exist_ok=True)

    simulator.simulate()
    simulator.save(out_simulation_dir)

    return out_simulation_dir / "sim_input.vcf" #TODO change hardcoded name later


def _run_sstar_score(vcf_path: pathlib.Path):

    # TODO change hardcoded names later
    ref_ind_file = vcf_path.parent / "sim_input.ref.ind.list"
    tgt_ind_file = vcf_path.parent / "sim_input.tgt.ind.list"
    out_file = vcf_path.parent / "sim_input.score.results" # TODO change hardcoded name later

    # NOTE the high-level `cal_s_star` function does it's own multiprocessing internally
    # NOTE so we can't call it here since we're already in a multiprocessing pool.
    # NOTE we get error `AssertionError: daemonic processes are not allowed to have children`
    # NOTE so instead we call the worker function directly here, once per target sample
    # NOTE and collect the score results

    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(str(vcf_path), ref_ind_file, tgt_ind_file, None, None)
    win_len = 50_000            # TODO set as parameter
    win_step = 10_000           # TODO set as parameter
    match_bonus = 5_000         # TODO set as parameter
    max_mismatch = 5            # TODO set as parameter
    mismatch_penalty = -10_000  # TODO set as parameter

    # TODO this is largely duplicate code to whats in cal_s_star, refactor into a shared function
    with open(out_file, 'w') as out:
        header = 'chrom\tstart\tend\tsample\tS*_score\tregion_ind_SNP_number\tS*_SNP_number\tS*_SNPs'
        out.write(header+'\n')

        for s, sample_name in enumerate(tgt_samples):
            score_lines = _cal_score_worker(s, sample_name, ref_data, tgt_data, win_len, win_step, match_bonus, max_mismatch, mismatch_penalty)
            for score in score_lines:
                out.write(score + '\n')

    return out_file