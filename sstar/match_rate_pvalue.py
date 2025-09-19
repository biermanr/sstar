import pathlib
import logging
from multiprocessing import Pool
import os
import numpy as np

from .simulate import MSPrimeSimulator
from .cal_s_star import _cal_score_worker
from .cal_match_rate import _cal_match_pct_ind, _read_score_file
from .utils import read_data, filter_data

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
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create an MSPrimeSimulator instance and define samples and parameters
    msprime_simulator = MSPrimeSimulator(demes_file=demes_file)
    msprime_simulator.define_sample(role="ref", population=ref_pop, num_samples=ref_size)
    msprime_simulator.define_sample(role="tgt", population=tgt_pop, num_samples=tgt_size)
    msprime_simulator.define_sample(role="src", population=src_pop, num_samples=src_size, time=src_sample_gen)
    msprime_simulator.define_params(mut_rate=mut_rate, recomb_rate=rec_rate, seq_length=seq_length)

    processes = max(1, min(os.cpu_count()-1, threads))

    # Run simulations with a multiprocessing pool
    with Pool(processes=processes) as p:
        sim_vcf_paths = p.starmap(
            _run_sstar_sim_vcf,
            [(msprime_simulator, output_dir, sim_num) for sim_num in range(num_sims)]
        )
        p.close()
        p.join()

    # Run score calculations with a multiprocessing pool
    with Pool(processes=processes) as p:
        score_paths = p.map(_run_sstar_score, sim_vcf_paths)
        p.close()
        p.join()

    # Run match rate calculations with a multiprocessing pool
    with Pool(processes=processes) as p:
        matchrate_paths = p.map(_run_sstar_matchrate, score_paths)
        p.close()
        p.join()

    # Read all matchrate results and calculate p-values for observed match rates
    sim_matchrates = []
    for matchrate_path in matchrate_paths:
        with open(matchrate_path, 'r') as f:
            next(f)  # Skip header
            for line in f:
                fields = line.strip().split('\t')
                _chrom, _start, _end, _sample, match_rate, _src_sample = fields
                sim_matchrates.append(float(match_rate))

    sim_matchrates = np.sort(sim_matchrates)
    num_matchrates = len(sim_matchrates)

    # Read observed match rates, calculate p-values, and write to output file
    out_path = output_dir / "observed_matchrate_pvalues.txt" # TODO change hardcoded name later
    with open(obs_matchrate_path, 'r') as f_in, open(out_path, 'w') as f_out:
        header = f_in.readline().strip()
        f_out.write(header + '\tp_value\n')
        for line in f_in:
            fields = line.strip().split('\t')
            _chrom, _start, _end, _sample, obs_match_rate, _src_sample = fields
            obs_match_rate = float(obs_match_rate)

            # using np.searchsorted for speed
            rank = np.searchsorted(sim_matchrates, obs_match_rate, side='right')
            p_value = 1.0 - (rank / num_matchrates)

            f_out.write(line.strip() + f'\t{p_value}\n')


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


def _run_sstar_score(vcf_path: pathlib.Path) -> pathlib.Path:

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

    chr_names = tgt_data.keys()
    for c in chr_names:
        # Remove variants observed in the reference populations
        # Assume 1 is the alt allele
        variants_not_in_ref = np.sum(ref_data[c]['GT'].is_hom_ref(),axis=1) == len(ref_samples)
        tgt_data = filter_data(tgt_data, c, variants_not_in_ref)

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


def _run_sstar_matchrate(score_path: pathlib.Path) -> pathlib.Path:

    # TODO change hardcoded names later
    vcf_path = score_path.parent / "sim_input.vcf"
    ref_ind_file = score_path.parent / "sim_input.ref.ind.list"
    tgt_ind_file = score_path.parent / "sim_input.tgt.ind.list"
    anc_ind_file = score_path.parent / "sim_input.src.ind.list"
    matchrate_path = score_path.parent / "sim_input.matchrate.results"

    # NOTE similar to scoring, the high-level `cal_match_pct` function does it's own multiprocessing internally
    # NOTE so we can't call it here since we're already in a multiprocessing pool. Instead we'll
    # NOTE call the `_cal_match_pct_ind` function directly here, once per target sample

    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(str(vcf_path), ref_ind_file, tgt_ind_file, anc_ind_file, None)

    chr_names = tgt_data.keys()

    mapped_intervals = None # note, normally uses `utils.read_mapped_region_file`
    data, windows, samples = _read_score_file(score_path, chr_names, tgt_samples)
    sample_size = len(samples)

    header = 'chrom\tstart\tend\tsample\tmatch_rate\tsrc_sample'
    with open(matchrate_path, 'w') as f_out:
        f_out.write(header+"\n")

        for t in samples:

            sample_out_lines = _cal_match_pct_ind(
                data=data[t],
                tgt_ind_index=tgt_samples.index(t),
                mapped_intervals=mapped_intervals,
                tgt_data=tgt_data,
                src_data=src_data,
                src_samples=src_samples,
                sample_size=sample_size,
            )

            for line in sample_out_lines:
                f_out.write(line + "\n")

    return matchrate_path
