from sstar import match_rate_pvalue
import pathlib
import hashlib
import collections

def md5(fname: pathlib.Path):
    """Calculate the MD5 checksum of a file."""
    with open(fname, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()

@pytest.mark.integration
def test_run_archaic_match_rate_pvalue(test_paths, tmp_path):
    output_dir = tmp_path
    num_sims = 2


    # NOTE not sure I set the correct ref/tgt/src populations here
    match_rate_pvalue.archaic_matchrate_pvalue(
        demes_file=test_paths.model_file,
        obs_matchrate_path=test_paths.expected_match_rate_file, #NOTE this is for Nean model. Need to make new test data
        output_dir=output_dir,
        num_sims=num_sims,
        ref_pop="Central",
        ref_size=2,
        tgt_pop="Bonobo",
        tgt_size=3,
        src_pop="Ghost",
        src_size=2,
        src_sample_gen=1000,
        mut_rate=1.0e-3,
        rec_rate=1.0e-8,
        seq_length=1000,
        threads=2,
    )

    # TODO currently hardcoded output path
    assert (output_dir / "observed_matchrate_pvalues.txt").exists()

    # Check that output files are created
    md5sums = collections.defaultdict(set)
    for i in range(num_sims):
        sim_dir = output_dir / f"sim{i}"
        assert sim_dir.exists()

        assert (sim_dir / "sim_input.vcf").exists()
        md5sums["vcfs"].add(md5(sim_dir / "sim_input.vcf"))
        
        assert (sim_dir / "sim_input.score.results").exists()
        md5sums["scores"].add(md5(sim_dir / "sim_input.score.results"))

        assert (sim_dir / "sim_input.matchrate.results").exists()
        md5sums["matchrates"].add(md5(sim_dir / "sim_input.matchrate.results"))

    # Check that files are unique across simulations for each type
    for ftype, sums in md5sums.items():
        assert len(sums) == num_sims, f"Not all {ftype} files are unique across simulations"