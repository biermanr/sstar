# Apache License Version 2.0
# Copyright 2022 Xin Huang
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest
from dataclasses import dataclass


@dataclass
class TestPaths:
    """Centralized test file paths and utilities."""
    
    # Base directories
    test_data_dir: str = "./tests/data"
    test_results_dir: str = "./tests/results"
    examples_dir: str = "./examples"
    
    # Individual list files
    ref_ind_file: str = "./examples/data/ind_list/ref.ind.list"
    tgt_ind_file: str = "./examples/data/ind_list/tgt.ind.list"
    src_ind_file: str = "./examples/data/ind_list/nean.ind.list"
    test_ref_ind_file: str = "./tests/data/test.ref.ind.list"
    test_tgt_ind_file: str = "./tests/data/test.tgt.ind.list"
    empty_ind_file: str = "./tests/data/test.empty.ind.list"
    
    # VCF files
    test_vcf_file: str = "./tests/data/test.score.data.vcf"
    match_rate_vcf_file: str = "./tests/data/test.match.rate.data.vcf"
    
    # Ancestral allele files
    anc_allele_file: str = "./tests/data/test.anc.allele.bed"
    empty_anc_allele_file: str = "./tests/data/test.empty.anc.allele.bed"
    
    # Model and configuration files
    model_file: str = "./examples/models/BonoboGhost_4K19_no_introgression.yaml"
    simulated_data_file: str = "./examples/data/simulated_data/gravel_asn_scale_60k.simulated.data"
    recomb_map_file: str = "./examples/data/real_data/hum.windows.50k.10k.recomb.map"
    mapped_regions_file: str = "./tests/data/test.mapped.region.bed"
    
    # Tract analysis files
    threshold_file: str = "./tests/data/test.tract.threshold"
    src1_match_pct_file: str = "./tests/data/test.tract.src1.match.rate"
    src2_match_pct_file: str = "./tests/data/test.tract.src2.match.rate"
    
    # Expected result files
    expected_quantile_file: str = "./tests/results/test.quantile.exp.summary"
    expected_match_rate_file: str = "./tests/results/test.match.rate.exp.results"
    expected_score_file: str = "./tests/results/test.score.exp.results"
    score_file: str = "./tests/results/test.match.rate.score.exp.results"
    expected_threshold_file: str = "./tests/results/test.threshold.exp.results"
    expected_tract_bed_file: str = "./tests/results/test.tract.exp.bed"
    expected_tract_bed_with_src_file: str = "./tests/results/test.tract.with.src.match.rate.exp.bed"
    expected_tract_src1_bed_file: str = "./tests/results/test.tract.exp.src1.bed"
    expected_tract_src2_bed_file: str = "./tests/results/test.tract.exp.src2.bed"
    
    # Output files
    output_match_rate_file: str = "./tests/results/test.match.rate.results"
    output_score_file: str = "./tests/results/test.score.results"
    output_threshold_file: str = "./tests/results/test.threshold.results"
    quantile_output_dir: str = "./tests/results/simulation"
    quantile_summary_file: str = "./tests/results/simulation/quantile.summary.txt"
    tract_output_prefix: str = "./tests/results/test.tract"
    tract_with_src_output_prefix: str = "./tests/results/test.tract.with.src.match.rate"


@pytest.fixture
def test_paths():
    """Fixture providing access to all test paths and utilities."""
    return TestPaths()
