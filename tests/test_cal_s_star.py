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

from sstar.cal_s_star import cal_s_star

def test_cal_s_star(test_paths, tmp_path):
    # Redirect output to a temporary file to avoid overwriting shared test data
    output_path = tmp_path / "cal_s_star_output.txt"
    cal_s_star(
        vcf=test_paths.test_vcf_file, 
        ref_ind_file=test_paths.test_ref_ind_file, 
        tgt_ind_file=test_paths.test_tgt_ind_file, 
        anc_allele_file=None, 
        output=str(output_path), 
        win_len=50000, 
        win_step=10000, 
        thread=1, 
        match_bonus=5000, 
        max_mismatch=5, 
        mismatch_penalty=-10000
    )
    
    with open(output_path, 'r') as f:
        result = f.read()
    
    with open(test_paths.expected_score_file, 'r') as f:
        expected_result = f.read()

    assert result == expected_result


def test_cal_s_star_multithread(test_paths, tmp_path):
    """Test that multithread execution produces the same results as single-thread execution."""
    single_output_path = tmp_path / "cal_s_star_single_thread.txt"
    multi_output_path = tmp_path / "cal_s_star_multi_thread.txt"

    # Run with single thread
    cal_s_star(
        vcf=test_paths.test_vcf_file, 
        ref_ind_file=test_paths.test_ref_ind_file, 
        tgt_ind_file=test_paths.test_tgt_ind_file, 
        anc_allele_file=None,
        output=single_output_path, 
        win_len=50000, 
        win_step=10000, 
        thread=1, 
        match_bonus=5000, 
        max_mismatch=5, 
        mismatch_penalty=-10000
    )
        
    # Run with multiple threads (use 2 threads for testing)
    cal_s_star(
        vcf=test_paths.test_vcf_file, 
        ref_ind_file=test_paths.test_ref_ind_file, 
        tgt_ind_file=test_paths.test_tgt_ind_file, 
        anc_allele_file=None, 
        output=multi_output_path, 
        win_len=50000, 
        win_step=10000, 
        thread=2, 
        match_bonus=5000, 
        max_mismatch=5, 
        mismatch_penalty=-10000
    )
        
    # Read results from both runs
    with open(single_output_path, 'r') as f:
        single_result = f.read()
    
    with open(multi_output_path, 'r') as f:
        multi_result = f.read()
    
    # Results should be identical
    assert single_result == multi_result, "Multithread results should match single thread results"