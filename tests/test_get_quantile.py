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
from sstar.get_quantile import get_quantile

@pytest.mark.integration
def test_get_quantile(test_paths):
    get_quantile(
        model=test_paths.model_file, 
        ms_dir='./ext/msdir', 
        seeds=[1,2,3], 
        N0=1000, 
        nsamp=22, 
        nreps=20000, 
        ref_index=4, 
        ref_size=20, 
        tgt_index=3, 
        tgt_size=2, 
        mut_rate=1.2e-8, 
        rec_rate=0.7e-8, 
        seq_len=40000, 
        snp_num_range=[25,30,5], 
        output_dir=test_paths.quantile_output_dir, 
        thread=2
    )
    
    with open(test_paths.quantile_summary_file, 'r') as f:
        result = f.read()
    
    with open(test_paths.expected_quantile_file, 'r') as f:
        expected_result = f.read()

    assert result == expected_result
