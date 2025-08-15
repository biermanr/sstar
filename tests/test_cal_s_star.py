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

import allel
import pytest
import numpy as np
from sstar.cal_s_star import cal_s_star


def test_cal_s_star(test_paths):
    cal_s_star(
        vcf=test_paths.test_vcf_file, 
        ref_ind_file=test_paths.test_ref_ind_file, 
        tgt_ind_file=test_paths.test_tgt_ind_file, 
        anc_allele_file=None, 
        output=test_paths.output_score_file, 
        win_len=50000, 
        win_step=10000, 
        thread=1, 
        match_bonus=5000, 
        max_mismatch=5, 
        mismatch_penalty=-10000
    )
    
    with open(test_paths.output_score_file, 'r') as f:
        result = f.read()
    
    with open(test_paths.expected_score_file, 'r') as f:
        expected_result = f.read()

    assert result == expected_result
