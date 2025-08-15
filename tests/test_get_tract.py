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
from sstar.get_tract import get_tract


def test_get_tract(test_paths):
    # Test 1: Basic tract extraction without match percentage files
    get_tract(threshold_file=test_paths.threshold_file, match_pct_files=None, output_prefix=test_paths.tract_output_prefix, diff=0)
    
    with open(f'{test_paths.tract_output_prefix}.bed', 'r') as f, open(test_paths.expected_tract_bed_file, 'r') as ef:
        assert f.read() == ef.read()

    # Test 2: Tract extraction with one source match percentage file
    get_tract(threshold_file=test_paths.threshold_file, match_pct_files=[test_paths.src1_match_pct_file], 
              output_prefix=test_paths.tract_with_src_output_prefix, diff=0)

    with open(f'{test_paths.tract_with_src_output_prefix}.bed', 'r') as f, open(test_paths.expected_tract_bed_with_src_file, 'r') as ef:
        assert f.read() == ef.read()

    # Test 3: Tract extraction with two source match percentage files
    get_tract(threshold_file=test_paths.threshold_file, match_pct_files=[test_paths.src1_match_pct_file, test_paths.src2_match_pct_file], 
              output_prefix=test_paths.tract_output_prefix, diff=0)
    
    with open(f'{test_paths.tract_output_prefix}.src1.bed', 'r') as f1, open(test_paths.expected_tract_src1_bed_file, 'r') as ef1:
        assert f1.read() == ef1.read()

    with open(f'{test_paths.tract_output_prefix}.src2.bed', 'r') as f2, open(test_paths.expected_tract_src2_bed_file, 'r') as ef2:
        assert f2.read() == ef2.read()
