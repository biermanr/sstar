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
import pandas as pd
import numpy as np
from sstar.cal_threshold import cal_threshold, _read_recomb_map

# WRITE SOME UNIT TESTS FOR cal_threshold!!!
# MAYBE START WITH _read_recomb_map
# AFTER I GET GOOD COVERAGE THEN ADD A FUNCTION FOR pyGAM AND TEST THAT AS A REPLACEMENT


def test_read_recomb_map(tmp_path):

    recomb_map = tmp_path / "recomb_map.txt"
    with open(recomb_map, 'w') as f:
        f.write("\t".join(["chr1",   "0", "1000", "0.5"]) + "\n")
        f.write("\t".join(["chr1", "500", "1500", "0.75"]) + "\n")
        f.write("\t".join(["chr1", "500", "1500", "0.85"]) + "\n")

    recomb_dict = _read_recomb_map(recomb_map) 

    assert len(recomb_dict) == 2
    assert recomb_dict["chr1:0-1000"] == 0.5
    assert recomb_dict["chr1:500-1500"] == 0.75

def test_read_recomb_map_incorrect_format(tmp_path):

    recomb_map = tmp_path / "recomb_map.txt"
    with open(recomb_map, 'w') as f:
        f.write("\t".join(["chr1",   "0", "1000"]) + "\n")
        f.write("\t".join(["chr1", "500", "1500", "0.75"]) + "\n")

    with pytest.warns(RuntimeWarning, match="Incorrect format"):
        _read_recomb_map(recomb_map)

def test_read_recomb_map_empty(tmp_path):

    recomb_map = tmp_path / "recomb_map.txt"
    with open(recomb_map, 'w') as f:
        pass

    with pytest.warns():
        recomb_dict = _read_recomb_map(recomb_map) 

    assert len(recomb_dict) == 0

@pytest.mark.integration
def test_cal_threshold(test_paths):
    cal_threshold(
        simulated_data=test_paths.simulated_data_file, 
        score_file=test_paths.expected_score_file, 
        recomb_rate=0, 
        recomb_map=test_paths.recomb_map_file, 
        quantile=0.99, 
        output=test_paths.output_threshold_file, 
        k=8
    )

    df1 = pd.read_csv(test_paths.output_threshold_file, sep="\t")
    df2 = pd.read_csv(test_paths.expected_threshold_file, sep="\t")

    assert df1.shape == df2.shape, "DataFrame shape mismatch"

    for col in df1.columns:
        assert col in df2.columns, f"Column '{col}' missing in expected output"

        if pd.api.types.is_float_dtype(df1[col]):
            assert np.allclose(
                df1[col], df2[col], rtol=1e-5, atol=1e-8, equal_nan=True
            ), f"Float column '{col}' differs"
        else:
            assert (df1[col].fillna("").astype(str).values == df2[col].fillna("").astype(str).values).all(), f"Column '{col}' differs"
