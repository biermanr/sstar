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
from sstar.cal_threshold import cal_threshold

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
