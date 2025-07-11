#!/bin/bash
set -eo pipefail

ml anaconda3/2024.10
conda activate sstar

sstar null-matchrate \
    --model examples/models/BonoboGhost_4K19_no_introgression.yaml \
    --ms-dir ext/msdir \
    --N0 1000 \
    --nsamp 22 \
    --nreps 20000 \
    --anc-index 4 \
    --anc-size 20 \
    --tgt-index 3 \
    --tgt-size 2 \
    --mut-rate 1.2e-8 \
    --rec-rate 0.7e-8 \
    --seq-len 40000 \
    --snp-num-range 25 30 5 \
    --output-dir null-matchrates \
    --thread 1

