#!/bin/bash
# A sample Bash script, by Ryan
module load anaconda
conda activate hgraph
python decode.py --rationale rationale.txt --model ckpt/chembl-h400beta0.3/model.20 --num_decode 10
