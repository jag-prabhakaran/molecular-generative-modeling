#!/bin/bash
# A sample Bash script, by Ryan
module load anaconda
conda activate hgraph
python decode.py --rationale rationale.txt --model ckpt/gsk3_jnk3_qed_sa/model.final --num_decode 10
