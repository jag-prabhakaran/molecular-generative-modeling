#!/bin/bash
# FILENAME:  translation

#SBATCH -A cis220051-gpu
#SBATCH -p gpu 
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=1     # Number of GPUs per node

#SBATCH --time=1:00:00  
#SBATCH -J translation          # Job name
#SBATCH -o jobs/translation.o          # Name of stdout output file
#SBATCH -e jobs/translation.e 

# Loads anaconda
module load anaconda
source activate hgraph
# python ./preprocess.py --train ./data/chembl/small.txt --vocab ./vocab.txt --ncpu 128 --mode single
#python ./preprocess.py --train ./data/qed/train_pairs.txt --vocab ./vocab_translation.txt --ncpu 128
# python ./translate.py --test ./data/qed/test.txt --vocab ./vocab_translation_logp06.txt --model ckpt/translation/model.6 --num_decode 5 > results_logp06.csv
python ./translate.py --test ./test.txt --vocab ./vocab_translation_logp06.txt --model ckpt/translation/model.11 --num_decode 1 > ./test_result.txt

