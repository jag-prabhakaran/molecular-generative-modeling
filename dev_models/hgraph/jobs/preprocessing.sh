#!/bin/bash
# FILENAME:  preprocessing

#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=00:45:00        # Total run time limit (hh:mm:ss)

#SBATCH -J preprocessing          # Job name
#SBATCH -o jobs/preprocessing.o          # Name of stdout output file
#SBATCH -e jobs/preprocessing.e 

# Loads anaconda
module load anaconda
source activate hgraph
# python ./preprocess.py --train ./data/chembl/small.txt --vocab ./vocab.txt --ncpu 128 --mode single
python ./preprocess.py --train ./data/logp06/train_pairs.txt --vocab ./vocab_translation_logp06.txt --ncpu 128

echo "done"