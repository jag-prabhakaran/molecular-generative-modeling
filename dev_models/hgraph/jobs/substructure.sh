#!/bin/bash
# FILENAME:  substructure

#SBATCH --nodes=1
#SBATCH --ntasks=128

#SBATCH --time=1:00:00
#SBATCH -J substructure          # Job name
#SBATCH -o jobs/substructure.o          # Name of stdout output file
#SBATCH -e jobs/substructure.e  


# Loads anaconda
module load anaconda
source activate hgraph
# python get_vocab.py --ncpu 128 < ./data/chembl/all.txt > ./vocab.txt
python ./get_vocab.py --ncpu 128 < ./data/logp06/mols.txt > ./vocab_translation_logp06.txt

echo "done"
