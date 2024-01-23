#!/bin/bash
# FILENAME:  training

#SBATCH -A cis220051-gpu
#SBATCH -p gpu 
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=1     # Number of GPUs per node

#SBATCH --time=48:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J training          # Job name
#SBATCH -o jobs/training.o          # Name of stdout output file
#SBATCH -e jobs/training.e   
 

# Loads anaconda
module load anaconda
source activate hgraph
# python ./train_generator.py --train ./train_processed/ --vocab ./vocab.txt  --load_model ./ckpt/chembl-pretrained/model.ckpt.430000 --save_dir ./ckpt/chembl-pretrained
python ./train_translator.py --train ./train_qed/ --vocab ./vocab_translation.txt --save_dir ckpt/translation

echo "done"

