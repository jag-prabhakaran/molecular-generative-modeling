#!/bin/bash
# FILENAME:  generate

#SBATCH -A cis220051-gpu
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1  # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=1     # Number of GPUs per node
#SBATCH --time=03:00:00        # Total run time limit (hh:mm:ss)
#SBATCH -J generate          # Job name
#SBATCH -o jobs/job.o          # Name of stdout output file
#SBATCH -e jobs/job.e   
#SBATCH -p gpu  

# Loads anaconda
module load anaconda
source activate hgraph
python ./finetune_generator.py --train ./data/chembl/short.txt --vocab ./vocab.txt --generative_model ./ckpt/chembl-pretrained/model.ckpt --chemprop_model ./data/chemprop --min_similarity 0.1 --max_similarity 0.5 --nsample 10000 --epoch 10 --threshold 0.0 --save_dir ./ckpt/finetune