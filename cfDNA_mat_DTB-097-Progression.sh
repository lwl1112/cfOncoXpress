#! /bin/bash -l
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --job-name=DTB-097-Progression-cfDNA_bigmat
##SBATCH --gres=gpu:a100:1 #a40
#SBATCH --time=96:00:00   # HH/MM/SS
#SBATCH --mem=96G   # memory requested, units available: K,M,G,T
#SBATCH --output=DTB-097-Progression-cfDNA_bigmat.out


date
cd
source .bashrc
python --version
python python/finalbigmat_all.py DTB-097-Progression-cfDNA
date
