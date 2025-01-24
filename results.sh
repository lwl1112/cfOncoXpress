#! /bin/bash -l
#SBATCH --partition=scu-cpu#,panda_physbio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=trainbycontrol01_res
##SBATCH --gres=gpu:a40:1
#SBATCH --time=120:00:00   # HH/MM/SS
#SBATCH --mem=128G   # memory requested, units available: K,M,G,T
#SBATCH --output=trainbycontrol01_res.out


cd
source .bashrc
# correlations for 13 samples with matched rna
python results.py
