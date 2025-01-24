#! /bin/bash -l
#SBATCH --partition=scu-cpu#,panda_physbio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=trainbycontrol01
##SBATCH --gres=gpu:a40:1
#SBATCH --time=120:00:00   # HH/MM/SS
#SBATCH --mem=128G   # memory requested, units available: K,M,G,T
#SBATCH --output=trainbycontrol01.out


cd
source .bashrc
# train by deep control01 and then test on all Herberts' 63 samples
python python/trainbycontrol01.py
