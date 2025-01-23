#! /bin/bash -l
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --job-name=DTB-097-Progression-cfDNA_simpson
##SBATCH --gres=gpu:a100:1 #a40
#SBATCH --time=96:00:00   # HH/MM/SS
#SBATCH --mem=128G   # memory requested, units available: K,M,G,T
#SBATCH --output=DTB-097-Progression-cfDNA_simpson.out

export TMPDIR=/scratch
echo $TMPDIR
date
cd
source .bashrc
python --version
python ../pythoncodes/cfDNA_compute_simpson_entropy_in_bins.py --sample-id DTB-097-Progression-cfDNA --bed-path ../csv_files/finalentropy1k.bed --out-dir csv_files/ --min-fragment-length 100 --max-fragment-length 300 --quality-threshold 20
date
