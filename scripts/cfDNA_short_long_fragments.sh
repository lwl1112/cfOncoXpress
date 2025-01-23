#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=cfDNA_sl
#SBATCH --time=96:00:00   # HH/MM/SS
#SBATCH --mem=150G   # memory requested, units available: K,M,G,T
#SBATCH --output=exp/slurm_out/frags/slurm-%x-%j.out

source ~/.bashrc

set -euxo pipefail
export TMPDIR=/scratch
echo $TMPDIR

cd $SLURM_SUBMIT_DIR

date 
hostname

for i in $(seq 1 ${SLURM_CPUS_PER_TASK}); do
    export LOG_SUBJOB=$i
    python -u ../pythoncodes/cfDNA_short_long_fragments.py &
done

exitcode=0
for job in $(jobs -p); do
    echo "waiting for job $job"
    if wait $job; then
        echo "job $job succeeded at $(date)"
    else
        exitcode=$?
        echo "job $job failed with exit code $exitcode at $(date)"
    fi
done

exit $exitcode
