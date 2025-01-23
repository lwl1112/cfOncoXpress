#! /bin/bash -l
 
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=wps_cfDNA
#SBATCH --time=96:00:00   # HH/MM/SS
#SBATCH --mem=128G   # memory requested, units available: K,M,G,T
#SBATCH --output=exp/slurm_out/wps/slurm-%x-%j.out


source ~/.bashrc

set -euxo pipefail
export TMPDIR=/scratch
echo $TMPDIR

if ! [ -z ${SLURM_SUBMIT_DIR+x} ]; then cd ${SLURM_SUBMIT_DIR}; fi

date 

wpsk=120
wpstau=1000
wpslambda=21
wpsalpha=51
relevantreadsize="120 180"
jobregionwindow=1000000

for i in $(seq 1 ${SLURM_CPUS_PER_TASK}); do
    export LOG_SUBJOB=$i
    python -u ../pythoncodes/cfDNA_wps.py --wps-k ${wpsk} --relevant-read-size ${relevantreadsize} --job-region-window ${jobregionwindow} --wps-tau ${wpstau} --wps-lambda ${wpslambda} --wps-alpha ${wpsalpha} &
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
