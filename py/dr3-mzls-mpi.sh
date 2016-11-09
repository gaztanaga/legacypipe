#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 100
#SBATCH -t 00:30:00
#SBATCH --account=desi
#SBATCH -J dr3-mzls-mpi
#SBATCH -o dr3-mzls-mpi.o%j
#SBATCH --mail-user=kburleigh@lbl.gov
#SBATCH --mail-type=END,FAIL
#SBATCH -L SCRATCH

source /scratch1/scratchdirs/desiproc/dr3-mzls/bashrc
#bcast
source /scratch1/scratchdirs/desiproc/code/yu-bcast/activate.sh
outdir=/scratch1/scratchdirs/desiproc/data-releases/dr3-mzls-mpi
bricklist=/scratch1/scratchdirs/desiproc/dr3-mzls/bricks-dr3-mzls2000.txt

if [ "$NERSC_HOST" == "cori" ]; then
    cores=32
elif [ "$NERSC_HOST" == "edison" ]; then
    cores=24
fi
threads=12
let tasks=${SLURM_JOB_NUM_NODES}*${cores}/${threads}

export OMP_NUM_THREADS=$threads
srun -n $tasks -N ${SLURM_JOB_NUM_NODES} -c $OMP_NUM_THREADS python runbrick_mpi4py.py \
     --brick_list $bricklist  --outdir $outdir \
     --jobid $SLURM_JOBID \
     --threads $OMP_NUM_THREADS 
#--zoom 1400 1600 1400 1600



#brick=$(head -n 1 kay_bricks2.txt)
#bri=$(echo $brick | head -c 3)
#srun -n 1 -c $OMP_NUM_THREADS python legacypipe/runbrick_mpi4py.py \
#     --skip \
#     --threads $OMP_NUM_THREADS \
#     --checkpoint checkpoints/${bri}/checkpoint-${brick}.pickle \
#     --pickle "pickles/${bri}/runbrick-%(brick)s-%%(stage)s.pickle" \
#     --brick $brick --outdir $outdir --nsigma 6 \
#     --no-wise --zoom 1400 1800 1400 1800 --force-all 
#     >> $log 2>&1

#    --force-all --no-write \
#    --skip-calibs \
#    --pipe \

echo DONE




