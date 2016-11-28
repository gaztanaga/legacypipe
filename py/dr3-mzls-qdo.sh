#!/bin/bash 

set -x
brick="$1"
bri="$(echo $brick | head -c 3)"

source /scratch1/scratchdirs/desiproc/dr3-mzls/bashrc
export PYTHONPATH=/scratch1/scratchdirs/desiproc/code/qdo:$PYTHONPATH
export PATH=/scratch1/scratchdirs/desiproc/code/qdo/bin:$PATH

outdir=/scratch1/scratchdirs/desiproc/data-releases/dr3-mzls-qdo
rundir=$SCRATCH/code/legacypipe/py

log="$outdir/logs/$bri/$brick/log.$SLURM_JOBID"
mkdir -p $(dirname $log)

echo Logging to: $log
echo Running on ${NERSC_HOST} $(hostname)

echo -e "\n\n\n" >> $log
echo "-----------------------------------------------------------------------------------------" >> $log
echo "PWD: $(pwd)" >> $log
echo "Modules:" >> $log
module list >> $log 2>&1
echo >> $log
echo "Environment:" >> $log
set | grep -v PASS >> $log
echo >> $log
ulimit -a >> $log
echo >> $log

echo -e "\nStarting on ${NERSC_HOST} $(hostname)\n" >> $log
echo "-----------------------------------------------------------------------------------------" >> $log

threads=8
export OMP_NUM_THREADS=$threads

echo outdir="$outdir", brick="$brick"
cd $rundir
srun -n 1 -c $OMP_NUM_THREADS python legacypipe/runbrick.py \
     --brick $brick \
     --skip \
     --threads $OMP_NUM_THREADS \
     --checkpoint $outdir/checkpoints/${bri}/${brick}.pickle \
     --pickle "$outdir/pickles/${bri}/runbrick-%(brick)s-%%(stage)s.pickle" \
     --outdir $outdir --nsigma 6 \
     --no-wise \
     >> $log 2>&1
#--zoom 1400 1600 1400 1600

#     --radec $ra $dec
#    --force-all --no-write \
#    --skip-calibs \
#
echo dr3-mzls-qdo DONE $SLURM_JOBID


