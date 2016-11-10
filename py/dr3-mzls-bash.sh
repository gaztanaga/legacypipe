#!/bin/bash -l

#SBATCH -p shared
#SBATCH -n 12
#SBATCH -t 01:00:00
#SBATCH --account=desi
#SBATCH -J dr3-mzls-bash
#SBATCH -o dr3-mzls-bash.o%j
#SBATCH --mail-user=kburleigh@lbl.gov
#SBATCH --mail-type=END,FAIL
#SBATCH -L SCRATCH

source /scratch1/scratchdirs/desiproc/dr3-mzls/bashrc
set -x

#outdir,statdir,brick
bri="$(echo $brick | head -c 3)"

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

threads=6
export OMP_NUM_THREADS=$threads

echo outdir="$outdir", brick="$brick"
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
rm $statdir/inq_$brick.txt

#     --radec $ra $dec
#    --force-all --no-write \
#    --skip-calibs \
#
echo dr3-mzls-bash DONE $SLURM_JOBID


