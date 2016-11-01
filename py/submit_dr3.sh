#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH -J dr3
#SBATCH -o dr3.o%j
#SBATCH --mail-user=kburleigh@lbl.gov
#SBATCH --mail-type=END,FAIL
#SBATCH -L SCRATCH

set -x
echo ########
echo set full_hpcp=special AND tractor_prod=yes in bashrc.ext
echo ########

#bcast
source /scratch2/scratchdirs/kaylanb/yu-bcase/activate.sh
outdir=testing_bcast_mpi
#outdir=$SCRATCH/dr3/legacypipe/py/runs/decam_dr3/bcast${did_bcast}

# From Aaron's cpy to edison scratch a few weeks ago
export UNWISE_COADDS_DIR=/scratch1/scratchdirs/desiproc/unwise-coadds/fulldepth:/scratch1/scratchdirs/desiproc/unwise-coadds
export UNWISE_COADDS_TIMERESOLVED_DIR=/scratch1/scratchdirs/desiproc/unwise-coadds/time_resolved_neo1
# from Dustin's email b4 Tucson workshopw
#export UNWISE_COADDS_DIR=/scratch1/scratchdirs/ameisner/unwise-coadds/fulldepth_sv:/scratch1/scratchdirs/desiproc/unwise-coadds
#export UNWISE_COADDS_TIMERESOLVED_DIR=/scratch1/scratchdirs/ameisner/unwise-coadds/time_resolved_dr3
#/scratch1/scratchdirs/desiproc/unwise-coadds/w3w4

#export LEGACY_SURVEY_DIR=/global/cscratch1/sd/desiproc/dr3
export LEGACY_SURVEY_DIR=/scratch2/scratchdirs/kaylanb/dr3/desiproc-dr3-template
export DUST_DIR=/project/projectdirs/cosmo/work/decam/modules/all/dust/v0_0
# above is equivalent to: module load dust/v0_0
#export DUST_DIR=/global/cscratch1/sd/desiproc/dust/v0_0
#export LEGACY_SURVEY_DIR=/scratch1/scratchdirs/desiproc/dr3/
#export DUST_DIR=/scratch1/scratchdirs/desiproc/dust/v0_0

##########
# https://github.com/legacysurvey/legacypipe/blob/master/bin/pipebrick-checkpoint.sh
#export PYTHONPATH=${PYTHONPATH}:.

# Force MKL single-threaded
# https://software.intel.com/en-us/articles/using-threaded-intel-mkl-in-multi-thread-application
# takes affect b/c bashrc has: module load intel 
export MKL_NUM_THREADS=1

# Try limiting memory to avoid killing the whole MPI job...
#ulimit -S -v 15000000
#ulimit -S -v 30000000
ulimit -a

# Make sure we're reading from Edison scratch
#module unload dust
#module load dust/scratch
### argh modules not seeming to work.

# Point to Legacypipe
# git checkout 6fad8727dc78  --> the latest version before dr4 refactoring

module unload tractor-hpcp
#export PYTHONPATH=~dstn/tractor:.:${PYTHONPATH}
# I installed this version tractor with
# bashrc to edison and tractor_procution == yes
# git clone tractor repo
# cd tractor; git checkout fe9babf -b dustins_version_dr3; make
export PYTHONPATH=.:/scratch2/scratchdirs/kaylanb/dr3/tractor:${PYTHONPATH}



#brick="$1"
#brick=2501p162
#outdir=$SCRATCH/dr3/runs/$brick

#mkdir -p $outdir/logs/$bri
#log="$outdir/logs/$bri/$brick.log"

#echo Logging to: $log
#echo Running on ${NERSC_HOST} $(hostname)
#
#echo -e "\n\n\n" >> $log
#echo "-----------------------------------------------------------------------------------------" >> $log
#echo "PWD: $(pwd)" >> $log
#echo "Modules:" >> $log
#module list >> $log 2>&1
#echo >> $log
#echo "Environment:" >> $log
#set | grep -v PASS >> $log
#echo >> $log
#ulimit -a >> $log
#echo >> $log
#
#echo -e "\nStarting on ${NERSC_HOST} $(hostname)\n" >> $log
#echo "-----------------------------------------------------------------------------------------" >> $log

#     --no-write \
##########################################

procs=6
export OMP_NUM_THREADS=$procs
module load mpi4py-hpcp
srun -n 4 -N 1 -c $OMP_NUM_THREADS python runbrick_mpi4py.py \
     --brick_list dr3_4bricks.txt --outdir $outdir \
     --jobid $SLURM_JOBID \
     --threads $OMP_NUM_THREADS --zoom 1400 1600 1400 1600



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




