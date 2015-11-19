#!/bin/bash -l
#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 06:00:00
#SBATCH -o cosmos.%j.out
#SBATCH -e cosmos.%j.err
#SBATCH -A desi

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=6

export SUBSET=0
export OUTDIR=/global/cscratch1/sd/kaylanb/desi/legacypipe/py/cosmos_s${SUBSET}
srun -n 1 -c $OMP_NUM_THREADS python legacypipe/runcosmos.py --threads $OMP_NUM_THREADS --radec 10:00:28.600 02:12:21.00 --decals-dir ../../legacypipe-dir/ --outdir $OUTDIR --splinesky --pixpsf --no-write --no-sdss --no-wise --force-all --subset $SUBSET

export SUBSET=1
export OUTDIR=/global/cscratch1/sd/kaylanb/desi/legacypipe/py/cosmos_s${SUBSET}
srun -n 1 -c $OMP_NUM_THREADS python legacypipe/runcosmos.py --threads $OMP_NUM_THREADS --radec 10:00:28.600 02:12:21.00 --decals-dir ../../legacypipe-dir/ --outdir $OUTDIR --splinesky --pixpsf --no-write --no-sdss --no-wise --force-all --subset $SUBSET

export SUBSET=2
export OUTDIR=/global/cscratch1/sd/kaylanb/desi/legacypipe/py/cosmos_s${SUBSET}
srun -n 1 -c $OMP_NUM_THREADS python legacypipe/runcosmos.py --threads $OMP_NUM_THREADS --radec 10:00:28.600 02:12:21.00 --decals-dir ../../legacypipe-dir/ --outdir $OUTDIR --splinesky --pixpsf --no-write --no-sdss --no-wise --force-all --subset $SUBSET

export SUBSET=3
export OUTDIR=/global/cscratch1/sd/kaylanb/desi/legacypipe/py/cosmos_s${SUBSET}
srun -n 1 -c $OMP_NUM_THREADS python legacypipe/runcosmos.py --threads $OMP_NUM_THREADS --radec 10:00:28.600 02:12:21.00 --decals-dir ../../legacypipe-dir/ --outdir $OUTDIR --splinesky --pixpsf --no-write --no-sdss --no-wise --force-all --subset $SUBSET

export SUBSET=4
export OUTDIR=/global/cscratch1/sd/kaylanb/desi/legacypipe/py/cosmos_s${SUBSET}
srun -n 1 -c $OMP_NUM_THREADS python legacypipe/runcosmos.py --threads $OMP_NUM_THREADS --radec 10:00:28.600 02:12:21.00 --decals-dir ../../legacypipe-dir/ --outdir $OUTDIR --splinesky --pixpsf --no-write --no-sdss --no-wise --force-all --subset $SUBSET

wait
