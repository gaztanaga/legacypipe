#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:20:00
#SBATCH -o mask2.%j.out
#SBATCH -e mask2.%j.err
#SBATCH -A desi
export ROOTDIR=$SCRATCH/desi/images/ptf/cosmos
export INDIR=mask
export OUTDIR=mask-2
#DW jobdw capacity=10GB access_mode=striped type=scratch
#DW stage_in source=$ROOTDIR/$INDIR destination=$DW_JOB_STRIPED/$INDIR type=directory
#DW stage_out source=$DW_JOB_STRIPED/$OUTDIR destination=$ROOTDIR/$OUTDIR type=directory

cd $SLURM_SUBMIT_DIR
srun -n 1 python rewrite_ptf_masks.py -indir $DW_JOB_STRIPED/$INDIR -search "PTF*.fits" -outdir $DW_JOB_STRIPED/$OUTDIR
