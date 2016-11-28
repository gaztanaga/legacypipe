#!/bin/bash

# Clean up
for i in `find . -maxdepth 1 -type f -name "dr3-mzls-bash.o*"|xargs grep "dr3-mzls-bash DONE"|cut -d ' ' -f 3|grep "^[0-9]*[0-9]$"`; do mv dr3-mzls-bash.o$i /scratch1/scratchdirs/desiproc/data-releases/dr3-mzls-bash/logs/;done

# Cancel ALL jobs and remove all inq.txt files
#for i in `squeue -u desiproc|grep dr3-mzls-b|cut -c10-18`; do scancel $i;done
#rm /scratch1/scratchdirs/desiproc/data-releases/dr3-mzls-bash/progress/*.txt

# Submit jobs
bricklist="$1"
if [ ! -e "$bricklist" ]; then
    echo file=$bricklist does not exist, quitting
    exit 999
fi

export outdir=/scratch1/scratchdirs/desiproc/data-releases/dr3-mzls-bash
export statdir="${outdir}/progress"
mkdir -p $statdir $outdir
for brick in `cat $bricklist`;do
    export brick="$brick"
    bri=$(echo $brick | head -c 3)
    tractor_fits=$outdir/tractor/$bri/tractor-$brick.fits
    if [ -e "$tractor_fits" ]; then
        echo skipping $brick, its done
        # Remove coadd, checkpoint, pickle files they are huge
        rm /scratch1/scratchdirs/desiproc/data-releases/dr3-mzls-bash/pickles/${bri}/runbrick-${brick}*.pickle
        rm /scratch1/scratchdirs/desiproc/data-releases/dr3-mzls-bash/checkpoints/${bri}/${brick}*.pickle
        rm /scratch1/scratchdirs/desiproc/data-releases/dr3-mzls-bash/coadd/${bri}/${brick}/*
    elif [ -e "$statdir/inq_$brick.txt" ]; then
        echo skipping $brick, its queued
    else
        echo submitting $brick
        touch $statdir/inq_$brick.txt
        sbatch dr3-mzls-bash.sh --export outdir,statdir,brick
    fi
done
