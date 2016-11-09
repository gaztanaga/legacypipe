#!/bin/bash

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
    elif [ -e "$statdir/inq_$brick.txt" ]; then
        echo skipping $brick, its queued
    else
        echo submitting $brick
        touch $statdir/inq_$brick.txt
        sbatch dr3-mzls-bash.sh --export outdir,statdir,brick
    fi
done
