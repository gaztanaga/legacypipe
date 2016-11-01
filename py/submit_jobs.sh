#!/bin/bash
export ver="2"
export outdir="/scratch2/scratchdirs/kaylanb/dr3/production/mzls_v${ver}"
#export outdir="testing"
export dr="${outdir}/in_progress"
mkdir -p $dr $outdir
for brick in `cat bricks_mzls_v2v3_50.txt`;do
    export brick="$brick"
    bri=$(echo $brick | head -c 3)
    tractor_fits=$outdir/tractor/$bri/tractor-$brick.fits
    if [ -e "$tractor_fits" ]; 
    then
        echo $brick >> $dr/finished.txt
        echo skipping $brick, its done
    elif [ -e "$dr/inq_$brick.txt" ] && [ ! -e "$dr/ran_$brick.txt" ]
    then
        echo skipping $brick, its queued or running
    elif [ ! -e "$dr/inq_$brick.txt" ] 
    then        
        echo submitting $brick
        sbatch submit_dr4.sh --export brick,ver,outdir,dr
        rm $dr/ran_$brick.txt
        touch $dr/inq_$brick.txt
    else
        echo AHHHHHHH should not reach this point!!!!!!!!!!!!!!!
    fi
done
